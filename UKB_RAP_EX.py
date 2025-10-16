#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
UK Biobank Phenotype Extraction Script
Optimized Version - by ChatGPT
"""
import dxpy
import dxdata
import pyspark
import pandas as pd
import re
from distutils.version import LooseVersion

# ============================================================
# 1️⃣ 初始化：加载 Dataset 
# ============================================================
def load_dataset_auto():
    ds_id = dxpy.find_one_data_object(typename='Dataset', name='app*.dataset', folder='/', name_mode='glob')['id']
    print(f"[INFO] Loaded Dataset: {ds_id}")
    return dxdata.load_dataset(id=ds_id)

def get_field_names(participant, field_ids):
    """自动补全字段名"""
    fields = []
    for fid in field_ids:
        fid = fid if fid.startswith("p") else "p" + fid
        fields += participant.find_fields(name_regex=rf'^{fid}(_i\d+)?(_a\d+)?$')
    return ['eid'] + sorted([f.name for f in fields], key=lambda n: LooseVersion(n))

# ============================================================
# 2️⃣ 公共函数：匹配、QC 及通用工具
# ============================================================
def _match_any(codes, targets, mode="prefix"):
    """前缀/精确匹配"""
    if not codes:
        return 0
    if mode == "exact":
        return int(any(str(c) in targets for c in codes))
    
    return int(any(any(str(c).startswith(t) for t in targets) for c in codes))

#- QC 函数
def qc_summary(df):
    """输出QC前的样本统计"""
    print("---- QC前样本统计 ----")
    print(f"Sample size: {len(df)}")
    print(f"Sex mismatch: {(df['p31'] != df['p22001']).sum()}")
    print(f"非白人: {(df['p22006'] != 1).sum()}")
    print(f"Sex chr 异常: {df['p22019'].notnull().sum()}")
    print(f"亲缘关系过多: {(df['p22021'] == 10).sum()}")
    print(f"杂合异常: {df['p22027'].notnull().sum()}")
    print("----------------------")

def perform_qc(df):
    """执行 QC 筛选"""
    return df[
        (df['p31'] == df['p22001']) &
        (df['p22006'] == 1) &
        (df['p22019'].isnull()) &
        (df['p22021'] != 10) &
        (df['p22027'].isnull())]

# ============================================================
# 3️⃣ 一次性加载数据：基础 + 协变量 + ALL_SOURCES + Pheno_SPEC
# ============================================================
def load_data_with_sources(dataset, sources, covariates, pheno_specs=None):
    """
    一次性加载字段：
    - 基础字段
    - 协变量字段
    - 所有 ALL_SOURCES 中的字段
    - 以及疾病定义中额外声明的字段 (extra_specs)
    """
    participant = dataset["participant"]
    print("[INFO] 开始一次性加载字段...")
    
    # ---- 基础字段 ----
    field_ids = ['31', '21022', '22001', '22006', '22009', '22019', '22021', '22027']
    # ---- 添加协变量字段 ----
    field_ids.extend([v.split('_')[0] for v in covariates.values()]) # 支持用户写 21001_i1 形式

    # ---- ALL_SOURCES 的字段 ----
    for src in sources.values():
        if isinstance(src, dict) and "fields" in src:
            src = src["fields"]
        field_ids.extend([f.replace("p", "") for f in src])
    # ---- 自动补充疾病定义中的字段 ----
    extra_fields = []
    if pheno_specs:
        for name, info in pheno_specs.items():
            extra_fields += info.get("fields", [])
        extra_fields = list(set(extra_fields))
    
    # ---- 检查哪些字段是“额外字段”（不在 ALL_SOURCES）----
    all_defined_fields = {f for src in sources.values() for f in (src["fields"] if isinstance(src, dict) and "fields" in src else src)}
    truly_extra_fields = [f for f in extra_fields if f not in all_defined_fields]

    if truly_extra_fields:
        print(f"[INFO] Pheno_SPEC 提供的额外字段有: {truly_extra_fields}")
        field_ids.extend([f.replace("p", "") for f in truly_extra_fields])
    else:
        print("[INFO] Pheno_SPEC 未提供额外字段。")

    # 去重
    field_ids = list(set(field_ids))
    field_names = get_field_names(participant, field_ids)

    # ---- Spark 数据提取 ----
    print(f"[INFO] 正在加载 {len(field_names)} 个字段 (包含基础、协变量、ALL_SOURCES 与 Pheno_SPEC 字段)...")
    engine = dxdata.connect(engine="spark")
    df = participant.retrieve_fields(names=field_names, engine=engine)
    pdf = df.toPandas()
    print(f"[INFO] 数据加载完成，共 {pdf.shape[0]} 行，{pdf.shape[1]} 列")

    return pdf

# ============================================================
# 4️⃣ 疾病匹配：在已加载数据上计算 has_trait, 并过滤无记录的
# ============================================================
def match_sources_on_data(pdf, disease_spec, control_spec, all_sources):
    """
    在已加载数据上匹配疾病来源：
        - 自动识别展开列（p41202_i0_a0/p41204_i1_a3...）
        - 支持 ICD10 前缀匹配与整数 code 精确匹配
        - 自动合并多来源 has_trait
        - 生成 has_trait (病例) 与 is_control (对照)
    """
    def _as_list(v):
        if v is None: return []
        if isinstance(v, (list, tuple, set)): return list(v)
        return [v]

    # ------------------------------
    # (1) 匹配病例 Pheno_SPEC
    # ------------------------------
    for name, info in disease_spec.items():
        base_fields = info.get("fields", all_sources.get(name, []))

        # ==== 收集所有匹配列（支持展开列名） ====
        available = []
        for fid in base_fields:
            cols = [c for c in pdf.columns if c == fid or c.startswith(fid + "_")]
            available.extend(cols)
        available = sorted(set(available))

        if not available:
            pdf[f"has_{name}"] = 0
            print(f"[WARN] {name}: 未找到任何字段，跳过。")
            continue

        # 无 codes → 二元变量
        codes = info.get("codes", set())
        mode = info.get("mode", "prefix")

        if not codes:
            pdf[f"has_{name}"] = 0
            for col in available:
                pdf[f"has_{name}"] |= pdf[col].apply(lambda x: int(x == 1))
            print(f"[INFO] {name}: 没有 codes，按二元字段 (==1) 处理。")
            continue

        # ==== 正常匹配 ====
        tmp_cols = []
        for col in available:
            col_flag = f"{name}_{col}_hit"
            pdf[col_flag] = pdf[col].apply(lambda x: _match_any(_as_list(x), targets=codes, mode=mode))
            tmp_cols.append(col_flag)

        pdf[f"has_{name}"] = pdf[tmp_cols].max(axis=1)
        print(f"[INFO] Finished matching {name}: {len(available)} fields used.")
    # ---- 合并多来源 ----
    src_cols = [f"has_{s}" for s in disease_spec.keys()]
    pdf["has_trait"] = pdf[src_cols].max(axis=1)

    print("--------------------------------------------------")
    print("[INFO] 各来源病例统计：")
    for s in disease_spec.keys():
        colname = f"has_{s}"
        if colname in pdf.columns:
            n_case = pdf[colname].sum()
            print(f"  - {s}: {n_case} 个病例")
        else:
            print(f"  - {s}: 未找到列（跳过）")
    print(f"[INFO] 已生成病例标签列 has_trait ({pdf['has_trait'].sum()} 个病例)")

    # ------------------------------
    # (2) 匹配对照 Control_SPEC
    # ------------------------------
    if control_spec:
        pdf["is_control"] = 1
        for name, info in control_spec.items():
            base_fields = info.get("fields", all_sources.get(name.replace("_exclude", ""), []))
            available = []
            for fid in base_fields:
                cols = [c for c in pdf.columns if c == fid or c.startswith(fid + "_")]
                available.extend(cols)
            available = sorted(set(available))
            if not available:
                print(f"[WARN] Control {name}: 未找到字段，跳过。")
                continue
            codes = info.get("codes", set())
            mode = info.get("mode", "prefix")
            tmp_hit = pd.Series(0, index=pdf.index)
            for col in available:
                tmp_hit |= pdf[col].apply(lambda x: _match_any(_as_list(x), targets=codes, mode=mode))
            pdf.loc[tmp_hit > 0, "is_control"] = 0
            print(f"[INFO] Control {name}: 过滤掉 {tmp_hit.sum()} 个样本")
        # 对照逻辑优先：非病例+未命中排除条件的才是对照
        pdf["is_control"] = ((pdf["has_trait"] == 0) & (pdf["is_control"] == 1)).astype(int)
    else:
        pdf["is_control"] = (pdf["has_trait"] == 0).astype(int)
        print(f"[INFO] 未提供 Control_SPEC，默认非病例为对照。")
    
    print(f"[INFO] 对照样本数: {pdf['is_control'].sum()}")

    # ============================================================
    # (3) 仅保留病例（has_trait==1）和对照（is_control==1）
    # ============================================================
    if "has_trait" in pdf.columns and "is_control" in pdf.columns:
        before_rows = len(pdf)
        pdf = pdf[(pdf["has_trait"] == 1) | (pdf["is_control"] == 1)]
        after_rows = len(pdf)
        print(f"[INFO] 保留病例(has_trait=1)和对照(is_control=1)样本，共 {after_rows} 个 (移除 {before_rows - after_rows} 个非病例/非对照样本)")

    return pdf

def filter_participants_by_sources(pdf, pheno_spec, all_sources):
    """
    根据给定表型定义 (Pheno_SPEC) 自动过滤出
    至少在一个来源(字段)中有记录的参与者。
    """
    def _has_data(v):
        if v is None: return False
        if isinstance(v, (list, tuple, set)): return len(v) > 0
        if isinstance(v, float) and pd.isna(v): return False
        return True
    
    keep_mask = pd.Series(False, index=pdf.index)
    
    for name, info in pheno_spec.items():
        base_fields = info.get("fields", all_sources.get(name, []))
        cols = []
        for fid in base_fields:
            cols += [c for c in pdf.columns if c == fid or c.startswith(fid + "_")]
        if not cols:
            print(f"[WARN] {name}: 未找到任何字段，跳过。")
            continue
        
        mask = pdf[cols].map(_has_data).any(axis=1)
        keep_mask |= mask
        print(f"[INFO] {name}: {mask.sum()} participants 有 {len(cols)} 列非空记录。")
    
    pdf_filtered = pdf[keep_mask].copy()
    print(f"[INFO] 共筛选出 {len(pdf_filtered)} participants (有至少一个来源数据)。")
    return pdf_filtered


# ============================================================
# 5️⃣ 生成表型文件：匹配 + QC + 重命名 + 输出 + 上传
# ============================================================
def generate_pheno(pdf, sources, covariates, pheno_spec, out_prefix, upload_path, upload=True, control_spec=None, n_pc=10):
    """
    在已加载数据上生成表型文件：
        1. 匹配疾病表型
        2. QC 筛选
        3. 允许使用额外字段定义对照样本
        3. 输出 .phe 文件
        4. 可选上传到 DNAnexus
    """
    print("="*60)
    print(f"[INFO] 开始生成表型：{out_prefix}")
    print("="*60)

    # ---- 匹配表型 ----
    pdf_matched = match_sources_on_data(pdf.copy(), pheno_spec, control_spec, sources)

    # 协变量命名映射：p{field} → covariate 名
    for cov_name, cov_field in covariates.items():
        # 拆分 base ID 和 instance
        m = re.match(r"^(\d+)(?:_i(\d+))?$", cov_field)
        if m:
            base_field = m.group(1)
            specified_instance = m.group(2)
        else:
            base_field = cov_field
            specified_instance = None

        pcol_prefix = f"p{base_field}"
        matched_cols = [c for c in pdf_matched.columns if c == pcol_prefix or c.startswith(pcol_prefix + "_")]

        if not matched_cols:
            print(f"[WARN] 协变量 {cov_name}: 未找到任何匹配列 ({pcol_prefix})")
            continue

        if specified_instance is not None:
            target_pat = f"{pcol_prefix}_i{specified_instance}"
            instance_cols = [c for c in matched_cols if c.startswith(target_pat)]
            if instance_cols:
                src = instance_cols[0]
                print(f"[INFO] 协变量 {cov_name}: 使用指定实例 {src}")
            else:
                src = matched_cols[0]
                print(f"[WARN] 协变量 {cov_name}: 指定期号 _i{specified_instance} 不存在，退回 {src}")
        else:
            i0_cols = [c for c in matched_cols if "_i0" in c]
            if i0_cols:
                src = i0_cols[0]
                print(f"[INFO] 协变量 {cov_name}: 未指定实例，默认使用 {src}")
            else:
                src = matched_cols[0]
                print(f"[INFO] 协变量 {cov_name}: 字段无 instance 展开，使用 {src}")

        pdf_matched = pdf_matched.assign(**{cov_name: pdf_matched[src]})
    # ---- QC ----
    qc_summary(pdf_matched)
    pdf_qced = perform_qc(pdf_matched)
    print(f"[INFO] QC后样本: {pdf_qced.shape[0]}")
    
    pdf_qced = pdf_qced.rename(columns=lambda x: re.sub('p22009_a','pc',x))
    base_map = {
        'eid': 'IID',
        'p31': 'sex',
        'p21022': 'age',
        'p22006': 'ethnic_group',
        'p22019': 'sex_chr_aneuploidy',
        'p22021': 'kinship_to_other',
        'p22027': 'outliers_for_heterozygosity_or_missing'}
    pdf_qced.rename(columns={k: v for k, v in base_map.items() if k in pdf_qced.columns}, inplace=True)

    # FID = IID
    if "IID" in pdf_qced.columns:
        pdf_qced["FID"] = pdf_qced["IID"]

    # ---- 选择输出列 ----
    base_cols = ['FID', 'IID', 'sex', 'age', 'has_trait']
    source_cols = [f"has_{name}" for name in pheno_spec.keys() if f"has_{name}" in pdf_qced.columns]
    cov_cols = [cov for cov in covariates.keys() if cov in pdf_qced.columns]
    pc_cols = [f'pc{i}' for i in range(1, n_pc + 1) if f'pc{i}' in pdf_qced.columns]

    out_cols = base_cols + source_cols + cov_cols + pc_cols
    missing = [c for c in out_cols if c not in pdf_qced.columns]
    if missing:
        print(f"[WARN] 以下输出列不存在，将被跳过: {missing}")
        out_cols = [c for c in out_cols if c in pdf_qced.columns]

    pdf_pheno = pdf_qced[out_cols]

    # ---- 清洗表型数据：删除含 NA 的样本（未被测评的个体）----
    rawrow = len(pdf_pheno)
    pdf_pheno = pdf_pheno.dropna()
    newrow = len(pdf_pheno)
    print(f"[INFO] 清洗数据：已删除含 NA 的样本 {rawrow - newrow} 个，剩余 {newrow} 个样本")
    
    # ---- 输出文件 ----
    outfile = f"{out_prefix}_combined.phe"
    pdf_pheno.to_csv(outfile, sep="\t", index=False)
    print(f"[INFO] 已保存表型文件: {outfile}")

    # ---- 统计病例数 ----
    affected = (pdf_pheno['has_trait'] > 0).sum()
    print(f"[INFO] {affected} participants "
          f"({affected / len(pdf_pheno) * 100:.1f}%) 为病例。")

    # ---- 上传 ----
    if upload:
        import subprocess
        try:
            subprocess.run(['dx', 'upload', outfile, '-p', '--path', upload_path, '--brief'], check=True)
            print(f"[UPLOAD OK] {outfile} → {upload_path}")
        except Exception as e:
            print(f"[UPLOAD FAIL] {outfile}: {e}")

    print("="*60)
    print(f"[INFO] 表型 {out_prefix} 生成完毕 ✅")
    print("="*60)
    return pdf_qced

# ============================================================
# 6️⃣ 用户配置与执行逻辑
# ============================================================
ALL_SOURCES = {
    # ICD10 主/次诊断 & 死亡记录
    "ICD10": ["p41202", "p41204", "p40001", "p40002"],
    # self reported non-cancer illness
    "SELF":  ["p20002"],
    # Mental health questionnaire
    "MENTAL": ["p20544"],
    # father / mother disease (family history)
    "PROXY": ["p20107", "p20110"]
}

# 可设置使用哪一 instance 期的数据
COVARIATES = {
    "BMI": "21001_i0",
    "alcohol_frequency": "1558_i0",
    "smoking_status": "20116_i0",
    "assessment_centre": "54_i0",
    "moderate_activity": "884_i0",
    "vigorous_activity": "904_i0"
}

# ---- 执行主流程 ----
dataset = load_dataset_auto()
pdf_all = load_data_with_sources(dataset, ALL_SOURCES, COVARIATES, pheno_specs=None)  # ✅ 一次加载


OUTPUT_PREFIX = "ADHD"  # 输出文件前缀
DX_UPLOAD_PATH = "/PRS_NPDs/"   # 上传路径
N_PC = 0

Pheno_SPEC = {
    "ICD10": {"codes": {"F90"}, "mode": "prefix"},
    "MENTAL": {"codes": {"18"}, "mode": "exact"} 
}

#- select control
# CONTROL_FILTERS = {
#     "ICD10_exclude": {"codes": {"F03"}, "mode": "prefix"},
#     "SELF_exclude":  {"codes": {"1263"}, "mode": "exact"},
#     "PROXY_exclude": {"codes": {"10"}, "mode": "exact"}
# }

# 根据表型来源筛选出有记录的人
pdf_filtered = filter_participants_by_sources(pdf_all, Pheno_SPEC, ALL_SOURCES)
generate_pheno(pdf_filtered, ALL_SOURCES, COVARIATES, Pheno_SPEC, OUTPUT_PREFIX, upload_path=DX_UPLOAD_PATH,  upload=True, n_pc=N_PC, control_spec=None)