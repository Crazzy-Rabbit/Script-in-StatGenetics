#! usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@File    :    GWASformatRSID.py
@Time    :    2025/02/13 11:35:39
@Author  :    Lulu Shi
@Mails   :    crazzy_rabbit@163.com
@line    :    https://github.com/Crazzy-Rabbit
'''
import re
import os
import math
import click
import logging
import numpy as np
import pandas as pd
import scipy.stats as stats

def get_logger(logger_name):
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG)
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter("[{asctime}] {levelname:.5s} | {name} - {message}", style="{"))
    logger.addHandler(handler)

    return logger

logger = get_logger(__name__)

@click.command('data_reformat')
@click.option('--gwas', help='gwas summary data', required=True)
@click.option('--chrom', help='Chromosome columns', default="CHR")
@click.option('--pos', help='SNP positions columns', default="POS")
@click.option('--snp', help='SNP id columns', default="SNP")
@click.option('--a1', help='Effect allele', default="A1")
@click.option('--a2', help='Non-effect allele', default="A2")
@click.option('--frq', help='Effect allele frequency', default="freq")
@click.option('--pval', help='Pvalue of the regression coefficient', default="p")
@click.option('--num', help='Sample size', default="N")
@click.option('--beta', help='Effect allele beta size', default=None)
@click.option('--se', help='Standard error of the regression coefficient', default=None)
@click.option('--OR', help='Odds ratio, will be transferred to linear scale', default=None)
@click.option('--zscore', help='zscore', default=None)
@click.option('--format', help='Format of output data', type=click.Choice(["SMR", "popcorn"]), case_sensitive=False, default="SMR")
@click.option('--outname', help='Outfile prefix', default=None)
def data_reformat(gwas, chrom, pos, snp, a1, a2, frq, pval, num, beta, se, OR, zscore, format, outname):
    if outname is None:
        outname = os.path.splitext(os.path.basename(gwas))[0] + "_reformat.txt"
        logger.info(f"Not provide outname, out prefix is set as {outname}.")
    gwas = pd.read_csv(gwas, sep="\t")
    old_name = gwas.columns

    name_update = {
        "CHR": chrom,
        "POS": pos,
        "SNP": snp,
        "A1": a1,
        "A2": a2,
        "b": beta,
        "se": se,
        "p": pval,
        "freq": frq,
        "N": num,
        "OR": OR,
        "z": zscore,
    }
    
    for key, value in name_update.items():
        if value is not None and value in gwas.columns:
            gwas = gwas.rename(columns={value: key})
        else:
            raise ValueError(f"Please check the columns you provide whether match your gwas columns!")
    new_name = gwas.columns  
    name_dict = {new_name[i]: old_name[i] for i in range(len(new_name))}
    
    # ensure A1, A2 is upper
    logger.info("\nEnsure A1, A2 is upper")
    gwas['A1'] = gwas['A1'].apply(lambda x: x.upper())
    gwas['A2'] = gwas['A2'].apply(lambda x: x.upper())
    # ensure SNP p <= 1, freq <= 1
    logger.info("\nEnsure SNP p <= 1, freq <= 1")
    gwas = gwas[gwas['p'] <= 1]
    gwas = gwas[gwas['freq'] <= 1]

    # gwas only contain OR p
    if "OR" in new_name and "p" in new_name:
        logger.info("\nTransfer OR to beta se, when beta = 0, set se as mean se of all SNP")
        gwas["b"] = gwas.OR.apply(lambda x: math.log(x) if x > 0 else None)     
        se_mean = (gwas[gwas['beta'] != 0]
                   .assign(se = np.sqrt(gwas['beta']**2 / stats.chi2.ppf(gwas['p'], 1)))
                   .mean()['se'])
        gwas['se'] = np.where(gwas['beta'] == 0, 
                              se_mean, 
                              np.sqrt(gwas['beta']**2 / stats.chi2.ppf(gwas['p'], 1)))
    
    interpreting = {
        "SNP": "Variant ID (e.g., rs number).",
        "A1": "Allele 1, interpreted as the effect allele for signed sumstat.",
        "A2": "Allele 2, interpreted as the non-effect allele for signed sumstat.",
        "b": "[linear/logistic] regression coefficient (0 → no effect; above 0 → A1 is trait/risk increasing).",
        "se": "Standard error of the regression coefficient.",
        "OR": "Odds ratio, will be transferred to linear scale.",
        "p": "P-Value.",
        "z": "Z-Value.",
        "N": "Sample size.",
        "freq": "Allele frequency of A1.",
        "CHR": "Chromsome.",
        "POS": "SNP positions.",
    }

    logger.info("\nIterpreting column names as follows:")
    for key, _value in interpreting.items():
        if key in new_name:
            logger.info(f"{name_dict[key]}: {interpreting[key]}")

    logger.info(f"Format the gwas data.")
    if format == "SMR":
        logger.info(f"Save gwas as format SMR: SNP\tA1\tA2\tfreq\tb\tse\tp\tN.")
        gwas[['SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N']].to_csv(outname, sep="\t", index=False, header=True)
    elif format == "popcorn":
        logger.info(f"Save gwas as format popcorn: SNP\tA1\tA2\taf\tbeta\tSE\tp\tN.")
        gwas = gwas.rename(columns={"freq": "af", "b": "beta", "se": "SE"})
        gwas[['SNP', 'A1', 'A2', 'af', 'beta', 'SE', 'p', 'N']].to_csv(outname, sep="\t", index=False, header=True)


def common_options(f):
    f = click.option('--gwas', help='gwas summary data, header as order: SNP A1 A2 freq b se p N', required=True)
    f = click.option('--format', help='Format of output data', type=click.Choice(["SMR", "popcorn"]), case_sensitive=False, default="SMR")
    f = click.option('--outname', help='Outfile prefix', default=None)
    return f

@click.command('variant_to_rsid')
@common_options()
@click.option('--dbsnp', help='Path to reference dnsnp file', default=None)
@click.option('--chunksize', help='Chunk size for loading dbsnp file', type=int, default=1000000)
def variant_to_rsid(gwas, dbsnp, chunksize, format, outname):
    """
    Convert variant id (Chr, Pos) to rsid using a dbSNP reference file.
    """
    logger.info("\nConverting the SNP position to rsid. This process may take some time.")
    if outname is None:
        outname = os.path.splitext(os.path.basename(gwas))[0] + "_rsid"
        logger.info(f"Not provide outname, out prefix is set as {outname}.")
        
    gwas = pd.read_csv(gwas, sep="\t")
    gwas.columns = ["SNP","A1","A2","freq","b","se","P","N"]
    # change gwas snpid to "chr_pos"
    gwas["id"] = gwas["CHR"].astype(str) + "_" + gwas["POS"].astype(str)
    gwas = gwas.drop_duplicates(subset="id").reset_index(drop=True)
    gwas.index = gwas.id

    unique_ids = set(gwas["id"])
    chr_format = gwas["Chr"].unique().astype(str)
    chr_format = [re.sub(r"\d+", "", value) for value in chr_format][1]

    # read dbSNP
    dtype = {"chr": str, "pos": str, "dbsnp": str}
    chunk_iter = pd.read_csv(
        dbsnp,
        chunksize=chunksize,
        sep="\t",
        header=None,
        dtype=dtype,
        names=["chr", "pos", "dbsnp"],
    )

    # Iterate over chunks
    matching_id = pd.DataFrame()
    for chunk in chunk_iter:
        # set dbsnp snpid as "chr_pos"
        chunk["id"] = chr_format + chunk["chr"] + "_" + chunk["pos"]
        matching_id = pd.concat(
            [matching_id, chunk[chunk["id"].isin(unique_ids)][["dbsnp", "id"]]]
        )

    # 避免重复匹配导致数据丢失
    matching_id = matching_id.drop_duplicates(subset=["dbsnp", "id"]).reset_index(drop=True)
    matching_id.index = matching_id.id
    
    # matching rsid
    gwas = gwas.merge(matching_id, on="id", how="left").reset_index(drop=True)
    num_fail = gwas["dbsnp"].isna().sum()
    logger.info(f"{num_fail} SNPs that did not convert to rsid, set these as CHR_POS_A1_A2")
    gwas["dbsnp"] = gwas["dbsnp"].fillna(gwas["id"] + "_" + gwas["A1"] + "_" + gwas["A2"])

    gwas["SNP"] = gwas.dbsnp

    # save the final gwas
    logger.info(f"\nWriting summary statistics for {len(gwas)} SNPs to {outname}.")
    if format == "SMR":
        logger.info(f"Save gwas as format SMR: SNP\tA1\tA2\tfreq\tb\tse\tp\tN.")
        gwas[['SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N']].to_csv(outname, sep="\t", index=False, header=True)
    elif format == "popcorn":
        logger.info(f"Save gwas as format popcorn: SNP\tA1\tA2\taf\tbeta\tSE\tp\tN.")
        gwas = gwas.rename(columns={"freq": "af", "b": "beta", "se": "SE"})
        gwas[['SNP', 'A1', 'A2', 'af', 'beta', 'SE', 'p', 'N']].to_csv(outname, sep="\t", index=False, header=True)


@click.command('effect_sample_size')
@common_options()
@click.option('--pop-prev', help='population prevalence of the trait', type=int, required=True)
@click.option('--sample-prev', help='sample prevalence of the trait, cal as case / (case + control)', type=int, required=True)
def effect_sample_size(gwas, pop_prev, sample_prev, outname):
    logger.info("\nCalculate the effect sample size and rescale the beta and se based on effect N.")
    if outname is None:
        outname = os.path.splitext(os.path.basename(gwas))[0] + "_eff_N.txt"
        logger.info(f"Not provide outname, out prefix is set as {outname}.")
    gwas = pd.read_csv(gwas, sep="\t")
    gwas.columns = ["SNP","A1","A2","freq","b","se","P","N"]

    K = pop_prev
    v = sample_prev
    i = stats.norm.pdf(stats.norm.ppf(1 - K)) / K
    n = gwas['N']
    n_eff = i**2 * v * (1 - v) / (1 - K)**2 * n # Yang et al. Gen Epi 2010
    z = gwas['b'] / gwas['se']
    p = gwas['freq']
    
    se = 1 / np.sqrt(2 * p * (1 - p) * (n_eff + z**2))
    b = z * se

    gwas['b'] = b
    gwas['se'] = se
    gwas['N'] = n_eff

    # save the final gwas
    if format == "SMR":
        logger.info(f"Save gwas as format SMR: SNP\tA1\tA2\tfreq\tb\tse\tp\tN.")
        gwas[['SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N']].to_csv(outname, sep="\t", index=False, header=True)
    elif format == "popcorn":
        logger.info(f"Save gwas as format popcorn: SNP\tA1\tA2\taf\tbeta\tSE\tp\tN.")
        gwas = gwas.rename(columns={"freq": "af", "b": "beta", "se": "SE"})
        gwas[['SNP', 'A1', 'A2', 'af', 'beta', 'SE', 'p', 'N']].to_csv(outname, sep="\t", index=False, header=True)
    else:
        raise ValueError(f"The progress only support the format of SMR and popcorn! ! !")


# Create the command group
@click.group()
def cli():
    "GWAS data reformat and variant to rsid" 
    pass
# Add the commands to the group
cli.add_command(data_reformat)
cli.add_command(variant_to_rsid)  
cli.add_command(effect_sample_size)

if __name__ == "__main__":
    cli()
