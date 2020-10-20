import pandas as pd
import numpy as np
import re
from ExpressionSearch import ExpressionSearch

def parse_probe_string(probe_string):
    re_match = re.search(r"\(\w+\)", probe_string)
    owner = re_match.group(0)[1:-1]
    owner_start = re_match.start(0) - 1
    gene_symbol = probe_string[:owner_start].strip().lower()
    if "*" in probe_string:
        homemade = True
    else:
        homemade = False

    return [gene_symbol, owner, homemade]

def load_hcr_spreadsheet(path):
    hcr_probes = pd.read_csv(path, index_col=None)
    del hcr_probes["Unnamed: 5"]
    hcr_probes = hcr_probes.T
    probe_data = []
    for idx in hcr_probes.index:
        initiator_probes = hcr_probes.loc[idx].dropna()
        for probe in initiator_probes:
            probe_data.append([idx] + parse_probe_string(probe))
    probe_data = pd.DataFrame(probe_data, columns=["Initiator", "Gene", "Owner", "Homemade"])
    return probe_data

searcher = ExpressionSearch.ExpressionSearcher()

hcr_probes = load_hcr_spreadsheet("hcr_probes.csv")

tissues = ["heart primordium", "heart rudiment", "heart tube"]

genes, genes_idx = searcher.filter_genes_by_tissue(hcr_probes["Gene"], tissues)
filtered_probes = hcr_probes.iloc[genes_idx, :]
genes["Initiator"] = list(filtered_probes["Initiator"])
genes["Owner"] = list(filtered_probes["Owner"])
genes["Homemade"] = list(filtered_probes["Homemade"])
genes = genes.fillna(0)
genes["PMID"] = genes["PMID"].astype(np.uint16)
genes["Year"] = genes["Year"].astype(np.uint16)
genes.to_csv("results.csv")