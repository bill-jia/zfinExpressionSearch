import pandas as pd
import os
import numpy as np

class ExpressionSearcher():
    def __init__(self, update_records=False):        
        dirname = os.path.dirname(__file__)
        self.expression = pd.read_csv(os.path.join(dirname, "wildtype-expression_fish.txt"), sep="\t", index_col=None)
        columns = [s.rstrip() for s in self.expression.columns]
        self.expression.columns = columns
        self.expression = self.expression.set_index("Super Structure Name")
        self.publications = pd.read_csv(os.path.join(dirname, "zfinpubs.txt"), sep="\t", index_col=None)
        columns = [s.rstrip() for s in self.publications.columns]
        self.publications.columns = columns
        self.publications = self.publications.set_index("Publication ID")
        self.anatomy = self._build_anatomy()
        if update_records:
            # pull data from zfin
            print("records updated")

    def filter_genes_by_tissue(self, genes, tissues):
        tissue_genes = []
        gene_idx = []
        for tissue in tissues:
            if tissue in self.anatomy:
                tissue_expression = self.expression.loc[self.anatomy[tissue]]
                for idx, gene in enumerate(genes):
                    gene_entry_indices = [i for i, g in enumerate(list(tissue_expression["Gene Symbol"])) if g == gene]
                    for i in gene_entry_indices:
                        try:
                            pub = tissue_expression.iloc[i]["Publication ID"]
                            pub_record = list(self.publications.loc[pub])[:-3]
                            pub_record[0] = int(pub_record[0])
                        except:
                            pub_record = [None]*5
                        tissue_genes.append([tissue, gene, tissue_expression.iloc[i]["Start Stage"], tissue_expression.iloc[i]["End Stage"]] + pub_record)
                        gene_idx.append(idx)
            else:
                print("%s does not exist in ZFIN" % tissue)
        
        tissue_genes = pd.DataFrame(tissue_genes, columns=["Tissue", "Gene Symbol", "Start Stage", "End Stage", "PMID", "Authors", "Title", "Journal", "Year"])
        return tissue_genes, gene_idx
    
    def _build_anatomy(self):
        dirname = os.path.dirname(__file__)
        anatomy_list =  pd.read_csv(os.path.join(dirname, "anatomy_synonyms.txt"), sep="\t", index_col=None)
        columns = [s.rstrip() for s in anatomy_list.columns]
        anatomy_list.columns = columns
        anatomy_lookup = {}
        for _, row in anatomy_list.iterrows():
            anatomy_lookup[row["Anatomy Synonym"]] = row["Anatomy Name"]
        anatomy_list =  pd.read_csv(os.path.join(dirname, "anatomy_item.txt"), sep="\t", index_col=None)
        columns = [s.rstrip() for s in anatomy_list.columns]
        anatomy_list.columns = columns
        for name in anatomy_list["Anatomy Name"]:
            anatomy_lookup[name] = name
        return anatomy_lookup