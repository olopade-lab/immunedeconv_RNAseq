import pandas as pd

base_dir = ("/Users/jbreynier/Desktop/Research/Olopade Lab/" +
            "Sheila Paper/immunedeconv_RNAseq/required_genes/")

# Immune genes:

dict_filename = {
    "cibersortx": "LM22.txt",
    "epic": "TRef_sigGenes.csv",
    "mcpcounter": "genes.txt",
    "quantiseq" : "TIL10_signature.txt",
    "timer": "geneMarker_BRCA.csv",
    "xcell": "genes_used_by_xCell.txt"
}

for method in dict_filename.keys():
    if method == "mcpcounter":
        df_required_genes = pd.read_csv(base_dir+method+"/"+dict_filename[method],
                                        header=0,
                                        sep="\t",
                                        usecols=[0, 3],
                                        names=["HGNC_symbol", "Ensembl_ID"],
                                        index_col=False)
    else:
        df_required_genes = pd.read_csv(base_dir+method+"/"+dict_filename[method],
                                        header=0,
                                        sep="\t",
                                        usecols=[0],
                                        names=["HGNC_symbol"],
                                        index_col=False)
    if method == "quantiseq":
        df_remove_genes = pd.read_csv(base_dir+method+"/TIL10_rmgenes.txt",
                                        header=0,
                                        sep="\t",
                                        usecols=[0],
                                        names=["HGNC_symbol"],
                                        index_col=False)
        df_required_genes = pd.concat([df_required_genes, df_remove_genes])
        df_required_genes = df_required_genes.drop_duplicates(keep=False)
    df_required_genes.to_csv(base_dir+f"immune_genes_{method}.csv", index=False)
    

# Example genes from TIMER2.0:

df_example_genes = pd.read_csv(base_dir+"timer2.0/exampleForLUAD.csv",
                                header=0,
                                usecols=[0],
                                names=["HGNC_symbol"],
                                index_col=False)
df_example_genes.to_csv(base_dir+f"example_genes_timer2.csv", index=False)