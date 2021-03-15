import pandas as pd
import numpy as np

list_methods = ["timer", "quantiseq", "xcell", "mcpcounter", "cibersortx", "epic"]
df_biomart = pd.read_csv("HGNC_biomart/biomart_database.csv")

for method in list_methods:
    file_name = f"required_genes/immune_genes_{method}.csv"
    df_immune_genes = pd.read_csv(file_name)
    df_biomart_specific = df_biomart
    if method == "mcpcounter":
        df_merged = df_immune_genes
    else:
        df_merged = pd.merge(df_immune_genes,
                                df_biomart_specific.drop_duplicates(subset=["HGNC_approved", "Ensembl_ID"]),
                                left_on="HGNC_symbol",  right_on="HGNC_approved", how="left")
        df_merged_missing_approved = df_merged[df_merged["HGNC_approved"].isnull()][["HGNC_symbol"]]
        df_merged.dropna(subset=["HGNC_approved"], inplace=True)
        df_biomart_specific = df_biomart_specific[~df_biomart_specific["Ensembl_ID"].isin(df_merged["Ensembl_ID"])]
        df_merged_previous = pd.merge(df_merged_missing_approved,
                                        df_biomart_specific.drop_duplicates(subset=["HGNC_previous", "Ensembl_ID"]),
                                        left_on="HGNC_symbol", right_on="HGNC_previous", how="left")
        df_merged_missing_previous = df_merged_previous[df_merged_previous["HGNC_previous"].isnull()][["HGNC_symbol"]]
        df_merged = pd.concat([df_merged, df_merged_previous])
        df_merged.dropna(subset=["Ensembl_ID"], inplace=True)
        df_biomart_specific = df_biomart_specific[~df_biomart_specific["Ensembl_ID"].isin(df_merged["Ensembl_ID"])]
        df_merged_alias = pd.merge(df_merged_missing_previous,
                                        df_biomart_specific.drop_duplicates(subset=["HGNC_alias", "Ensembl_ID"]),
                                        left_on="HGNC_symbol", right_on="HGNC_alias", how="left")
        df_merged_missing_alias = df_merged_alias[df_merged_alias["HGNC_alias"].isnull()][["HGNC_symbol"]]
        df_merged = pd.concat([df_merged, df_merged_alias])
        if method == "xcell":
            df_merged = df_merged[(df_merged["Ensembl_ID"] != "ENSG00000136628") & 
                                            (df_merged["HGNC_symbol"] != "QARS")]
    list_missing = df_merged[df_merged["Ensembl_ID"].isna()]["HGNC_symbol"].to_list()
    if len(list_missing) > 0:
        print(f"Missing genes for {method}:")
        print(*list_missing, sep = ", ")
        print("")
    df_merged.dropna(subset=["Ensembl_ID"], inplace=True)
    df_biomart_specific = df_biomart_specific[~df_biomart_specific["Ensembl_ID"].isin(df_merged["Ensembl_ID"])]
    df_remaining = df_biomart_specific.loc[:,["Ensembl_ID", "HGNC_approved"]]
    df_remaining.drop_duplicates(inplace=True)
    df_remaining.drop_duplicates(subset=["Ensembl_ID"], inplace=True, keep=False)
    df_remaining.columns = ["Ensembl_ID", "HGNC_symbol"]
    df_merged = pd.concat([df_merged.loc[:, ["Ensembl_ID", "HGNC_symbol"]], df_remaining])
    df_merged.to_csv(f"conversion_tables/{method}_conversion.csv", index=False)