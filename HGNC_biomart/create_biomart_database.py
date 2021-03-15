import pandas as pd

df_biomart = pd.read_csv("results.txt", sep="\t",
                            header=0, names=["Ensembl_ID", "HGNC_previous", "HGNC_alias", "HGNC_approved"])
df_biomart.dropna(subset=["Ensembl_ID"], inplace=True)
df_biomart = df_biomart[["Ensembl_ID", "HGNC_approved", "HGNC_previous", "HGNC_alias"]]
df_biomart.to_csv("biomart_database.csv", index=False)