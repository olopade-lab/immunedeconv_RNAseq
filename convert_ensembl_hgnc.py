import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='''Ensembl ID to HGNC symbol conversion script for 
                                                cell-type deconvolution''')
parser.add_argument('-i', '--input', type=str, required=True, help='''Path of the input csv expression file
                                                                    (1st column contains gene names)''')
parser.add_argument('-o', '--output_prefix', type=str, required=True, help='''Path/prefix of the output file''')
parser.add_argument('-p', '--program', type=str, nargs='+', required=False, 
                    default=["timer", "quantiseq", "xcell", "mcpcounter", "cibersortx", "epic"],
                    choices=["timer", "quantiseq", "xcell", "mcpcounter", "cibersortx", "epic"],
                    help='''List of programs for which to do conversion''')
args = parser.parse_args()

df_input = pd.read_csv(args.input)
ensembl_colname = df_input.columns[0]
for method in args.program:
    conversion_table = pd.read_csv(f"conversion_tables/{method}_conversion.csv")
    df_output = pd.merge(conversion_table, df_input, left_on="Ensembl_ID", right_on=ensembl_colname, how="left")
    df_output.drop(["Ensembl_ID", ensembl_colname], axis=1, inplace=True)
    df_output.rename(columns = {'HGNC_symbol':ensembl_colname}, inplace=True)
    df_output.dropna(inplace=True)
    if method == "cibersortx":
        df_output.to_csv(f"{args.output_prefix}_{method}.txt", sep="\t", index=False)
    else:
        df_output.to_csv(f"{args.output_prefix}_{method}.csv", index=False)
    df_required_genes = pd.read_csv(f"required_genes/immune_genes_{method}.csv")
    df_required_genes = pd.merge(df_output, df_required_genes, left_on=ensembl_colname, right_on="HGNC_symbol", how="right")
    list_missing = df_required_genes[df_required_genes[ensembl_colname].isna()]["HGNC_symbol"].to_list()
    if len(list_missing) > 0:
        print(f"The following immune HGNC symbols for {method} "+
                "do not have an Ensembl ID equivalent in this input file:")
        print(*list_missing, sep = ", ")
        print("")
