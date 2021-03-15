import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='''Reformatting script for CIBERSORTx results''')
parser.add_argument('-i', '--input', type=str, required=True, help='''Path of the CIBERSORTx results csv file''')
parser.add_argument('-o', '--output_prefix', type=str, required=True, help='''Path/prefix of the output file''')
args = parser.parse_args()
df_map = pd.read_csv("cell_type_mapping/cibersortx_celltype_map.csv")
df_cibersortx = pd.read_csv(args.input)
df_output = df_cibersortx.transpose()
df_output.rename(columns=df_output.iloc[0], inplace=True)
df_output["cibersortx_names"] = df_output.index
df_output.drop("Mixture", inplace=True)
df_output.reset_index(drop=True, inplace=True)
df_output = pd.merge(df_map, df_output)
df_output.drop("cibersortx_names", axis=1, inplace=True)
df_output.rename(columns={"immunedeconv_names": "cell_type"}, inplace=True)
df_output.to_csv(args.output_prefix+"_reformatted.csv", index=False)
