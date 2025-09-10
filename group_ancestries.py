#!/usr/bin/env python3

import pandas as pd

# Read in combined data table
df = pd.read_csv("lrrk2_data.tsv", low_memory=False, sep="\t")

# Aggregate hom, het, and total counts across ancestries for each variant
cols_to_sum = ["het_PD", "hom_PD", "total_PD", "het_HC", "hom_HC", "total_HC"]
agg_dict = {col: "sum" for col in cols_to_sum}
for col in df.columns:
    if col not in cols_to_sum + ["variant"]:
        agg_dict[col] = "first"
df_grouped = df.groupby("variant", as_index=False).agg(agg_dict)

# Save aggregated table
df_grouped.to_csv("lrrk2_grouped.tsv", index=False, sep="\t")
