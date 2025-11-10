#!/usr/bin/env python3

import pandas as pd

# Read in combined data table
df_orig = pd.read_csv("lrrk2_data.tsv", low_memory=False, sep="\t")
ancestry_groups = dict(tuple(df_orig.groupby("Ancestry")))

# Find all variants present across all ancestries
df_variants = df_orig[["variant"]].drop_duplicates(subset=["variant"]).copy()
df_variants["pos"] = df_variants["variant"].str.split(":").str[1].astype(int)
df_variants.sort_values(by="pos", ascending=True, inplace=True)
df_variants = df_variants[["variant"]].reset_index(drop=True)

# Include all variants, filling in those not present in a given ancestry
count_cols = [
    "het_PD", "hom_PD", "total_PD", "het_HC", "hom_HC", "total_HC", 
    "FH_neg_PD", "FH_pos_PD", "FH_unk_PD", "FH_neg_HC", "FH_pos_HC", "FH_unk_HC",
    "Age_10", "Age_20", "Age_30", "Age_40", "Age_50", "Age_60",
    "Age_70", "Age_80", "Age_90", "Age_100",
]
df_final = []
for anc, df_ancestry in ancestry_groups.items():
    df_counts = df_ancestry[count_cols + ["variant", "Age_min", "Age_max", "Age_median"]].copy()
    df_merged = df_variants.merge(df_counts, on="variant", how="left")
    df_merged[["total_PD", "total_HC"]] = df_merged[["total_PD", "total_HC"]].apply(
        lambda col: col.fillna(col.max())
    )
    df_merged[count_cols] = df_merged[count_cols].fillna(0)
    df_merged["Ancestry"] = anc
    df_final.append(df_merged)
df_final = pd.concat(df_final)

# Aggregate counts and age min/max/median across ancestries for each variant
agg_dict = {col: "sum" for col in count_cols}
agg_dict["Age_min"] = "min"
agg_dict["Age_max"] = "max"
### TODO: Change this
agg_dict["Age_median"] = "mean"
for col in df_final.columns:
    if col not in count_cols + ["variant", "Age_min", "Age_max", "Age_median"]:
        agg_dict[col] = "first"
df_combined = df_final.groupby("variant", as_index=False).agg(agg_dict)
df_combined["Ancestry"] = "Combined"

# Combine individual ancestry data with combined data
df_final = [df_combined, df_final]
df_final = pd.concat(df_final)

# Merge annotations back in with the counts
df_annotations = df_orig.drop(
    columns=count_cols+["SNP", "Age_min", "Age_max", "Age_median", "Ancestry", "Otherinfo11"], 
    errors="ignore",
)
df_annotations = df_annotations.drop_duplicates()

# Merge counts for each ancestry
df_final = pd.merge(df_final, df_annotations, on="variant", how="left")
df_no_counts = df_final[(df_final["Ancestry"] == "Combined") & (df_final[["het_PD","hom_PD","het_HC","hom_HC"]].sum(axis=1) == 0)]
df_final = df_final.loc[~df_final["variant"].isin(df_no_counts["variant"]), :]
df_final.to_csv("lrrk2_combined.tsv", sep="\t", index=False)
