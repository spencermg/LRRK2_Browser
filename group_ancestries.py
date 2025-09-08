#!/usr/bin/env python3

import pandas as pd

df = pd.read_csv("/Users/grantsm/Desktop/projects/lrrk2_browser/lrrk2_data.csv", low_memory=False)
ancestries = list(df["Ancestry"].unique())

# df_counts = []
# for ancestry in ancestries:
#     df_ancestry = df[df["Ancestry"] == ancestry]
#     df_ancestry.rename(columns={"het_PD":f"het_PD_{ancestry}", "hom_PD":f"hom_PD_{ancestry}", "het_HC":f"het_HC_{ancestry}", "hom_HC":f"hom_HC_{ancestry}", "total_PD":f"total_PD_{ancestry}", "total_HC":f"total_HC_{ancestry}"}, inplace=True)
#     df_ancestry = df_ancestry.loc[:, [f"het_PD_{ancestry}", f"hom_PD_{ancestry}", f"het_HC_{ancestry}", f"hom_HC_{ancestry}", f"total_PD_{ancestry}", f"total_HC_{ancestry}"]]
#     df_counts.append(df_ancestry)

cols_to_sum = ["het_PD", "hom_PD", "total_PD", "het_HC", "hom_HC", "total_HC"]

agg_dict = {col: "sum" for col in cols_to_sum}
for col in df.columns:
    if col not in cols_to_sum + ["variant"]:
        agg_dict[col] = "first"
df_grouped = df.groupby("variant", as_index=False).agg(agg_dict)

df_grouped.to_csv("lrrk2_grouped.csv", index=False)
