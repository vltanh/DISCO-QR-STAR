import argparse
from pathlib import Path

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description="Plot results from multiple runs.")
    parser.add_argument("--input", type=str, required=True, help="Input CSV file path.")
    parser.add_argument("--output", type=str, required=True, help="Output directory for plots.")
    return parser.parse_args()

args = parse_args()
input_file = Path(args.input)
output_dir = Path(args.output)

# num_species,dup_rate,loss_rate_indicator,hILS,g_type,run_id,method,ncd
df = pd.read_csv(input_file)

output_dir.mkdir(parents=True, exist_ok=True)

# Experiment 1: Varying the duplication rate
# num_species = 20
# g_type \in [50, 100]
# hILS \in [False, True]

df_0 = df[(df["num_species"] == 20) & (df["g_type"].isin(["50", "100"]))]

fig, axes = plt.subplots(1, 2, figsize=(10, 4), dpi=300, tight_layout=True)

# First plot is for hILS = False
df_0_hILS_false = df_0[df_0["hILS"] == False]
sns.boxplot(
    x="dup_rate", y="ncd", hue="method",
    data=df_0_hILS_false, 
    showmeans=True,
    palette="Set3",
    ax=axes[0],
)
axes[0].grid(True, axis='y', linestyle='--', alpha=0.7)
axes[0].set_title("Low ILS")
axes[0].set_ylabel("Error (nCD)")
axes[0].set_xlabel("Duplication Rate")
axes[0].set_ylim(0, 1.0)
axes[0].legend_.remove()

# Second plot is for hILS = True
df_0_hILS_true = df_0[df_0["hILS"] == True]
sns.boxplot(
    x="dup_rate", y="ncd", hue="method",
    data=df_0_hILS_true, 
    showmeans=True,
    palette="Set3",
    ax=axes[1],
)
axes[1].grid(True, axis='y', linestyle='--', alpha=0.7)
axes[1].set_title("High ILS")
axes[1].set_ylabel("Error (nCD)")
axes[1].set_xlabel("Duplication Rate")
axes[1].set_ylim(0, 1.0)
axes[1].legend_.remove()

handles, labels = axes.flatten()[0].get_legend_handles_labels()
fig.legend(
    handles,
    labels,
    loc='upper center',
    ncol=3,
    bbox_to_anchor=(0.5, 1.05),
    fancybox=True,
)

# Save the figure
fig.savefig(output_dir / "exp1.pdf", bbox_inches='tight')

# Experiment 2: Varying GTEE
# num_species = 20
# dup_rate = 1e-12
# g_type \in [50, 100]
# hILS \in [False, True]

df_1 = df[(df["num_species"] == 20) & (df["dup_rate"] == "1e-12") & (df["g_type"].isin(["50", "100"]))]
fig, axes = plt.subplots(1, 2, figsize=(10, 4), dpi=300, tight_layout=True)

# First plot is for hILS = False
df_1_hILS_false = df_1[df_1["hILS"] == False]
sns.boxplot(
    x="g_type", y="ncd", hue="method",
    data=df_1_hILS_false, 
    showmeans=True,
    palette="Set3",
    ax=axes[0],
)
axes[0].grid(True, axis='y', linestyle='--', alpha=0.7)
axes[0].set_title("Low ILS")
axes[0].set_ylabel("Error (nCD)")
axes[0].set_xlabel("GTEE")
axes[0].set_ylim(0, 1.0)
