import argparse
from pathlib import Path

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd


def plot_experiment_2x2(df, x_axis, x_axis_name, ylim=(0.0, 0.65)):
    methods = ["DISCO+QR", "DISCO+QR-STAR"]
    hILS_values = [False, True]
    titles = [
        ["Low ILS / DISCO+QR", "Low ILS / DISCO+QR-STAR"],
        ["High ILS / DISCO+QR", "High ILS / DISCO+QR-STAR"],
    ]
    y_axis = "ncd"
    hue = "sampling_mult"

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), dpi=300, tight_layout=True)

    for i, hILS in enumerate(hILS_values):
        for j, method in enumerate(methods):
            df_filtered = df[(df["hILS"] == hILS) & (df["method"] == method)]
            sns.boxplot(
                x=x_axis,
                y=y_axis,
                hue=hue,
                data=df_filtered,
                showmeans=True,
                fliersize=2,
                palette="Set3",
                ax=axes[i, j],
            )
            axes[i, j].grid(True, axis="y", linestyle="--", alpha=0.7)
            axes[i, j].set_title(titles[i][j])
            axes[i, j].set_ylabel("Error (nCD)")
            axes[i, j].set_xlabel(x_axis_name)
            axes[i, j].set_ylim(*ylim)
            try:
                axes[i, j].legend_.remove()
            except AttributeError:
                pass

    handles, labels = axes.flatten()[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        ncol=4,
        bbox_to_anchor=(0.5, 1.05),
        fancybox=True,
    )
    return fig


def parse_args():
    parser = argparse.ArgumentParser(description="Plot results from multiple runs.")
    parser.add_argument("--input", type=str, required=True, help="Input CSV file path.")
    parser.add_argument(
        "--output", type=str, required=True, help="Output directory for plots."
    )
    return parser.parse_args()


args = parse_args()
input_file = Path(args.input)
output_dir_raw = Path(args.output)

# num_species,unrooted_s_tree,n_genes,dup_rate,loss_rate_indicator,hILS,g_type,run_id,method,sampling_method,sampling_mult,ncd
df = pd.read_csv(input_file)

output_dir_raw.mkdir(parents=True, exist_ok=True)

for species_tree in df["unrooted_s_tree"].unique():
    output_dir = output_dir_raw / species_tree
    output_dir.mkdir(parents=True, exist_ok=True)

    # Experiment 1: Varying the duplication rate
    fig = plot_experiment_2x2(
        df[
            (df["num_species"] == 20)
            & (df["unrooted_s_tree"] == species_tree)
            & (df["g_type"].isin(["true"]))
            & (df["n_genes"] == 1000)
            & (df["dup_rate"].isin([1e-9, 5e-10, 1e-10, 1e-11, 1e-12, 1e-13]))
            & (df["loss_rate_indicator"] == 1)
            & (df["sampling_mult"].isin([1, 5, 10, 50]))
        ],
        x_axis="dup_rate",
        x_axis_name="Duplication Rate",
    )
    fig.savefig(output_dir / "exp1.pdf", bbox_inches="tight")

    # Experiment 2: Varying GTEE
    fig = plot_experiment_2x2(
        df[
            (df["num_species"] == 20)
            & (df["unrooted_s_tree"] == species_tree)
            & (df["g_type"].isin(["true", "50", "100", "500"]))
            & (df["n_genes"] == 100)
            & (df["dup_rate"] == 1e-12)
            & (df["loss_rate_indicator"] == 1)
            & (df["sampling_mult"].isin([1, 5, 10, 50]))
        ],
        x_axis="g_type",
        x_axis_name="Sequence Length (bp)",
    )
    fig.savefig(output_dir / "exp2.pdf", bbox_inches="tight")

    # Experiment 3: Varying the number of gene trees
    fig = plot_experiment_2x2(
        df[
            (df["num_species"] == 20)
            & (df["unrooted_s_tree"] == species_tree)
            & (df["g_type"].isin(["true"]))
            & (df["n_genes"].isin([50, 100, 500, 1000]))
            & (df["dup_rate"] == 1e-12)
            & (df["loss_rate_indicator"] == 1)
            & (df["sampling_mult"].isin([1, 5, 10, 50]))
        ],
        x_axis="n_genes",
        x_axis_name="Number of Gene Trees",
    )
    fig.savefig(output_dir / "exp3.pdf", bbox_inches="tight")

    # Experiment 4: Varying the number of species
    fig = plot_experiment_2x2(
        df[
            (df["num_species"].isin([20, 50, 100]))
            & (df["unrooted_s_tree"] == species_tree)
            & (df["g_type"].isin(["true"]))
            & (df["n_genes"] == 100)
            & (df["dup_rate"] == 1e-12)
            & (df["loss_rate_indicator"] == 1)
            & (df["sampling_mult"].isin([1, 5, 10, 50]))
        ],
        x_axis="num_species",
        x_axis_name="Number of Species",
    )
    fig.savefig(output_dir / "exp4.pdf", bbox_inches="tight")

# Compare between ASTRID and ASTRAL
fig = plot_experiment_2x2(
    df[
        (df["num_species"] == 20)
        & (df["unrooted_s_tree"].isin(["trues", "astrid", "astral"]))
        & (df["g_type"].isin(["true"]))
        & (df["n_genes"] == 1000)
        & (df["dup_rate"] == 1e-12)
        & (df["loss_rate_indicator"] == 1)
        & (df["sampling_mult"].isin([1, 5, 10, 50]))
    ],
    x_axis="unrooted_s_tree",
    x_axis_name="Unrooted Species Tree Estimation Method",
    ylim=(0.0, 0.65),
)
fig.savefig(output_dir_raw / "exp_u_s.pdf", bbox_inches="tight")
