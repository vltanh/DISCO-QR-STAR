import argparse
from pathlib import Path

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd


def plot_comparison(
    df,
    filters,
    hue,
    x,
    xlabel,
    y="ncd",
    ylabel="Error (nCD)",
    order=None,
    ylim=(0.0, 0.65),
):
    df_filtered = df.copy()
    for key, value in filters.items():
        if isinstance(value, list):
            # Handle None explicitly: keep rows where key is in value or (None in value and pd.isna)
            mask = df_filtered[key].isin([v for v in value if v is not None])
            if None in value:
                mask |= df_filtered[key].isna()
            df_filtered = df_filtered[mask]
        else:
            if value is None:
                df_filtered = df_filtered[df_filtered[key].isna()]
            else:
                df_filtered = df_filtered[df_filtered[key] == value]

    df_filtered["ILS Level"] = df_filtered["hILS"].map(
        {False: "Low ILS", True: "High ILS"}
    )

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), dpi=300, tight_layout=True)

    for i, ils_level in enumerate(["Low ILS", "High ILS"]):
        sns.boxplot(
            data=df_filtered[df_filtered["ILS Level"] == ils_level],
            x=x,
            y=y,
            hue=hue,
            palette="Set3",
            order=order,
            ax=axes[i],
            showmeans=True,
            fliersize=2,
            medianprops={"color": "red", "linewidth": 2, "alpha": 0.7},
        )
        axes[i].grid(True, axis="y", linestyle="--", alpha=0.7)
        axes[i].set_title(ils_level)
        axes[i].set_xlabel(xlabel)
        axes[i].set_ylabel(ylabel if i == 0 else "")
        axes[i].set_ylim(*ylim)
        if i == 1:
            axes[i].tick_params(axis="y", left=False, labelleft=False)
        axes[i].legend_.remove()

    handles, labels = axes.flatten()[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.08))

    return fig


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot focused results from multiple runs."
    )
    parser.add_argument("--input", type=str, required=True, help="Input CSV file path.")
    parser.add_argument(
        "--output", type=str, required=True, help="Output directory for plots."
    )
    return parser.parse_args()


args = parse_args()
input_file = Path(args.input)
output_dir = Path(args.output)
output_dir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(input_file)

df["unrooted_s_tree"] = df["unrooted_s_tree"].map(
    {"trues": "True", "astrid": "ASTRID", "astral": "ASTRAL"}
)
df["num_species"] = df["num_species"].map({20: 21, 50: 51, 100: 101})

# Experiment 1: Multiplicity of quintet sampling
fig = plot_comparison(
    df,
    filters={
        "num_species": 21,
        "n_genes": 1000,
        "g_type": "true",
        "dup_rate": 1e-12,
        "loss_rate_indicator": 1,
        "unrooted_s_tree": "True",
        "sampling_mult": [1, 5, 10, 50],
    },
    hue="method",
    x="sampling_mult",
    xlabel="Quintet Sampling Multiplicity",
)
fig.savefig(output_dir / "exp_multiplicity_true.pdf", bbox_inches="tight")
plt.close(fig)

fig = plot_comparison(
    df,
    filters={
        "num_species": 21,
        "n_genes": 1000,
        "g_type": "true",
        "dup_rate": 1e-12,
        "loss_rate_indicator": 1,
        "unrooted_s_tree": "ASTRAL",
        "sampling_mult": [1, 5, 10, 50],
    },
    hue="method",
    x="sampling_mult",
    xlabel="Quintet Sampling Multiplicity",
)
fig.savefig(output_dir / "exp_multiplicity_astral.pdf", bbox_inches="tight")
plt.close(fig)

# Experiment 2: Method for unrooted species tree estimation
fig = plot_comparison(
    df,
    filters={
        "num_species": 21,
        "n_genes": 1000,
        "g_type": "true",
        "dup_rate": 1e-12,
        "loss_rate_indicator": 1,
        "sampling_mult": 50,
        "unrooted_s_tree": ["True", "ASTRID", "ASTRAL"],
    },
    hue="method",
    x="unrooted_s_tree",
    xlabel="Method for Unrooted Species Tree Estimation",
    order=["True", "ASTRAL", "ASTRID"],
)
fig.savefig(output_dir / "exp_unrooted-species-tree.pdf", bbox_inches="tight")
plt.close(fig)

# With STRIDE
fig = plot_comparison(
    df,
    filters={
        "num_species": 21,
        "n_genes": 1000,
        "g_type": "true",
        "dup_rate": 1e-12,
        "loss_rate_indicator": 1,
        "sampling_mult": [None, 50],
        "sampling_method": [None, "le"],
        "unrooted_s_tree": ["True", "ASTRID", "ASTRAL"],
    },
    hue="method",
    x="unrooted_s_tree",
    xlabel="Method for Unrooted Species Tree Estimation",
    order=["True", "ASTRAL", "ASTRID"],
)
fig.savefig(output_dir / "exp_unrooted-species-tree_wSTRIDE.pdf", bbox_inches="tight")

for unrooted_s_tree in df["unrooted_s_tree"].unique():
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / unrooted_s_tree).mkdir(parents=True, exist_ok=True)

    # Experiment 3: Varying the duplication rate
    fig = plot_comparison(
        df,
        filters={
            "num_species": 21,
            "n_genes": 1000,
            "g_type": "true",
            "dup_rate": [1e-9, 5e-10, 1e-10, 1e-11, 1e-12, 1e-13],
            "loss_rate_indicator": 1,
            "sampling_mult": 50,
            "unrooted_s_tree": unrooted_s_tree,
        },
        hue="method",
        x="dup_rate",
        xlabel="Duplication Rate",
    )
    fig.savefig(
        output_dir / unrooted_s_tree / "exp_duplication-rate.pdf", bbox_inches="tight"
    )
    plt.close(fig)

    # With STRIDE
    fig = plot_comparison(
        df,
        filters={
            "num_species": 21,
            "n_genes": 1000,
            "g_type": "true",
            "dup_rate": [1e-9, 5e-10, 1e-10, 1e-11, 1e-12, 1e-13],
            "loss_rate_indicator": 1,
            "sampling_mult": [None, 50],
            "sampling_method": [None, "le"],
            "unrooted_s_tree": unrooted_s_tree,
        },
        hue="method",
        x="dup_rate",
        xlabel="Duplication Rate",
    )
    fig.savefig(
        output_dir / unrooted_s_tree / "exp_duplication-rate_wSTRIDE.pdf",
        bbox_inches="tight",
    )
    plt.close(fig)

    # Experiment 4: Varying the GTEE
    fig = plot_comparison(
        df,
        filters={
            "num_species": 21,
            "n_genes": 1000,
            "g_type": ["true", "50", "100", "500"],
            "dup_rate": 1e-12,
            "loss_rate_indicator": 1,
            "sampling_mult": 50,
            "unrooted_s_tree": unrooted_s_tree,
        },
        hue="method",
        x="g_type",
        xlabel="Sequence Length (bp)",
    )
    fig.savefig(
        output_dir / unrooted_s_tree / "exp_sequence-length.pdf", bbox_inches="tight"
    )
    plt.close(fig)

    # With STRIDE
    fig = plot_comparison(
        df,
        filters={
            "num_species": 21,
            "n_genes": 1000,
            "g_type": ["true", "50", "100", "500"],
            "dup_rate": 1e-12,
            "loss_rate_indicator": 1,
            "sampling_mult": [None, 50],
            "sampling_method": [None, "le"],
            "unrooted_s_tree": unrooted_s_tree,
        },
        hue="method",
        x="g_type",
        xlabel="Sequence Length (bp)",
    )
    fig.savefig(
        output_dir / unrooted_s_tree / "exp_sequence-length_wSTRIDE.pdf",
        bbox_inches="tight",
    )
    plt.close(fig)

    # Experiment 5: Varying the number of gene trees
    fig = plot_comparison(
        df,
        filters={
            "num_species": 21,
            "n_genes": [50, 100, 500, 1000],
            "g_type": "true",
            "dup_rate": 1e-12,
            "loss_rate_indicator": 1,
            "sampling_mult": 50,
            "unrooted_s_tree": unrooted_s_tree,
        },
        hue="method",
        x="n_genes",
        xlabel="Number of Gene Trees",
    )
    fig.savefig(
        output_dir / unrooted_s_tree / "exp_number-of-gene-trees.pdf",
        bbox_inches="tight",
    )
    plt.close(fig)

    # With STRIDE
    fig = plot_comparison(
        df,
        filters={
            "num_species": 21,
            "n_genes": [50, 100, 500, 1000],
            "g_type": "true",
            "dup_rate": 1e-12,
            "loss_rate_indicator": 1,
            "sampling_mult": [None, 50],
            "sampling_method": [None, "le"],
            "unrooted_s_tree": unrooted_s_tree,
        },
        hue="method",
        x="n_genes",
        xlabel="Number of Gene Trees",
    )
    fig.savefig(
        output_dir / unrooted_s_tree / "exp_number-of-gene-trees_wSTRIDE.pdf",
        bbox_inches="tight",
    )
    plt.close(fig)

    # Experiment 6: Varying the number of species
    fig = plot_comparison(
        df,
        filters={
            "num_species": [21, 51, 101],
            "n_genes": 1000,
            "g_type": "true",
            "dup_rate": 1e-12,
            "loss_rate_indicator": 1,
            "sampling_mult": 50,
            "unrooted_s_tree": unrooted_s_tree,
        },
        hue="method",
        x="num_species",
        xlabel="Number of Species",
    )
    fig.savefig(
        output_dir / unrooted_s_tree / "exp_number-of-species.pdf", bbox_inches="tight"
    )
    plt.close(fig)

    # With STRIDE
    fig = plot_comparison(
        df,
        filters={
            "num_species": [21, 51, 101],
            "n_genes": 1000,
            "g_type": "true",
            "dup_rate": 1e-12,
            "loss_rate_indicator": 1,
            "sampling_mult": [None, 50],
            "sampling_method": [None, "le"],
            "unrooted_s_tree": unrooted_s_tree,
        },
        hue="method",
        x="num_species",
        xlabel="Number of Species",
    )
    fig.savefig(
        output_dir / unrooted_s_tree / "exp_number-of-species_wSTRIDE.pdf",
        bbox_inches="tight",
    )
    plt.close(fig)
