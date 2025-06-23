import argparse
from pathlib import Path

import pandas as pd


def save_sorted_results(output_file, results):
    df = pd.DataFrame(results)
    df.sort_values(
        by=[
            "num_species",
            "unrooted_s_tree",
            "n_genes",
            "dup_rate",
            "loss_rate_indicator",
            "hILS",
            "g_type",
            "run_id",
            "method",
            "sampling_method",
            "sampling_mult",
        ],
        inplace=True,
    )
    df.reset_index(drop=True, inplace=True)
    df.to_csv(output_file, index=False)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Aggregate results from multiple runs."
    )
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Input directory containing the results.",
    )
    parser.add_argument(
        "--reference",
        type=str,
        help="Reference directory",
        default=None,
    )
    parser.add_argument(
        "--output", type=str, required=True, help="Output CSV file path."
    )
    return parser.parse_args()


args = parse_args()
input_dir = Path(args.input)
output_file = Path(args.output)
if args.reference:
    reference_dir = Path(args.reference)
    assert (
        reference_dir.is_dir()
    ), f"Reference directory {reference_dir} does not exist."
else:
    reference_dir = None

assert input_dir.is_dir(), f"Input directory {input_dir} does not exist."

assert (
    output_file.parent.is_dir()
), f"Output directory {output_file.parent} does not exist."

hILS_list = [False, True]
dup_rate_list = ["1e-9", "1e-10", "5e-10", "1e-11", "1e-12", "1e-13"]
loss_rate_indicator_list = [1, 0]
num_species_list = [20, 50, 100]
run_id_list = range(1, 11)
n_genes_list = [50, 100, 500, 1000]
g_type_list = ["true", 50, 100, 500]
mult_list = [1, 5, 10, 50]

results = []

for loss_rate_indicator in loss_rate_indicator_list:
    for hILS in hILS_list:
        for dup_rate in dup_rate_list:
            for num_species in num_species_list:
                for run_id in run_id_list:
                    id = f"{num_species}_gdl_{dup_rate}_{loss_rate_indicator}"
                    if hILS:
                        id = f"{id}_hILS"

                    if reference_dir:
                        subref_dir = Path(f"{reference_dir}/{id}/{run_id:02d}")
                        if not subref_dir.exists():
                            continue

                    print(f"Running id {id} {run_id}")

                    subroot_dir = Path(f"{input_dir}/{id}/{run_id:02d}")

                    if not subroot_dir.is_dir():
                        # print(
                        #     f"Output directory {subroot_dir} does not exist. Skipping {id} {run_id}."
                        # )
                        continue

                    for g_type in g_type_list:
                        if reference_dir:
                            g_fp = subref_dir / f"g_{g_type}.trees"
                            if not g_fp.is_file():
                                continue

                        for n_genes in n_genes_list:
                            for unrooted_s_tree in ["trues", "astrid", "astral"]:
                                ncd_stride_fp = (
                                    subroot_dir
                                    / f"{g_type}g"
                                    / f"{n_genes}"
                                    / "stride"
                                    / f"{unrooted_s_tree}"
                                    / "s_rooted_est.score"
                                )

                                if ncd_stride_fp.is_file():
                                    ncd_stride = float(
                                        ncd_stride_fp.read_text().strip()
                                    )

                                    results.append(
                                        {
                                            "num_species": num_species,
                                            "unrooted_s_tree": unrooted_s_tree,
                                            "n_genes": n_genes,
                                            "dup_rate": dup_rate,
                                            "loss_rate_indicator": loss_rate_indicator,
                                            "hILS": hILS,
                                            "g_type": g_type,
                                            "run_id": run_id,
                                            "method": "STRIDE",
                                            "sampling_method": None,
                                            "sampling_mult": None,
                                            "ncd": ncd_stride,
                                        }
                                    )

                                for mult in mult_list:
                                    ncd_qr_fp = (
                                        subroot_dir
                                        / f"{g_type}g"
                                        / f"{n_genes}"
                                        / "disco"
                                        / f"{unrooted_s_tree}"
                                        / "qr"
                                        / "le"
                                        / f"{mult}"
                                        / f"s_rooted_est.score"
                                    )

                                    if ncd_qr_fp.is_file():
                                        ncd_qr = float(ncd_qr_fp.read_text().strip())

                                        results.append(
                                            {
                                                "num_species": num_species,
                                                "unrooted_s_tree": unrooted_s_tree,
                                                "n_genes": n_genes,
                                                "dup_rate": dup_rate,
                                                "loss_rate_indicator": loss_rate_indicator,
                                                "hILS": hILS,
                                                "g_type": g_type,
                                                "run_id": run_id,
                                                "method": f"DISCO+QR",
                                                "sampling_method": "le",
                                                "sampling_mult": mult,
                                                "ncd": ncd_qr,
                                            }
                                        )

                                    ncd_qrstar_fp = (
                                        subroot_dir
                                        / f"{g_type}g"
                                        / f"{n_genes}"
                                        / "disco"
                                        / f"{unrooted_s_tree}"
                                        / "qrstar"
                                        / "le"
                                        / f"{mult}"
                                        / f"s_rooted_est.score"
                                    )

                                    if ncd_qrstar_fp.is_file():
                                        ncd_qrstar = float(
                                            ncd_qrstar_fp.read_text().strip()
                                        )

                                        results.append(
                                            {
                                                "num_species": num_species,
                                                "unrooted_s_tree": unrooted_s_tree,
                                                "n_genes": n_genes,
                                                "dup_rate": dup_rate,
                                                "loss_rate_indicator": loss_rate_indicator,
                                                "hILS": hILS,
                                                "g_type": g_type,
                                                "run_id": run_id,
                                                "method": f"DISCO+QR-STAR",
                                                "sampling_method": "le",
                                                "sampling_mult": mult,
                                                "ncd": ncd_qrstar,
                                            }
                                        )

                save_sorted_results(output_file, results)
save_sorted_results(output_file, results)
