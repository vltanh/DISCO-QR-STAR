import argparse
from pathlib import Path

import pandas as pd


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
        "--output", type=str, required=True, help="Output CSV file path."
    )
    return parser.parse_args()


args = parse_args()
input_dir = Path(args.input)
output_file = Path(args.output)

assert input_dir.is_dir(), f"Input directory {input_dir} does not exist."
assert input_dir.exists(), f"Input directory {input_dir} does not exist."

assert (
    output_file.parent.is_dir()
), f"Output directory {output_file.parent} does not exist."
assert (
    output_file.parent.exists()
), f"Output directory {output_file.parent} does not exist."

hILS_list = [False, True]
dup_rate_list = ["1e-9", "1e-10", "5e-10", "1e-11", "1e-12", "1e-13"]
loss_rate_indicator_list = [0, 1]
num_species_list = [20, 50, 100]
run_id_list = range(1, 11)
n_genes_list = [50, 100, 500, 1000]
g_type_list = ["true", 50, 100, 500]
mult_list = [1, 5, 10, 50]

results = []
for hILS in hILS_list:
    for loss_rate_indicator in loss_rate_indicator_list:
        for dup_rate in dup_rate_list:
            for num_species in num_species_list:
                for run_id in run_id_list:
                    print(f"==========================")
                    print(f"Running id {run_id}")

                    id = f"{num_species}_gdl_{dup_rate}_{loss_rate_indicator}"
                    if hILS:
                        id = f"{id}_hILS"

                    subroot_dir = Path(f"{input_dir}/{id}/{run_id:02d}")

                    if not subroot_dir.is_dir():
                        print(
                            f"Output directory {subroot_dir} does not exist. Skipping {id} {run_id}."
                        )
                        continue

                    for n_genes in n_genes_list:
                        for g_type in g_type_list:
                            for mult in mult_list:
                                ncd_qr_fp = (
                                    subroot_dir
                                    / f"{g_type}g"
                                    / f"{n_genes}"
                                    / "disco"
                                    / "astrid"
                                    / "qr"
                                    / "le"
                                    / f"{mult}"
                                    / f"s_rooted_est.score"
                                )

                                if not ncd_qr_fp.is_file():
                                    print(f"File {ncd_qr_fp} does not exist.")
                                else:
                                    ncd_qr = float(ncd_qr_fp.read_text().strip())

                                    results.append(
                                        {
                                            "num_species": num_species,
                                            "true_species_tree": False,
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
                                    / "astrid"
                                    / "qrstar"
                                    / "le"
                                    / f"{mult}"
                                    / f"s_rooted_est.score"
                                )

                                if not ncd_qrstar_fp.is_file():
                                    print(f"File {ncd_qrstar_fp} does not exist.")
                                else:
                                    ncd_qrstar = float(
                                        ncd_qrstar_fp.read_text().strip()
                                    )

                                    results.append(
                                        {
                                            "num_species": num_species,
                                            "true_species_tree": False,
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

                                ncd_qr_trues_fp = (
                                    subroot_dir
                                    / f"{g_type}g"
                                    / f"{n_genes}"
                                    / "disco"
                                    / "trues"
                                    / "qr"
                                    / "le"
                                    / f"{mult}"
                                    / f"s_rooted_est.score"
                                )

                                if not ncd_qr_trues_fp.is_file():
                                    print(f"File {ncd_qr_trues_fp} does not exist.")
                                else:
                                    ncd_qr_trues = float(
                                        ncd_qr_trues_fp.read_text().strip()
                                    )

                                    results.append(
                                        {
                                            "num_species": num_species,
                                            "true_species_tree": True,
                                            "n_genes": n_genes,
                                            "dup_rate": dup_rate,
                                            "loss_rate_indicator": loss_rate_indicator,
                                            "hILS": hILS,
                                            "g_type": g_type,
                                            "run_id": run_id,
                                            "method": f"DISCO+QR",
                                            "sampling_method": "le",
                                            "sampling_mult": mult,
                                            "ncd": ncd_qr_trues,
                                        }
                                    )

                                ncd_qrstar_trues_fp = (
                                    subroot_dir
                                    / f"{g_type}g"
                                    / f"{n_genes}"
                                    / "disco"
                                    / "trues"
                                    / "qrstar"
                                    / "le"
                                    / f"{mult}"
                                    / f"s_rooted_est.score"
                                )

                                if not ncd_qrstar_trues_fp.is_file():
                                    print(f"File {ncd_qrstar_trues_fp} does not exist.")
                                else:
                                    ncd_qrstar_trues = float(
                                        ncd_qrstar_trues_fp.read_text().strip()
                                    )

                                    results.append(
                                        {
                                            "num_species": num_species,
                                            "true_species_tree": True,
                                            "n_genes": n_genes,
                                            "dup_rate": dup_rate,
                                            "loss_rate_indicator": loss_rate_indicator,
                                            "hILS": hILS,
                                            "g_type": g_type,
                                            "run_id": run_id,
                                            "method": f"DISCO+QR-STAR",
                                            "sampling_method": "le",
                                            "sampling_mult": mult,
                                            "ncd": ncd_qrstar_trues,
                                        }
                                    )

df = pd.DataFrame(results)
df.sort_values(
    by=[
        "num_species",
        "true_species_tree",
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
