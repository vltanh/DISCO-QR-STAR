import argparse
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Aggregate results from multiple runs."
    )
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Input directory containing the results",
    )
    parser.add_argument(
        "--reference",
        type=str,
        required=True,
        help="Reference directory",
    )
    return parser.parse_args()


args = parse_args()
inp_dir = Path(args.input)
ref_dir = Path(args.reference)

assert inp_dir.is_dir(), f"Input directory {inp_dir} does not exist."
assert inp_dir.exists(), f"Input directory {inp_dir} does not exist."

assert ref_dir.is_dir(), f"Reference directory {ref_dir} does not exist."
assert ref_dir.exists(), f"Reference directory {ref_dir} does not exist."

hILS_list = [False, True]
dup_rate_list = ["1e-9", "1e-10", "5e-10", "1e-11", "1e-12", "1e-13"]
loss_rate_indicator_list = [0, 1]
num_species_list = [20, 50, 100]
run_id_list = range(1, 11)
n_genes_list = [50, 100, 500, 1000]
g_type_list = ["true", 50, 100, 500]
mult_list = [1, 5, 10, 50]  # [1, 5, 10, 50]

for hILS in hILS_list:
    for loss_rate_indicator in loss_rate_indicator_list:
        for dup_rate in dup_rate_list:
            for num_species in num_species_list:
                for run_id in run_id_list:
                    id = f"{num_species}_gdl_{dup_rate}_{loss_rate_indicator}"
                    if hILS:
                        id = f"{id}_hILS"

                    subref_dir = Path(f"{ref_dir}/{id}/{run_id:02d}")
                    if not subref_dir.exists():
                        continue

                    subinp_dir = Path(f"{inp_dir}/{id}/{run_id:02d}")
                    if not subinp_dir.exists():
                        print(f"Missing input directory: {subinp_dir}")
                        continue

                    for g_type in g_type_list:
                        g_fp = subref_dir / f"g_{g_type}.trees"
                        if not g_fp.is_file():
                            continue

                        for n_genes in n_genes_list:
                            for us_est_method in ["astral", "astrid", "trues"]:
                                stride_ncd_fp = (
                                    subinp_dir
                                    / f"{g_type}g"
                                    / f"{n_genes}"
                                    / "stride"
                                    / us_est_method
                                    / f"s_rooted_est.score"
                                )

                                if not stride_ncd_fp.is_file():
                                    print(f"Missing file: {stride_ncd_fp}")
                                    continue

                                for rs_est_method in ["qr", "qrstar"]:
                                    for mult in mult_list:
                                        ncd_fp = (
                                            subinp_dir
                                            / f"{g_type}g"
                                            / f"{n_genes}"
                                            / "disco"
                                            / us_est_method
                                            / rs_est_method
                                            / "le"
                                            / f"{mult}"
                                            / f"s_rooted_est.score"
                                        )

                                        if not ncd_fp.is_file():
                                            print(f"Missing file: {ncd_fp}")
                                            continue
