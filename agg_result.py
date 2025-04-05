import argparse
from pathlib import Path

import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description="Aggregate results from multiple runs.")
    parser.add_argument("--input", type=str, required=True, help="Input directory containing the results.")
    parser.add_argument("--output", type=str, required=True, help="Output CSV file path.")
    return parser.parse_args()

args = parse_args()
input_dir = Path(args.input)
output_file = Path(args.output)

assert input_dir.is_dir(), f"Input directory {input_dir} does not exist."
assert input_dir.exists(), f"Input directory {input_dir} does not exist."

assert output_file.parent.is_dir(), f"Output directory {output_file.parent} does not exist."
assert output_file.parent.exists(), f"Output directory {output_file.parent} does not exist."

is_complete = []
results = []
for hILS in [False, True]:
    for loss_rate_indicator in [0, 1]:
        for dup_rate in ['1e-9', '1e-10', '5e-10', '1e-11', '1e-12', '1e-13']:
            for num_species in [20, 50, 100]:
                for run_id in range(1, 11):
                    print(f"==========================")
                    print(f"Running id {run_id}")

                    id = f"{num_species}_gdl_{dup_rate}_{loss_rate_indicator}"
                    if hILS:
                        id = f"{id}_hILS"

                    output_dir = Path(f"{input_dir}/{id}/{run_id:02d}")
                    if not output_dir.is_dir():
                        print(f"Output directory {output_dir} does not exist. Skipping {id} {run_id}.")
                        for g_type in ["true", 50, 100, 500]:
                            is_complete.append({
                                "num_species": num_species,
                                "dup_rate": dup_rate,
                                "loss_rate_indicator": loss_rate_indicator,
                                "hILS": hILS,
                                "g_type": g_type,
                                "run_id": run_id,
                                "method": "DISCO+QR",
                                "is_complete": False,
                            })
                            is_complete.append({
                                "num_species": num_species,
                                "dup_rate": dup_rate,
                                "loss_rate_indicator": loss_rate_indicator,
                                "hILS": hILS,
                                "g_type": g_type,
                                "run_id": run_id,
                                "method": "DISCO+QR-STAR",
                                "is_complete": False,
                            })
                        continue

                    for g_type in ["true", 50, 100, 500]:
                        # score_qr=$(cat ${output_dir}/disco_qr_le_${g_type}g_ncd.score)
                        # score_qrstar=$(cat ${output_dir}/disco_qrstar_le_${g_type}g_ncd.score)
                        ncd_qr_fp = output_dir / f"disco_qr_le_{g_type}g_ncd.score"
                        ncd_qrstar_fp = output_dir / f"disco_qrstar_le_{g_type}g_ncd.score"
                        
                        is_continue = True
                        
                        if not ncd_qr_fp.is_file():
                            print(f"File {ncd_qr_fp} does not exist. Skipping {id} {run_id}.")
                            is_complete.append({
                                "num_species": num_species,
                                "dup_rate": dup_rate,
                                "loss_rate_indicator": loss_rate_indicator,
                                "hILS": hILS,
                                "g_type": g_type,
                                "run_id": run_id,
                                "method": "DISCO+QR",
                                "is_complete": False,
                            })
                            is_continue = False
                        
                        if not ncd_qrstar_fp.is_file():
                            print(f"File {ncd_qrstar_fp} does not exist. Skipping {id} {run_id}.")
                            is_complete.append({
                                "num_species": num_species,
                                "dup_rate": dup_rate,
                                "loss_rate_indicator": loss_rate_indicator,
                                "hILS": hILS,
                                "g_type": g_type,
                                "run_id": run_id,
                                "method": "DISCO+QR-STAR",
                                "is_complete": False,
                            })
                            is_continue = False
                            
                        if not is_continue:
                            continue
                        
                        ncd_qr = float(ncd_qr_fp.read_text().strip())
                        ncd_qrstar = float(ncd_qrstar_fp.read_text().strip())
                        
                        results.append({
                            "num_species": num_species,
                            "dup_rate": dup_rate,
                            "loss_rate_indicator": loss_rate_indicator,
                            "hILS": hILS,
                            "g_type": g_type,
                            "run_id": run_id,
                            "method": "DISCO+QR",
                            "ncd": ncd_qr,
                        })
                        results.append({
                            "num_species": num_species,
                            "dup_rate": dup_rate,
                            "loss_rate_indicator": loss_rate_indicator,
                            "hILS": hILS,
                            "g_type": g_type,
                            "run_id": run_id,
                            "method": "DISCO+QR-STAR",
                            "ncd": ncd_qrstar,
                        })
                        
df = pd.DataFrame(results)
df.sort_values(by=["num_species", "dup_rate", "loss_rate_indicator", "hILS", "g_type", "run_id"], inplace=True)
df.reset_index(drop=True, inplace=True)
df.to_csv(output_file, index=False)

df = pd.DataFrame(is_complete)
df.sort_values(by=["num_species", "dup_rate", "loss_rate_indicator", "hILS", "g_type", "run_id"], inplace=True)
df.reset_index(drop=True, inplace=True)
df.to_csv(output_file.parent / "is_complete.csv", index=False)