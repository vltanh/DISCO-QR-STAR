#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --output=slurm_output/run/slurm-%j.out
#SBATCH --job-name="project"
#SBATCH --partition=eng-instruction
#SBATCH --account=25sp-cs581a-eng
#SBATCH --mem=64G

root_dir="output/trees"

for hILS in false true
do
    for loss_rate_indicator in 0 1
    do
        for dup_rate in 1e-9 1e-10 5e-10 1e-11 1e-12 1e-13
        do
            for num_species in 20 50 100
            do
                for run_id in {01..10}
                do
                    id=${num_species}_gdl_${dup_rate}_${loss_rate_indicator}
                    if [ $hILS == "true" ]; then
                        id=${id}_hILS
                    fi

                    subroot_dir=${root_dir}/${id}/${run_id}

                    for g_type in true 50 100 500
                    do
                        for n_genes in 50 100 500 1000
                        do
                            if [ ! -f data/trees/${id}/${run_id}/s_tree.trees ]; then
                                continue
                            fi

                            mkdir -p ${subroot_dir}/${g_type}g/${n_genes}/

                            mv ${subroot_dir}/g_${g_type}_${n_genes}.trees ${subroot_dir}/${g_type}g/${n_genes}/g_multi.trees

                            mkdir -p ${subroot_dir}/${g_type}g/${n_genes}/disco/
                            mv ${subroot_dir}/disco_${g_type}g_${n_genes}.trees ${subroot_dir}/${g_type}g/${n_genes}/disco/g_single.trees
                            mv ${subroot_dir}/disco_${g_type}g_${n_genes}_done ${subroot_dir}/${g_type}g/${n_genes}/disco/done
                            mv ${subroot_dir}/disco_${g_type}g_${n_genes}.out ${subroot_dir}/${g_type}g/${n_genes}/disco/run.out
                            mv ${subroot_dir}/disco_${g_type}g_${n_genes}.err ${subroot_dir}/${g_type}g/${n_genes}/disco/run.err

                            mkdir -p ${subroot_dir}/${g_type}g/${n_genes}/disco/astrid/
                            mv ${subroot_dir}/astriddisco_${g_type}g_${n_genes}_s.tree ${subroot_dir}/${g_type}g/${n_genes}/disco/astrid/s_unrooted_est.tree
                            mv ${subroot_dir}/astriddisco_${g_type}g_${n_genes}_s.tree.1 ${subroot_dir}/${g_type}g/${n_genes}/disco/astrid/s_unrooted_est.tree.1
                            mv ${subroot_dir}/astriddisco_${g_type}g_${n_genes}_s.tree.2 ${subroot_dir}/${g_type}g/${n_genes}/disco/astrid/s_unrooted_est.tree.2
                            mv ${subroot_dir}/astriddisco_${g_type}g_${n_genes}_s.tree.3 ${subroot_dir}/${g_type}g/${n_genes}/disco/astrid/s_unrooted_est.tree.3
                            mv ${subroot_dir}/astriddisco_${g_type}g_${n_genes}_s_done ${subroot_dir}/${g_type}g/${n_genes}/disco/astrid/done
                            mv ${subroot_dir}/astriddisco_${g_type}g_${n_genes}_s.out ${subroot_dir}/${g_type}g/${n_genes}/disco/astrid/run.out
                            mv ${subroot_dir}/astriddisco_${g_type}g_${n_genes}_s.err ${subroot_dir}/${g_type}g/${n_genes}/disco/astrid/run.err

                            mkdir -p ${subroot_dir}/${g_type}g/${n_genes}/disco/astrid/qr/le/
                            mv ${subroot_dir}/disco_qr_le_${g_type}g_${n_genes}_s.tree ${subroot_dir}/${g_type}g/${n_genes}/disco/astrid/qr/le/s_rooted_est.tree
                            mv ${subroot_dir}/disco_qr_le_${g_type}g_${n_genes}_s_done ${subroot_dir}/${g_type}g/${n_genes}/disco/astrid/qr/done
                            mv ${subroot_dir}/disco_qr_le_${g_type}g_${n_genes}_s.out ${subroot_dir}/${g_type}g/${n_genes}/disco/astrid/qr/run.out
                            mv ${subroot_dir}/disco_qr_le_${g_type}g_${n_genes}_s.err ${subroot_dir}/${g_type}g/${n_genes}/disco/astrid/qr/run.err
                            mv ${subroot_dir}/disco_qr_le_${g_type}g_${n_genes}_ncd.score ${subroot_dir}/${g_type}g/${n_genes}/disco/astrid/qr/s_rooted_est.score
                            mv ${subroot_dir}/${g_type}g/${n_genes}/disco/astrid/qr/done ${subroot_dir}/${g_type}g/${n_genes}/disco/astrid/qr/le/done
                            mv ${subroot_dir}/${g_type}g/${n_genes}/disco/astrid/qr/run.out ${subroot_dir}/${g_type}g/${n_genes}/disco/astrid/qr/le/run.out
                            mv ${subroot_dir}/${g_type}g/${n_genes}/disco/astrid/qr/run.err ${subroot_dir}/${g_type}g/${n_genes}/disco/astrid/qr/le/run.err
                            mv ${subroot_dir}/${g_type}g/${n_genes}/disco/astrid/qr/s_rooted_est.score ${subroot_dir}/${g_type}g/${n_genes}/disco/astrid/qr/le/s_rooted_est.score

                            mkdir -p ${subroot_dir}/${g_type}g/${n_genes}/disco/astrid/qrstar/le/
                            mv ${subroot_dir}/disco_qrstar_le_${g_type}g_${n_genes}_s.tree ${subroot_dir}/${g_type}g/${n_genes}/disco/astrid/qrstar/le/s_rooted_est.tree
                            mv ${subroot_dir}/disco_qrstar_le_${g_type}g_${n_genes}_s_done ${subroot_dir}/${g_type}g/${n_genes}/disco/astrid/qrstar/le/done
                            mv ${subroot_dir}/disco_qrstar_le_${g_type}g_${n_genes}_s.out ${subroot_dir}/${g_type}g/${n_genes}/disco/astrid/qrstar/le/run.out
                            mv ${subroot_dir}/disco_qrstar_le_${g_type}g_${n_genes}_s.err ${subroot_dir}/${g_type}g/${n_genes}/disco/astrid/qrstar/le/run.err
                            mv ${subroot_dir}/disco_qrstar_le_${g_type}g_${n_genes}_ncd.score ${subroot_dir}/${g_type}g/${n_genes}/disco/astrid/qrstar/le/s_rooted_est.score

                            mkdir -p ${subroot_dir}/${g_type}g/${n_genes}/disco/trues/qr/le/
                            mv ${subroot_dir}/disco_qr_le_${g_type}g_${n_genes}_trues.tree ${subroot_dir}/${g_type}g/${n_genes}/disco/trues/qr/le/s_rooted_est.tree
                            mv ${subroot_dir}/disco_qr_le_${g_type}g_${n_genes}_trues_done ${subroot_dir}/${g_type}g/${n_genes}/disco/trues/qr/le/done
                            mv ${subroot_dir}/disco_qr_le_${g_type}g_${n_genes}_trues.out ${subroot_dir}/${g_type}g/${n_genes}/disco/trues/qr/le/run.out
                            mv ${subroot_dir}/disco_qr_le_${g_type}g_${n_genes}_trues.err ${subroot_dir}/${g_type}g/${n_genes}/disco/trues/qr/le/run.err
                            mv ${subroot_dir}/disco_qr_le_${g_type}g_${n_genes}_trues_ncd.score ${subroot_dir}/${g_type}g/${n_genes}/disco/trues/qr/le/s_rooted_est.score

                            mkdir -p ${subroot_dir}/${g_type}g/${n_genes}/disco/trues/qrstar/le/
                            mv ${subroot_dir}/disco_qrstar_le_${g_type}g_${n_genes}_trues.tree ${subroot_dir}/${g_type}g/${n_genes}/disco/trues/qrstar/le/s_rooted_est.tree
                            mv ${subroot_dir}/disco_qrstar_le_${g_type}g_${n_genes}_trues_done ${subroot_dir}/${g_type}g/${n_genes}/disco/trues/qrstar/le/done
                            mv ${subroot_dir}/disco_qrstar_le_${g_type}g_${n_genes}_trues.out ${subroot_dir}/${g_type}g/${n_genes}/disco/trues/qrstar/le/run.out
                            mv ${subroot_dir}/disco_qrstar_le_${g_type}g_${n_genes}_trues.err ${subroot_dir}/${g_type}g/${n_genes}/disco/trues/qrstar/le/run.err
                            mv ${subroot_dir}/disco_qrstar_le_${g_type}g_${n_genes}_trues_ncd.score ${subroot_dir}/${g_type}g/${n_genes}/disco/trues/qrstar/le/s_rooted_est.score
                        done
                    done
                done
            done
        done
    done
done