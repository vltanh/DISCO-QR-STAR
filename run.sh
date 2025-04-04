#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --output=slurm_output/run/slurm-%j.out
#SBATCH --job-name="project"
#SBATCH --partition=eng-instruction
#SBATCH --account=25sp-cs581a-eng
#SBATCH --mem=64G

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
                    echo "=========================="
                    echo "Running id ${run_id}"

                    id=${num_species}_gdl_${dup_rate}_${loss_rate_indicator}
                    if [ $hILS == "true" ]; then
                        id=${id}_hILS
                    fi

                    input_dir=data/trees/${id}/${run_id}
                    if [ ! -d ${input_dir} ]; then
                        echo "Input directory ${input_dir} does not exist. Skipping ${id} ${run_id}."
                        continue
                    fi
                    echo "Input directory ${input_dir} exists. Proceeding with ${id} ${run_id}."

                    output_dir=output/trees/${id}/${run_id}
                    mkdir -p $output_dir

                    # == True tree ==

                    for g_type in true 50 100 500
                    do
                        input_tree=${input_dir}/g_${g_type}.trees

                        if [ -f ${input_tree} ]; then
                            echo "Gene tree ${id} ${run_id} ${g_type} exists."
                        else
                            echo "Gene tree does not exist. Skipping ${id} ${run_id} ${g_type}."
                            continue
                        fi

                        if [ -f ${output_dir}/disco_${g_type}g_done ]; then
                            echo "DISCO already done"
                        else
                            echo "Running DISCO on gene trees"

                            python DISCO/disco.py \
                                -i ${input_tree} \
                                -o ${output_dir}/disco_${g_type}g.trees \
                                -d _ \
                            1>${output_dir}/disco_${g_type}g.log 2>${output_dir}/disco_${g_type}g.err

                            if [ -f ${output_dir}/disco_${g_type}g.trees ]; then
                                touch ${output_dir}/disco_${g_type}g_done
                            fi
                        fi

                        if [ -f ${output_dir}/disco_qr_le_${g_type}g_s_done ]; then
                            echo "QR already done"
                        else
                            echo "Running QR on gene trees"

                            python Quintet-Rooting/quintet_rooting.py \
                                -t ${input_dir}/s_tree.trees  \
                                -g ${output_dir}/disco_${g_type}g.trees \
                                -o ${output_dir}/disco_qr_le_${g_type}g_s.tree \
                                -sm LE \
                                -rs 0 \
                            1>${output_dir}/disco_qr_le_${g_type}g_s.log 2>${output_dir}/disco_qr_le_${g_type}g_s.err

                            if [ -f ${output_dir}/disco_qr_le_${g_type}g_s.tree ]; then
                                touch ${output_dir}/disco_qr_le_${g_type}g_s_done
                            fi
                        fi

                        if [ -f ${output_dir}/disco_qr_le_${g_type}g_s.tree ]; then
                            echo "Computing nCD on DISCO+QR tree"

                            python ncd.py \
                                -t1 ${input_dir}/s_tree.trees \
                                -t2 ${output_dir}/disco_qr_le_${g_type}g_s.tree \
                            >${output_dir}/disco_qr_le_${g_type}g_ncd.score
                        fi

                        if [ -f ${output_dir}/disco_qrstar_le_${g_type}g_s_done ]; then
                            echo "QR-STAR already done"
                        else
                            echo "Running QR-STAR on gene trees"

                            python Quintet-Rooting/quintet_rooting.py \
                                -t ${input_dir}/s_tree.trees  \
                                -g ${output_dir}/disco_${g_type}g.trees \
                                -o ${output_dir}/disco_qrstar_le_${g_type}g_s.tree \
                                -sm LE \
                                -c STAR \
                                -rs 0 \
                            1>${output_dir}/disco_qrstar_le_${g_type}g_s.log 2>${output_dir}/disco_qrstar_le_${g_type}g_s.err

                            if [ -f ${output_dir}/disco_qrstar_le_${g_type}g_s.tree ]; then
                                touch ${output_dir}/disco_qrstar_le_${g_type}g_s_done
                            fi
                        fi

                        if [ -f ${output_dir}/disco_qrstar_le_${g_type}g_s.tree ]; then
                            echo "Computing nCD on DISCO+QR-STAR tree"

                            python ncd.py \
                                -t1 ${input_dir}/s_tree.trees \
                                -t2 ${output_dir}/disco_qrstar_le_${g_type}g_s.tree \
                            >${output_dir}/disco_qrstar_le_${g_type}g_ncd.score
                        fi

                        if [ -f ${output_dir}/disco_qr_le_${g_type}g_ncd.score ] && [ -f ${output_dir}/disco_qrstar_le_${g_type}g_ncd.score ]; then
                            echo "nCD scores computed"
                            score_qr=$(cat ${output_dir}/disco_qr_le_${g_type}g_ncd.score)
                            score_qrstar=$(cat ${output_dir}/disco_qrstar_le_${g_type}g_ncd.score)
                            echo "nCD(QR): ${score_qr}"
                            echo "nCD(QR-STAR): ${score_qrstar}"
                            if [ $(echo "$score_qr < $score_qrstar" | bc) -eq 1 ]; then
                                echo "QR is better than QR-STAR"
                                # read -p "Press enter to continue"
                            elif [ $(echo "$score_qr > $score_qrstar" | bc) -eq 1 ]; then
                                echo "QR-STAR is better than QR"
                                # read -p "Press enter to continue"
                            else
                                echo "QR and QR-STAR are equal"
                            fi
                        fi
                    done
                done
            done
        done
    done
done