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

                    for g_type in true 50 100 500
                    do
                        input_tree_raw=${input_dir}/g_${g_type}.trees

                        if [ ! -f ${input_tree_raw} ]; then
                            echo "Gene tree does not exist. Skipping ${id} ${run_id} ${g_type}."
                            continue
                        fi

                        for n_genes in 50 100 500 1000
                        do
                            input_tree=${output_dir}/${g_type}g/${n_genes}/g_multi.trees

                            if [ ! -f ${input_tree} ]; then
                                mkdir -p ${output_dir}/${g_type}g/${n_genes}/
                                head -n ${n_genes} ${input_tree_raw} > ${input_tree}
                            else
                                echo "Gene tree file with ${n_genes} genes already exists."
                            fi

                            if [ ! -f ${input_tree} ]; then
                                echo "Gene tree does not exist. Skipping ${id} ${run_id} ${g_type} ${n_genes}."
                                continue
                            fi

                            echo "Running on ${id} ${run_id} ${g_type} ${n_genes}"

                            if [ -f ${output_dir}/${g_type}g/${n_genes}/disco/done ]; then
                                echo "DISCO already done"
                            else
                                echo "Running DISCO on gene trees"

                                mkdir -p ${output_dir}/${g_type}g/${n_genes}/disco/

                                /usr/bin/time -v python DISCO/disco.py \
                                    -i ${input_tree} \
                                    -o ${output_dir}/${g_type}g/${n_genes}/disco/g_single.trees \
                                    -d _ \
                                1>${output_dir}/${g_type}g/${n_genes}/disco/run.out 2>${output_dir}/${g_type}g/${n_genes}/disco/run.err

                                if [ -f ${output_dir}/${g_type}g/${n_genes}/disco/g_single.trees ]; then
                                    touch ${output_dir}/${g_type}g/${n_genes}/disco/done
                                fi
                            fi

                            if [ ! -f ${output_dir}/${g_type}g/${n_genes}/disco/g_single.trees ]; then
                                echo "DISCO tree does not exist. Skipping ${id} ${run_id} ${g_type} ${n_genes}."
                                continue
                            fi

                            for s_est_method in trues astrid
                            do
                                if [ -f ${output_dir}/${g_type}g/${n_genes}/disco/${s_est_method}/done ]; then
                                    echo "${s_est_method} already done"
                                else
                                    mkdir -p ${output_dir}/${g_type}g/${n_genes}/disco/${s_est_method}/

                                    if [ $s_est_method == "astrid" ]; then
                                        /usr/bin/time -v ./ASTRID/bazel-bin/src/ASTRID \
                                            -i ${output_dir}/${g_type}g/${n_genes}/disco/g_single.trees \
                                            -o ${output_dir}/${g_type}g/${n_genes}/disco/${s_est_method}/s_unrooted_est.tree \
                                        1>${output_dir}/${g_type}g/${n_genes}/disco/${s_est_method}/run.out 2>${output_dir}/${g_type}g/${n_genes}/disco/${s_est_method}/run.err
                                    elif [ $s_est_method == "trues" ]; then
                                        cp ${input_dir}/s_tree.trees ${output_dir}/${g_type}g/${n_genes}/disco/${s_est_method}/s_unrooted_est.tree
                                    fi

                                    if [ -f ${output_dir}/${g_type}g/${n_genes}/disco/${s_est_method}/s_unrooted_est.tree ]; then
                                        touch ${output_dir}/${g_type}g/${n_genes}/disco/${s_est_method}/done
                                    fi
                                fi

                                if [ ! -f ${output_dir}/${g_type}g/${n_genes}/disco/${s_est_method}/s_unrooted_est.tree ]; then
                                    echo "${s_est_method} tree does not exist. Skipping ${id} ${run_id} ${g_type} ${n_genes}."
                                    continue
                                fi

                                for mult in 1 5 10 50
                                do
                                    if [ -f ${output_dir}/${g_type}g/${n_genes}/disco/${s_est_method}/qr/le/${mult}/done ]; then
                                        echo "QR already done"
                                    else
                                        echo "Running QR on gene trees"

                                        mkdir -p ${output_dir}/${g_type}g/${n_genes}/disco/${s_est_method}/qr/le/${mult}

                                        /usr/bin/time -v python Quintet-Rooting/quintet_rooting.py \
                                            -t ${output_dir}/${g_type}g/${n_genes}/disco/${s_est_method}/s_unrooted_est.tree  \
                                            -g ${output_dir}/${g_type}g/${n_genes}/disco/g_single.trees \
                                            -o ${output_dir}/${g_type}g/${n_genes}/disco/${s_est_method}/qr/le/${mult}/s_rooted_est.tree \
                                            -sm LE \
                                            -rs 0 \
                                            -mult ${mult} \
                                        1>${output_dir}/${g_type}g/${n_genes}/disco/${s_est_method}/qr/le/${mult}/run.out 2>${output_dir}/${g_type}g/${n_genes}/disco/${s_est_method}/qr/le/${mult}/run.err

                                        if [ -f ${output_dir}/${g_type}g/${n_genes}/disco/${s_est_method}/qr/le/${mult}/s_rooted_est.tree ]; then
                                            touch ${output_dir}/${g_type}g/${n_genes}/disco/${s_est_method}/qr/le/${mult}/done
                                        fi
                                    fi

                                    if [ -f ${output_dir}/${g_type}g/${n_genes}/disco/${s_est_method}/qr/le/${mult}/done ]; then
                                        echo "Computing nCD on DISCO+QR tree"

                                        python ncd.py \
                                            -t1 ${input_dir}/s_tree.trees \
                                            -t2 ${output_dir}/${g_type}g/${n_genes}/disco/${s_est_method}/qr/le/${mult}/s_rooted_est.tree \
                                        >${output_dir}/${g_type}g/${n_genes}/disco/${s_est_method}/qr/le/${mult}/s_rooted_est.score
                                    fi

                                    if [ -f ${output_dir}/${g_type}g/${n_genes}/disco/${s_est_method}/qrstar/le/${mult}/done ]; then
                                        echo "QR-STAR already done"
                                    else
                                        echo "Running QR-STAR on gene trees"

                                        mkdir -p ${output_dir}/${g_type}g/${n_genes}/disco/${s_est_method}/qrstar/le/${mult}

                                        /usr/bin/time -v python Quintet-Rooting/quintet_rooting.py \
                                            -t ${output_dir}/${g_type}g/${n_genes}/disco/${s_est_method}/s_unrooted_est.tree  \
                                            -g ${output_dir}/${g_type}g/${n_genes}/disco/g_single.trees \
                                            -o ${output_dir}/${g_type}g/${n_genes}/disco/${s_est_method}/qrstar/le/${mult}/s_rooted_est.tree \
                                            -sm LE \
                                            -c STAR \
                                            -rs 0 \
                                            -mult ${mult} \
                                        1>${output_dir}/${g_type}g/${n_genes}/disco/${s_est_method}/qrstar/le/${mult}/run.out 2>${output_dir}/${g_type}g/${n_genes}/disco/${s_est_method}/qrstar/le/${mult}/run.err

                                        if [ -f ${output_dir}/${g_type}g/${n_genes}/disco/${s_est_method}/qrstar/le/${mult}/s_rooted_est.tree ]; then
                                            touch ${output_dir}/${g_type}g/${n_genes}/disco/${s_est_method}/qrstar/le/${mult}/done
                                        fi
                                    fi

                                    if [ -f ${output_dir}/${g_type}g/${n_genes}/disco/${s_est_method}/qrstar/le/${mult}/done ]; then
                                        echo "Computing nCD on DISCO+QR-STAR tree"

                                        python ncd.py \
                                            -t1 ${input_dir}/s_tree.trees \
                                            -t2 ${output_dir}/${g_type}g/${n_genes}/disco/${s_est_method}/qrstar/le/${mult}/s_rooted_est.tree \
                                        >${output_dir}/${g_type}g/${n_genes}/disco/${s_est_method}/qrstar/le/${mult}/s_rooted_est.score
                                    fi

                                    if [ -f ${output_dir}/${g_type}g/${n_genes}/disco/${s_est_method}/qr/le/${mult}/s_rooted_est.score ] && [ -f ${output_dir}/${g_type}g/${n_genes}/disco/${s_est_method}/qrstar/le/${mult}/s_rooted_est.score ]; then
                                        echo "nCD scores computed"
                                        score_qr=$(cat ${output_dir}/${g_type}g/${n_genes}/disco/${s_est_method}/qr/le/${mult}/s_rooted_est.score)
                                        score_qrstar=$(cat ${output_dir}/${g_type}g/${n_genes}/disco/${s_est_method}/qrstar/le/${mult}/s_rooted_est.score)
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
        done
    done
done