#!/bin/bash
#SBATCH --time=05-00:00:00
#SBATCH --nodes=1
#SBATCH --output=slurm_output/stride/slurm-%j.out
#SBATCH --job-name="disco+qrstar"
#SBATCH --partition=tallis
#SBATCH --mem=32G

for loss_rate_indicator in 1 # 0 1
do
    for hILS in true # true false
    do
        for dup_rate in 1e-9 # 1e-9 1e-10 5e-10 1e-11 1e-12 1e-13
        do
            for num_species in 20 50 100 # 20 50 100
            do
                for run_id in {01..10} # {01..10}
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

                    for g_type in true 50 100 500 # true 50 100 500
                    do
                        input_tree_raw=${input_dir}/g_${g_type}.trees

                        if [ ! -f ${input_tree_raw} ]; then
                            echo "Gene tree does not exist. Skipping ${id} ${run_id} ${g_type}."
                            continue
                        fi

                        for n_genes in 50 100 500 1000 # 50 100 500 1000
                        do
                            input_tree=${output_dir}/${g_type}g/${n_genes}/g_multi.trees
                            input_gtrees_dir=${output_dir}/${g_type}g/${n_genes}/split/

                            echo ${input_gtrees_dir}

                            if [ ! -d ${input_gtrees_dir} ]; then
                                mkdir -p ${output_dir}/${g_type}g/${n_genes}/
                                head -n ${n_genes} ${input_tree_raw} > ${input_tree}
                                
                                mkdir -p ${output_dir}/${g_type}g/${n_genes}/split/
                                split ${input_tree} -l 1 --additional-suffix=.tree ${output_dir}/${g_type}g/${n_genes}/split/g_multi_
                            else
                                echo "Gene trees with ${n_genes} genes already exists."
                            fi

                            if [ ! -d ${input_gtrees_dir} ] || [ ! -f ${input_tree} ]; then
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

                            for s_est_method in trues astrid astral # trues astrid astral
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
                                    elif [ $s_est_method == "astral" ]; then
                                        /usr/bin/time -v ./ASTER/bin/astral4 \
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

                                if [ -f ${output_dir}/${g_type}g/${n_genes}/stride/${s_est_method}/done ]; then
                                    echo "STRIDE already done"
                                else
                                    echo "Running STRIDE"

                                    mkdir -p ${output_dir}/${g_type}g/${n_genes}/stride/${s_est_method}/

                                    /usr/bin/time -v python STRIDE/stride/stride.py \
                                        -d ${output_dir}/${g_type}g/${n_genes}/split/ \
                                        -s dash \
                                        -S ${output_dir}/${g_type}g/${n_genes}/disco/${s_est_method}/s_unrooted_est.tree  \
                                        -o ${output_dir}/${g_type}g/${n_genes}/stride/${s_est_method}/s_rooted_est \
                                    1>${output_dir}/${g_type}g/${n_genes}/stride/${s_est_method}/run.out 2>${output_dir}/${g_type}g/${n_genes}/stride/${s_est_method}/run.err

                                    if [ -f ${output_dir}/${g_type}g/${n_genes}/stride/${s_est_method}/s_rooted_est.tree ]; then
                                        touch ${output_dir}/${g_type}g/${n_genes}/stride/${s_est_method}/done
                                    fi
                                fi

                                if [ -f ${output_dir}/${g_type}g/${n_genes}/stride/${s_est_method}/done ]; then
                                    echo "Computing nCD on STRIDE tree"

                                    python ncd.py \
                                        -t1 ${input_dir}/s_tree.trees \
                                        -t2 ${output_dir}/${g_type}g/${n_genes}/stride/${s_est_method}/s_rooted_est.tree \
                                    >${output_dir}/${g_type}g/${n_genes}/stride/${s_est_method}/s_rooted_est.score
                                fi

                                echo "Done with ${id} ${run_id} ${g_type} ${n_genes} ${s_est_method}"
                            done
                        done
                    done
                done
            done
        done
    done
done
