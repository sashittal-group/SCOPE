configfile: "config.yaml"

seeds = [ i for i in range(config["nseeds"])]

import itertools

# def filter_combinator(combinator):
#    def filtered_combinator(*args, **kwargs):
#        for wc_comb in combinator(*args, **kwargs):
#            #if int(wc_comb[1][-1]) >= int(wc_comb[2][-1]):
#            if int(wc_comb[1][-1]) == int(wc_comb[2][-1]):
#                yield wc_comb
#    return filtered_combinator
#
# filtered_product = filter_combinator(itertools.product)

K = range(7, 26)
T = [0.60, 0.65, 0.70, 0.72, 0.75, 0.80]
F = [0.1, 0.2, 0.25, 0.3, 0.4]

rule all:
    input:
        # # simulation
        # expand('data/simulation/ground_truth/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_multi_state_tree_node_character_matrix.parquet', \
        #     itertools.product, seed=seeds, ncells=config['ncells'], n_mutation_groups=config['n_mutation_groups'], mutation_group_sizes=config['mutation_group_sizes'], \
        #     nclusters=config['nclusters'], cov=config['coverage']),
        # scope input
        # expand("data/simulation/scope_input/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/F_plus.csv", \
        #     itertools.product, seed=seeds, ncells=config['ncells'], n_mutation_groups=config['n_mutation_groups'], mutation_group_sizes=config['mutation_group_sizes'], \
        #     nclusters=config['nclusters'], cov=config['coverage']),
        # # scope input kmeans
        # expand("data/simulation/scope_input_kmeans/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/kmeans_clones.csv", \
        #     itertools.product, seed=seeds, ncells=config['ncells'], n_mutation_groups=config['n_mutation_groups'], mutation_group_sizes=config['mutation_group_sizes'], \
        #     nclusters=config['nclusters'], cov=config['coverage']),
        # # # scope input kmeans known k
        # expand("data/simulation/scope_input_kmeans_known_k/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/kmeans_clones.csv", \
        #     itertools.product, seed=seeds, ncells=config['ncells'], n_mutation_groups=config['n_mutation_groups'], mutation_group_sizes=config['mutation_group_sizes'], \
        #     nclusters=config['nclusters'], cov=config['coverage']),
        # # phertilizer input
        # # expand("data/simulation/phertilizer_input/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/snv_counts.tsv", \
        # #     itertools.product, seed=seeds, ncells=config['ncells'], n_mutation_groups=config['n_mutation_groups'], mutation_group_sizes=config['mutation_group_sizes'], \
        # #     nclusters=config['nclusters'], cov=config['coverage']),
        # # scope on simulation
        # expand("data/simulation/scope_output/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/summary.txt", \
        #     itertools.product, seed=seeds, ncells=config['ncells'], n_mutation_groups=config['n_mutation_groups'], mutation_group_sizes=config['mutation_group_sizes'], \
        #     nclusters=config['nclusters'], cov=config['coverage']),
        # # scope kmeans on simulation
        # expand("data/simulation/scope_output_kmeans/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/summary.txt", \
        #     itertools.product, seed=seeds, ncells=config['ncells'], n_mutation_groups=config['n_mutation_groups'], mutation_group_sizes=config['mutation_group_sizes'], \
        #     nclusters=config['nclusters'], cov=config['coverage']),
        # # scope kmeans known k on simulation
        # expand("data/simulation/scope_output_kmeans_known_k/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/summary.txt", \
        #     itertools.product, seed=seeds, ncells=config['ncells'], n_mutation_groups=config['n_mutation_groups'], mutation_group_sizes=config['mutation_group_sizes'], \
        #     nclusters=config['nclusters'], cov=config['coverage']),
        # # phertilizer on simulation
        # expand("data/simulation/phertilizer_output/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/cell_clusters.csv", \
        #     itertools.product, seed=seeds, ncells=config['ncells'], n_mutation_groups=config['n_mutation_groups'], mutation_group_sizes=config['mutation_group_sizes'], \
        #     nclusters=config['nclusters'], cov=config['coverage']),
        # # post process clone
        # expand("data/simulation/scope_post_kmeans_known_k/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/kmeans_cleaned_clones.csv", \
        #     itertools.product, seed=seeds, ncells=config['ncells'], n_mutation_groups=config['n_mutation_groups'], mutation_group_sizes=config['mutation_group_sizes'], \
        #     nclusters=config['nclusters'], cov=config['coverage']),
        # post process clone 2
        # expand("data/simulation/scope_post_kmeans/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/kmeans_cleaned_clones.csv", \
        #     itertools.product, seed=seeds, ncells=config['ncells'], n_mutation_groups=config['n_mutation_groups'], mutation_group_sizes=config['mutation_group_sizes'], \
        #     nclusters=config['nclusters'], cov=config['coverage']),
        # post williams
        # expand("data/williams/scratch/{sample_ids}/scope_mut/summary.txt", \
        #     itertools.product, sample_ids=config['sample_ids']),
        # scope on laks
        expand("data/laks/scope/outputs/k_{k}_t_{t}/summary.txt", k=K, t=T),
        # scope on laks filtered
        expand("data/laks/scope/outputs/filtered/k_{k}_t_{t}_f_{f}/summary.txt", k=K, t=T, f=F)



rule simulate:
    output:
        readcount_matrix="data/simulation/ground_truth/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_read_count.parquet",
        variantcount_matrix="data/simulation/ground_truth/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_variant_count.parquet",
        character_matrix="data/simulation/ground_truth/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_character_matrix_without_noise.parquet",
        copy_number_table="data/simulation/ground_truth/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_copy_numbers.parquet",
        mutation_group_table="data/simulation/ground_truth/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_mutation_group.parquet",
        mutation_to_bin_table="data/simulation/ground_truth/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_mutation_to_bin_mapping.parquet",
        B_truth="data/simulation/ground_truth/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_multi_state_tree_node_character_matrix.parquet",
    params:
        prefix="data/simulation/ground_truth/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim",
        fn=config['fn'],
        fp=config['fp'],
        mutation_rate=config['lambda'],
        vaf_threshold=config['vaf_threshold'],
        read_threshold=config['read_threshold'],
        missing=config['missing_rate'],
        k=config['k']
    log:
        std="data/simulation/ground_truth/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim.log", 
        err="data/simulation/ground_truth/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim.err.log",
    benchmark: "data/simulation/ground_truth/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim.benchmark", 
    shell:
        "python src/simulate_data.py -s {wildcards.seed} -o {params.prefix} "
        " -n {wildcards.ncells} -m {wildcards.n_mutation_groups} -d {params.missing} --size {wildcards.mutation_group_sizes}" 
        " -p {wildcards.nclusters} -a {params.fp} -b {params.fn} -k {params.k} -l {params.mutation_rate} --cov {wildcards.cov} --maxcn 6"
        " > {log.std} 2> {log.err}"

rule phertilizer_input:
    input:
        readcount_matrix="data/simulation/ground_truth/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_read_count.parquet",
        variantcount_matrix="data/simulation/ground_truth/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_variant_count.parquet",
        mutation_to_bin_table="data/simulation/ground_truth/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_mutation_to_bin_mapping.parquet",
    output:
        snv_counts="data/simulation/phertilizer_input/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/snv_counts.tsv",
        binned_counts="data/simulation/phertilizer_input/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/binned_read_counts.csv",
    resources:
        mem_mb=10240
    log:
        std="data/simulation/phertilizer_input/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/log",
        err="data/simulation/phertilizer_input/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/err.log",
    benchmark: "data/simulation/phertilizer_input/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/benchmark",
    shell:
        "python src/data_processors/process_simulated_data_for_phertilizer.py -i data/simulation/ground_truth/n{wildcards.ncells}_m{wildcards.n_mutation_groups}_size{wildcards.mutation_group_sizes}_cov{wildcards.cov}_p{wildcards.nclusters}_s{wildcards.seed}/sim "
        " -o data/simulation/phertilizer_input/ "
        " > {log.std} 2> {log.err}"

rule scope_input:
    input:
        character_matrix="data/simulation/ground_truth/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_character_matrix_without_noise.parquet",
        readcount_matrix="data/simulation/ground_truth/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_read_count.parquet",
        variantcount_matrix="data/simulation/ground_truth/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_variant_count.parquet",
        copy_number_table="data/simulation/ground_truth/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_copy_numbers.parquet",
        mutation_group_table="data/simulation/ground_truth/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_mutation_group.parquet",
    output:
        F_plus="data/simulation/scope_input/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/F_plus.csv",
        F_minus="data/simulation/scope_input/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/F_minus.csv",
        F_bar="data/simulation/scope_input/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/F_bar.csv",
        clone_sizes="data/simulation/scope_input/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/clone_sizes.csv",
    log:
        std="data/simulation/scope_input/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/log",
        err="data/simulation/scope_input/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/err.log",
    benchmark: "data/simulation/scope_input/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/benchmark",
    shell:
        "python src/data_processors/process_simulated_data_for_scope.py -i data/simulation/ground_truth/n{wildcards.ncells}_m{wildcards.n_mutation_groups}_size{wildcards.mutation_group_sizes}_cov{wildcards.cov}_p{wildcards.nclusters}_s{wildcards.seed}/sim "
        " -o data/simulation/scope_input/ --clone-table _mutation_group.parquet"
        " > {log.std} 2> {log.err}"

rule scope_input_kmeans:
    input:
        character_matrix="data/simulation/ground_truth/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_character_matrix_without_noise.parquet",
        readcount_matrix="data/simulation/ground_truth/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_read_count.parquet",
        variantcount_matrix="data/simulation/ground_truth/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_variant_count.parquet",
        copy_number_table="data/simulation/ground_truth/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_copy_numbers.parquet",
    output:
        F_plus="data/simulation/scope_input_kmeans/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/F_plus.csv",
        F_minus="data/simulation/scope_input_kmeans/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/F_minus.csv",
        F_bar="data/simulation/scope_input_kmeans/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/F_bar.csv",
        clone_sizes="data/simulation/scope_input_kmeans/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/clone_sizes.csv",
        kmeans_clones="data/simulation/scope_input_kmeans/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/kmeans_clones.csv",
    log:
        std="data/simulation/scope_input_kmeans/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/log",
        err="data/simulation/scope_input_kmeans/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/err.log",
    benchmark: "data/simulation/scope_input_kmeans/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/benchmark",
    shell:
        "python src/data_processors/kmeans_clustering_for_scope.py -i data/simulation/ground_truth/n{wildcards.ncells}_m{wildcards.n_mutation_groups}_size{wildcards.mutation_group_sizes}_cov{wildcards.cov}_p{wildcards.nclusters}_s{wildcards.seed}/sim "
        " -o data/simulation/scope_input_kmeans/ "
        " > {log.std} 2> {log.err}"

rule scope_input_kmeans_known_k:
    input:
        character_matrix="data/simulation/ground_truth/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_character_matrix_without_noise.parquet",
        readcount_matrix="data/simulation/ground_truth/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_read_count.parquet",
        variantcount_matrix="data/simulation/ground_truth/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_variant_count.parquet",
        copy_number_table="data/simulation/ground_truth/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_copy_numbers.parquet",
    output:
        F_plus="data/simulation/scope_input_kmeans_known_k/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/F_plus.csv",
        F_minus="data/simulation/scope_input_kmeans_known_k/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/F_minus.csv",
        F_bar="data/simulation/scope_input_kmeans_known_k/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/F_bar.csv",
        clone_sizes="data/simulation/scope_input_kmeans_known_k/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/clone_sizes.csv",
        kmeans_clones="data/simulation/scope_input_kmeans_known_k/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/kmeans_clones.csv",
    log:
        std="data/simulation/scope_input_kmeans_known_k/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/log",
        err="data/simulation/scope_input_kmeans_known_k/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/err.log",
    benchmark: "data/simulation/scope_input_kmeans_known_k/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/benchmark",
    shell:
        "python src/data_processors/kmeans_clustering_for_scope.py -i data/simulation/ground_truth/n{wildcards.ncells}_m{wildcards.n_mutation_groups}_size{wildcards.mutation_group_sizes}_cov{wildcards.cov}_p{wildcards.nclusters}_s{wildcards.seed}/sim "
        " -o data/simulation/scope_input_kmeans_known_k/ --start_k {wildcards.n_mutation_groups} --end_k {wildcards.n_mutation_groups} "
        " > {log.std} 2> {log.err}"

rule scope_on_simulation:
    input:
        F_plus="data/simulation/scope_input/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/F_plus.csv",
        F_minus="data/simulation/scope_input/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/F_minus.csv",
        clone_sizes="data/simulation/scope_input/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/clone_sizes.csv",
    output:
        solution_summary="data/simulation/scope_output/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/summary.txt",
    log:
        std="data/simulation/scope_output/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/log",
        err="data/simulation/scope_output/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/err.log",
    benchmark: "data/simulation/scope_output/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/benchmark",
    shell:
        "python -m src.benchmark.run_scope_on_simulated_data -i data/simulation/scope_input/n{wildcards.ncells}_m{wildcards.n_mutation_groups}_size{wildcards.mutation_group_sizes}_cov{wildcards.cov}_p{wildcards.nclusters}_s{wildcards.seed} "
        " -o data/simulation/scope_output "
        " > {log.std} 2> {log.err}"

rule scope_kmeans_on_simulation:
    input:
        F_plus="data/simulation/scope_input_kmeans/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/F_plus.csv",
        F_minus="data/simulation/scope_input_kmeans/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/F_minus.csv",
        clone_sizes="data/simulation/scope_input_kmeans/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/clone_sizes.csv",
    output:
        solution_summary="data/simulation/scope_output_kmeans/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/summary.txt",
    log:
        std="data/simulation/scope_output_kmeans/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/log",
        err="data/simulation/scope_output_kmeans/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/err.log",
    benchmark: "data/simulation/scope_output_kmeans/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/benchmark",
    shell:
        "python -m src.benchmark.run_scope_on_simulated_data -i data/simulation/scope_input_kmeans/n{wildcards.ncells}_m{wildcards.n_mutation_groups}_size{wildcards.mutation_group_sizes}_cov{wildcards.cov}_p{wildcards.nclusters}_s{wildcards.seed} "
        " -o data/simulation/scope_output_kmeans "
        " > {log.std} 2> {log.err}"

rule scope_kmeans_known_k_on_simulation:
    input:
        F_plus="data/simulation/scope_input_kmeans_known_k/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/F_plus.csv",
        F_minus="data/simulation/scope_input_kmeans_known_k/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/F_minus.csv",
        clone_sizes="data/simulation/scope_input_kmeans_known_k/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/clone_sizes.csv",
    output:
        solution_summary="data/simulation/scope_output_kmeans_known_k/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/summary.txt",
    log:
        std="data/simulation/scope_output_kmeans_known_k/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/log",
        err="data/simulation/scope_output_kmeans_known_k/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/err.log",
    benchmark: "data/simulation/scope_output_kmeans_known_k/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/benchmark",
    shell:
        "python -m src.benchmark.run_scope_on_simulated_data -i data/simulation/scope_input_kmeans_known_k/n{wildcards.ncells}_m{wildcards.n_mutation_groups}_size{wildcards.mutation_group_sizes}_cov{wildcards.cov}_p{wildcards.nclusters}_s{wildcards.seed} "
        " -o data/simulation/scope_output_kmeans_known_k "
        " > {log.std} 2> {log.err}"

rule post_process_clones:
    input:
        F_bar="data/simulation/scope_input_kmeans_known_k/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/F_bar.csv",
    output:
        corrected_clones="data/simulation/scope_post_kmeans_known_k/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/kmeans_cleaned_clones.csv",
    log:
        std="data/simulation/scope_post_kmeans_known_k/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/log",
        err="data/simulation/scope_post_kmeans_known_k/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/err.log",
    benchmark: "data/simulation/scope_post_kmeans_known_k/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/benchmark",
    shell:
        "python src/data_processors/post_process_mutation_groups.py -s n{wildcards.ncells}_m{wildcards.n_mutation_groups}_size{wildcards.mutation_group_sizes}_cov{wildcards.cov}_p{wildcards.nclusters}_s{wildcards.seed} "
        " > {log.std} 2> {log.err}"

rule post_process_clones_2:
    input:
        F_bar="data/simulation/scope_input_kmeans/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/F_bar.csv",
    output:
        corrected_clones="data/simulation/scope_post_kmeans/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/kmeans_cleaned_clones.csv",
    log:
        std="data/simulation/scope_post_kmeans/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/log",
        err="data/simulation/scope_post_kmeans/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/err.log",
    benchmark: "data/simulation/scope_post_kmeans/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/benchmark",
    shell:
        "python src/data_processors/post_process_mutation_groups_2.py -s n{wildcards.ncells}_m{wildcards.n_mutation_groups}_size{wildcards.mutation_group_sizes}_cov{wildcards.cov}_p{wildcards.nclusters}_s{wildcards.seed} "
        " > {log.std} 2> {log.err}"

rule phertilizer_on_simulation:
    input:
        snv_counts="data/simulation/phertilizer_input/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/snv_counts.tsv",
        binned_counts="data/simulation/phertilizer_input/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/binned_read_counts.csv",
    output:
        phertilizer_cell_clusters="data/simulation/phertilizer_output/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/cell_clusters.csv",
        phertilizer_snv_clusters="data/simulation/phertilizer_output/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/snv_clusters.csv",
        phertilizer_tree_txt="data/simulation/phertilizer_output/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/tree.txt",
        phertilizer_tree_image="data/simulation/phertilizer_output/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/tree.png",
    conda:
        "phertilizer"
    log:
        std="data/simulation/phertilizer_output/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/log",
        err="data/simulation/phertilizer_output/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/err.log",
    benchmark: "data/simulation/phertilizer_output/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/benchmark",
    shell:
        "phertilizer -f data/simulation/phertilizer_input/n{wildcards.ncells}_m{wildcards.n_mutation_groups}_size{wildcards.mutation_group_sizes}_cov{wildcards.cov}_p{wildcards.nclusters}_s{wildcards.seed}/snv_counts.tsv "
        " --bin_count_data data/simulation/phertilizer_input/n{wildcards.ncells}_m{wildcards.n_mutation_groups}_size{wildcards.mutation_group_sizes}_cov{wildcards.cov}_p{wildcards.nclusters}_s{wildcards.seed}/binned_read_counts.csv "
        " --tree data/simulation/phertilizer_output/n{wildcards.ncells}_m{wildcards.n_mutation_groups}_size{wildcards.mutation_group_sizes}_cov{wildcards.cov}_p{wildcards.nclusters}_s{wildcards.seed}/tree.png "
        " --tree_text data/simulation/phertilizer_output/n{wildcards.ncells}_m{wildcards.n_mutation_groups}_size{wildcards.mutation_group_sizes}_cov{wildcards.cov}_p{wildcards.nclusters}_s{wildcards.seed}/tree.txt "
        " -n data/simulation/phertilizer_output/n{wildcards.ncells}_m{wildcards.n_mutation_groups}_size{wildcards.mutation_group_sizes}_cov{wildcards.cov}_p{wildcards.nclusters}_s{wildcards.seed}/cell_clusters.csv "
        " -m data/simulation/phertilizer_output/n{wildcards.ncells}_m{wildcards.n_mutation_groups}_size{wildcards.mutation_group_sizes}_cov{wildcards.cov}_p{wildcards.nclusters}_s{wildcards.seed}/snv_clusters.csv "
        " -j 10 -s 4 "
        " > {log.std} 2> {log.err}"

rule williams:
    output:
        williams_scope="data/williams/scratch/{sample_ids}/scope_mut/summary.txt",
    log:
        std="data/williams/scratch/{sample_ids}/scope_mut/log",
        err="data/williams/scratch/{sample_ids}/scope_mut/err.log",
    benchmark: "data/williams/scratch/{sample_ids}/scope_mut/benchmark",
    shell:
        "python -m src.run_williams --sample {wildcards.sample_ids} "
        " > {log.std} 2> {log.err}"

rule scope_laks:
    input:
        loh="data/laks/scope/loh-counts.csv",
        F_plus="data/laks/scope/cell_fractions/k_{k}_t_{t}/F_plus.csv",
        F_minus="data/laks/scope/cell_fractions/k_{k}_t_{t}/F_minus.csv",
        labels="data/laks/scope/kmeans_labels/k_{k}.csv",
    output:
        summary_file_laks="data/laks/scope/outputs/k_{k}_t_{t}/summary.txt"
    log:
        std="data/laks/scope/outputs/k_{k}_t_{t}/log",
        err="data/laks/scope/outputs/k_{k}_t_{t}/err.log",
    benchmark: "data/laks/scope/outputs/k_{k}_t_{t}/benchmark",
    shell:
        """
        python -m src.laks.run_scope \
            --loh_conflicts_filepath {input.loh} \
            --F_plus_filepath {input.F_plus} \
            --F_minus_filepath {input.F_minus} \
            --mutation_labels_filepath {input.labels} \
            --output_directory data/laks/scope/outputs/k_{wildcards.k}_t_{wildcards.t} \
            > {log.std} 2> {log.err}
        """

rule scope_laks_filtered:
    input:
        loh="data/laks/scope/loh-counts.csv",
        F_plus="data/laks/scope/cell_fractions/filtered/k_{k}_t_{t}_f_{f}/F_plus.csv",
        F_minus="data/laks/scope/cell_fractions/filtered/k_{k}_t_{t}_f_{f}/F_minus.csv",
        labels="data/laks/scope/kmeans_labels/filtered/k_{k}_t_{t}_f_{f}.csv",
    output:
        summary_file_laks_filtered="data/laks/scope/outputs/filtered/k_{k}_t_{t}_f_{f}/summary.txt"
    log:
        std="data/laks/scope/outputs/filtered/k_{k}_t_{t}_f_{f}/log",
        err="data/laks/scope/outputs/filtered/k_{k}_t_{t}_f_{f}/err.log",
    benchmark: "data/laks/scope/outputs/filtered/k_{k}_t_{t}_f_{f}/benchmark",
    shell:
        """
        python -m src.laks.run_scope \
            --loh_conflicts_filepath {input.loh} \
            --F_plus_filepath {input.F_plus} \
            --F_minus_filepath {input.F_minus} \
            --mutation_labels_filepath {input.labels} \
            --output_directory data/laks/scope/outputs/filtered/k_{wildcards.k}_t_{wildcards.t}_f_{wildcards.f} \
            > {log.std} 2> {log.err}
        """

