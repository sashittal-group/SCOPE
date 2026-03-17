configfile: "config.yaml"

seeds = [ i for i in range(config["nseeds"])]
Ks = [i for i in range(10, 21)]
lambdas = [f"{int(l*100):02d}" for l in config['lambda']]

import itertools


rule all:
    input:
        # # simulation
        # expand('data/simulation/ground_truth/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_multi_state_tree_node_character_matrix.parquet', \
        #     itertools.product, seed=seeds, ncells=config['ncells'], n_mutation_groups=config['n_mutation_groups'], mutation_group_sizes=config['mutation_group_sizes'], \
        #     nclusters=config['nclusters'], cov=config['coverage']),
        # # simulation with loss
        # expand('data/simulation/ground_truth_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_multi_state_tree_node_character_matrix.parquet', \
        #     itertools.product, seed=seeds, ncells=config['ncells'], n_mutation_groups=config['n_mutation_groups'], mutation_group_sizes=config['mutation_group_sizes'], \
        #     nclusters=config['nclusters'], cov=config['coverage'], lambda_pct=lambdas),
        # # scope input
        # expand("data/simulation/scope_input/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/F_plus.csv", \
        #     itertools.product, seed=seeds, ncells=config['ncells'], n_mutation_groups=config['n_mutation_groups'], mutation_group_sizes=config['mutation_group_sizes'], \
        #     nclusters=config['nclusters'], cov=config['coverage']),
        # # scope input kmeans
        # expand("data/simulation/scope_input_kmeans/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/kmeans_clones.csv", \
        #     itertools.product, seed=seeds, ncells=config['ncells'], n_mutation_groups=config['n_mutation_groups'], mutation_group_sizes=config['mutation_group_sizes'], \
        #     nclusters=config['nclusters'], cov=config['coverage']),
        # scope input kmeans with loss
        expand("data/simulation/scope_input_kmeans_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/kmeans_clones.csv", \
            itertools.product, seed=seeds, ncells=config['ncells'], n_mutation_groups=config['n_mutation_groups'], mutation_group_sizes=config['mutation_group_sizes'], \
            nclusters=config['nclusters'], cov=config['coverage'], lambda_pct=lambdas),
        # # scope input kmeans known k
        # expand("data/simulation/scope_input_kmeans_known_k/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/kmeans_clones.csv", \
        #     itertools.product, seed=seeds, ncells=config['ncells'], n_mutation_groups=config['n_mutation_groups'], mutation_group_sizes=config['mutation_group_sizes'], \
        #     nclusters=config['nclusters'], cov=config['coverage']),
        # # phertilizer input
        # expand("data/simulation/phertilizer_input/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/snv_counts.tsv", \
        #     itertools.product, seed=seeds, ncells=config['ncells'], n_mutation_groups=config['n_mutation_groups'], mutation_group_sizes=config['mutation_group_sizes'], \
        #     nclusters=config['nclusters'], cov=config['coverage']),
        # # phertilizer input with loss
        # expand("data/simulation/phertilizer_input_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/snv_counts.tsv", \
        #     itertools.product, seed=seeds, ncells=config['ncells'], n_mutation_groups=config['n_mutation_groups'], mutation_group_sizes=config['mutation_group_sizes'], \
        #     nclusters=config['nclusters'], cov=config['coverage'], lambda_pct=lambdas),
        # # pharming input
        # expand("data/simulation/pharming_input/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/read_counts.tsv", \
        #     itertools.product, seed=seeds, ncells=config['ncells'], n_mutation_groups=config['n_mutation_groups'], mutation_group_sizes=config['mutation_group_sizes'], \
        #     nclusters=config['nclusters'], cov=config['coverage']),
        # # pharming input with loss
        # expand("data/simulation/pharming_input_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/read_counts.tsv", \
        #     itertools.product, seed=seeds, ncells=config['ncells'], n_mutation_groups=config['n_mutation_groups'], mutation_group_sizes=config['mutation_group_sizes'], \
        #     nclusters=config['nclusters'], cov=config['coverage'], lambda_pct=lambdas),
        # # sbmclone input
        # expand("data/simulation/sbmclone_input/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/matrix.csv", \
        #     itertools.product, seed=seeds, ncells=config['ncells'], n_mutation_groups=config['n_mutation_groups'], mutation_group_sizes=config['mutation_group_sizes'], \
        #     nclusters=config['nclusters'], cov=config['coverage']),
        # # scope on simulation
        # expand("data/simulation/scope_output/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/summary.txt", \
        #     itertools.product, seed=seeds, ncells=config['ncells'], n_mutation_groups=config['n_mutation_groups'], mutation_group_sizes=config['mutation_group_sizes'], \
        #     nclusters=config['nclusters'], cov=config['coverage']),
        # # scope kmeans on simulation
        # expand("data/simulation/scope_output_kmeans/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/summary.txt", \
        #     itertools.product, seed=seeds, ncells=config['ncells'], n_mutation_groups=config['n_mutation_groups'], mutation_group_sizes=config['mutation_group_sizes'], \
        #     nclusters=config['nclusters'], cov=config['coverage']),
        # scope kmeans on simulation with loss
        expand("data/simulation/scope_output_kmeans_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/summary.txt", \
            itertools.product, seed=seeds, ncells=config['ncells'], n_mutation_groups=config['n_mutation_groups'], mutation_group_sizes=config['mutation_group_sizes'], \
            nclusters=config['nclusters'], cov=config['coverage'], lambda_pct=lambdas),
        # # scope kmeans known k on simulation
        # expand("data/simulation/scope_output_kmeans_known_k/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/summary.txt", \
        #     itertools.product, seed=seeds, ncells=config['ncells'], n_mutation_groups=config['n_mutation_groups'], mutation_group_sizes=config['mutation_group_sizes'], \
        #     nclusters=config['nclusters'], cov=config['coverage']),
        # # phertilizer on simulation
        # expand("data/simulation/phertilizer_output/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/cell_clusters.csv", \
        #     itertools.product, seed=seeds, ncells=config['ncells'], n_mutation_groups=config['n_mutation_groups'], mutation_group_sizes=config['mutation_group_sizes'], \
        #     nclusters=config['nclusters'], cov=config['coverage']),
        # # phertilizer on simulation with loss
        # expand("data/simulation/phertilizer_output_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/cell_clusters.csv", \
        #     itertools.product, seed=seeds, ncells=config['ncells'], n_mutation_groups=config['n_mutation_groups'], mutation_group_sizes=config['mutation_group_sizes'], \
        #     nclusters=config['nclusters'], cov=config['coverage'], lambda_pct=lambdas),
        # pharming on simulation
        # expand("data/simulation/pharming_output/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/solutions.pkl", \
        #     itertools.product, seed=seeds, ncells=config['ncells'], n_mutation_groups=config['n_mutation_groups'], mutation_group_sizes=config['mutation_group_sizes'], \
        #     nclusters=config['nclusters'], cov=config['coverage']),
        # pharming on simulation with loss
        # expand("data/simulation/pharming_output_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/solutions.pkl", \
        #     itertools.product, seed=seeds, ncells=config['ncells'], n_mutation_groups=config['n_mutation_groups'], mutation_group_sizes=config['mutation_group_sizes'], \
        #     nclusters=config['nclusters'], cov=config['coverage'], lambda_pct=lambdas),
        # # sbmclone on simulation
        # expand("data/simulation/sbmclone_output/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/blockmatrix.csv", \
        #     itertools.product, seed=seeds, ncells=config['ncells'], n_mutation_groups=config['n_mutation_groups'], mutation_group_sizes=config['mutation_group_sizes'], \
        #     nclusters=config['nclusters'], cov=config['coverage']),
        # # post process clone
        # expand("data/simulation/scope_post_kmeans_known_k/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/kmeans_cleaned_clones.csv", \
        #     itertools.product, seed=seeds, ncells=config['ncells'], n_mutation_groups=config['n_mutation_groups'], mutation_group_sizes=config['mutation_group_sizes'], \
        #     nclusters=config['nclusters'], cov=config['coverage']),
        # # post process clone 2
        # expand("data/simulation/scope_post_kmeans/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/kmeans_cleaned_clones.csv", \
        #     itertools.product, seed=seeds, ncells=config['ncells'], n_mutation_groups=config['n_mutation_groups'], mutation_group_sizes=config['mutation_group_sizes'], \
        #     nclusters=config['nclusters'], cov=config['coverage']),
        # post process clone 2 with loss
        expand("data/simulation/scope_post_kmeans_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/kmeans_cleaned_clones.csv", \
            itertools.product, seed=seeds, ncells=config['ncells'], n_mutation_groups=config['n_mutation_groups'], mutation_group_sizes=config['mutation_group_sizes'], \
            nclusters=config['nclusters'], cov=config['coverage'], lambda_pct=lambdas),
        # # split inputs williams
        # expand("data/williams/scratch/separate.txt"),
        # # make input williams
        # expand("outputs/scope/williams/{sample_ids}/mutation_clusters/k_10/F_plus.csv", \
        #     itertools.product, sample_ids=config['sample_ids']),
        # # post williams
        # expand("outputs/scope/williams/{sample_ids}/solutions_mut/k_{Ks}/summary.txt", \
        #     itertools.product, sample_ids=config['sample_ids'], Ks=Ks),
        # split inputs williams
        expand("data/williams/scratch/separate.txt"),
        # make input williams
        expand("outputs/scope/williams/{sample_ids}/mutation_clusters/k_10/F_plus.csv", \
            itertools.product, sample_ids=config['sample_ids']),
        # post williams
        expand("outputs/scope/williams/{sample_ids}/solutions_mut/k_{Ks}/summary.txt", \
            itertools.product, sample_ids=config['sample_ids'], Ks=Ks),


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
        fn=config["fn"],
        fp=config["fp"],
        mutation_rate=0.50,
        vaf_threshold=config["vaf_threshold"],
        read_threshold=config["read_threshold"],
        missing=config["missing_rate"],
    log:
        std="data/simulation/ground_truth/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim.log",
        err="data/simulation/ground_truth/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim.err.log",
    benchmark:
        "data/simulation/ground_truth/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim.benchmark",
    shell:
        (
            "python -m src.simulation.simulate_data "
            "-s {wildcards.seed} "
            "-o {params.prefix} "
            "-n {wildcards.ncells} "
            "-m {wildcards.n_mutation_groups} "
            "--size {wildcards.mutation_group_sizes} "
            "-d {params.missing} "
            "-p {wildcards.nclusters} "
            "-a {params.fp} "
            "-b {params.fn} "
            "-l {params.mutation_rate} "
            "--cov {wildcards.cov} "
            "--maxcn 6 "
            "> {log.std} 2> {log.err}"
        )

rule simulate_with_loss:
    output:
        readcount_matrix="data/simulation/ground_truth_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_read_count.parquet",
        variantcount_matrix="data/simulation/ground_truth_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_variant_count.parquet",
        character_matrix="data/simulation/ground_truth_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_character_matrix_without_noise.parquet",
        copy_number_table="data/simulation/ground_truth_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_copy_numbers.parquet",
        mutation_group_table="data/simulation/ground_truth_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_mutation_group.parquet",
        mutation_to_bin_table="data/simulation/ground_truth_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_mutation_to_bin_mapping.parquet",
        B_truth="data/simulation/ground_truth_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_multi_state_tree_node_character_matrix.parquet",

    params:
        prefix="data/simulation/ground_truth_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim",
        fn=config["fn"],
        fp=config["fp"],
        vaf_threshold=config["vaf_threshold"],
        read_threshold=config["read_threshold"],
        missing=config["missing_rate"],
        lambda_val=lambda wc: float(wc.lambda_pct) / 100,

    log:
        std="data/simulation/ground_truth_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim.log",
        err="data/simulation/ground_truth_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim.err.log",

    benchmark:
        "data/simulation/ground_truth_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim.benchmark",

    shell:
        (
            "python -m src.simulation.simulate_data "
            "-s {wildcards.seed} "
            "-o {params.prefix} "
            "-n {wildcards.ncells} "
            "-m {wildcards.n_mutation_groups} "
            "--size {wildcards.mutation_group_sizes} "
            "-d {params.missing} "
            "-p {wildcards.nclusters} "
            "-a {params.fp} "
            "-b {params.fn} "
            "-l {params.lambda_val} "
            "--cov {wildcards.cov} "
            "--maxcn 6 "
            "--mutation-loss "
            "> {log.std} 2> {log.err}"
        )


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
        "python -m src.simulation.process_simulated_data_for_phertilizer -i data/simulation/ground_truth/n{wildcards.ncells}_m{wildcards.n_mutation_groups}_size{wildcards.mutation_group_sizes}_cov{wildcards.cov}_p{wildcards.nclusters}_s{wildcards.seed}/sim "
        " -o data/simulation/phertilizer_input/ "
        " > {log.std} 2> {log.err}"


rule phertilizer_input_with_loss:
    input:
        readcount_matrix="data/simulation/ground_truth_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_read_count.parquet",
        variantcount_matrix="data/simulation/ground_truth_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_variant_count.parquet",
        mutation_to_bin_table="data/simulation/ground_truth_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_mutation_to_bin_mapping.parquet",
    output:
        snv_counts="data/simulation/phertilizer_input_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/snv_counts.tsv",
        binned_counts="data/simulation/phertilizer_input_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/binned_read_counts.csv",
    resources:
        mem_mb=10240
    log:
        std="data/simulation/phertilizer_input_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/stdout.log",
        err="data/simulation/phertilizer_input_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/stderr.log",
    benchmark:
        "data/simulation/phertilizer_input_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/benchmark.txt"
    shell:
        """
        python -m src.simulation.process_simulated_data_for_phertilizer \
            -i data/simulation/ground_truth_with_loss_{wildcards.lambda_pct}/n{wildcards.ncells}_m{wildcards.n_mutation_groups}_size{wildcards.mutation_group_sizes}_cov{wildcards.cov}_p{wildcards.nclusters}_s{wildcards.seed}/sim \
            -o data/simulation/phertilizer_input_with_loss_{wildcards.lambda_pct}/ \
            > {log.std} 2> {log.err}
        """


rule pharming_input:
    input:
        readcount_matrix="data/simulation/ground_truth/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_read_count.parquet",
        variantcount_matrix="data/simulation/ground_truth/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_variant_count.parquet",
        mutation_to_bin_table="data/simulation/ground_truth/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_mutation_to_bin_mapping.parquet",
    output:
        read_counts="data/simulation/pharming_input/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/read_counts.tsv",
        copy_numbers="data/simulation/pharming_input/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/copy_numbers.csv",
        segments="data/simulation/pharming_input/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/segments.txt",
    resources:
        mem_mb=10240
    log:
        std="data/simulation/pharming_input/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/log",
        err="data/simulation/pharming_input/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/err.log",
    benchmark: "data/simulation/pharming_input/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/benchmark",
    shell:
        "python -m src.simulation.process_simulated_data_for_pharming -i data/simulation/ground_truth/n{wildcards.ncells}_m{wildcards.n_mutation_groups}_size{wildcards.mutation_group_sizes}_cov{wildcards.cov}_p{wildcards.nclusters}_s{wildcards.seed}/sim "
        " -o data/simulation/pharming_input/ "
        " > {log.std} 2> {log.err}"


rule pharming_input_with_loss:
    input:
        readcount_matrix="data/simulation/ground_truth_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_read_count.parquet",
        variantcount_matrix="data/simulation/ground_truth_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_variant_count.parquet",
        mutation_to_bin_table="data/simulation/ground_truth_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_mutation_to_bin_mapping.parquet",

    output:
        read_counts="data/simulation/pharming_input_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/read_counts.tsv",
        copy_numbers="data/simulation/pharming_input_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/copy_numbers.csv",
        segments="data/simulation/pharming_input_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/segments.txt",

    params:
        prefix="data/simulation/ground_truth_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim",
        outdir="data/simulation/pharming_input_with_loss_{lambda_pct}/",

    resources:
        mem_mb=10240

    log:
        std="data/simulation/pharming_input_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/pharming_input.log",
        err="data/simulation/pharming_input_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/pharming_input.err.log",

    benchmark:
        "data/simulation/pharming_input_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/pharming_input.benchmark",

    shell:
        (
            "python -m src.simulation.process_simulated_data_for_pharming "
            "-i {params.prefix} "
            "-o {params.outdir} "
            "> {log.std} 2> {log.err}"
        )


rule sbmclone_input:
    input:
        readcount_matrix="data/simulation/ground_truth/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_read_count.parquet",
        variantcount_matrix="data/simulation/ground_truth/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_variant_count.parquet",
    output:
        matrix="data/simulation/sbmclone_input/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/matrix.csv",
    params:
        threshold=0.1
    log:
        std="data/simulation/sbmclone_input/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/log",
        err="data/simulation/sbmclone_input/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/err.log",
    benchmark: "data/simulation/sbmclone_input/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/benchmark",
    shell:
        "python -m src.simulation.process_simulated_data_for_sbmclone -i data/simulation/ground_truth/n{wildcards.ncells}_m{wildcards.n_mutation_groups}_size{wildcards.mutation_group_sizes}_cov{wildcards.cov}_p{wildcards.nclusters}_s{wildcards.seed}/sim "
        " -o data/simulation/sbmclone_input/ --threshold {params.threshold} "
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
        "python -m src.simulation.process_simulated_data_for_scope -i data/simulation/ground_truth/n{wildcards.ncells}_m{wildcards.n_mutation_groups}_size{wildcards.mutation_group_sizes}_cov{wildcards.cov}_p{wildcards.nclusters}_s{wildcards.seed}/sim "
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
        "python -m src.simulation.kmeans_clustering_for_scope -i data/simulation/ground_truth/n{wildcards.ncells}_m{wildcards.n_mutation_groups}_size{wildcards.mutation_group_sizes}_cov{wildcards.cov}_p{wildcards.nclusters}_s{wildcards.seed}/sim "
        " -o data/simulation/scope_input_kmeans/ "
        " > {log.std} 2> {log.err}"


rule scope_input_kmeans_with_loss:
    input:
        character_matrix="data/simulation/ground_truth_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_character_matrix_without_noise.parquet",
        readcount_matrix="data/simulation/ground_truth_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_read_count.parquet",
        variantcount_matrix="data/simulation/ground_truth_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_variant_count.parquet",
        copy_number_table="data/simulation/ground_truth_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_copy_numbers.parquet",
    output:
        F_plus="data/simulation/scope_input_kmeans_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/F_plus.csv",
        F_minus="data/simulation/scope_input_kmeans_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/F_minus.csv",
        F_bar="data/simulation/scope_input_kmeans_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/F_bar.csv",
        clone_sizes="data/simulation/scope_input_kmeans_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/clone_sizes.csv",
        kmeans_clones="data/simulation/scope_input_kmeans_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/kmeans_clones.csv",
    log:
        std="data/simulation/scope_input_kmeans_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/stdout.log",
        err="data/simulation/scope_input_kmeans_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/stderr.log",
    benchmark:
        "data/simulation/scope_input_kmeans_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/benchmark.txt",
    shell:
        """
        python -m src.simulation.kmeans_clustering_for_scope \
            -i data/simulation/ground_truth_with_loss_{wildcards.lambda_pct}/n{wildcards.ncells}_m{wildcards.n_mutation_groups}_size{wildcards.mutation_group_sizes}_cov{wildcards.cov}_p{wildcards.nclusters}_s{wildcards.seed}/sim \
            -o data/simulation/scope_input_kmeans_with_loss_{wildcards.lambda_pct}/ \
            > {log.std} 2> {log.err}
        """

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
        "python -m src.simulation.kmeans_clustering_for_scope -i data/simulation/ground_truth/n{wildcards.ncells}_m{wildcards.n_mutation_groups}_size{wildcards.mutation_group_sizes}_cov{wildcards.cov}_p{wildcards.nclusters}_s{wildcards.seed}/sim "
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
        "python -m src.simulation.run_scope_on_simulated_data -i data/simulation/scope_input/n{wildcards.ncells}_m{wildcards.n_mutation_groups}_size{wildcards.mutation_group_sizes}_cov{wildcards.cov}_p{wildcards.nclusters}_s{wildcards.seed} "
        " -o data/simulation/scope_output "
        " > {log.std} 2> {log.err}"

rule scope_kmeans_on_simulation:
    input:
        F_plus="data/simulation/scope_input_kmeans/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/F_plus.csv",
        F_minus="data/simulation/scope_input_kmeans/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/F_minus.csv",
        clone_sizes="data/simulation/scope_input_kmeans/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/clone_sizes.csv",
    output:
        solution_summary="data/simulation/scope_output_kmeans/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/summary.txt",
        X_kmeans="data/simulation/scope_output_kmeans/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/solution_0/X.csv",
    log:
        std="data/simulation/scope_output_kmeans/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/log",
        err="data/simulation/scope_output_kmeans/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/err.log",
    benchmark:
        "data/simulation/scope_output_kmeans/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/benchmark",
    shell:
        (
            "python -m src.simulation.run_scope_on_simulated_data "
            "-i data/simulation/scope_input_kmeans/n{wildcards.ncells}_m{wildcards.n_mutation_groups}_size{wildcards.mutation_group_sizes}_cov{wildcards.cov}_p{wildcards.nclusters}_s{wildcards.seed} "
            " -o data/simulation/scope_output_kmeans/ "
            "> {log.std} 2> {log.err}"
        )


rule scope_kmeans_on_simulation_with_loss:
    input:
        F_plus="data/simulation/scope_input_kmeans_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/F_plus.csv",
        F_minus="data/simulation/scope_input_kmeans_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/F_minus.csv",
        clone_sizes="data/simulation/scope_input_kmeans_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/clone_sizes.csv",
    output:
        solution_summary="data/simulation/scope_output_kmeans_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/summary.txt",
        X_kmeans="data/simulation/scope_output_kmeans_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/solution_0/X.csv",
    log:
        std="data/simulation/scope_output_kmeans_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/stdout.log",
        err="data/simulation/scope_output_kmeans_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/stderr.log",
    benchmark:
        "data/simulation/scope_output_kmeans_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/benchmark.txt",
    shell:
        """
        python -m src.simulation.run_scope_on_simulated_data \
            -i data/simulation/scope_input_kmeans_with_loss_{wildcards.lambda_pct}/n{wildcards.ncells}_m{wildcards.n_mutation_groups}_size{wildcards.mutation_group_sizes}_cov{wildcards.cov}_p{wildcards.nclusters}_s{wildcards.seed} \
            -o data/simulation/scope_output_kmeans_with_loss_{wildcards.lambda_pct}/ \
            > {log.std} 2> {log.err}
        """


rule scope_kmeans_known_k_on_simulation:
    input:
        F_plus="data/simulation/scope_input_kmeans_known_k/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/F_plus.csv",
        F_minus="data/simulation/scope_input_kmeans_known_k/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/F_minus.csv",
        clone_sizes="data/simulation/scope_input_kmeans_known_k/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/clone_sizes.csv",
    output:
        solution_summary="data/simulation/scope_output_kmeans_known_k/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/summary.txt",
        X_kmeans_known_K="data/simulation/scope_output_kmeans_known_k/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/solution_0/X.csv",
    log:
        std="data/simulation/scope_output_kmeans_known_k/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/log",
        err="data/simulation/scope_output_kmeans_known_k/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/err.log",
    benchmark: "data/simulation/scope_output_kmeans_known_k/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/benchmark",
    shell:
        "python -m src.simulation.run_scope_on_simulated_data -i data/simulation/scope_input_kmeans_known_k/n{wildcards.ncells}_m{wildcards.n_mutation_groups}_size{wildcards.mutation_group_sizes}_cov{wildcards.cov}_p{wildcards.nclusters}_s{wildcards.seed} "
        " -o data/simulation/scope_output_kmeans_known_k "
        " > {log.std} 2> {log.err}"

rule post_process_clones:
    input:
        F_bar="data/simulation/scope_input_kmeans_known_k/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/F_bar.csv",
        X_kmeans_known_K="data/simulation/scope_output_kmeans_known_k/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/solution_0/X.csv",
    output:
        corrected_clones="data/simulation/scope_post_kmeans_known_k/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/kmeans_cleaned_clones.csv",
    log:
        std="data/simulation/scope_post_kmeans_known_k/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/log",
        err="data/simulation/scope_post_kmeans_known_k/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/err.log",
    benchmark: "data/simulation/scope_post_kmeans_known_k/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/benchmark",
    shell:
        "python -m src.simulation.post_process_mutation_groups -s n{wildcards.ncells}_m{wildcards.n_mutation_groups}_size{wildcards.mutation_group_sizes}_cov{wildcards.cov}_p{wildcards.nclusters}_s{wildcards.seed} "
        " > {log.std} 2> {log.err}"

rule post_process_clones_2:
    input:
        F_bar="data/simulation/scope_input_kmeans/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/F_bar.csv",
        X_kmeans="data/simulation/scope_output_kmeans/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/solution_0/X.csv",
    output:
        corrected_clones="data/simulation/scope_post_kmeans/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/kmeans_cleaned_clones.csv",
    log:
        std="data/simulation/scope_post_kmeans/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/log",
        err="data/simulation/scope_post_kmeans/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/err.log",
    benchmark: "data/simulation/scope_post_kmeans/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/benchmark",
    shell:
        "python -m src.simulation.post_process_mutation_groups_2 -s n{wildcards.ncells}_m{wildcards.n_mutation_groups}_size{wildcards.mutation_group_sizes}_cov{wildcards.cov}_p{wildcards.nclusters}_s{wildcards.seed} "
        " > {log.std} 2> {log.err}"

rule post_process_clones_2_with_loss:
    input:
        F_bar="data/simulation/scope_input_kmeans_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/F_bar.csv",
        X_kmeans="data/simulation/scope_output_kmeans_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/solution_0/X.csv",
    output:
        corrected_clones="data/simulation/scope_post_kmeans_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/kmeans_cleaned_clones.csv",
    log:
        std="data/simulation/scope_post_kmeans_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/stdout.log",
        err="data/simulation/scope_post_kmeans_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/stderr.log",
    benchmark:
        "data/simulation/scope_post_kmeans_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/benchmark.txt",
    shell:
        """
        python -m src.simulation.post_process_mutation_groups_2 \
            -s n{wildcards.ncells}_m{wildcards.n_mutation_groups}_size{wildcards.mutation_group_sizes}_cov{wildcards.cov}_p{wildcards.nclusters}_s{wildcards.seed} \
            -i scope_input_kmeans_with_loss_{wildcards.lambda_pct} \
            -o scope_output_kmeans_with_loss_{wildcards.lambda_pct} \
            -p scope_post_kmeans_with_loss_{wildcards.lambda_pct} \
            -g ground_truth_with_loss_{wildcards.lambda_pct} \
            > {log.std} 2> {log.err}
        """

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


rule phertilizer_on_simulation_with_loss:
    input:
        snv_counts="data/simulation/phertilizer_input_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/snv_counts.tsv",
        binned_counts="data/simulation/phertilizer_input_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/binned_read_counts.csv",
    output:
        phertilizer_cell_clusters="data/simulation/phertilizer_output_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/cell_clusters.csv",
        phertilizer_snv_clusters="data/simulation/phertilizer_output_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/snv_clusters.csv",
        phertilizer_tree_txt="data/simulation/phertilizer_output_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/tree.txt",
        phertilizer_tree_image="data/simulation/phertilizer_output_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/tree.png",
    conda:
        "phertilizer"
    log:
        std="data/simulation/phertilizer_output_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/stdout.log",
        err="data/simulation/phertilizer_output_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/stderr.log",
    benchmark:
        "data/simulation/phertilizer_output_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/benchmark.txt",
    shell:
        """
        phertilizer \
            -f {input.snv_counts} \
            --bin_count_data {input.binned_counts} \
            --tree {output.phertilizer_tree_image} \
            --tree_text {output.phertilizer_tree_txt} \
            -n {output.phertilizer_cell_clusters} \
            -m {output.phertilizer_snv_clusters} \
            -j 10 -s 4 \
            > {log.std} 2> {log.err}
        """


rule pharming_on_simulation:
    input:
        read_counts="data/simulation/pharming_input/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/read_counts.tsv",
        copy_numbers="data/simulation/pharming_input/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/copy_numbers.csv",
        segments="data/simulation/pharming_input/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/segments.txt",
    output:
        pharming_solution_pkl="data/simulation/pharming_output/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/solutions.pkl",
        pharming_tree_png="data/simulation/pharming_output/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/tree.png",
        pharming_labels="data/simulation/pharming_output/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/labels.csv",
    conda:
        "pharming"
    log:
        std="data/simulation/pharming_output/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/log",
        err="data/simulation/pharming_output/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/err.log",
    benchmark: "data/simulation/pharming_output/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/benchmark",
    shell:
        "pharming -f data/simulation/pharming_input/n{wildcards.ncells}_m{wildcards.n_mutation_groups}_size{wildcards.mutation_group_sizes}_cov{wildcards.cov}_p{wildcards.nclusters}_s{wildcards.seed}/read_counts.tsv "
        " -c data/simulation/pharming_input/n{wildcards.ncells}_m{wildcards.n_mutation_groups}_size{wildcards.mutation_group_sizes}_cov{wildcards.cov}_p{wildcards.nclusters}_s{wildcards.seed}/copy_numbers.csv "
        " -s 0 -l 1000 --top_n 1 --dcfs data/simulation/pharming_input/dcfs.txt "
        " --root_x 2 --root_y 0 "
        " --ninit-segs 10 "
        " --thresh-prop 0.01 "
        " --sum-condition --collapse --cell-threshold 10 "
        " -P data/simulation/pharming_output/n{wildcards.ncells}_m{wildcards.n_mutation_groups}_size{wildcards.mutation_group_sizes}_cov{wildcards.cov}_p{wildcards.nclusters}_s{wildcards.seed}/solutions.pkl "
        " -O data/simulation/pharming_output/n{wildcards.ncells}_m{wildcards.n_mutation_groups}_size{wildcards.mutation_group_sizes}_cov{wildcards.cov}_p{wildcards.nclusters}_s{wildcards.seed} "
        " --tree data/simulation/pharming_output/n{wildcards.ncells}_m{wildcards.n_mutation_groups}_size{wildcards.mutation_group_sizes}_cov{wildcards.cov}_p{wildcards.nclusters}_s{wildcards.seed}/tree.png "
        " --labels data/simulation/pharming_output/n{wildcards.ncells}_m{wildcards.n_mutation_groups}_size{wildcards.mutation_group_sizes}_cov{wildcards.cov}_p{wildcards.nclusters}_s{wildcards.seed}/labels.csv "
        " > {log.std} 2> {log.err}"



rule pharming_on_simulation_with_loss:
    input:
        read_counts="data/simulation/pharming_input_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/read_counts.tsv",
        copy_numbers="data/simulation/pharming_input_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/copy_numbers.csv",
        segments="data/simulation/pharming_input_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/segments.txt",

    output:
        pharming_solution_pkl="data/simulation/pharming_output_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/solutions.pkl",
        pharming_tree_png="data/simulation/pharming_output_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/tree.png",
        pharming_labels="data/simulation/pharming_output_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/labels.csv",

    params:
        outdir="data/simulation/pharming_output_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}",

    conda:
        "pharming"

    log:
        std="data/simulation/pharming_output_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/pharming.log",
        err="data/simulation/pharming_output_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/pharming.err.log",

    benchmark:
        "data/simulation/pharming_output_with_loss_{lambda_pct}/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/pharming.benchmark",

    shell:
        (
            "pharming "
            "-f {input.read_counts} "
            "-c {input.copy_numbers} "
            "-s 0 -l 1000 --top_n 1 "
            "--dcfs data/simulation/pharming_input/dcfs.txt "
            "--root_x 2 --root_y 0 "
            "--ninit-segs 10 "
            "--thresh-prop 0.01 "
            "--sum-condition --collapse --cell-threshold 10 "
            "-P {output.pharming_solution_pkl} "
            "-O {params.outdir} "
            "--tree {output.pharming_tree_png} "
            "--labels {output.pharming_labels} "
            "> {log.std} 2> {log.err}"
        )


rule sbmclone_on_simulation:
    input:
        matrix="data/simulation/sbmclone_input/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/matrix.csv",
    output:
        blockmatrix="data/simulation/sbmclone_output/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/blockmatrix.csv",
    conda:
        "sbmclone"
    log:
        std="data/simulation/sbmclone_output/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/log",
        err="data/simulation/sbmclone_output/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/err.log",
    benchmark: "data/simulation/sbmclone_output/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/benchmark",
    shell:
        "python ../SBMClone/sbmclone.py data/simulation/sbmclone_input/n{wildcards.ncells}_m{wildcards.n_mutation_groups}_size{wildcards.mutation_group_sizes}_cov{wildcards.cov}_p{wildcards.nclusters}_s{wildcards.seed}/matrix.csv "
        "-o data/simulation/sbmclone_output/n{wildcards.ncells}_m{wildcards.n_mutation_groups}_size{wildcards.mutation_group_sizes}_cov{wildcards.cov}_p{wildcards.nclusters}_s{wildcards.seed}/"
        " > {log.std} 2> {log.err}"

rule split_inputs_williams:
    output:
        williams_separate="data/williams/scratch/separate.txt",
    log:
        std="data/williams/scratch/separate.log",
        err="data/williams/scratch/separate.err.log",
    benchmark: "data/williams/scratch/separate.benchmark",
    shell:
        "python -m src.williams.split_input_data "
        " > {log.std} 2> {log.err}"


rule make_inputs_williams:
    output:
        williams_F_plus="outputs/scope/williams/{sample_ids}/mutation_clusters/k_10/F_plus.csv",
        williams_F_minus="outputs/scope/williams/{sample_ids}/mutation_clusters/k_10/F_minus.csv",
        williams_F_loh_counts="outputs/scope/williams/{sample_ids}/loh_conflicts.csv"
    log:
        std="outputs/scope/williams/{sample_ids}/log",
        err="outputs/scope/williams/{sample_ids}/err.log",
    benchmark: "outputs/scope/williams/{sample_ids}/benchmark",
    shell:
        "python -m src.williams.make_input_data --sample {wildcards.sample_ids} "
        " > {log.std} 2> {log.err}"


rule scope_williams:
    input:
        williams_F_plus="outputs/scope/williams/{sample_ids}/mutation_clusters/k_{Ks}/F_plus.csv",
        williams_F_minus="outputs/scope/williams/{sample_ids}/mutation_clusters/k_{Ks}/F_minus.csv",
        williams_F_loh_counts="outputs/scope/williams/{sample_ids}/loh_conflicts.csv"
    output:
        williams_scope="outputs/scope/williams/{sample_ids}/solutions_mut/k_{Ks}/summary.txt",
    log:
        std="outputs/scope/williams/{sample_ids}/solutions_mut/k_{Ks}/log",
        err="outputs/scope/williams/{sample_ids}/solutions_mut/k_{Ks}/err.log",
    benchmark: "outputs/scope/williams/{sample_ids}/solutions_mut/k_{Ks}/benchmark",
    shell:
        "python -m src.williams.run_scope_on_williams --sample {wildcards.sample_ids} -k {wildcards.Ks} "
        " > {log.std} 2> {log.err}"
