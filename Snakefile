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

rule all:
    input:
        # simulation
        expand('data/simulation/ground_truth/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_read_count.csv', \
            itertools.product, seed=seeds, ncells=config['ncells'], n_mutation_groups=config['n_mutation_groups'], mutation_group_sizes=config['mutation_group_sizes'], \
            nclusters=config['nclusters'], cov=config['coverage']),

rule simulate:
    output:
        readcount_matrix="data/simulation/ground_truth/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_read_count.csv",
        variantcount_matrix="data/simulation/ground_truth/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_variant_count.csv",
        character_matrix="data/simulation/ground_truth/n{ncells}_m{n_mutation_groups}_size{mutation_group_sizes}_cov{cov}_p{nclusters}_s{seed}/sim_character_matrix_without_noise.csv",
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
