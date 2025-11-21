
def get_approx_F_error(solutions, F_bar):

    errors = []

    for solution in solutions:

        X, B, U, F, G = solution

        cn_clusters = G.index.to_list()
        mutations_selected = X[X[0] > 0.5].index.to_list()

        error = 0
        for cn_cluster in cn_clusters:
            for mutation in mutations_selected:
                if G.loc[cn_cluster, mutation] > 0.5:
                    e = 0
                elif B.loc[cn_cluster, mutation] > 0.5:
                    e = 1 - F_bar.loc[cn_cluster, mutation]
                    error += e
                else:
                    e = F_bar.loc[cn_cluster, mutation]
                    error += e

        errors.append(error)
    
    return errors


def rank_phylogenies(solutions, F_bar):
    errors = get_approx_F_error(solutions, F_bar)
    indexed = list(enumerate(errors))
    ranked = sorted(indexed, key=lambda x: x[1])
    return ranked
