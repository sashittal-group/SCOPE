import pandas as pd
import numpy as np
import gurobipy as gp
from gurobipy import Model, GRB, quicksum


def solve_cncff(
    F_plus: pd.DataFrame,
    F_minus: pd.DataFrame,
    n_clones: int,
    cluster_weights: list = None,
    time_limit: int = 600,
    cluster_parent_restriction_pairs = None,
    cluster_not_at_root: bool= False,
    found_Bs: list[np.array]= None,
    Xs: np.array=None,
    only_mutation_tree_variant=False,
    only_cn_tree_variant=False,
):
    params = {
        "WLSACCESSID": '85dfeec1-65a1-402f-9425-7465d0f3a229',
        "WLSSECRET": 'f755a3a4-dba0-444d-8bc8-d8e781d8a05c',
        "LICENSEID": 2687964,
    }
    env = gp.Env(params=params)

    # TODO: Check if F_plus and F_minus have the same clones and mutations
    clusters = F_plus.index.tolist()
    mutations = F_plus.columns.tolist()
    clones = np.arange(n_clones).tolist()

    n_mutations = len(mutations)
    n_clusters = len(clusters)
    
    if cluster_weights is None:
        cluster_weights = [1] * n_clusters

    model = Model("CNCFF", env=env)
    # Define variables
    x = model.addVars(n_mutations, vtype=GRB.BINARY, name="x")
    x.start = 0
    b = model.addVars(n_clones + n_clusters, n_mutations, vtype=GRB.BINARY, name="b")
    u = model.addVars(n_clusters, n_clones, lb=0, ub=1, vtype=GRB.CONTINUOUS, name="u")
    g = model.addVars(n_clusters, n_mutations, vtype=GRB.BINARY, name="g")
    a = model.addVars(n_clusters, n_clones, n_mutations, lb=0, vtype=GRB.CONTINUOUS, name="a")
    f = model.addVars(n_clusters, n_mutations, vtype=GRB.CONTINUOUS, name="f")
    y = model.addVars(n_mutations, n_mutations, vtype=GRB.BINARY, name="y")
    z = model.addVars(n_mutations, n_mutations, vtype=GRB.BINARY, name="z")
    
    # Ignore all variables in the solution pool except b
    u.PoolIgnore = 1
    g.PoolIgnore = 1
    a.PoolIgnore = 1
    f.PoolIgnore = 1
    y.PoolIgnore = 1
    z.PoolIgnore = 1

    # Constraints
    # (1)
    for l in range(n_clusters):
        for k in range(n_clones):
            for j in range(n_mutations):
                model.addConstr(a[l, k, j] <= u[l, k], name=f"a[{l}][{k}][{j}] <= u[{l}][{k}]")
                model.addConstr(a[l, k, j] <= b[k, j], name=f"a[{l}][{k}][{j}] <= b[{k}][{j}]")
                model.addConstr(a[l, k, j] >= u[l, k] + b[k, j] - 1, name=f"a[{l}][{k}][{j}] >= u[{l}][{k}] + b[{k}][{j}]")

    for l in range(n_clusters):
        for j in range(n_mutations):
            model.addConstr(f[l, j] == quicksum(a[l, k, j] for k in range(n_clones)), name=f"f[{l}][{j}] = sum(a[{l}][k][{j}] for k in range(n_clones))") 

    # (2)
    for l in range(n_clusters):
        for j in range(n_mutations):
            model.addConstr(f[l, j] >= F_minus.iloc[l, j] * x[j], name=f"f[{l}][{j}] >= f_minus[{l}][{j}] * x[{j}]")
            model.addConstr(f[l, j] <= F_plus.iloc[l, j] * x[j], name=f"f[{l}][{j}] <= f_plus[{l}][{j}] * x[{j}]")
    
    # (3)
    for l in range(n_clusters):
        model.addConstr(quicksum(u[l, j] for j in range(n_clones)) == 1)

    # (4)
    for j in range(n_mutations):
        for l in range(n_clones + n_clusters):
            model.addConstr(b[l, j] <= x[j], name=f"b[{l}][{j}] <= x[{j}]")
    
    for j in range(n_mutations):
        for i in range(n_clusters):
            model.addConstr(f[i,j] <= x[j], name=f"f[{l}][{j}] = x[{j}]")
    
    for j in range(n_mutations):
        model.addConstr(quicksum(b[i, j] for i in range(n_clones)) >= x[j], name=f"sum_i b[{i}][{j}] >= x[{j}]")

    # (5)
    for j in range(n_mutations):
        model.addConstr(quicksum(g[i, j] for i in range(n_clusters)) <= n_clusters - (n_clusters - 1) * x[j], name=f"g[{l}][{j}] <= p - (p - 1) * x[{j}]")


    # (6)
    for i in range(n_clusters):
        for j in range(n_mutations):
            if F_plus.iloc[i, j] == 1:
                model.addConstr(f[i, j] >= 1 - g[i, j], name=f"f[{i}][{j}] >= 1 - g[{i}][{j}]")
            elif F_minus.iloc[i, j] == 0:
                model.addConstr(f[i, j] <= g[i, j], name=f"f[{i}][{j}] <= g[{i}][{j}]")
            else:
                model.addConstr(g[i, j] == 1, name=f"g[{i}][{j}] == 1")

    # (7)
    for i in range(n_clusters):
        l = n_clones + i
        for j in range(n_mutations):
            model.addConstr(b[l, j] <= f[i, j], name=f"b[{l}][{j}] <= f[{i}][{j}]")
            model.addConstr(b[l, j] <= 1 - g[i, j], name=f"b[{l}][{j}] <= 1 - g[{i}][{j}]")
            model.addConstr(b[l, j] >= f[i, j] - g[i, j], name=f"b[{l}][{j}] >= f[{i}][{j}] - g[{i}][{j}]")

    # (8)
    for j in range(n_clones + n_clusters):
        for c in range(n_mutations):
            for d in range(n_mutations):
                if c != d:
                    # b_{jc} <= b_{jd} + (1 - y_{cd})
                    model.addConstr(b[j, c] <= b[j, d] + (1 - y[c,d]))

                    # y_{cd} + y_{dc} + z_{cd} >= 1
                    model.addConstr(y[c,d] + y[d,c] + z[c,d] >= 1)

                    # if c < d: TODO: Should this if condition be here ?
                    # b_{jc} + b_{jd} <= 2 - z_{cd}
                    model.addConstr(b[j, c] + b[j, d] <= 2 - z[c,d])
    

    # (9) If clone gained in cluster then in Tree cluster -> clone
    for j1 in range(n_mutations):
        for j2 in range(n_mutations):
            for i in range(n_clusters):
                l = i + n_clones
                model.addConstr(y[j2, j1] >= b[l, j1] + g[i, j2] - 1)

    # (10)    
    for i in range(n_clusters):
        for j in range(n_mutations):
            model.addConstr(f[i, j] * 20 >= g[i, j] * x[j])
            model.addConstr(f[i, j] * 20 <= 20 - g[i, j] * x[j])
    
    # (10) ? At least one cluster should have mutation frequency of 0.05 to avoid all-zero solution
    for j in range(n_mutations):
        model.addConstr(quicksum(f[i, j] for i in range(n_clusters)) * 20 >= x[j], name=f"sum_i f[{i}][{j}] >= 0.05 * x[{j}]")
    

    # (*) ? Fixing clones to mutation for restricting redundant mutation trees
    # for i in range(n_clones):
    #     model.addConstr(b[i, i] >= x[i], name=f"b[{i}][{i}] >= x[{i}]")
    
    # for i in range(n_clones):
    #     for j in range(n_mutations):
    #         model.addConstr(b[i, j] <= x[i])
    #         model.addConstr(b[i, j] <= x[j])
    
    # (*) ? bit-encoded ordering
    for i in range(n_clones-1):
        lhs = sum((2**j) * b[i,j]     for j in range(n_mutations))
        rhs = sum((2**j) * b[i+1,j]   for j in range(n_mutations))
        model.addConstr(lhs <= rhs - 1 , name=f"order_{i}")
        if   i == 0: model.addConstr(lhs == 0)
        elif i == 1: model.addConstr(lhs >= 1)
    
    # (*) B rows cannot be same
    # for i in range(n_clones):
    #     for j in range(n_clones):
    #         if i != j:
    #             for k in range(n_mutations):
    #                 model.addConstr(d_vars[i,j,k] >= b[i,k] * x[k] - b[j,k] * x[k])
    #                 model.addConstr(d_vars[i,j,k] >= b[j,k] * x[k] - b[i,k] * x[k])
    #                 model.addConstr(d_vars[i,j,k] * x[k] <= b[i,k] + b[j,k])
    #                 model.addConstr(d_vars[i,j,k] * x[k] <= 2 - (b[i,k] + b[j,k]))
                
    #             model.addConstr(quicksum(d_vars[i, j, k] for k in range(n_mutations)) >= x[i] * x[j], name=f"diff_{i}_{j} >= 1")
    

    # (12) If in Tree cluster -> clone then clone gained in cluster 
    # share = model.addVars(n_clusters, n_clones, vtype=GRB.BINARY, name="share")
    # more  = model.addVars(n_clusters, n_clones, vtype=GRB.BINARY, name="more")
    # child = model.addVars(n_clusters, n_clones, vtype=GRB.BINARY, name="child")
    # M = n_mutations

    # for i in range(n_clusters):
    #     l = i + n_clones
    #     for k in range(n_clones):
    #         for j in range(n_mutations):
    #             model.addConstr(b[l, j] + b[k, j] - 1 >= share[i, k])
    
    # for i in range(n_clusters):
    #     l = i + n_clones
    #     for k in range(n_clones):
    #         diff = quicksum(b[l,j] - b[k,j] for j in range(n_mutations))
    #         model.addConstr(diff >= 0 - M*(1-more[i,k]))
    #         model.addConstr(diff <= M*more[i,k])
    
    # for i in range(n_clusters):
    #     for k in range(n_clones):
    #         model.addConstr(child[i,k] <= share[i,k])
    #         model.addConstr(child[i,k] <= more[i,k])
    #         model.addConstr(child[i,k] >= share[i,k] + more[i,k] - 1)
    #         model.addConstr(child[i,k] <= g[i, k])

    # (#) Enforce cluster parent restrictions
    if cluster_parent_restriction_pairs is not None:
        for c, d in cluster_parent_restriction_pairs:
            for j in range(n_mutations):
                model.addConstr(g[c, j] + b[d + n_clones, j] <= 1)
        # p01 = model.addVars(n_clusters, n_clusters, n_mutations, vtype=GRB.BINARY, name="p01")
        # p10 = model.addVars(n_clusters, n_clusters, n_mutations, vtype=GRB.BINARY, name="p10")
        # p11 = model.addVars(n_clusters, n_clusters, n_mutations, vtype=GRB.BINARY, name="p11")
        # p00 = model.addVars(n_clusters, n_clusters, n_mutations, vtype=GRB.BINARY, name="p00")

        # for c in range(n_clusters):
        #     for d in range(n_clusters):
        #         if c != d:
        #             for j in range(n_mutations):
        #                 cl = c + n_clones
        #                 dl = d + n_clones
        #                 model.addConstr(p01[c, d, j] >= b[cl, j] - b[dl, j])
        #                 model.addConstr(p10[c, d, j] >= b[dl, j] - b[cl, j])
        #                 model.addConstr(p11[c, d, j] >= b[cl, j] + b[dl, j] - 1)
        #                 model.addConstr(p00[c, d, j] >= 1 - b[cl, j] - b[dl, j])
                    
        #                 model.addConstr(p01[c, d, j] + p10[c, d, j] + p11[c, d, j] + p00[c, d, j] <= 1)    

        # for c, d in cluster_parent_restriction_pairs:
        #     model.addConstr(quicksum(p10[d, c, j] for j in range(n_mutations)) >= 1)



    # Cluster cannot be at root
    if cluster_not_at_root:
        for i in range(n_clusters):
            model.addConstr(quicksum(b[i + n_clones, j] for j in range(n_mutations)) >= 1)


    # Do not allow previously found solutions
    if found_Bs is not None:
        check_rows = n_clones + n_clusters
        if only_mutation_tree_variant: check_rows = n_clones
        for found_B in found_Bs:
            expr = gp.LinExpr()
            for i in range(check_rows):
                for j in range(n_mutations):
                    if found_B[i, j] > 0.5: 
                        expr += (1 - b[i, j])
                    else:
                        expr += b[i, j]
            
            model.addConstr(expr >= 1, "new solution")
        
        if only_cn_tree_variant:
            for found_B in found_Bs:
                expr = gp.LinExpr()
                for i in range(n_clusters):
                    l = i + n_clones
                    for j in range(n_mutations):
                        if found_B[l, j] > 0.5: 
                            expr += (1 - b[l, j])
                        else:
                            expr += b[l, j]
                
                model.addConstr(expr >= 1, "new CN tree")
    
    # Fix X if provided
    if Xs is not None:
        print(Xs)
        for j in range(n_mutations):
            if Xs[j] > 0.5: model.addConstr(x[j] == 1)
            else: model.addConstr(x[j] == 0)

    model.update()

    for v in model.getVars():
        if not (v.VarName.startswith('b') or v.VarName.startswith('x')):
            v.PoolIgnore = 1
    
    model.update()
    model.setObjective(quicksum((x[j] * cluster_weights[j]) for j in range(n_mutations)), GRB.MAXIMIZE)

    model.Params.TimeLimit = time_limit
    
    model.optimize()

    solution = None
    if model.Status == gp.GRB.OPTIMAL:
        X_df = pd.DataFrame([x[j].X for j in range(n_mutations)], index=mutations)

        b_values = [
            [b[i, j].X for j in range(n_mutations)]
            for i in range(n_clones + n_clusters)
        ]
        B_df = pd.DataFrame(b_values, index=clones + clusters, columns=mutations)

        u_values = [
            [u[i, j].X for j in range(n_clones)]
            for i in range(n_clusters)
        ]
        U_df = pd.DataFrame(u_values, index=clusters, columns=clones)

        f_values = [
            [f[i, j].X for j in range(n_mutations)]
            for i in range(n_clusters)
        ]
        F_df = pd.DataFrame(f_values, index=clusters, columns=mutations)

        g_values = [
            [g[i, j].X for j in range(n_mutations)]
            for i in range(n_clusters)
        ]
        G_df = pd.DataFrame(g_values, index=clusters, columns=mutations)

        solution = (X_df, B_df, U_df, F_df, G_df)
    
    return solution, model



if __name__ == "__main__":
    # Example usage
    F_plus = pd.DataFrame(
        np.array([[0.2, 1], [0.4, 1]]),
        index=["clone1", "clone2"],
        columns=["mut1", "mut2"],
    )
    F_minus = pd.DataFrame(
        np.array([[0.1, 0.2], [0.3, 0.4]]),
        index=["clone1", "clone2"],
        columns=["mut1", "mut2"],
    )
    X, B, U, F, G = solve_cncff(F_plus, F_minus, 2)
    print(X)
    print(B)
    print(U)
    print(F)
    print(G)
