import pandas as pd
import numpy as np
import gurobipy as gp
from gurobipy import Model, GRB, quicksum


def solve_cppme_with_fixed_clones(
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
                model.addConstr(a[l, k, j] <= u[l, k])
                model.addConstr(a[l, k, j] <= b[k, j])
                model.addConstr(a[l, k, j] >= u[l, k] + b[k, j] - 1)

    for l in range(n_clusters):
        for j in range(n_mutations):
            model.addConstr(f[l, j] == quicksum(a[l, k, j] for k in range(n_clones))) 

    # (2)
    for l in range(n_clusters):
        for j in range(n_mutations):
            model.addConstr(f[l, j] >= F_minus.iloc[l, j] * x[j])
            model.addConstr(f[l, j] <= F_plus.iloc[l, j] * x[j])
    
    # (3)
    for l in range(n_clusters):
        model.addConstr(quicksum(u[l, j] for j in range(n_clones)) == 1)

    # (4)
    for j in range(n_mutations):
        for l in range(n_clones + n_clusters):
            model.addConstr(b[l, j] <= x[j])
    
    for j in range(n_mutations):
        for i in range(n_clusters):
            model.addConstr(f[i,j] <= x[j])
    
    for j in range(n_mutations):
        model.addConstr(quicksum(b[i, j] for i in range(n_clones)) >= x[j])

    # (5)
    for j in range(n_mutations):
        model.addConstr(quicksum(g[i, j] for i in range(n_clusters)) <= n_clusters - (n_clusters - 1) * x[j])


    # (6)
    for i in range(n_clusters):
        for j in range(n_mutations):
            if F_plus.iloc[i, j] == 1:
                model.addConstr(f[i, j] >= 1 - g[i, j])
            elif F_minus.iloc[i, j] == 0:
                model.addConstr(f[i, j] <= g[i, j])
            else:
                model.addConstr(g[i, j] == 1)

    # (7)
    for i in range(n_clusters):
        l = n_clones + i
        for j in range(n_mutations):
            model.addConstr(b[l, j] <= f[i, j])
            model.addConstr(b[l, j] <= 1 - g[i, j])
            model.addConstr(b[l, j] >= f[i, j] - g[i, j])

    # (8)
    for j in range(n_clones + n_clusters):
        for c in range(n_mutations):
            for d in range(n_mutations):
                if c != d:
                    # b_{jc} <= b_{jd} + (1 - y_{cd})
                    model.addConstr(b[j, c] <= b[j, d] + (1 - y[c,d]))

                    # y_{cd} + y_{dc} + z_{cd} >= 1
                    model.addConstr(y[c,d] + y[d,c] + z[c,d] >= 1)

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
    
    # (10) At least one cluster should have mutation frequency of 0.05 to avoid all-zero solution
    for j in range(n_mutations):
        model.addConstr(quicksum(f[i, j] for i in range(n_clusters)) * 20 >= x[j])
    
    # (*) bit-encoded ordering
    for i in range(n_clones-1):
        lhs = sum((2**j) * b[i,j]     for j in range(n_mutations))
        rhs = sum((2**j) * b[i+1,j]   for j in range(n_mutations))
        model.addConstr(lhs <= rhs - 1)
        if   i == 0: model.addConstr(lhs == 0)
        elif i == 1: model.addConstr(lhs >= 1)
    
    # (#) Enforce cluster parent restrictions
    if cluster_parent_restriction_pairs is not None:
        for c, d in cluster_parent_restriction_pairs:
            for j in range(n_mutations):
                model.addConstr(g[c, j] + b[d + n_clones, j] <= 1)

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
            model.addConstr(expr >= 1)
        
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
                
                model.addConstr(expr >= 1)
    
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


def solve_cppme(F_plus, F_minus, mutation_cluster_labels=None, 
                cluster_weights=None, forbidden_pairs=None, cluster_not_at_root=False,
                only_mutation_tree_variant=False):
    if mutation_cluster_labels:
        cluster_weights = mutation_cluster_labels.groupby("mutation_group").size().tolist()
    elif cluster_weights is None:
        raise ValueError("both mutation cluster labels and cluster weights cannot be None")

    total_clones = F_plus.shape[1] + 1

    n_clones = total_clones
    for group in F_plus.columns:
        mask = (F_plus[group] < 1) & (F_minus[group] > 0)
        subclonal_clusters = F_plus.index[mask].tolist()
        if len(subclonal_clusters) > 1: n_clones -= 1
    
    while n_clones >= 1:
        print(f"TRYING {n_clones} CLONES")
        solution, model = solve_cppme_with_fixed_clones(F_plus, F_minus, n_clones=n_clones,
                                        cluster_weights=cluster_weights, time_limit=60 * 60 * 24, 
                                        cluster_parent_restriction_pairs=forbidden_pairs, 
                                        cluster_not_at_root=cluster_not_at_root, 
                                        only_mutation_tree_variant=only_mutation_tree_variant)
        if solution is None: n_clones -= 1
        else: break

    print(f"THERE IS A SOLUTION WITH {n_clones} CLONES")

    X, _, _, _, _ = solution
    X = X.to_numpy().T[0]

    print(f"THERE IS A SOLUTION WITH {n_clones} CLONES")

    found_Bs = []
    solutions = []
    best_val = 0
    for i in range(1000):    
        solution, model = solve_cppme_with_fixed_clones(F_plus, F_minus, n_clones=n_clones,
                                                cluster_weights=cluster_weights, time_limit=24* 60 * 10, 
                                                cluster_parent_restriction_pairs=forbidden_pairs, 
                                                cluster_not_at_root=cluster_not_at_root,
                                                found_Bs=found_Bs, Xs=X, 
                                                only_mutation_tree_variant=only_mutation_tree_variant)
        if solution is None:
            break
        
        print("SOLUTION", i)
        found_B = solution[1].astype(int).to_numpy()
        found_Bs.append(found_B)
        solutions.append(solution)
        best_val = max(best_val, model.ObjVal)
        
    return solutions, best_val
