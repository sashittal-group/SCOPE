import pandas as pd
import numpy as np
import gurobipy as gp
from gurobipy import Model, GRB, quicksum


def solve_cncff(
    F_plus: pd.DataFrame,
    F_minus: pd.DataFrame,
    n_clones: int,
    cluster_weights: list = None,
    n_solutions: int = 10,
    time_limit: int = 600,
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
    clones = list(np.arange(n_clones))

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
    d_vars = model.addVars(n_clones, n_clones, n_mutations, vtype=GRB.BINARY, name="d")
    c_vars = model.addVars(n_clusters, n_clones, vtype=GRB.BINARY, name="c")
    
    # Ignore all variables in the solution pool except b
    u.PoolIgnore = 1
    g.PoolIgnore = 1
    a.PoolIgnore = 1
    f.PoolIgnore = 1
    y.PoolIgnore = 1
    z.PoolIgnore = 1
    d_vars.PoolIgnore = 1

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
        model.addConstr(quicksum(u[l, j] for j in range(n_mutations)) == 1)

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
    
    # (9) ? Fixing clones to mutation for restricting redundant mutation trees
    for i in range(n_clones):
        model.addConstr(b[i, i] >= x[i], name=f"b[{i}][{i}] >= x[{i}]")
    
    # (9) ? bit-encoded ordering
    # for i in range(n_clones-1):
    #     lhs = sum((2**j) * b[i,j]     for j in range(n_mutations))
    #     rhs = sum((2**j) * b[i+1,j]   for j in range(n_mutations)) - 1
    #     model.addConstr(lhs <= rhs, name=f"order_{i}")

    # for j in range(n_mutations):
    #     model.addConstr(b[n_clones, j] >= b[n_clones + 1, j], name=f"order_0_{j}")

    # (9) B rows cannot be same
    for i in range(n_clones):
        for j in range(n_clones):
            if i != j:
                for k in range(n_mutations):
                    model.addConstr(d_vars[i,j,k] >= b[i,k] - b[j,k])
                    model.addConstr(d_vars[i,j,k] >= b[j,k] - b[i,k])
                    model.addConstr(d_vars[i,j,k] <= b[i,k] + b[j,k])
                    model.addConstr(d_vars[i,j,k] <= 2 - (b[i,k] + b[j,k]))
                
                model.addConstr(quicksum(d_vars[i, j, k] for k in range(n_mutations)) >= 1, name=f"diff_{i}_{j} >= 1")

    for j in range(n_mutations):
        model.addConstr(quicksum(b[i, j] for i in range(n_clones)) >= x[i])
    
    # (10) ? At least one cluster should have mutation frequency of 0.05 to avoid all-zero solution
    for j in range(n_mutations):
        model.addConstr(quicksum(f[i, j] for i in range(n_clusters)) * 20 >= x[j], name=f"sum_i f[{i}][{j}] >= 0.05 * x[{j}]")
    

    # (11) If clone gained in cluster then in Tree cluster -> clone
    for j1 in range(n_mutations):
        for j2 in range(n_mutations):
            for i in range(n_clusters):
                l = i + n_clones
                model.addConstr(y[j2, j1] >= b[l, j1] + g[i, j2] - 1)
    
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

    model.update()

    for v in model.getVars():
        if not (v.VarName.startswith('b') or v.VarName.startswith('x')):
            v.PoolIgnore = 1

    model.update()
    model.setObjective(quicksum((x[j] * cluster_weights[j]) for j in range(n_mutations)), GRB.MAXIMIZE)

    model.Params.TimeLimit = time_limit
    model.Params.PoolSearchMode = 2
    model.Params.PoolSolutions = n_solutions
    model.Params.PoolGap = 0

    model.optimize()

    # Form Dataframes with solution
    n_solutions_found = model.SolCount

    print(f"Number of solutions found: {n_solutions_found}")

    solutions = []
    for i in range(n_solutions_found):
        model.setParam(GRB.Param.SolutionNumber, i)
        print(f"Solution {i+1}:")
        print("Objective value:", model.PoolObjVal)

        X_df = pd.DataFrame([x[j].Xn for j in range(n_mutations)], index=mutations)

        b_values = [
            [b[i, j].Xn for j in range(n_mutations)]
            for i in range(n_clones + n_clusters)
        ]
        B_df = pd.DataFrame(b_values, index=clones + clusters, columns=mutations)

        u_values = [
            [u[i, j].Xn for j in range(n_mutations)]
            for i in range(n_clusters)
        ]
        U_df = pd.DataFrame(u_values, index=clusters, columns=clones)

        f_values = [
            [f[i, j].Xn for j in range(n_mutations)]
            for i in range(n_clusters)
        ]
        F_df = pd.DataFrame(f_values, index=clusters, columns=mutations)

        g_values = [
            [g[i, j].Xn for j in range(n_mutations)]
            for i in range(n_clusters)
        ]
        G_df = pd.DataFrame(g_values, index=clusters, columns=mutations)

        solutions.append((X_df, B_df, U_df, F_df, G_df))

    return solutions



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
