import pandas as pd
import numpy as np
import gurobipy as gp
from gurobipy import Model, GRB, quicksum

params = {
    "WLSACCESSID": '85dfeec1-65a1-402f-9425-7465d0f3a229',
    "WLSSECRET": 'f755a3a4-dba0-444d-8bc8-d8e781d8a05c',
    "LICENSEID": 2687964,
}
env = gp.Env(params=params)

def solve_cncff(
    F_plus: pd.DataFrame,
    F_minus: pd.DataFrame,
    n_clones: int,
    cluster_weights: list = None,
    n_solutions: int = 10,
    time_limit: int = 600,
):
    # TODO: Check if F_plus and F_minus have the same clones and mutations
    clusters = F_plus.index.tolist()
    mutations = F_plus.columns.tolist()

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
        model.addConstr(quicksum(u[l, j] for j in range(n_mutations)) == 1)

    # (4)
    for j in range(n_mutations):
        for l in range(n_clones + n_clusters):
            model.addConstr(b[l, j] <= x[j], name=f"b[{l}][{j}] <= x[{j}]")
    
    for j in range(n_mutations):
        for i in range(n_clusters):
            model.addConstr(f[i,j] <= x[j], name=f"f[{l}][{j}] = x[{j}]")

    # (5)
    for j in range(n_mutations):
        model.addConstr(quicksum(g[i, j] for i in range(n_clusters)) <= n_clusters - (n_clusters - 1) * x[j], name=f"g[{l}][{j}] <= p - (p - 1) * x[{j}]")


    # (6)
    for i in range(n_clusters):
        for j in range(n_mutations):
            if F_minus.iloc[i, j] == 0:
                model.addConstr(f[i, j] <= g[i, j], name=f"f[{i}][{j}] <= g[{i}][{j}]")
            elif F_plus.iloc[i, j] == 1:
                model.addConstr(f[i, j] >= 1 - g[i, j], name=f"f[{i}][{j}] >= 1 - g[{i}][{j}]")
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

                    if c < d:
                        # b_{jc} + b_{jd} <= 2 - z_{cd}
                        model.addConstr(b[j, c] + b[j, d] <= 2 - z[c,d])
    
    model.update()
    model.setObjective(quicksum((x[j] * cluster_weights[j]) for j in range(n_mutations)), GRB.MAXIMIZE)

    model.Params.TimeLimit = time_limit
    model.Params.PoolSearchMode = 2
    model.Params.PoolSolutions = n_solutions

    model.optimize()

    n_solutions_found = model.SolCount

    print(f"Number of solutions found: {n_solutions_found}")

    solutions = []
    for i in range(n_solutions_found):
        model.setParam(GRB.Param.SolutionNumber, i)
        print(f"Solution {i+1}:")
        print("Objective value:", model.PoolObjVal)

        X_df = pd.DataFrame([x[j].Xn for j in range(n_mutations)])

        b_values = [
            [b[i, j].Xn for j in range(n_mutations)]
            for i in range(n_clones + n_clusters)
        ]
        B_df = pd.DataFrame(b_values)

        u_values = [
            [u[i, j].Xn for j in range(n_mutations)]
            for i in range(n_clusters)
        ]
        U_df = pd.DataFrame(u_values)

        f_values = [
            [f[i, j].Xn for j in range(n_mutations)]
            for i in range(n_clusters)
        ]
        F_df = pd.DataFrame(f_values)

        g_values = [
            [g[i, j].Xn for j in range(n_mutations)]
            for i in range(n_clusters)
        ]
        G_df = pd.DataFrame(g_values)

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
