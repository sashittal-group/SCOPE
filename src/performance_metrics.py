import numpy as np

def normalized_mutation_matrix_error(A: np.array, A_hat: np.array):
    return np.mean(np.abs(A - A_hat))
