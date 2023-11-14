import numpy as np
import math
import code_utils as cu

def genetic_algorithm_for_permutation(s, H, t, population_size=10, generations=100):
    m, n = H.shape
    current_weight = math.inf
    best_permutation = np.identity(n, dtype=int)
    
    for _ in range(generations):
        population = [cu.random_permutation_matrix(n) for _ in range(population_size)]

        for perm in population:
            # Esta es la parte que seria la funcion de fitness
            H_hat = cu.multiply_matrices(H, perm)
            U, H_hat = cu.gaussian_elimination(H_hat, start_column=n-m)
            if np.array_equal(H_hat[:, n-m:], np.identity(m, dtype=int)):
                s_bar = cu.apply_transforms(U, s)
                e_hat = np.zeros(n-m, dtype=int)
                e_hat = np.hstack((e_hat, s_bar.transpose()[0]))
                current_weight_candidate = np.sum(e_hat)

                if current_weight_candidate < current_weight:
                    current_weight = current_weight_candidate
                    best_permutation = perm

    return cu.multiply_matrices(np.array([e_hat]), best_permutation)

def concat_matrix(matrix):
    concatenated_string = ''.join(map(str, matrix.flatten()))
    return concatenated_string

if __name__ == "__main__":
    print("Prange's algorithm with Genetic Algorithm for Permutation")
    n = 12
    k = 4
    t = 2
    H = np.array([[0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1],
                  [0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0],
                  [0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1],
                  [1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0],
                  [1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1],
                  [0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0],
                  [0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1],
                  [0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1]])

    r = np.array([[1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])  # 110000000000
    s = cu.find_syndrome(H, r)

    e = genetic_algorithm_for_permutation(s, H, t)
    computed_s = cu.find_syndrome(H, e)

    # In theory, H * e^t = s
    print("Error vector: \n", concat_matrix(e))
    print("Syndrome: \n", concat_matrix(computed_s))
    print("Actual syndrome: \n ", concat_matrix(s))
    print("Error weight: ", np.sum(e))

    successes = 0
    for i in range(1000):
        e = genetic_algorithm_for_permutation(s, H, t)
        computed_s = cu.find_syndrome(H, e)
        if np.array_equal(computed_s, s):
            successes += 1
        if i % 30 == 0 or np.array_equal(computed_s, s):
            print("SUCCESS:", np.array_equal(computed_s, s))
            print("Error vector: \n", concat_matrix(e))
            print("Syndrome: \n", concat_matrix(computed_s))
            print("Actual syndrome: \n ", concat_matrix(s))
            print("Error weight: ", np.sum(e))
            print("---------------------")

    print("Successes: ", successes)
    print("Success rate: ", successes / 1000)
