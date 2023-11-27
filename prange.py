import numpy as np
import math
import code_utils as cu
import time

def prange(s, H, t):
    """
    s: syndrome
    H: parity check matrix
    t: error correction capability
    """
    m,n = H.shape
    current_weight = math.inf
    while current_weight > t:
        permutation_chosen = False
        while not permutation_chosen:
            rand_permutation = cu.random_permutation_matrix(n)
            H_hat = cu.multiply_matrices(H, rand_permutation)
            #U son las transformaciones
            U, H_hat = cu.gaussian_elimination(H_hat, start_column=n-m)
            if np.array_equal(H_hat[:,n-m:], np.identity(m, dtype=int)):
                # permutacion es valida si h_hat es full rank
                permutation_chosen = True
        
        s_bar = cu.apply_transforms(U, s)
        e_hat = np.zeros(n-m, dtype=int)
        e_hat = np.hstack((e_hat, s_bar.transpose()[0]))

        current_weight = np.sum(e_hat)
    return cu.multiply_matrices(np.array([e_hat]), rand_permutation)

def concat_matrix(matrix):
    # Convierte la matriz NumPy a una cadena
    concatenated_string = ''.join(map(str, matrix.flatten()))

    return concatenated_string


if __name__ == "__main__":
    print("Prange's algorithm")
    h_file = open("./test/goppa_h.txt", "r")
    H = []
    for line in h_file:
        H.append([int(x) for x in list(line.strip())])
    h_file.close()
    H = np.array(H)

    test_file = open("./test/goppa_test.txt", "r")
    test = []
    for line in test_file:
        test.append([int(x) for x in list(line.strip())])   
    test_file.close()

    
    n = H.shape[1]
    t=2
    print("Test ID \t Outcome \t Expected \t\t Result \t iterations \t time")
    test_id = 0 
    total_iters = 0
    total_time = 0
    for r in test:
        iterations = 1
        recieved = np.array([r])

        start_time = time.time()
        while iterations < 1000:
            s = cu.find_syndrome(H, recieved)
            e = prange(s,H,t)
            computed_s = cu.find_syndrome(H, e)
            if np.array_equal(s, computed_s):
                break
            iterations += 1
        end_time = time.time()
        elapsed_time = end_time - start_time

        print(f"{test_id}\t \t {np.array_equal(s, computed_s)} \t {s.transpose()[0]} \t {computed_s.transpose()[0]} \t {iterations} \t \t {elapsed_time}s")
        test_id += 1
        total_iters += iterations
        total_time += elapsed_time
    print(f"Average iterations: {total_iters/test_id}")
    print(f"Average time: {total_time/test_id}s")
