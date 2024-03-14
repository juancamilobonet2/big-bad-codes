import math
import code_utils as cu
import time
from sage.all import *

def prange(s, H, t):
    """
    s: syndrome
    H: parity check matrix
    t: error correction capability
    """
    m,n = H.dimensions()
    current_weight = math.inf
    while current_weight > t:
        permutation_chosen = False
        while not permutation_chosen:
            rand_permutation = cu.random_permutation_matrix(n)
            H_hat =H*rand_permutation
            #U son las transformaciones
            U, H_hat = cu.gaussian_elimination(H_hat, start_column=n-m)
            
            if H_hat.rank() == m:
                # permutacion es valida si h_hat es full rank
                permutation_chosen = True
        s_bar = U*s
        e_hat = zero_vector(n-m)
        e_hat = vector(e_hat.list() + s_bar.list())

        current_weight = e_hat.hamming_weight()
    return cu.multiply_matrices(e_hat, rand_permutation.transpose())

def concat_matrix(matrix):
    # Convierte la matriz NumPy a una cadena
    concatenated_string = ''.join(map(str, matrix.flatten()))

    return concatenated_string


if __name__ == "__main__":
    print("Prange's algorithm")
    # H= cu.file_to_matrix("./data/simple_h.txt")
    n,k = 10,5
    F = GF(11)
    C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
    H = C.parity_check_matrix()

    # test_file = open("./data/goppa_test.txt", "r")
    # test = []
    # for line in test_file:
    #     test.append([int(x) for x in list(line.strip())])   
    # test_file.close()   
    codeword = C.random_element()
    t= (C.minimum_distance()-1)//2
    print(C.minimum_distance())
    print(t)

    Chan = channels.StaticErrorRateChannel(C.ambient_space(), t)
    test = [Chan(codeword)]
    print(codeword)
    print(test[0])
    n = H.dimensions()[1]
    print("Test ID \t Outcome \t Expected \t\t Result \t iterations \t time")
    test_id = 0 
    total_iters = 0
    total_time = 0
    for r in test:
        iterations = 1
        recieved = r

        start_time = time.time()
        while iterations < 1000:
            s = cu.find_syndrome(H, recieved)
            e = prange(s,H,t)
            computed_s = cu.find_syndrome(H, e)
            if s == computed_s:
                break
            iterations += 1
        end_time = time.time()
        elapsed_time = end_time - start_time

        print(f"{test_id}\t \t {s == computed_s} \t {s} \t {computed_s} \t {iterations} \t \t {elapsed_time}s")
        print(e)
        test_id += 1
        total_iters += iterations
        total_time += elapsed_time
    print(f"Average iterations: {total_iters/test_id}")
    print(f"Average time: {total_time/test_id}s")
