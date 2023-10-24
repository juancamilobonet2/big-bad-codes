import numpy as np
import math
import random
import code_utils as cu

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



if __name__ == "__main__":
    print("Prange's algorithm")
    n = 6
    k = 3
    t = 1
    H = np.array([[1,1,0,1,0,0],[0,1,1,0,1,0],[1,1,1,0,0,1]])
    r = np.array([[1,1,0,1,1,0]])
    s = cu.find_syndrome(H,r)

    e = prange(s,H,t)
    computed_s = cu.find_syndrome(H, e)

    # en teoria, H*e^t = s
    print("Error vector: \n", e)
    print("Syndrome: \n", computed_s)
    print("actual syndrome: \n ", s)
    print("Error weight: ", np.sum(e))

    successes = 0
    for i in range(1000):
        e = prange(s,H,t)
        computed_s = cu.find_syndrome(H, e)
        if np.array_equal(computed_s, s):
            successes += 1
    
    print("Successes: ", successes)
    print("Success rate: ", successes/1000)


