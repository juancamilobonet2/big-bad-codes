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
    n = H.shape[1]
    current_weight = math.inf
    while current_weight > t:
        rand_permutation = np.random.randint(2,size=(n,n))
        H_hat = cu.multiply_matrices(H, rand_permutation)
        print(H_hat)


    




if __name__ == "__main__":
    print("Prange's algorithm")
    n = 6
    k = 3
    t = 1
    H = np.array([[1,1,0,1,0,0],[0,1,1,0,1,0],[1,1,1,0,0,1]])
    r = np.array([[1,1,0,1,1,0]])
    s = cu.find_syndrome(H,r)
    print(s)

    rand_permutation = np.random.randint(2,size=(n,n))
    H_hat = cu.multiply_matrices(H, rand_permutation)
    print(H_hat)
