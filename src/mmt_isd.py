import math
import code_utils as cu
import time
from sage.all import *

def mmt_isd(H, s, w, l, p, num_iterations=1000):
    # Implementation of May-Meurer-Thomae ISD algorithm
    # Input: Parity check matrix H (n-k x n), Syndrome s (length n-k), Natural number w.
    #PARAMETERS
    # l must be less than k-l
    # p no se que es p
    p1 = p2 = p//2

    l1 = cu.choose(p, p1)

    for _ in range(num_iterations):
        rand_permutation = cu.random_permutation_matrix(H.ncols())


    pass