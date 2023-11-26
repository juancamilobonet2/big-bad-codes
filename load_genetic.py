import numpy as np
import genetic_prange as gp
import code_utils as cu
import time

# Generates a matrix of the form (I, I, I, ...) as much times as 
# times_identity. row_num determines the size of the identity.
def generate_superI_matrix(row_num, times_identity):
    return np.concatenate(
        [np.identity(row_num, dtype=int) for _ in range(times_identity)], 
        axis=1
    )

def load(number_of_Is, row_num):
    t = (number_of_Is-1)//2
    HyG = generate_superI_matrix(row_num, number_of_Is)
    information = np.zeros((1, row_num), dtype=int)
    information[0][0]=1
    error = np.zeros((1, row_num*number_of_Is), dtype=int)
    error[0] = [1 if i<t else 0 for i in range(row_num*number_of_Is)]
    received = np.bitwise_xor(cu.multiply_matrices(information, HyG), error)
    print(received)
    s = cu.find_syndrome(HyG, received)
    gstart_time = time.time()
    e = gp.genetic_prange(20, 4, 0, s, HyG, t)
    computed_s = cu.find_syndrome(HyG, e)
    gelapsed_time = time.time() - gstart_time
    print(f"{np.array_equal(s, computed_s)} \n Received {s.transpose()[0]} \n Computed {computed_s.transpose()[0].astype(int)} \n Time {gelapsed_time} ")
        
for i in range(101, 261, 20):
    print(f"-------- Iteration {i} --------")
    load(i, 100)