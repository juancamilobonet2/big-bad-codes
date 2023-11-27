import numpy as np
import genetic_prange as gp
import code_utils as cu
import time
import prange as pr

# Generates a matrix of the form (I, I, I, ...) as much times as 
# times_identity. row_num determines the size of the identity.
def generate_superI_matrix(row_num, times_identity):
    return np.concatenate(
        [np.identity(row_num, dtype=int) for _ in range(times_identity)], 
        axis=1
    )

def load(number_of_Is, row_num):
    # Maximum number of errors that can be corrected.
    t = (number_of_Is-1)//2
    # The matrix that will be used for the problem.
    G = generate_superI_matrix(row_num, number_of_Is)
    # If G is written as (I, A) of size kxn, H is (A^T, I), 
    # where the second I has dim n-k.
    H = np.concatenate(
        [generate_superI_matrix(row_num, number_of_Is-1).transpose(), 
         np.identity((number_of_Is-1)*row_num)], 
        axis=1
    )
    # The information we will decode (its easier to generate the message that way).
    information = np.zeros((1, row_num), dtype=int)
    information[0][0]=1
    # The arbitrary error we will use in the received message.
    error = np.zeros((1, row_num*number_of_Is), dtype=int)
    error[0] = [1 if i<t//2 else 0 for i in range(row_num*number_of_Is)]
    # Received message and its syndrome.
    received = np.bitwise_xor(cu.multiply_matrices(information, G), error)
    s = cu.find_syndrome(H, received)
    # Experiments time.
    rstart_time = time.time()
    re = pr.prange(s, H, t)
    relapsed_time = time.time() - rstart_time
    print("Ya por lo menos")
    gstart_time = time.time()
    ge = gp.genetic_prange(20, 4, 0, s, H, t)
    gelapsed_time = time.time() - gstart_time
    # Syndrome of the original error induced.
    computed_gs = cu.find_syndrome(H, ge)
    computed_rs = cu.find_syndrome(H, re)
    print(f'''
Original syndrome: {s.transpose()[0]}    
*** Prange
Computed: {computed_rs.transpose()[0].astype(int)}
Equal: {np.array_equal(s, computed_rs)}
Time: {relapsed_time}s
*** Genetic Prange
Computed: {computed_gs.transpose()[0].astype(int)}
Equal: {np.array_equal(s, computed_gs)} 
Time: {gelapsed_time}
''')
        
for i in range(41, 42, 20):
    print(f"-------- Iteration {i} --------")
    load(i, 7) 