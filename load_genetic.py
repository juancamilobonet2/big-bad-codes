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

def load(k,n, load=False):
    if load:
        G = cu.file_to_matrix(f'./test/goppa_g.txt')
        H = cu.file_to_matrix(f'./test/goppa_h.txt')
        m = "Not calculated"
        gelapsed_time = "Not calculated"
        t=2
    else:
        rstart_time = time.time()
        G, H = cu.gen_g_h(40,30)
        m = cu.calculate_min_distance(G)
        t= (m-1)//2
        gelapsed_time = time.time() - rstart_time
    print("Generation of G: \n", G, "\n and H: \n", H)
    print("Minimum distance: ", m)
    print("time: ", gelapsed_time)
    print("-----------------------------------")

    word = cu.gen_word(G)
    error = cu.gen_error(t, G.shape[1])

    # Received message and its syndrome.
    received = cu.add_matrices(word, error)
    s = cu.find_syndrome(H, received)
    # Experiments time.
    print(f"Original syndrome: {s.transpose()[0]}")
    print(f"Original word: {word.transpose()}")
    print(f"error: {error}")
    print(f"received: {received}")

    rstart_time = time.time()
    re = pr.prange(s, H, t)
    relapsed_time = time.time() - rstart_time
    computed_rs = cu.find_syndrome(H, re)
    print(f'''*** Prange
                Computed: {computed_rs.transpose()[0].astype(int)}
                Equal: {np.array_equal(s, computed_rs)}
                Time: {relapsed_time}s''')
    gstart_time = time.time()
    ge = gp.genetic_prange(20, 4, 0, s, H, t)
    gelapsed_time = time.time() - gstart_time
    computed_gs = cu.find_syndrome(H, ge)
    print(f'''*** Genetic Prange
                    Computed: {computed_gs.transpose()[0].astype(int)}
                    Equal: {np.array_equal(s, computed_gs)} 
                    Time: {gelapsed_time}
                    ''')

for i in range(15, 42, 40):
    print(f"-------- Iteration {i} --------")
    load(i, 6, load=True) 