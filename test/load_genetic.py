from sage.all import *
import time
import sys
sys.path.append("./src/")

import genetic_prange as gp
import code_utils as cu
import prange as pr

def load(load=False):
    # if load:
    #     G = cu.file_to_matrix(f'./data/goppa_g.txt')
    #     H = cu.file_to_matrix(f'./data/goppa_h.txt')
    #     m = "Not calculated"
    #     gelapsed_time = "Not calculated"
    #     t=2
    # else:
    #     rstart_time = time.time()
    #     G, H = cu.gen_g_h(40,30)
    #     m = cu.calculate_min_distance(G)
    #     gelapsed_time = time.time() - rstart_time
    # print("Generation of G: \n", G, "\n and H: \n", H)
    # print("Minimum distance: ", m)
    # print("time: ", gelapsed_time)
    # print("-----------------------------------")
    #GOPPA
    # R = F['x']; (x,) = R._first_ngens(1) 
    # g = x**Integer(3) +x+ Integer(1)
    # L = [a for a in F.list() if g(a) != 0]
    # C = codes.GoppaCode(g, L)
    # # E = codes.encoders.GoppaCodeEncoder(C)
    # H = C.parity_check_matrix()
    # G = C.generator_matrix()
    # codeword = C.random_element()
    # t= (C.minimum_distance()-1)//2

    #REED SOLOMON
    # n,k = 10,5
    # F = GF(11)
    # C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
    # H = C.parity_check_matrix()
    # codeword = C.random_element()
    # t= (C.minimum_distance()-1)//2
    n, k, w, H_transpose, s_transpose = cu.read_dc_file('./data/challenge_goppa_mc_80.txt')
    H = H_transpose.transpose()
    H = H.change_ring(GF(2))
    H = block_matrix([[H, identity_matrix(GF(2), n - k)]], subdivide=False)
    

    print(n)
    print(k)
    print(w)
    # print(H)
    print(s_transpose)
    print(H.dimensions())

    # Chan = channels.StaticErrorRateChannel(C.ambient_space(), t)
    # test = [Chan(codeword)]

    # word = codeword

    # Received message and its syndrome.
    # received = 
    # s = cu.find_syndrome(H, received)
    s = s_transpose
    # Experiments time.
    print(f"Original syndrome: {s}")
    # print(f"Original word: {word}")
    # print(f"received: {received}")

    rstart_time = time.time()
    re = pr.prange(s, H, w)
    relapsed_time = time.time() - rstart_time
    computed_rs = cu.find_syndrome(H, re)
    print(f'''*** Prange
                Computed: {computed_rs}
                Equal: {(s == computed_rs) and (re.hamming_weight() == w)}
                Time: {relapsed_time}s''')
    gstart_time = time.time()
    ge = gp.genetic_prange(math.inf, 1000, 0.5, s, H, w)
    gelapsed_time = time.time() - gstart_time
    computed_gs = cu.find_syndrome(H, ge)
    print(f'''*** Genetic Prange
                    Computed: {computed_gs}
                    Equal: {(s == computed_gs) and (ge.hamming_weight() == w)}
                    Time: {gelapsed_time}
                    ''')

    # error_real = "000000000000000000010000000000000000000000000000"
    # error_real = vector(GF(2), [int(i) for i in error_real])
    # print(ge)
    # print(error_real)
    # print(H*error_real)
    # print(s)
    # print(f"right? {cu.find_syndrome(H, error_real) == s}")

if __name__ == '__main__':
    print(f"-------- Iteration 1 --------")
    load(load=True) 