from sage.all import *
import time
import sys
sys.path.append("./src/")

import genetic_prange as gp
import code_utils as cu
import prange as pr







def run_goppa_mc():
    n, k, w, H_transpose, s_transpose = cu.read_dc_file('./data/challenge_goppa_mc_48.txt')
    H = H_transpose.transpose()
    H = H.change_ring(GF(2))
    H = block_matrix([[H, identity_matrix(GF(2), n - k)]], subdivide=False)
    s = s_transpose

    # Experiments time.
    print("starting experiments")
    prange_start_time = time.time()
    prange_result = run_prange(s,H,w)
    prange_elapsed_time = time.time() - prange_start_time
    print("prange done")

    genetic_prange_start_time = time.time()
    genetic_prange_result = run_genetic_prange(s,H,w)
    genetic_prange_elapsed_time = time.time() - genetic_prange_start_time

    # print_results(H, s, prange_result, prange_elapsed_time, w)
    # print_results(H, s, genetic_prange_result, genetic_prange_elapsed_time, w)
    return H, s, prange_result, prange_elapsed_time, genetic_prange_result, genetic_prange_elapsed_time



def run_reed_solomon(n, k, q):
    #REED SOLOMON
    F = GF(q)
    C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
    H = C.parity_check_matrix()
    codeword = C.random_element()
    t= (C.minimum_distance()-1)//2

    Chan = channels.StaticErrorRateChannel(C.ambient_space(), t)
    received = Chan(codeword)
    s = cu.find_syndrome(H, received)

    # Experiments time.
    print("starting experiments")
    prange_start_time = time.time()
    prange_result = run_prange(s,H,t)
    prange_elapsed_time = time.time() - prange_start_time
    print("prange done")

    genetic_prange_start_time = time.time()
    genetic_prange_result = run_genetic_prange(s,H,t)
    genetic_prange_elapsed_time = time.time() - genetic_prange_start_time
    
    # print_results_with_original(H, s, codeword, received, prange_result, prange_elapsed_time, t)
    # print_results_with_original(H, s, codeword, received, genetic_prange_result, genetic_prange_elapsed_time, t)

    return H, s, prange_result, prange_elapsed_time, genetic_prange_result, genetic_prange_elapsed_time

def run_golay(q, extended):
    C = codes.GolayCode(GF(Integer(q)), extended)
    H = C.parity_check_matrix()
    codeword = C.random_element()
    t= (C.minimum_distance()-1)//2
    n,k = H.dimensions()
    print(f"dimensions: {n}x{k}")

    Chan = channels.StaticErrorRateChannel(C.ambient_space(), t)
    received = Chan(codeword)
    s = cu.find_syndrome(H, received)

    # Experiments time.\
    print("starting experiments")
    prange_start_time = time.time()
    prange_result = run_prange(s,H,t)
    prange_elapsed_time = time.time() - prange_start_time
    print("prange done")

    genetic_prange_start_time = time.time()
    genetic_prange_result = run_genetic_prange(s,H,t)
    genetic_prange_elapsed_time = time.time() - genetic_prange_start_time

    # print_results_with_original(H, s, codeword, received, prange_result, prange_elapsed_time, t)

    # print_results_with_original(H, s, codeword, received, genetic_prange_result, genetic_prange_elapsed_time, t)
    return H, s, prange_result, prange_elapsed_time, genetic_prange_result, genetic_prange_elapsed_time

def run_reed_muller(order, variables, q):
    C = codes.ReedMullerCode(GF(q), order, variables)
    H = C.parity_check_matrix()
    codeword = C.random_element()
    t= (C.minimum_distance()-1)//2

    n,k = H.dimensions()
    print(f"dimensions: {n}x{k}")

    Chan = channels.StaticErrorRateChannel(C.ambient_space(), t)
    received = Chan(codeword)
    s = cu.find_syndrome(H, received)

    # Experiments time.
    print("starting experiments")
    prange_start_time = time.time()
    prange_result = run_prange(s,H,t)
    prange_elapsed_time = time.time() - prange_start_time
    print("prange done")

    genetic_prange_start_time = time.time()
    genetic_prange_result = run_genetic_prange(s,H,t)
    genetic_prange_elapsed_time = time.time() - genetic_prange_start_time
    
    # print_results_with_original(H, s, codeword, received, prange_result, prange_elapsed_time, t)
    # print_results_with_original(H, s, codeword, received, genetic_prange_result, genetic_prange_elapsed_time, t)

    return H, s, prange_result, prange_elapsed_time, genetic_prange_result, genetic_prange_elapsed_time

def run_BCH(length, designed_distance, q):
    C = codes.BCHCode(GF(q), length, designed_distance)
    H = C.parity_check_matrix()
    codeword = C.random_element()
    t= (C.minimum_distance()-1)//2

    n,k = H.dimensions()
    print(f"dimensions: {n}x{k}")

    Chan = channels.StaticErrorRateChannel(C.ambient_space(), t)
    received = Chan(codeword)
    s = cu.find_syndrome(H, received)

    # Experiments time.
    print("starting experiments")
    prange_start_time = time.time()
    prange_result = run_prange(s,H,t)
    prange_elapsed_time = time.time() - prange_start_time
    print("prange done")

    genetic_prange_start_time = time.time()
    genetic_prange_result = run_genetic_prange(s,H,t)
    genetic_prange_elapsed_time = time.time() - genetic_prange_start_time
    
    # print_results_with_original(H, s, codeword, received, prange_result, prange_elapsed_time, t)
    # print_results_with_original(H, s, codeword, received, genetic_prange_result, genetic_prange_elapsed_time, t)
    return H, s, prange_result, prange_elapsed_time, genetic_prange_result, genetic_prange_elapsed_time

def run_prange(s,H, num_errors):
    computed_error = pr.prange(s, H, num_errors)
    return computed_error

def run_genetic_prange(s,H, num_errors):
    computed_error = gp.genetic_prange(1_000, 10, 0.9, s, H, num_errors)
    return computed_error

def print_results(H, s, computed_error, elapsed_time, t):
    print("---------------------------------------------------------------------")
    print(f"Original syndrome: {s}")
    print(f"Computed error: {computed_error}")

    print(f"Correct syndrome?: {s == cu.find_syndrome(H, computed_error)}")
    print(f'Correct error weight?: {computed_error.hamming_weight() == t}')
    print(f"Time: {elapsed_time}s")

def print_results_with_original(H, s, original_word, received_word, computed_error, elapsed_time, t):
    print("---------------------------------------------------------------------")
    print(f"Original syndrome: {s}")
    print(f"Original word: {original_word}")
    print(f"Error: {received_word- original_word}")
    print(f"Computed error: {computed_error}")

    print(f"Correct syndrome?: {s == cu.find_syndrome(H, computed_error)}")
    print(f'Correct error weight?: {computed_error.hamming_weight() == t}')
    print(f"Correct error?: {received_word-original_word == computed_error}")
    print(f"Time: {elapsed_time}s")

def run_many(n):
    results = []
    for i in range(n):
        print(f"Running experiment {i}")
        results.append(run_goppa_mc())

    return results

def process_results(results):
    prange_times = []
    genetic_prange_times = []
    prange_correct = 0
    genetic_prange_correct = 0
    for result in results:
        prange_times.append(result[3])
        genetic_prange_times.append(result[5])
        if result[1] == cu.find_syndrome(result[0], result[2]):
            prange_correct += 1

        if result[1] == cu.find_syndrome(result[0], result[4]):
            genetic_prange_correct += 1

    print(f"Average prange time: {sum(prange_times)/len(prange_times)}")
    print(f"Average genetic prange time: {sum(genetic_prange_times)/len(genetic_prange_times)}")

    print(f"Prange correct: {prange_correct}/{len(results)}")
    print(f"Genetic prange correct: {genetic_prange_correct}/{len(results)}")

    return prange_times, genetic_prange_times


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

    n, k, w, H_transpose, s_transpose = cu.read_dc_file('./data/challenge_goppa_mc_48.txt')
    H = H_transpose.transpose()
    H = H.change_ring(GF(2))
    H = block_matrix([[H, identity_matrix(GF(2), n - k)]], subdivide=False)
    

    print(n)
    print(k)
    print(w)
    # print(H)
    print(s_transpose)
    print(H.dimensions())


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

if __name__ == '__main__':
    # run_reed_solomon(40, 15, 17)
    # run_goppa_mc()
    # run_golay(2, True)
    # run_reed_muller(2, 4, 2)
    # run_BCH(15, 7, 2)
    results = run_many(100)
    process_results(results)