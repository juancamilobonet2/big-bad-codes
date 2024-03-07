import random as rand
import copy as copy_lib
# import numpy as np
import math
import code_utils as cu
import time
from sage.all import *

def rand_base_individual(column_num: int) -> list[tuple[int]]:
    return [(rand.randint(0, column_num-1), rand.randint(0, column_num-1))]

# Receives two individuals and generates an offspring from them.
# Does not modify the original ones.
def crossover(ind1: list[tuple[int]], ind2: list[tuple[int]]) -> (list[tuple[int]], list[tuple[int]]):
    # Decide which genes will be passed to the first individual form
    # original individual 1 and original individual 2.
    genes_from_1 = [rand.choice([True, False]) for position in ind1]
    genes_from_2 = [rand.choice([True, False]) for position in ind2]
    # New_ind1 receives all True positions
    new_ind1 = [perm for index, perm in enumerate(ind1) if genes_from_1[index]]\
        + [perm for index, perm in enumerate(ind2) if genes_from_2[index]]
    new_ind2 = [perm for index, perm in enumerate(ind1) if not genes_from_1[index]]\
        + [perm for index, perm in enumerate(ind2) if not genes_from_2[index]]
    return (new_ind1, new_ind2)

# Mutates a single individual.
# It's a not an in place operation.
# column_num is the number of columns of the H matrix.
def mutation(individual: list[tuple[int]], column_num: int, mutation_rate: float):
    if len(individual) <= 0:
        return []

    # Chooses the random permutation in the individual to mutate.
    random_position = rand.randint(0, len(individual)-1)
    # Chooses between the origin column of the permutation and its destiny column.
    # If it is 0, the origin was chosen.
    origin_or_destiny = rand.randint(0,1)
    # Generates the new column to be put, must be part of H.
    random_new_column = rand.randint(0, column_num-1)
    new_ind = copy_lib.deepcopy(individual)
    old_tuple = new_ind[random_position]
    if origin_or_destiny == 0:
        new_ind[random_position] = (random_new_column, old_tuple[1])
    else:
        new_ind[random_position] = (old_tuple[0], random_new_column)
    return new_ind

# Basically makes all values a curve that is 1 when individual_weight = t.
# Greater than 1 when w < t and less than one when w > t. 
def fitness(s, e, parity_check_matrix, t):
    fitness = 0
    new_syndrome = cu.find_syndrome(parity_check_matrix, e)
    if s==new_syndrome:
        fitness = len(e)
    for i in range(len(e)):
        if e[i] == 0:
            fitness += 1
    return fitness

# Returns a pair, the first is the fitness and the second the error iff
# fitness greater than 0.
def modified_prange(s, H, t, ind):
    """
    s: syndrome
    H: parity check matrix
    t: error correction capability
    """
    m,n = H.dimensions()
    P = cu.permutation_matrix(n, ind)
    H_hat = cu.multiply_matrices(H, P)
        #U son las transformaciones
    U, H_hat = cu.gaussian_elimination(H_hat, start_column=n-m)
    # if H_hat.rank() == m:
    #     # permutacion es valida si h_hat es full rank
    #     return (0, zero_vector(n-m), ind)
    s_bar = U*s
    e_hat = zero_vector(n-m)
    e_hat = vector(e_hat.list() + s_bar.list())

    e = e_hat*P.transpose()
    return (fitness(s, e, H, t), e, ind)

def genetic_prange(max_iters, number_of_inds, mutation_rate, s, H, t):
    # print("-------------------------------------------------")
    # print("Genetic Prange")
    # print(f"{s}, {H}, {t}")
    m,n = H.dimensions()
    inds = [rand_base_individual(n) for _ in range(number_of_inds)]
    results = [modified_prange(s, H, t, ind) for ind in inds]
    best_fit = max(results, key=lambda item: item[0])
    best_weight = best_fit[1].hamming_weight()
    contador = 0
    while best_weight > t and contador < max_iters:
        contador+=1
        inds = next_gen(results, n, mutation_rate)
        results = [modified_prange(s, H, t, ind) for ind in inds]
        best_fit = max(results, key=lambda item: item[0])
        best_weight = best_fit[1].hamming_weight()
    return best_fit[1]

def next_gen(results, column_num, mutation_rate=1):
    survivors_list = survivors(results)
    next_gen = copy_lib.deepcopy(survivors_list)

    while len(next_gen) < len(results):
        parent1 = rand.choice(survivors_list)
        parent2 = rand.choice(survivors_list)
        child1, child2 = crossover(parent1, parent2)
        child1 = mutation(child1, column_num, mutation_rate)
        child2 = mutation(child2, column_num, mutation_rate)
        next_gen.append(child1)
        next_gen.append(child2)

    next_gen = next_gen[:len(results)]
    return next_gen

def survivors(results):
    fitness_sum = sum([result[0] for result in results])
    survivors_list = []
    ordered = sorted(results, key=lambda item: item[0])
    survivor_number = rand.randint(0, fitness_sum)
    cumulative_fitness = 0

    for ind in ordered:
        # print(len(survivors_list))
        # print(ind)
        survivors_list.append(ind[2])
        cumulative_fitness += ind[0]
        if cumulative_fitness >= survivor_number:
            break
    
    return survivors_list

def concat_matrix(matrix):
    # Convierte la matriz NumPy a una cadena
    concatenated_string = ''.join(map(str, matrix.flatten()))
    return concatenated_string

if __name__ == "__main__":
    print("Prange's algorithm")
    # H = cu.file_to_matrix("./data/goppa_h.txt")
    # G = cu.file_to_matrix("./test/goppa_g.txt")
    # H = []
    # for line in h_file:
    #     H.append([int(x) for x in list(line.strip())])
    # h_file.close()
    # H = np.array(H)

    # test_file = open("./data/goppa_test.txt", "r")
    # test = []
    # for line in test_file:
    #     test.append(vector([int(x) for x in list(line.strip())]))   
    # test_file.close()
    n,k = 10,5
    F = GF(11)
    C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
    H = C.parity_check_matrix()
    codeword = C.random_element()
    t= (C.minimum_distance()-1)//2
    
    n = H.dimensions()[1]
    # t=2
    Chan = channels.StaticErrorRateChannel(C.ambient_space(), t)
    test = [Chan(codeword)]
    print("Test ID \t Outcome \t Expected \t\t Result \t iterations \t time")
    test_id = 0 
    total_iters = 0
    total_time = 0
    for r in test:
        recieved = r

        start_time = time.time()
        s = cu.find_syndrome(H, recieved)
        e = genetic_prange(10, 10, 0, s,H,t)
        computed_s = cu.find_syndrome(H, e)
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"{test_id}\t \t {s==computed_s} \t {s} \t {computed_s} \t TODO \t \t {elapsed_time}s")
        test_id += 1
        total_time += elapsed_time
    print(f"Average iterations: {total_iters/test_id}")
    print(f"Average time: {total_time/test_id}s")