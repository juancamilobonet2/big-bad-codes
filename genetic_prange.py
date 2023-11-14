import random
import copy

def rand_base_individual(column_num):
    return [[random.randint(0, column_num-1), random.randint(0, column_num-1)]]

def crossover(ind1, ind2):
    new_ind1, new_ind2 = [], []
    # Decide which genes will be passed to the first individual form
    # original individual 1 and original individual 2.
    genes_from_1 = [random.choice([True, False]) for position in ind1]
    genes_from_2 = [random.choice([True, False]) for position in ind2]
    # Append to new_ind 1 the original gen of ind1 if gen from 1 was present. 
    # Append to new_ind2 if not.
    for index, gen_presence in enumerate(genes_from_1):
        if gen_presence:
            new_ind1.append(ind1[index])
        else:
            new_ind2.append(ind1[index])
    # The same but with individual 2.
    for index, gen_presence in enumerate(genes_from_2):
        if gen_presence:
            new_ind1.append(ind2[index])
        else:
            new_ind2.append(ind2[index])
    return (new_ind1, new_ind2)

def mutation(individual, column_num):
    random_position = random.randint(0, len(individual)-1)
    origin_or_destiny = random.randint(0,1)
    random_new_column = random.randint(0, column_num-1)
    new_ind = copy.deepcopy(individual)
    new_ind[random_position][origin_or_destiny] = random_new_column
    return new_ind

# Basically makes all values a curve that is 1 when individual_weight = t.
# Greater than 1 when w < t and less than one when w > t. 
def fitness(individual_weight, column_num, t):
    return round((column_num-(individual_weight-t))/column_num, 4)

import numpy as np
import math
import code_utils as cu
import time

# Returns a pair, the first is the fitness and the second the error iff
# fitness greater than 0.
def modified_prange(s, H, t, P):
    """
    s: syndrome
    H: parity check matrix
    t: error correction capability
    """
    m,n = H.shape
    H_hat = cu.multiply_matrices(H, P)
        #U son las transformaciones
    U, H_hat = cu.gaussian_elimination(H_hat, start_column=n-m)
    if not np.array_equal(H_hat[:,n-m:], np.identity(m, dtype=int)):
        # permutacion es valida si h_hat es full rank
        return (0, None)
    s_bar = cu.apply_transforms(U, s)
    e_hat = np.zeros(n-m, dtype=int)
    e_hat = np.hstack((e_hat, s_bar.transpose()[0]))
    current_weight = np.sum(e_hat)
    return (fitness(current_weight, n, t), 
            cu.multiply_matrices(np.array([e_hat]), P))

def genetic_prange(number_of_inds, mutation_rate, s, H, t):
    m,n = H.shape
    inds = [rand_base_individual(n) for _ in range(number_of_inds)]
    results = [modified_prange(s, H, t, cu.permutation_matrix(n, ind)) for ind in inds]
    best_fit = max(results, key=lambda item: item[0])
    while best_fit[0] < 1:
        print("Best fit", best_fit)
        inds = next_gen(inds, results)
        results = [modified_prange(s, H, t, cu.permutation_matrix(n, ind)) for ind in inds]
        best_fit = max(results, key=lambda item: item[0])
    return best_fit[1]


def next_gen(inds, results):
    next_gen = []
    survivors_list = survivors(results)
    for index in survivors_list:
        next_gen.append(inds[index])
    indexed_results = [[index, *ind] for index, ind in enumerate(results)]
    indexed_results.sort(reverse=True, key=lambda item: item[1])
    for index in range(0, len(inds)-len(survivors_list)-1, 2):
        ind1_index = indexed_results[index][0]
        ind2_index = indexed_results[index+1][0]
        childs = crossover(inds[ind1_index], inds[ind2_index])
        next_gen.append(childs[0])
        next_gen.append(childs[1])
    return next_gen

def survivors(results):
    fitnesses = [perm[0] for perm in results]
    ponderated_fitnesses = fitnesses/np.sqrt(np.sum(fitnesses)**2) 
    survivors_list = []
    for index, individual in enumerate(ponderated_fitnesses):
        # Puede ser cambiado entre > o < dependiendo de si queremos que pasen más
        # de la gen anterior o menos.
        if individual > random.random():
            survivors_list.append(index)
    return survivors_list


def concat_matrix(matrix):
    # Convierte la matriz NumPy a una cadena
    concatenated_string = ''.join(map(str, matrix.flatten()))
    return concatenated_string

if __name__ == "__main__":
    print("Prange's algorithm")
    h_file = open("./test/goppa_h.txt", "r")
    H = []
    for line in h_file:
        H.append([int(x) for x in list(line.strip())])
    h_file.close()
    H = np.array(H)

    test_file = open("./test/goppa_test.txt", "r")
    test = []
    for line in test_file:
        test.append([int(x) for x in list(line.strip())])   
    test_file.close()
    
    n = H.shape[1]
    t=2
    print("Test ID \t Outcome \t Expected \t\t Result \t iterations \t time")
    test_id = 0 
    total_iters = 0
    total_time = 0
    for r in test:
        iterations = 1
        recieved = np.array([r])

        start_time = time.time()
        while iterations < 1000:
            s = cu.find_syndrome(H, recieved)
            e = genetic_prange(10, 0, s,H,t)
            computed_s = cu.find_syndrome(H, e)
            if np.array_equal(s, computed_s):
                break
            iterations += 1
        end_time = time.time()
        elapsed_time = end_time - start_time

        print(f"{test_id}\t \t {np.array_equal(s, computed_s)} \t {s.transpose()[0]} \t {computed_s.transpose()[0]} \t {iterations} \t \t {elapsed_time}s")
        test_id += 1
        total_iters += iterations
        total_time += elapsed_time
    print(f"Average iterations: {total_iters/test_id}")
    print(f"Average time: {total_time/test_id}s")