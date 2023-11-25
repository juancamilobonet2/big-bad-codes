import random
import copy
import numpy as np
import math
import code_utils as cu
import time

def rand_base_individual(column_num: int) -> list[tuple[int]]:
    return [[random.randint(0, column_num-1), random.randint(0, column_num-1)]]

# Receives two individuals and generates an offspring from them.
# Does not modify the original ones.
def crossover(ind1: list[tuple[int]], ind2: list[tuple[int]]) -> (list[tuple[int]], list[tuple[int]]):
    # Decide which genes will be passed to the first individual form
    # original individual 1 and original individual 2.
    genes_from_1 = [random.choice([True, False]) for position in ind1]
    genes_from_2 = [random.choice([True, False]) for position in ind2]
    # New_ind1 receives all True positions
    new_ind1 = [perm for index, perm in enumerate(ind1) if genes_from_1[index]]\
        + [perm for index, perm in enumerate(ind2) if genes_from_2[index]]
    new_ind2 = [perm for index, perm in enumerate(ind1) if not genes_from_1[index]]\
        + [perm for index, perm in enumerate(ind2) if not genes_from_2[index]]
    return (new_ind1, new_ind2)

# Mutates a single individual.
# It's a not an in place operation.
# column_num is the number of columns of the H matrix.
def mutation(individual: list[tuple[int]], column_num: int):
    if len(individual) <= 0:
        return []

    # Chooses the random permutation in the individual to mutate.
    random_position = random.randint(0, len(individual)-1)
    # Chooses between the origin column of the permutation and its destiny column.
    # If it is 0, the origin was chosen.
    origin_or_destiny = random.randint(0,1)
    # Generates the new column to be put, must be part of H.
    random_new_column = random.randint(0, column_num-1)
    new_ind = copy.deepcopy(individual)
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
    if np.array_equal(s, new_syndrome):
        fitness = len(e[0])
    
    for i in range(len(e[0])):
        if e[0][i] == 0:
            fitness += 1

    return fitness

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

    e = cu.multiply_matrices(np.array([e_hat]), P)
    return (fitness(s, e, H, t), e)

def genetic_prange(max_iters, number_of_inds, mutation_rate, s, H, t):
    m,n = H.shape
    inds = [rand_base_individual(n) for _ in range(number_of_inds)]
    results = [modified_prange(s, H, t, cu.permutation_matrix(n, ind)) for ind in inds]
    best_fit = max(results, key=lambda item: item[0])
    best_weight = np.sum(best_fit[1])
    contador = 0
    while best_weight > t and contador < max_iters:
        contador+=1
        inds = next_gen(inds, results, n, mutation_rate)
        results = [modified_prange(s, H, t, cu.permutation_matrix(n, ind)) for ind in inds]
        best_fit = max(results, key=lambda item: item[0])
        best_weight = np.sum(best_fit[1])
    return best_fit[1]

def next_gen(num_inds, results, column_num, mutation_rate=1):
    survivors = survivors(results)
    next_gen = copy.deepcopy(survivors)

    while len(next_gen) < num_inds:
        parent1 = random.choice(survivors)
        parent2 = random.choice(survivors)
        child1, child2 = crossover(parent1, parent2)
        child1 = mutation(child1, column_num, mutation_rate)
        child2 = mutation(child2, column_num, mutation_rate)
        next_gen.append(child1)
        next_gen.append(child2)

    next_gen = [next_gen[:num_inds]]
    return next_gen

def survivors(results):
    fitness_sum = sum([result[0] for result in results])
    survivors_list = []
    ordered = sorted(results, key=lambda item: item[0])
    survivor_number = random.randint(0, fitness_sum)
    cumulative_fitness = 0

    for ind in ordered:
        survivors_list.append(ind[1])
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
        recieved = np.array([r])

        start_time = time.time()
        s = cu.find_syndrome(H, recieved)
        e = genetic_prange(1000, 10, 0, s,H,t)
        computed_s = cu.find_syndrome(H, e)
        end_time = time.time()
        elapsed_time = end_time - start_time
        print("hola")
        #print(f"{test_id}\t \t {np.array_equal(s, computed_s)} \t {s.transpose()[0]} \t {computed_s.transpose()[0]} \t {100} \t \t {elapsed_time}s")
        test_id += 1
        total_time += elapsed_time
    print(f"Average iterations: {total_iters/test_id}")
    print(f"Average time: {total_time/test_id}s")