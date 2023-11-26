import random
import numpy as np

# cosas generales de coding theory

# creo que en vd no hace falta
# class Binary_Matrix:
#     def __init__(self,n,k, matrix=[]):
#         self.n = n
#         self.k = k
#         if matrix == []:
#             self.matrix = np.zeros((k,n),dtype=int)
#         else:
#             self.matrix = matrix
    
#     def __str__(self):
#         return str(self.matrix)

#     def get_matrix(self):
#         return self.matrix
    
# operaciones modulo 2
def add_matrices(matrix1, matrix2):
    return (matrix1 + matrix2) % 2

def multiply_matrices(matrix1, matrix2):
    return (matrix1 @ matrix2) % 2

def find_syndrome(parity_check_matrix, vector):
    return (multiply_matrices(parity_check_matrix, vector.transpose()))

def random_permutation_matrix(n):
    """
    Returns a random permutation matrix of size n
    """
    permutation = [i for i in range(1,n+1)]
    permutation = np.random.permutation(permutation)
    permutation_matrix = np.zeros((n,n),dtype=int)
    for i in range(n):
        permutation_matrix[i,permutation[i]-1] = 1
    return permutation_matrix

def gaussian_elimination(matrix, start_column=0):
    """
    Gaussian elimination modulo 2
    """
    m, n = matrix.shape
    # U son las transformaciones que se le hacen a la matriz
    U = []
    for i in range(m):
        column = i + start_column
        # find row with 1 in the ith place
        for j in range(i,m):
            if matrix[j,column] == 1:
                matrix[[i,j]] = matrix[[j,i]]
                transform = np.identity(m, dtype=int)
                transform[[i,j]] = transform[[j,i]]
                U.append(transform)
                break
        for j in range(m):
            if matrix[j,column] == 1 and j != i:
                matrix[j] = add_matrices(matrix[j], matrix[i])
                transform = np.identity(m, dtype=int)
                transform[j] = add_matrices(transform[j], transform[i])
                U.append(transform)
    return (U,matrix)

def apply_transforms(U, matrix):
    """
    Applies the transformations in U to matrix
    """
    for i in range(len(U)):
        matrix = multiply_matrices(U[i], matrix)
    return matrix

# Number of columns.
def permutation_matrix(n, list_of_perms):
    identity = np.identity(n)
    for perm in list_of_perms:
        origin, end = perm
        for row_num in range(len(identity)):
            identity[row_num][origin], identity[row_num][end] = \
                identity[row_num][end], identity[row_num][origin]
    return identity

def file_to_matrix(file_name):
    file = open(file_name, "r")
    matrix = []
    for line in file:
        matrix.append([int(x) for x in list(line.strip())])
    file.close()
    return np.array(matrix)

def matrix_to_file(matrix, file_name):
    file = open(file_name, "w")
    for row in matrix:
        for element in row:
            file.write(str(element))
        file.write("\n")
    file.close()

def gen_g_h(n, k):
    """
    Generates a random generator matrix of size k x n
    """
    I = np.identity(k, dtype=int)
    A = np.random.randint(2, size=(k, n-k))

    G = np.concatenate((I, A), axis=1)

    I_h = np.identity(n-k, dtype=int)
    H = np.concatenate((A.transpose(), I_h), axis=1)

    return (G,H)

def calculate_min_distance(G):
    """
    Calculates the minimum distance of a code given its G
    """
    vectors = []
    for i in range(1,(2**G.shape[0])):
        vectors.append(np.array([[int(x) for x in list(np.binary_repr(i, width=G.shape[0]))]]))
    
    min_distance = G.shape[1]
    for vector in vectors:
        word = multiply_matrices(vector, G)
        weight = np.sum(word)
        if weight < min_distance and weight != 0:
            min_distance = weight

    return min_distance

def gen_error(t,n):
    """
    Generates a random error vector of size n with t ones
    """
    error = np.zeros((1,n),dtype=int)
    error[0,:t] = 1
    np.random.shuffle(error[0])
    return error
    

if __name__ == "__main__":
    G, H = gen_g_h(20,8)
    print(G)
    print(H)
    print(calculate_min_distance(G))