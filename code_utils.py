import numpy as np
from sage.all import *

# cosas generales de coding theory
    
def create_matrix(p, values):
    return matrix(GF(Integer(p)), values)

def create_vector(p, values):
    return vector(GF(Integer(p)), values)

# DEPRECATED: just add them normally
def add_matrices(matrix1, matrix2):
    return (matrix1 + matrix2)

# DEPRECATED: just multiply them normally
def multiply_matrices(matrix1, matrix2):
    return (matrix1 * matrix2)

def find_syndrome(parity_check_matrix, vector):
    return (parity_check_matrix * vector)

def permutation_matrix(n, two_perms):
    p = matrix.identity(n)
    for perm in two_perms:
        origin, end = perm
        p.swap_columns(origin, end)
    return p

def random_permutation_matrix(n):
    """
    Returns a random permutation matrix of size n
    """
    permutation = [i for i in range(1,n+1)]
    permutation = np.random.permutation(permutation)
    permutation_matrix = matrix.zero(n)
    for i in range(n):
        permutation_matrix[i,permutation[i]-1] = 1
    return matrix(permutation_matrix)

def gaussian_elimination(matrix_in, start_column=0):
    """
    Gaussian elimination modulo 2
    """
    
    m = matrix_in.nrows()
    mat = copy(matrix_in)
    mat.subdivide(col_lines=start_column)
    submat = mat.subdivision(0,1)
    determinant = submat.determinant()
    if determinant == 0:
        return (matrix.identity(m), matrix_in)
    else:
        U = submat.inverse()
        return (U, U*matrix_in)
    

    # U son las transformaciones que se le hacen a la matriz
    # U = matrix.identity(m)
    # for i in range(m):
    #     column = i + start_column
    #     # find row with 1 in the ith place
    #     for j in range(i,m):
    #         if matrix_in[j,column] == 1:
    #             matrix_in[[i,j]] = matrix_in[[j,i]]
    #             transform = matrix.identity(m)
    #             transform[[i,j]] = transform[[j,i]]
    #             # U.append(transform)
    #             U = transform * U
    #             break
    #     for j in range(m):
    #         if matrix_in[j,column] == 1 and j != i:
    #             matrix_in[j] = add_matrices(matrix_in[j], matrix_in[i])
    #             transform = matrix.identity(m)
    #             transform[j] = add_matrices(transform[j], transform[i])
    #             # U.append(transform)
    #             U = transform * U
    # p = permutation_matrix(matrix.ncolumns(), [(0+i, start_column+i) for i in range(matrix.ncolumns()-start_column)])
    # m = matrix * p
    # m, u = m.echelon_form(transformation=True)

    return (U,matrix_in)

def apply_transforms(U, matrix):
    """
    Applies the transformations in U to matrix
    """
    for i in range(len(U)):
        matrix = multiply_matrices(U[i], matrix)
    return matrix

# Number of columns.
def permutation_matrix(n, list_of_perms):
    identity = matrix.identity(n)
    for perm in list_of_perms:
        origin, end = perm
        for row_num in range(n):
            identity[row_num,origin], identity[row_num,end] = \
                identity[row_num,end], identity[row_num,origin]
    return identity


def file_to_matrix(file_name):
    file = open(file_name, "r")
    p = 0
    matrix = []
    lines = file.readlines()
    p = int(lines.pop(0))
    for line in lines:
        matrix.append([int(x) for x in list(line.strip())])
    file.close()
    return create_matrix(p, matrix)

def matrix_to_file(matrix, file_name):
    file = open(file_name, "w")
    file.write(str(matrix.base_ring().order()))
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

def gen_word(G):
    """
    Generates a random word of size k with G
    """
    vector = G.linear_combination_of_rows() 
    word = multiply_matrices(vector[0], G)
    return word

    

if __name__ == "__main__":
    G, H = gen_g_h(20,8)
    print(G)
    print(H)
    print(calculate_min_distance(G))