import code_utils as cu
from sage.all import *

def test_matrix_creation():
    try:
        matrix = cu.create_matrix(2, [[1,1],[1,0]])
        return True
    except Exception as err:
        print(err)
        return False

def test_matrix_addition():
    try:
        matrix1 = cu.create_matrix(2, [[1,1],[1,0]])
        matrix2 = cu.create_matrix(2, [[1,1],[0,1]])
        matrix3 = cu.add_matrices(matrix1, matrix2)
        return matrix1 + matrix2 == matrix3
    except Exception as err:
        print(err)
        return False

def test_matrix_multiplication():
    try:
        matrix1 = cu.create_matrix(2, [[1,1],[1,0]])
        matrix2 = cu.create_matrix(2, [[1,1],[0,1]])
        matrix3 = cu.multiply_matrices(matrix1, matrix2)
        return matrix1 * matrix2 == matrix3
    except Exception as err:
        print(err)
        return False

if __name__ == "__main__":
    print("Testing matrix creation")
    print("Test 1: ", "Passed" if test_matrix_creation() else "Failed")
    print("Testing matrix addition")
    print("Test 1: ", "Passed" if test_matrix_addition() else "Failed")
    print("Testing matrix multiplication")
    print("Test 1: ", "Passed" if test_matrix_multiplication() else "Failed")

    print(type(cu.random_permutation_matrix(5)))
   

    n,k = 10,5
    F = GF(11)
    C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
    H = C.parity_check_matrix()

    U, m2 = cu.gaussian_elimination(H, start_column=n-k)
    print(m2)
    print(m2.rank())