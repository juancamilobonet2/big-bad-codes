import numpy as np

# cosas generales de coding theory

class Binary_Matrix:
    def __init__(self,n,k, matrix=[]):
        self.n = n
        self.k = k
        if matrix == []:
            self.matrix = np.zeros((k,n),dtype=int)
        else:
            self.matrix = matrix
    
    def __str__(self):
        return str(self.matrix)

    def get_matrix(self):
        return self.matrix
    
# operaciones modulo 2
def add_matrices(matrix1, matrix2):
    return (matrix1 + matrix2) % 2

def multiply_matrices(matrix1, matrix2):
    return (matrix1 @ matrix2) % 2

def find_syndrome(parity_check_matrix, vector):
    return (multiply_matrices(parity_check_matrix, vector.transpose()))

#test
if __name__ == "__main__":
    parity_check=np.array([[1,1,0,1,0,0],[0,1,1,0,1,0],[1,1,1,0,0,1]])
    r = np.array([[1,1,0,1,1,0]])
    print(find_syndrome(parity_check,r))
