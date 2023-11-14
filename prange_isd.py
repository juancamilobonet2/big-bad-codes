import random
from typing import List

N = 12  # Assuming a default size for N (adjust as needed)
K = 4  # Assuming a default size for K (adjust as needed)

def fisher_yates_shuffle(src: List[int]) -> List[int]:
    dst = src.copy()
    for i in range(len(dst) - 1, 0, -1):
        idx = random.randint(0, i)
        dst[idx], dst[i] = dst[i], dst[idx]
    return dst

def permute_matrix_columns(src: List[List[int]], new_pos: List[int]) -> List[List[int]]:
    dst = []
    tmp_row = []

    start_pos = list(range(N))
    new_pos = fisher_yates_shuffle(start_pos)

    for i in range(N - K):
        tmp_row = [src[i][new_pos[j]] for j in range(N)]
        dst.append(tmp_row.copy())
    
    return dst

def gaussian_elimination(h: List[List[int]], s: List[int], sigma: List[int]) -> bool:
    tmp_h = [row.copy() for row in h]
    sigma[:] = s.copy()
    tmp_b = 0

    for col in range(K):
        row = 0
        while row < N - K:
            if tmp_h[row][col]:
                tmp_h[row], tmp_h[row] = tmp_h[row], tmp_h[row]
                tmp_b = sigma[row]
                sigma[row] = sigma[row]
                sigma[row] = tmp_b
                break
            row += 1
        
        # print("row", row)
        # print("tmp_h", tmp_h, len(tmp_h), len(tmp_h[0]))
        if not tmp_h[row-1][col]:
            continue

        for i in range(N - K):
            if i != row and tmp_h[i][col]:
                tmp_h[i] = [a ^ b for a, b in zip(tmp_h[i], tmp_h[row])]
                if sigma[row]:
                    sigma[i] = 1 - sigma[i]
        
        row += 1

    for row in tmp_h:
        if not row[K:]:
            return False
    
    return True

def randomize(h: List[List[int]], s: List[int], p: List[int], sigma: List[int]) -> None:
    hp = []
    full_rank = False
    while not full_rank:
        hp = permute_matrix_columns(h, p)
        full_rank = gaussian_elimination(hp, s, sigma)

def search_prange(sigma: List[int], error_weight: int, epsilon: List[int]) -> bool:
    sigma_weight = sum(sigma)
    success = (sigma_weight == error_weight)
    epsilon[:] = [0] * K
    if success:
        for i in range(N - K):
            epsilon[K + i] = sigma[i]
    return success

def isd(h: List[List[int]], y: List[int], w: int, e: List[int]) -> None:
    # print("h:", h)
    # print("N:", N)
    # print("K:", K)
    # print("y:", y)
    # for i in range(len(h)):
    #     for j in range(len(h[i])):
    #         print("i", i)
    #         print("j", j)
    #         print("h", h)
    #         print("h[i][j]", h[i][j])
    #         print("y", y)
    #         print("y[j]", y[j])
    #         print("h[i][j] & y[j]", h[i][j] & y[j])
    #         print(sum([h[i][j] & y[j]]))
    #         print("---------")
    s = [sum([(h[i][j] & y[j]) for j in range(len(h[i]))]) % 2 for i in range(len(h))]
    sig = [0] * (N - K)
    eps = [0] * N
    p = list(range(N))
    success = False

    random.seed()

    while not success:
        randomize(h, s, p, sig)
        success = search_prange(sig, w, eps)

    e[:] = eps.copy()

# Example usage
h = [
    [0,0,1,0,1,0,0,0,0,0,0,1],
    [0,1,1,1,0,0,0,0,1,0,0,0],
    [0,0,0,0,1,1,0,1,0,1,1,1],
    [1,1,0,0,0,1,1,1,0,0,0,0],
    [1,1,1,1,0,0,1,0,1,1,1,1],
    [0,0,1,0,0,0,0,0,1,1,1,0],
    [0,1,0,0,0,0,0,1,0,1,0,1],
    [0,1,1,1,1,0,1,1,1,0,0,1]
]  # Replace with your actual matrix
y = [0,0,0,0,1,0,0,0,0,0,0,0]  
# Replace with your actual vector
w = 1  # Replace with your actual weight

e = [0] * N
isd(h, y, w, e)
print("Resulting error vector:", e)
