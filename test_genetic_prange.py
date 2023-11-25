import pytest
import genetic_prange as gp

def test_crossover():
    ind1 = [(0, 1), (1, 2), (3, 2)]
    ind2 = [(0, 2), (1, 1)]
    for _ in range(1000):
        ind3, ind4 = gp.crossover(ind1, ind2)
        # Assert no new elements were created.
        assert len(ind1)+len(ind2) == len(ind3)+len(ind4)
        # Assert all elements of new inds existed in the old ones.
        assert all([perm in ind1 or perm in ind2 for perm in ind3])
        assert all([perm in ind1 or perm in ind2 for perm in ind4])
        ind1, ind2 = ind3, ind4

def test_mutation():
    ind1 = [(0, 1)]
    ind1_2 = gp.mutation(ind1, 2)
    assert ind1_2 == [(0, 1)] or ind1_2 == [(1, 1)] or ind1_2 == [(0, 0)]
    ind2 = ind1
    ind2_2 = ind1_2
    for i in range(10):
        ind2 = ind2 + ind1_2
        ind2_2 = gp.mutation(ind2, 10)
        # The length is always maintained after the operation.
        assert len(ind2_2) == 2 ** (i+1)

