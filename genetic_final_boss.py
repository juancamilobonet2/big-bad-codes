from pyeasyga import pyeasyga
import random

def rand_base_individual(column_num):
    return [[random.randint(0, column_num-1), random.randint(0, column_num-1)]]
