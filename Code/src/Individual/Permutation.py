'''
Created on 24 Jan 2020

@author: ja9508
'''

import random


def create_new_permutation(permu_size):   
    """Generates a u.a.r. permutation
    
        @param permu_size : size of the permutation wanted to generate.
    """
    # Fill the vector with the values 1, 2, 3, ..., n  
    aux_list = list(range(permu_size))
    new_perm = []
  
    # While vector has elements get a random number from the vector and print      
    while aux_list != []: 
        aux_size = len(aux_list)
        # Generate a random number within the index range 
        index = random.randint(0, aux_size - 1) 
      
        # Get random number from the vector  
        new_perm.append(aux_list[index])  
      
        # Remove the number from the vector (swap the last element and the chosen index element, and pop it)
        aux_list[index], aux_list[aux_size - 1] = aux_list[aux_size - 1], aux_list[index]
        aux_list.pop()  
    return new_perm

def is_permutation(permutation):
    """Check whether the given list is a permutation or not
    
        @param permutation : list of integer values between 0 and len(list) - 1
    """
    repeated = []
    for i in range(len(permutation)):
        if i not in permutation or i in repeated:
            return False
        repeated.append(i)
    return True

""" Generates a random permutation using the package numpy
"""
def sample_random_permutation(permu_size):
    new_perm = random.sample(range(permu_size),permu_size)
    return new_perm
#    return list(numpy.random.permutation(permu_size))

def inverse(sigma):
    permutation = list(range(len(sigma)))
    return [x for _,x in sorted(zip(sigma,permutation))]

def reverse(sigma):
    return list(reversed(sigma))

def compose(permutation1,permutation2):
    if len(permutation1) != len(permutation2): return False
    permutation = []
    for i in range(len(permutation1)):
        permutation.append(permutation1[permutation2[i]])
    return permutation