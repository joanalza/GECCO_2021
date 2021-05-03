# -*- coding: utf-8 -*-

'''
Created on 23 Nov 2020

@author: ja9508
'''

########## qap.py ##########
import random
import  sys

Infinity = 1.e10000
LOG = True	# whether or not to print intermediate solutions

class DQAP():

#     Initialisation for the class TSP given the number of cities and the distance matrix.
#     It is possible to use either symmetric or asymmetric TSP instances. 
    def __init__(self, siize =None, distance_matrix=None, flow_matrix=None, rotation_list = None):
        self.size = siize
        self.d_matrix = distance_matrix
        self.f_matrix = flow_matrix
        if rotation_list == None:
            self.masks = [list(range(self.size))]
        else:
            self.masks = rotation_list
    
    def read_rotation_mask(self, filename):
        try:
            with open(filename) as f:
                self.masks = []
                for i, line in enumerate(f):
                    if "," in line:
                        self.masks.append([int(n) for n in line.strip().split(',')] )
        except:
            print("Dynamic environment unknown. Running static version...") 
    
    def evaluate(self, genes, change_index = 0):
        if len(genes) != self.size: return False
        cost = 0
        for i in range(0, len(genes)):
            for j in range(0, len(genes) ):
                cost += self.f_matrix[i,j] * self.d_matrix[self.masks[change_index][genes[i]],self.masks[change_index][genes[j]]]
        return cost      

    def evaluate__(self, genes, change_index = 0):
        """Evaluate solution 'genes' and 
        create additional cost information for incremental evaluation."""
        if len(genes) != self.size: return False
        delta = {}
        for i in range(0, len(genes) ):
            for j in range(0, len(genes) ):
                delta[i,j] = 0
                for k in range(0, len(genes) ):
                    delta[i,j] += self.f_matrix[i,k] * self.d_matrix[j,genes[k]]
        cost = 0
        for i in range(0,len(genes)):
            cost += delta[i,genes[i]]
        return cost,delta


def generate_random_sQAP(problem_size):
    """Make data for a random problem of size 'n'."""
    f = {}      # for holding n x n flow matrix
    d = {}      # for holding n x n distance matrix
    
    for i in range(problem_size):
        f[i,i] = 0
        d[i,i] = 0
    for i in range(problem_size-1):
        for j in range(i+1,problem_size):
            f[i,j] = random.randint(1,99)
            f[j,i] = f[i,j]
            d[i,j] = random.randint(1,99)
            d[j,i] = d[i,j]

    return DQAP(siize =problem_size, distance_matrix=d, flow_matrix=f)

def return_best(solution1, solution2):
    if solution1[1] < solution2[1]:
        return solution1
    else:
        return solution2

def generate_random_aQAP(problem_size):
    """Make data for a random problem of size 'n'."""
    f = {}      # for holding n x n flow matrix
    d = {}      # for holding n x n distance matrix
    
    for i in range(problem_size):
        f[i,i] = 0
        d[i,i] = 0
    for i in range(problem_size-1):
        for j in range(i+1,problem_size):
            f[i,j] = random.randint(0,10)
            f[j,i] = random.randint(0,10)
            d[i,j] = random.randint(1,10)
            d[j,i] = random.randint(1,10)

    return DQAP(siize =problem_size, distance_matrix=d, flow_matrix=f)

def read_qaplib_instance(filename):
    """Read data for a QAP problem from file in QAPLIB format."""
    try:
        if len(filename)>3 and filename[-3:] == ".gz":  # file compressed with gzip
            import gzip
            f = gzip.open(filename, "rb")
        else:   # usual, uncompressed file
            f = open(filename)
    except IOError:
        print("could not open file", filename)
        exit(-1)

    data = f.read()
    f.close()

    try:
        pass
        data = data.split()
        n = int(data.pop(0))
        f = {}  # for n times n flow matrix
        d = {}  # for n times n distance matrix
        for i in range(n):
            for j in range(n):
                f[i,j] = int(data.pop(0))
        for i in range(n):
            for j in range(n):
                d[i,j] = int(data.pop(0))
    except IndexError:
        print("inconsistent data on QAP file", filename)
        exit(-1)
    return DQAP(siize =n, distance_matrix=d, flow_matrix=f)


def is_assymmetric(problem):
    a = set([])
    for i in range(problem.size):
        for j in range(problem.size):
            if problem.d_matrix[i,j] != problem.d_matrix[j,i]:
                a.add("d")
            elif problem.f_matrix[i,j] != problem.f_matrix[j,i]:
                a.add("f")
                
    if a == [] : return True
    else: return a 

if __name__ == "__main__":
    rndseed = 1
    print("1")
    if len(sys.argv) == 1:
        print("1.1")
        problem = read_qaplib_instance("Problem_instances/QAP/lipa40a.qap")
    else:
        print("1.2")
        instance = sys.argv[1]
        problem = read_qaplib_instance(instance)  # create the distance matrix
    print("Assymmetric matrix: " + str(is_assymmetric(problem)))
    print("done")
        