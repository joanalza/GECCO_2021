"""
Created on 24 Nov 2020
@author: ja9508
"""
    
import random
import math
import sys
import timeit

import Permutation

## number of perms at each dist
def num_perms_at_K_dist(n):
    sk = [[0] * int(n*(n-1)/2+1) for i in range(n+1)]
    for i in range(n+1):
        sk[i][0] = 1
    for i in range(1,1+n):
        for j in range(1,int(i*(i-1)/2+1)):
            if j - i >= 0 :
                sk[i][j] = sk[i][j-1]+ sk[i-1][j] - sk[i-1][j-i]
            else:
                sk[i][j] = sk[i][j-1]+ sk[i-1][j]
    return sk

def cumulative_sum(lists): 
    return [sum(lists[0:x + 1]) for x in range(0, len(lists))] 

def v2ranking(v, n):
    rem = list(range(n))
    rank = [0] * n# np.zeros(n,dtype=np.int)
    for i in range(len(v)):
        rank[i] = rem[v[i]]
        rem.pop(v[i])
    return rank

## random permutations at Cayley distance
def Cayley_by_X(n, dist):  
    # Generate a random binary decomposition vector
    x = [0] * (n-1)
    for j in range(dist): x[j] = 1
    random.shuffle(x)  
    
    # Generate perm from X vector
    sigma = [i for i in range(n)]
    for k in range(len(x)):
        if x[k] == 1:
            j = random.randint(k+1,n-1)
            sigma[k],sigma[j] = sigma[j],sigma[k]
    return tuple(sigma)

def random_perm_at_C_dist1(n, dist):
    def stirling(n,k):
        n1=n
        k1=k
        if n<=0:
            return 1
         
        elif k<=0:
            return 0
         
        elif (n==0 and k==0):
            return -1
         
        elif n!=0 and n==k:
            return 1
         
        elif n<k:
            return 0
     
        else:
            temp1=stirling(n1-1,k1)
            temp1=k1*temp1
            return (k1*(stirling(n1-1,k1)))+stirling(n1-1,k1-1)
        
    long_cycle = [False] * n
    sigma = [i for i in range(n)]
    k = n - dist
    n1 = n
    while k > 1 :
        ran1 =random.random()
        if ran1 < (stirling(n-1, k-1) / stirling(n, k)):
            long_cycle[ n - 1 ] = False;
            k = k - 1
        else:
            long_cycle[ n - 1 ] = True;
        n = n - 1;
        
    sigma_inv_ = Permutation.create_new_permutation(n1)
    for i in range(n-1):
        sigma[sigma_inv_[ i ]] = sigma_inv_[i+1];
    sigma[sigma_inv_[n-1]] = sigma_inv_[0];
    
    for i in range(n,n1):
        if long_cycle[i]:
            ran2 = random.randint(0,i-1)
            sigma[ i ] = sigma[ ran2 ];
            sigma[ ran2 ] = i;
    return tuple(sigma)

# Adapted random partition, return a list of columns
def random_inv_FS1(n,d):
    lis = n-d
    if lis<0: return False
    results = [1 for i in range(lis)]
    for value in range(d):
        x = random.randrange(lis)
        results[x] = results[x] + 1
    return sorted(results)

def random_inv_FS2(n,d):
    lis = n-d
    if lis<0: return False
    something = list(range(n))
    results = []
    for p in subsets_k(something, lis):
        fs = sorted([len(i) for i in p])
        if fs not in results:
            results.append(fs)
    return results[random.randrange(len(results))]

def subsets_k(collection, k): yield partition_k(collection, k, k)
def partition_k(collection, min, k):
  if len(collection) == 1:
    yield [ collection ]
    return

  first = collection[0]
  for smaller in partition_k(collection[1:], min - 1, k):
    if len(smaller) > k: continue
    if len(smaller) >= min:
      for n, subset in enumerate(smaller):
        yield smaller[:n] + [[ first ] + subset]  + smaller[n+1:] 
    if len(smaller) < k: yield [ [ first ] ] + smaller

def generate_SYT(n,invFS):
    def corner_box(fs):
        corners =[]
        elems = []
        for i in range(len(fs)):
            if fs[i] not in elems and fs[i] != 0:
                corners.append(i)
                elems.append(fs[i])      
        return corners
    
    def hook_length_formula(fs):
        f = math.factorial
        p=j=-1
        d={}
        for col in fs:
            j=j+1
            i=0
            while i<col:  
                a=d[i]=d.get(i,j)
                p=p*(col-i+j-a)
                i+=1
        return f(sum(fs))/(-p)
    
    # Generate Ferrers diagram (by rows)
    invFS_copy = [w for w in invFS]
    SYT = [[0]*i for i in invFS]
    
    for r in range(n-1,-1,-1):
        # Get corners
        corners = corner_box(invFS_copy)
        probs = []
        
        # For each corner
        for cb in corners:
            # Delete box
            sub_diagram = [invFS_copy[icol] - 1 if icol == cb else invFS_copy[icol] for icol in range(len(invFS_copy))]
            
            # Calculate hook length product
            # Add it to probabilities
            probs.append(hook_length_formula(sub_diagram))
        probs = [float(i)/sum(probs) for i in probs]
        col = weighted_choice(corners, probs)
        
        last = max([i for i, e in enumerate(SYT[col]) if e == 0])
        SYT[col][last] = r
        invFS_copy[col] = invFS_copy[col] - 1
    
    SYT.reverse()
    return SYT

def weighted_choice(choices, weights):
    r = random.random()
    upto = 0
    for i in range(len(choices)):
        if upto + weights[i] >= r:
            return choices[i]
        upto += weights[i]
    # Shouldnt get here
    return random.sample(choices,1)[0]

def choose(n,r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)

## random permutations at Kendalls-t distance
def random_perm_at_K_dist(n, dist):
    # param sk is the results of the function num_perms_at_dist(n)
    sk = num_perms_at_K_dist(n)
    i = 0
    probs = [0.0] * (n)
    v = [0] * n
    while i<n and dist > 0 :
        rest_max_dist = (n - i - 1 ) * ( n - i - 2 ) / 2
        if rest_max_dist  >= dist:
            probs[0] = sk[n-i-1][dist]
        else:
            probs[0] = 0
        mi = min(dist + 1 , n - i )
        for j in range(1,mi):
            if rest_max_dist + j >= dist: probs[j] = float(sk[n-i-1][ dist-j])
            else: probs[ j ] = 0.0
        prob_sum = sum(probs[:mi])
        probabilities = [0.0] * mi
        for k in range(mi):
            probabilities[k] = probs[k] / prob_sum
        v[i] = weighted_choice(range(mi),probabilities)
        dist -= v[i]
        i += 1
    return tuple(v2ranking(v,n))

def perms_at_dist_tab(n):
   s = [[1]]
   for i in range(n):
       s.append([1])
       for d in range(i):
           s[i+1].append(i * s[i][d] + s[i][d+1])
       s[i+1].append(0)
   return s

#Then, the following code does uniform random sampling of the permutations at distance d of the identity permutation:

def random_perm_at_C_dist(n, d):
   Snd = perms_at_dist_tab(n)
   # "Random permutation at a given distance of the identity permutation"
   p = list(range(n)) # This could also be any other "centre" permutation
   for i in range(n-1, 0, -1):
       if random.random() >= Snd[i][d] / Snd[i+1][d]:
           j = random.randrange(i)
           p[i], p[j] = p[j], p[i]
           d -= 1
   return tuple(p)

# This is basically just the first loop of the other implementations, where instead
# of recording that a cycle must be introduced at a certain point, we just perform
# a swap effectively introducing such a new cycle. This is also nearly identical
# to the classic so-called Knuth's shuffle algorithm, which is simply:

def shuffle(n):
   p = list(range(n))
   for i in range(n-1, 0, -1):
       j = random.randrange(i+1)
       p[i], p[j] = p[j], p[i]
   return p



def random_perm_at_H_dist(n,dist):
    # Keep shuffling the array
    # If the new permutation is not valid (v[i]==i), we break and start from scratch.
    if dist<2: return False
    identity = [i for i in range(n)]
    derange = random.sample(identity,dist)
    while True:
        v = [i for i in range(len(derange))]
        for j in range(len(v) - 1, -1, -1):
            p = random.randint(0, j)
            if v[p] == j:
                break
            else:
                v[j], v[p] = v[p], v[j]
        else:
            for i in range(len(derange)):
                identity[derange[i]] = derange[v[i]]
            return tuple(identity) #Used for counter
            # return identity
    return False

def random_perm_at_U_dist(n,dist):
    def search_on_list(list_lists, elem):
        for i in range(len(list_lists)):
            if elem in list_lists[i]:
                return i,list_lists[i].index(elem)
        return False
                
        
    # Generate a Ferrers diagram (shape)
    fd = random_inv_FS2(n,dist)
    # fd = [1,1,3]
    
    # Generate two standard (young) tableaux of a given shape
    P = generate_SYT(n,fd)
    Q = generate_SYT(n,fd)
    
    # Construct permu from standard tableux (Schensted Correspondence)
    sigma = [-1] * n
    for r in range(n-1,-1,-1):
        i,j = search_on_list(Q, r)
        u = P[i][j] 
        del(P[i][j])
        del(Q[i][j])
        while i != 0:
            i = i-1
            largest = []
            for y in P[i]: 
                if y < u: largest.append(y)
            largest = P[i].index(max(largest))
            u, P[i][largest] = P[i][largest],u
        sigma[n-r-1] = u
    return sigma
    
    
if __name__ == "__main__":   
    
    start_time = timeit.default_timer()
    
    N = 10
    D = 4

    # enumerate all derangements for testing
    from collections import Counter
    counter = []
    
    # make M probes for each derangement
    M = 10000
    # M=10
    for _ in range(M):
        # generate a random derangement
        p = random_perm_at_U_dist(N,D)
        # is it really?
        # assert p in counter
        # ok, record it
        counter.append(p)
    
    end_time = timeit.default_timer()
    
    print(len(Counter(counter)))
    print(Counter(counter)) # Needs perms as tuple()
    
    print("Time spent for ",M," iterations: ", end_time - start_time)
        