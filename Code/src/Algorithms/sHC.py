# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 16:14:10 2020
@author: ja9508
"""

import csv
import sys
import random
import timeit
import os.path
import math

sys.path.insert(1, "src/Individual")
import Permutation
sys.path.insert(1, "src/Algorithms")
import Population
import Distances
sys.path.insert(1, "src/Problems")
import DQAP


"""
Main
"""
def main():
    start_time = timeit.default_timer()
    wait_time = random.uniform(0,20)
    random.seed(1)
    
    # Parametrization of the required variables
    output_name = summary_name = local_optima_name = None
    problem = DQAP.generate_random_sQAP(10)
    create_new_member = Permutation.create_new_permutation
    maxGens = 1000
    
    # Algorithm´s parametrisation
    mutation_strategy = random_neighbour # swap_random
    
    # Rotation algorithms
    stuck_index = 0
    stuck = problem.size
    update_cooling = linear_cooling
    cooling_type = "linear"
    metric = "C"
    max_perms = problem.size - 1
    min_perms = 1
    
    # Dynamics initialisation
    alg_version = "restart" #standard, rotation1, rotation2, restart
    
    # Update variables if there are inputs
    if len(sys.argv) != 1:
        if "?" in sys.argv:
            print_usage()
            return 0
        else:
            parameters = dict([ (sys.argv[i], sys.argv[i+1]) for i in range(1,len(sys.argv),2)])
            if 'instance' in parameters and "qap" in parameters['instance']:
                problem = DQAP.read_qaplib_instance(parameters['instance'])
                create_new_member = Permutation.create_new_permutation
                stuck = problem.size
                max_perms = problem.size - 1
            if 'dynamic' in parameters: 
                problem.read_rotation_mask(parameters['dynamic'])
            if 'result' in parameters: 
                output_name = parameters['result']
            if 'summary' in parameters: 
                summary_name = parameters['summary']
            if 'stop' in parameters: 
                maxGens = int(parameters['stop'])
                stuck = int(maxGens * 0.01)
            if 'seed' in parameters: 
                random.seed(int(parameters['seed']))
            if 'algorithm' in parameters: 
                alg_version = parameters['algorithm']
            if 'metric' in parameters and "K" in parameters['metric']:
                metric = parameters['metric']
                mask_generator = Distances.random_perm_at_K_dist
                max_perms = Distances.choose(problem.size, 2)
                min_perms = 1
            elif 'metric' in parameters and "C" in parameters['metric']:
                metric = parameters['metric']
                mask_generator = Distances.random_perm_at_C_dist
                max_perms = problem.size - 1
                min_perms = 1
            elif 'metric' in parameters and "H" in parameters['metric']:
                metric = parameters['metric']
                mask_generator = Distances.random_perm_at_H_dist
                max_perms = problem.size
                min_perms = 2
            elif 'metric' in parameters and "U" in parameters['metric']:
                metric = parameters['metric']
                mask_generator = Distances.random_perm_at_U_dist
                max_perms = problem.size - 1
                min_perms = 1
            elif 'metric' not in parameters:
                max_perms = problem.size - 1
                min_perms = 1
                mask_generator = Distances.random_perm_at_C_dist
            if 'maxRate'in parameters: 
                max_perms = int(max_perms * float(parameters['maxRate']))
            if 'trials' in parameters: 
                stuck = int(parameters['trials'])
            if 'waiting' in parameters: 
                wait_time = float(parameters['waiting'])
            if 'cooling' in parameters:
                try:
                    cooling_type = parameters['cooling']
                    if type(int(parameters['cooling'])) == int:
                        update_cooling = linear_cooling
                        max_perms = int(parameters['cooling'])  
                        min_perms = int(parameters['cooling'])
                except:
                    if "linear" in parameters['cooling']:
                        update_cooling = linear_cooling
                    elif "exponential" in parameters['cooling']:
                        update_cooling = exponential_cooling
                        min_perms = math.log(float(min_perms)/max_perms) / maxGens # Stored in 'min_perms' just to use the same updating function structure

    # Initialise all the necessary parameters
    rotated_gens = generations = evaluations = 0
    change = 0
    max_rotation_gens = stuck
    changeDetected = rotated = False
    convergence = [["Generation","Best","Pop","Original","nRotations","Rotated","Distance","Stucked","Algorithm"]]
    summary = [["Algorithm","Metric","Best","nRotations","Trials","Distance","Time","Permutation"]]
    local_optimas = [["Metric","Distance","Algorithm","Mask","Encoding","Fitness"]]
    
        
    # Initialise population generating a random population (the solutions are already evaluated and the average is calculated too)      
    population, evaluations = Population.create_random_individual(problem.size, create_new_member, problem.evaluate, evaluations, change)
    best = population
    neighbourhood = all_swap_combinations(best[0])
    
    # Iteration: terminate?
    while (generations < maxGens):
        generations += 1
        
        # Get the distance by the generation
        distance = update_cooling(max_perms, min_perms, generations, maxGens)
        
        # Check if there has been an environmental change
        changeDetected = detect_change(population, problem.evaluate, change)
        
        # Apply dynamic approach
        if "restart" in alg_version and changeDetected:
            # Restart the entire population and sort it
            population, evaluations = Population.create_random_individual(problem.size, create_new_member, problem.evaluate, evaluations, change)
            neighbourhood = all_swap_combinations(population[0])
        elif "restart" in alg_version and stuck_index == stuck:
            # Store the information of the LO
            if len(neighbourhood) == 0:
                local_optimas.append([metric, distance, problem.masks[len(problem.masks)-1], alg_version, population[0], population[1]])
            
            # Restart the population
            problem.masks.append(problem.masks[0])
            population, evaluations = Population.create_random_individual(problem.size, create_new_member, problem.evaluate, evaluations, change)
            stuck_index = 0 
            neighbourhood = all_swap_combinations(population[0])
        elif "rotation2" in alg_version and rotated_gens == max_rotation_gens and rotated:
            # Return to initial instance
            undo_rotation = Permutation.compose((problem.masks[change]),population[0])
            change = 0
            rotated = False
            stuck_index = rotated_gens = 0 
            
            # Re-evaluate population
            population = (undo_rotation, problem.evaluate(undo_rotation, change))
            evaluations += 1
            
            neighbourhood = all_swap_combinations(best[0])
            
            # Decrease number of trials
            max_rotation_gens = int(math.ceil(float(max_rotation_gens)/2))
        elif "rotation2" in alg_version and stuck_index == stuck and not rotated: 
            # Store information of the LO
            if len(neighbourhood) == 0:
                local_optimas.append([metric, distance, problem.masks[len(problem.masks)-1], alg_version, population[0], population[1]])
            
            # Generate a random permu considering the input
            mask = mask_generator(problem.size, distance)
            
            # Append it to the problem masks and increment change index
            problem.masks.append(mask)
            change = len(problem.masks) - 1
            rotated = True
            stuck_index = 0 
            rotated_gens = 0
            
            # Re-evaluate population
            population = (population[0], problem.evaluate(population[0], change))
            evaluations += 1
            neighbourhood = all_swap_combinations(population[0])
            
        elif "rotation1" in alg_version and stuck_index == stuck and not rotated:
            # Store information of the LO
            if len(neighbourhood) == 0:
                local_optimas.append([metric, distance,problem.masks[len(problem.masks)-1],alg_version,population[0],population[1]])
            
            # Generate a random permu considering the input
            mask = mask_generator(problem.size, distance)
            
            # Append it to the problem masks and increment change index
            problem.masks.append(mask)
            change = len(problem.masks) - 1
            rotated = True
            stuck_index = 0 
            
            # Re-evaluate population
            population = (population[0], problem.evaluate(population[0], change))
            evaluations += 1
            
            neighbourhood = all_swap_combinations(population[0])
        
        elif "rotation1" in alg_version and stuck_index == stuck and rotated:
            # Return to initial instance
            undo_rotation = Permutation.compose((problem.masks[change]),population[0])
            change = 0
            rotated = False
            
            # Re-evaluate population
            population = (undo_rotation, problem.evaluate(undo_rotation, change))
            evaluations += 1
            
            # Store information of the LO
            if len(neighbourhood) == 0:
                local_optimas.append([metric, distance,problem.masks[len(problem.masks)-1],alg_version,population[0],population[1]])
            
            # If new pop after a rotation is better, update best solution
            if population[1]<best[1]: 
                best = population
                neighbourhood = all_swap_combinations(best[0])
            else: 
                population = best
            
            # Store info and go to begining of the iteration
            convergence.append([generations,best[1],population[1],problem.evaluate(population[0]),len(problem.masks)-1,rotated,distance,stuck_index,alg_version])
            continue
        else:
            if changeDetected:
                # Re-evaluate the population and sort it
                population, evaluations = Population.evaluate_population(population, problem.evaluate, change, evaluations)
        
        # Store the information of each generation
        convergence.append([generations,best[1],population[1],problem.evaluate(population[0]),len(problem.masks)-1,rotated,distance,stuck_index,alg_version])
        
        # Generation of the new sample
        neighbour_chromosome, neighbourhood = mutation_strategy(neighbourhood) if len(neighbourhood) > 0 else swap_random(population)
        neighbour = (neighbour_chromosome, problem.evaluate(neighbour_chromosome, change))
        
        # Choose best, check if the search is stuck and update population
        new_population = DQAP.return_best(population,neighbour)
        if population == new_population:
            stuck_index += 1
        else:
            stuck_index = 0
            neighbourhood = all_swap_combinations(new_population[0])
        population = new_population
        
        # Increase rotated generations (for 'rotation2' algorithm)
        rotated_gens += 1
        
        # Update best from initial instance
        if change == 0 and population[1]<best[1]: 
            best = population
            neighbourhood = all_swap_combinations(best[0])
    
    # Write the results on a csv file
    print(output_name," execution time: ", timeit.default_timer() - start_time)
    if output_name == None:
        pass
    else:
        random.seed()
        print("Execution paused ", wait_time, " seg. to avoid data loosing.")
        timeit.time.sleep(wait_time)
        if os.path.isfile(output_name): convergence.pop(0) # Remove headers
        with open(output_name, 'a+') as file:
            writer = csv.writer(file)
            writer.writerows(convergence)
            
    if summary_name != None:
        try:
            summary.append([alg_version,metric,best[1],len(problem.masks)-1,stuck,cooling_type,timeit.default_timer() - start_time,best[0]])
        except:
            summary.append([alg_version,"",best[1],len(problem.masks)-1,stuck,cooling_type,timeit.default_timer() - start_time,best[0]])
        random.seed(wait_time)
        wait_time = random.uniform(0,10)
        print("Execution paused ", wait_time, " seg. to avoid data loosing.")
        timeit.time.sleep(wait_time)
        if os.path.isfile(summary_name): summary.pop(0) # Remove headers
        with open(summary_name, 'a+') as file:
            writer1 = csv.writer(file)
            writer1.writerows(summary)
            
    if local_optima_name != None:
#        wait_time = random.uniform(0,10)
        print("Execution paused ", wait_time, " seg. to avoid data loosing.")
        timeit.time.sleep(wait_time)
        if os.path.isfile(local_optima_name): local_optimas.pop(0) # Remove headers
        with open(local_optima_name, 'a+') as file:
            writer1 = csv.writer(file)
            writer1.writerows(local_optimas)
            
    return convergence

"""
COMPLEMENTARY FUNCTIONS
"""        
def print_usage():
    my_list = ["HC.py","instance","[intance file]","dynamic","[dynamic file]","result","[results file]",
               "stop","[maximum generations]","seed","[seed]",
               "algorithm","{standard, rotation1, rotation2, restart}*",
               "metric","{K,C,U,H}","trials","[generations looking for inprovement]",
               "cooling","{1:maximum distance,linear,exponential}",
               "waiting","[seconds]"]
    print(" ".join(my_list))
    
def detect_change(population, evaluation_function, i_change):
    if population[1] != evaluation_function(population[0], i_change):
        return True
    return False

def print_output(result):
    print("\n".join(map(str,result)))
    
def staggered_cooling(max_dist, min_dist, generation, maximum_generations):
    if max_dist != min_dist:
        return((int)(max_dist - (float(max_dist) * generation / maximum_generations) + 1))
    return(max_dist)
    
def linear_cooling(max_dist, min_dist, generation, maximum_generations):
    return(int(round((min_dist - max_dist) * (float(generation - 1)/(maximum_generations - 1))) + max_dist))
    
def exponential_cooling(max_dist, min_dist, generation, maximum_generations):
    return(int(round(max_dist * math.exp(min_dist * generation))))

def swap_random(seq):
    copy = [i for i in seq[0]]
    idx = range(len(copy))
    i1, i2 = random.sample(idx, 2)
    copy[i1], copy[i2] = copy[i2], copy[i1]
    return(copy,[])

def random_neighbour(neighbourhood):
    idx = random.randrange(len(neighbourhood))
    sample = neighbourhood[idx]
    neighbourhood.pop(idx)
    return sample, neighbourhood  

def all_swap_combinations(s):
    def combinations(iterable, r):
        pool = tuple(iterable)
        n = len(pool)
        for indices in all_permutations(range(n), r):
            if sorted(indices) == list(indices):
                yield tuple(pool[i] for i in indices)
                
    def all_permutations(iterable, r = None):
        pool = tuple(iterable)
        n = len(pool)
        r = n if r == None else r
        if r > n:
            return
        indices = list(range(n))
        cycles = list(range(n, n-r, -1))
        yield tuple(pool[i] for i in indices[:r])
        while n:
            for i in reversed(range(r)):
                cycles[i] -= 1
                if cycles[i] == 0:
                    indices[i:] = indices[i+1:] + indices[i:i+1]
                    cycles[i] = n - i
                else:
                    j = cycles[i]
                    indices[i], indices[-j] = indices[-j], indices[i]
                    yield tuple(pool[i] for i in indices[:r])
                    break
            else:
                return
    result = []
    for idx1, idx2 in combinations(range(len(s)),2):
        swapped_s = list(s)
        swapped_s[idx1], swapped_s[idx2] = swapped_s[idx2], swapped_s[idx1]
        result.append((swapped_s))
    return result
    
if __name__ == "__main__":
    result = main()
