'''
Created on 27 Jan 2020

@author: ja9508
'''
        
def create_random_population(pop_size, sol_size, sol_generator, evaluation_function, fitness_evaluations, change):
    """Initialises a random populaiton composed by randomly created individuals.

    Parameters:
    @param pop_size : size of the wanted population
    @param sol_size : the size of the solutions of the population
    @param sol_generator -- function 'sol_generator' to generate a random solutions
    """
    population = [None] * pop_size
    avg = 0
    for i in range(pop_size):
        chromosome = sol_generator(sol_size)
        fitness = evaluation_function(chromosome, change)
        fitness_evaluations += 1
        avg += fitness
        population[i] = (chromosome, fitness)
    return population, avg/pop_size, fitness_evaluations # Is it necessary to normalise the fitness values?

def create_random_individual(sol_size, sol_generator, evaluation_function, fitness_evaluations, change):
    chromosome = sol_generator(sol_size)
    fitness = evaluation_function(chromosome, change)
    return (chromosome, fitness), fitness_evaluations + 1

def print_population(population):
    """Print the individuals and their quality of given population

    Parameters:
    @param population : population wanted to print
    """
    for individual in range(len(population)):
        print(individual)
        
def population_calculate_average(population):
    """Calculates the average of a given population

    Parameters:
    @param population : population wanted to print. The population must be a list of tuples containing the 
            chromosome and its quality (in that order).
    """
    return sum(individual[1] for individual in population) / len(population)

def evaluate_population(population, evaluation_function, i_change, evaluations):
    evaluations += len(population)
    return [(sol, evaluation_function(sol, i_change)) for (sol, fit) in population], evaluations
        
