#!/usr/bin/python
# Author: Qijin Chau
#
# HOW TO RUN: All you have to do is run backpack.py :)     
#
# Questions + Answers:
#
# 1. Define the problem as a genetic algorithm
#       The individuals/people of a population will be represented as
#       the possible backpacks that can be filled using the given boxes.
#       
#       These individuals will have a discrete representation as
#       their genetic encoding. See full description in question 2.
#       
#       A fitness level will be given to each of these individuals
#       and the generation will evolve with the fittest individuals
#       from the previous children along with children of these
#       individuals.
#       
#       In my genetic algorithm, I use Elitism and keep the top half 
#       of the old generation and have that generation reproduce and
#       regenerate the rest of the population. This as a whole is the
#       culling technique of 50% called for in the spec.
#       
#       Thus 50% of the new generation will be the elite half from
#       the old generation and the remaining 50% will be the children
#       of these elite individuals. 
# 
#
# 2. Provide the genome for the problem
#       As stated above, individuals will have a discrete representation
#       as their genetic encoding. This discrete representation will be
#       a 12-element binary array because there are also 12 possible
#       boxes that can go inside each backpack. 
#       
#       Each index of the array will corresspond to the BOXES array 
#       (defined after the box class definition) which will hold all the boxes 
#       that can be in the backpacks along with their value and weight.
#       
#       Thus, for example, if the genetic encoding has a 1 in the first
#       index, then the individual or the backpack will contain Box 1.
#       Having these boxes will represent the phenotypes resulting from
#       having the gene to contain these boxes. 
# 
# 
# 3. Define all the fringe operations
#
#       Fringe operations used in this algorithm are single-point/multi-point 
#       mutations and one-section/multi-section crossovers. 
#           
#           - A single-point mutation will randomly modify a single gene of a 
#             given individual. 
#           - A multi-point mutation will randomly modify 1-4 genes of a 
#             given individual. 
#
#           - A one-section crossover will take two individuals and randomly 
#             select a point on each of the two individuals and swapping the 
#             genetic material around this point. 
#           - A multi-section crossover will take two individuals and randomly
#             select 1-4 sections on each of the two individuals and swapping 
#             the genetic material around this point.
#       
#       My genetic algorithm uses a 8% probability for single point/section 
#       fringe operations and 3% probability for multi point/section fringe 
#       operations to occur on a child that was produced.
# 
#
# 4. Cull your population by 50% at every generation
#
#       A variation of Rank-based selection that consists of taking only 
#       the top 50% of individuals in the ranked list of the current 
#       population and reproducing with these elite 50% to regenerate the 
#       rest of the population.
#
#
# NOTES:
#   Using the Genetic Algorithm, I found that there are two possible
#   options for the most valuable backpack. 
#   
#   These backpacks contain either
#       Boxes: 1,2,5,6,7,8,10  or  Boxes: 1,2,3,6,8,10,11
#   
#   This is because Box 5 + Box 7 = Box 3 + Box 11 in terms of both
#   value and weight.


import random
import numpy as np


# Define box object that can go in the backpack
# Has a value and a weight
class box:
    def __init__(self, value, weight):
        # The value of the box which we want to maximize overall
        self.value = value
        
        # The weight of the box
        # Sum of all the box weights has to be less than or equal to 250
        self.weight = weight



# Set all the boxes available with their value and weight respectively.
# Each element and its index in this BOXES array will correspond to the 
# indices of the individuals in each population because these 
# individuals will also be represented as an array. 
BOXES = [box(6, 20), box(5, 30), box(8, 60), box(7, 90), 
        box(6, 50), box(9, 70), box(4, 30), box(5, 30),
        box(4, 70), box(9, 20), box(2, 20), box(1, 60)]


# Weight capacity that the backpack can hold 
CAPACITY = 250

# Size of initial population
# Each individual will be represented as a 12-element
# array filled with some permutation of 0s and 1s
INITIAL_POPULATION_SIZE = 250

# Maximum number of generations that the algorithm will go through
MAX_NUMBER_OF_GENERATION = 20

# Optimum fitness score for the backpack
OPTIMUM = 44

# Worst fitness score possible for a backpack
MINIMUM_FITNESS = 0

# 8% change for a single-point fringe operation to occur
SINGLE_FRINGE_OP_THRESHOLD = 8

# 3% change for a multi-point fringe operation to occur
MULTI_FRINGE_OP_THRESHOLD = 3



# Initialize the starting population based on the max population passed in
def InitializePopulation(max_population):
    return [generateIndividual() for individual in range (0, max_population)]


# Creates an array of the length of the number of boxes we have so 12.
# Then, in this array put in a random permutation of 0's and 1's.
# A "1" means we have the box in the backpack, a "0" means we don't have it.
# The indices correspond to the BOXES indices declared at the top.
def generateIndividual():
    return [random.randint(0, 1) for count in range (0, len(BOXES))]


# Create a child using a random crossover point from the parents
def reproduce(mother, father):
    cutoff = random.randint(0, len(mother))
    child = mother[:cutoff] + father[cutoff:]

    return child


# Returns the reproduction probability of the individual proportionate
# to the rank of the individual in its population based on its fitness
def getReproductionProbability(individual, population):
    total_fitness = 0
    for individual in population:
        total_fitness = total_fitness + fitness(individual)

    return fitness(individual)/total_fitness


# Returns the fitness value of the individual/permutation passed in.
# A higher fitness score means a more valuable backpack and 
# a more fit individual overall. 
def fitness(permutation):
    total_value = 0
    total_weight = 0
    index = 0
    for box in permutation:        
        if index >= len(BOXES):
            break

        # If the index is 1, then the current box is in the backpack 
        # for this permutation so account for this box in the calculation
        if (box == 1):
            total_value += BOXES[index].value
            total_weight += BOXES[index].weight
        
        index = index + 1
        
    # If total weight exceeds the weight capacity, then we
    # can't consider this permutation so return 0 for fitness
    if total_weight > CAPACITY:
        return 0
    else:
        # total_value is the same as the fitness of the backpack/individual
        return total_value


# Fringe operation that changes a random element/gene in the permutation 
# array from a 0 -> 1 (not having the box in the backpack to having it)
# or 1 -> 0 (having the box in the backpack to not having it)
def mutate(permutation):
    index_to_mutate = random.randint(0, len(permutation)-1)
    if permutation[index_to_mutate] == 1:
        permutation[index_to_mutate] = 0
    else:
        permutation[index_to_mutate] = 1


# Fringe operation that performs the mutation fringe operation above
# on multiple points
def multi_mutate(permutation, num_points):
    for i in range(0, num_points):
        mutate(permutation)


# Fringe operation that crosses over two random sections of two parent 
# permutation arrays which represent their genes in order to produce 
# two children with these mixed sections/genes
def crossover(first_permutation, second_permutation):
    threshold = random.randint(1, len(first_permutation)-1)
    
    temp_1 = first_permutation[threshold:]
    temp_2 = second_permutation[threshold:]
    
    first_permutation = first_permutation[:threshold]
    second_permutation = second_permutation[:threshold]
    first_permutation.extend(temp_2)
    second_permutation.extend(temp_1)


# Fringe operation that performs the crossover fringe operation above
# on multiple sections
def multi_crossover(first_permutation, second_permutation, num_sections):
    for i in range(0, num_sections):
        crossover(first_permutation, second_permutation)



# Use Elitism of top 50% to reproduce and regenerate the rest of the population 
def evolvePopulation(population):
    # 8% chance for a mutation or crossover fringe operation
    single_mutation_threshold = SINGLE_FRINGE_OP_THRESHOLD
    single_crossover_threshold = SINGLE_FRINGE_OP_THRESHOLD

    # 3% chance for a multi-mutation or multi-crossover fringe operation
    multi_mutation_threshold = MULTI_FRINGE_OP_THRESHOLD
    multi_crossover_threshold = MULTI_FRINGE_OP_THRESHOLD

    # Note: "parents" represent the current generation or population
    # that we are going to use to create the new generation.

    # Sort the parents based on fitness in ascending order
    all_parents = sorted(population, key=lambda x: fitness(x))

    # Size of the current population passed in which will be the
    # same size as the new generation that we want
    new_generation_size = int(len(all_parents))

    # Get the elite 50% of the current population
    elite_parents = all_parents[int((new_generation_size / 2)):]

    # Initialize the new generation 
    new_generation = elite_parents
    

    # Half this new generation will be the elite 50% of the old generation and
    # half this new generation will be the children made from random parents
    # of this elite 50% reproducing
    while len(new_generation) < new_generation_size:
        # Pick a random mother and father from the top half of the
        # current generation
        mother = elite_parents[random.randint(0, len(elite_parents)-1)]
        father = elite_parents[random.randint(0, len(elite_parents)-1)]   
        child = reproduce(mother, father)

        # Single-point Mutate based on the a small probability
        if random.randint(1, 100) <= single_mutation_threshold:
            mutate(child)
        # Multi-point Mutate based on the a small probability
        if random.randint(1, 100) <= multi_mutation_threshold:
            # Random number of points
            num_points = random.randint(1, 4)
            multi_mutate(child, num_points)

        # Single-section Crossover on two individuals in this new 
        # generation based on the a small probability
        if random.randint(1, 100) <= single_crossover_threshold and len(new_generation) >= 2:
            
            individual_1 = new_generation[random.randint(0, len(new_generation)-1)]
            individual_2 = new_generation[random.randint(0, len(new_generation)-1)]
            crossover(individual_1, individual_2)
        # Multi-section Crossover on two individuals in this new 
        # generation based on the a small probability
        if random.randint(1, 100) <= multi_crossover_threshold and len(new_generation) >= 2:
            
            individual_1 = new_generation[random.randint(0, len(new_generation)-1)]
            individual_2 = new_generation[random.randint(0, len(new_generation)-1)]
            # Random number of sections
            num_sections = random.randint(1, 4)
            multi_crossover(individual_1, individual_2, num_sections)


        # Add the child to the new generation
        new_generation.append(child)

    return new_generation


# Returns a string with that will be used as a print statement
# to print out all the boxes that are inside the most valuable
# backpacks
def getBestSolution(boxes, permutations):
    string = "Using our Genetic Algorithm, after " + str(MAX_NUMBER_OF_GENERATION)
    string = string + " generations,\nwe have found that the most valuable" 
    string = string + " backpacks contain:\n\n"
    for backpack in permutations:
        index = 0
        for box in boxes:
            if backpack[index] == 1:
                string = string + "Box" + str(index + 1) + ", "
            index = index + 1

        string = string[:-2]
        string = string + "\nWith a fitness value of: "
        string = string + str(fitness(backpack)) + "\n\n"

    return string



# Runs the Genetic Algorithm
def GeneticAlgorithm_Backpack():
    generation = 1
    population = InitializePopulation(INITIAL_POPULATION_SIZE)
    best_individual = []
    
    # Runs through each generation until we reached the maximum 
    # number of generations allowed
    for gen in range(0, MAX_NUMBER_OF_GENERATION):

        print ("GENERATION %d" % generation)
        
        # Declare variables for the total and best fitness of the generation
        total_fitness = 0
        best_fitness = MINIMUM_FITNESS
        
        # Runs through each individual in the current population
        for individual in population:
            # If the fitness is optimum and we don't already have this
            # solution, then add it to the list of possible solutions
            if fitness(individual) == OPTIMUM and individual not in best_individual:
                best_individual.append(individual)

            # If we found an individual in the generation with a better fitness 
            if fitness(individual) > best_fitness:
                best_fitness = fitness(individual)

            # Add the fitness of the current individual to the total
            total_fitness = total_fitness + fitness(individual)
            

        print ("Average fitness = %d" % (total_fitness/len(population)))
        if gen not in range(0, 5) and gen != (MAX_NUMBER_OF_GENERATION - 1):
            print()
        
        # Print the best fitness for only the first 5 generations and
        # the last generation because these generation will have the
        # most chance for vary in best fitness.
        if gen in range(0, 5) or gen == (MAX_NUMBER_OF_GENERATION - 1):
            print ("Best fitness in this generation: %d\n" % (best_fitness))     
        
        # Evolve the population to get the next generation that
        # is hopefully more fit overall than the current generation
        population = evolvePopulation(population)
        generation = generation + 1


    # Prints the most fit indiviudals that we found in the generations
    # These individuals represent the most valuable backpacks possible
    print()
    print(getBestSolution(BOXES, best_individual))


# Main function
if __name__ == "__main__":
    print ("The population of each generation is %d individuals\n" % (INITIAL_POPULATION_SIZE))
    GeneticAlgorithm_Backpack()