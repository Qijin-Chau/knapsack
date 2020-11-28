README

Artificial Intelligence

1. Define the problem as a genetic algorithm
       
       The individuals/people of a population will be represented as
       the possible backpacks that can be filled using the given boxes.
       
       These individuals will have a discrete representation as
       their genetic encoding. See full description in question 2.
       
       A fitness level will be given to each of these individuals
       and the generation will evolve with the fittest individuals
       from the previous children along with children of these
       individuals.
       
       In my genetic algorithm, I use Elitism and keep the top half 
       of the old generation and have that generation reproduce and
       regenerate the rest of the population. This as a whole is the
       culling technique of 50% called for in the spec.
       
       Thus 50% of the new generation will be the elite half from
       the old generation and the remaining 50% will be the children
       of these elite individuals. 
 

 2. Provide the genome for the problem
 
       As stated above, individuals will have a discrete representation
       as their genetic encoding. This discrete representation will be
       a 12-element binary array because there are also 12 possible
       boxes that can go inside each backpack. 
       
       Each index of the array will corresspond to the BOXES array 
       (defined after the box class definition) which will hold all the boxes 
       that can be in the backpacks along with their value and weight.
       
       Thus, for example, if the genetic encoding has a 1 in the first
       index, then the individual or the backpack will contain Box 1.
       Having these boxes will represent the phenotypes resulting from
       having the gene to contain these boxes. 
 
 
 3. Define all the fringe operations

       Fringe operations used in this algorithm are single-point/multi-point 
       mutations and one-section/multi-section crossovers. 
           
           - A single-point mutation will randomly modify a single gene of a 
             given individual. 
           - A multi-point mutation will randomly modify 1-4 genes of a 
             given individual. 

           - A one-section crossover will take two individuals and randomly 
             select a point on each of the two individuals and swapping the 
             genetic material around this point. 
           - A multi-section crossover will take two individuals and randomly
             select 1-4 sections on each of the two individuals and swapping 
             the genetic material around this point.
       
       My genetic algorithm uses a 8% probability for single point/section 
       fringe operations and 3% probability for multi point/section fringe 
       operations to occur on a child that was produced.
 

 4. Cull your population by 50% at every generation

       A variation of Rank-based selection that consists of taking only 
       the top 50% of individuals in the ranked list of the current 
       population and reproducing with these elite 50% to regenerate the 
       rest of the population.
