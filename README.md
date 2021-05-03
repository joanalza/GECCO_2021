# TOWARDS THE LANDSCAPE ROTATION AS A PERTURBATION STRATEGY ON THE QUADRATIC ASSIGNMENT PROBLEM

## Joan Alza, Mark Bartlett, Josu Ceberio, John McCall

This repository contains the code and results of the paper "Towards the Landscape Rotation As a Perturbation Strategy on the Quadratic Assignment Problem", published on the workshop "Evolutionary Computation for Permutation Problems" on GECCO'21.

This work has been coded in Python, and it can be found on the './Code/' folder. It is structured as follows:

- **Algorithms**: this folder contains the stochastic hill-climbing algorithm (sHC.py), some functions about to work with permutation distances (Distances.py) and a file containing population based functions (Population.py).
- **Individual**: this folder contains some functions to work with permutations (Permutation.py). 
- **Problems**: it contains the QAP able to work with rotated landscapes (DQAP.py).

The code can bee easy launched by executing `python sHC.py` and the input parameters wanted. For extra help, type the following code on  './Code/' directory:
```bash
python src/Algorithms/sHC.py '?'
```
which will return the possible arguments that can be used and the variables that they may take. In short, these are the most interesting parameters that may be passed as input to the launcher:

- Instance: instance of the 
- Results: path and filename of the output file. It will contain information about the generation, best solution found, candidate solution, the cooling parameter,etc.
- Seed: used for testing purposes
- Algorithm: version of the sHC, "standard", "rotation1", "rotation2", "restart".
- Stop: maximum budget (total number of iterations).
- Metric: permutation distance metric used, "K" for Kendall's-t, "C" for Cayley and "H" for Hamming distance.
- Trials: number of iterations to compare before considering the algorithm trapped.
- Cooling: the scheme used for the landscape rotation. It may be an integer within the possible distances for the given metric or linear or exponential cooling schedules.
