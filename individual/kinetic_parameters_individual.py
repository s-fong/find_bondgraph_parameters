# Template of a specific module with known kinetic rate constants from literature
# This is called from find_BG_parameters.py
# example reaction       [A][B] <-> [C]
# Replace 'individual' in the filename with the name of the module you are solving for

# / INPUTS /
#   dims: dict containing number of rows and number of columns of stoichiometric matrix N
#   V: dict containing cell volume values

# / OUTPUTS /
#   k_kinetic: list of forward and reverse rate constants
#   N_cT: vector of any constraints
#   K_C: RHS corresponding to constraints vector
#   W: vector of volumes relevant to K number of species
 
import numpy as np 

def kinetic_parameters(dims, V):


    # Set the kinetic rate constants.
    # convert all units to be consistent with each other - uM or mM are fine, but theuy ALL have to be of the same order. 
    # if reactions are irreversible, made reversible by letting kr ~= 0 and kr > 0

    num_cols = dims['num_cols']
    num_rows = dims['num_rows']

    fastKineticConstant = 1e6
    smallReverse = 1e-2

    k1p = 16                 # 1/s
    k1m = smallReverse           # 1/s            
    # k2p = 16                 # 1/s
    # k2m = smallReverse           # 1/s            
    # k3p = fastKineticConstant
    # k3m = smallReverse

    # Any detailed balance is obeyed for closed loops of known equilibrium constant (specified in literature)
    if False:
        k2m = k1m * k2p / k1p
    

    k_kinetic = [
        k1p, # all forward rates
        k1m, # all reverse rates
        ]

    # CONSTRAINTS
    # each row of N_cT corresponds to a new constraint to the system, if known
    N_cT = []
    
    if False:           
        id = [1,2,3] # example vector of indices of species A,B,C in total stoichiometry 
        N_cT = np.zeros(num_rows+num_cols)
        N_cT[num_cols + id[2]] = 1   # C
        N_cT[num_cols + id[0]] = -1  # A
        N_cT[num_cols + id[2]] = -1  # B

    K_C = [1]

    # volume vector
    W = list(np.append([1] * num_cols, [V['V_myo']] * num_rows))

    return (k_kinetic, [N_cT], K_C, W) 











