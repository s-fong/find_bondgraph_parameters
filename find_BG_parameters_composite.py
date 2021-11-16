# find bond-graph parameters for a system with multiple modules

# require separate folders within this directory containing each module's kinetic_parameters.py file and data files

# write out cellml file in text form

# prints out error between given kinetic parameters, and parameters found back-transforming the bond-graph parameters

import os
import sys
import importlib
import json
import csv
import math
import numpy as np
from scipy.linalg import null_space
import sympy
from sympy import Matrix, S, nsimplify
from fractions import Fraction

def read_IDs(path):
    data = []
    with open(path,'r') as f:
        reader = csv.reader(f)
        for row in reader:
            data.append(row[0])
        f.close()
    return data

def load_matrix(stoich_path):
    matrix = []
    with open(stoich_path,'r') as f:
        reader = csv.reader(f,delimiter=',')
        for row in reader:
            matrix.append([int(r) for r in row])
        f.close()
    return matrix

# def rational_nullspace(A, max_denom = 10):
    # v = null_space(A)
    # vFrac = [[Fraction(num).limit_denominator(max_denominator=max_denom) for num in row] for row in v]

    # vRat = [] #np.zeros([len(vFrac),len(vFrac[0])])
    # if not v.any():
        # return []
    # for row in vFrac:
        # largest_denom = max([res.denominator for res in row])
        # vRat.append( [vi.numerator for vi in row] )
    # return vRat

def calcT(I_vec,num_rows):
    num_cols = len(I_vec)
    T = np.zeros([num_rows,num_cols])
    for i in range(num_cols):
        T[I_vec[i]][i] = 1
    
    return T

if __name__ == "__main__":

    # Set directories
    current_dir = os.getcwd()
    main_dir = os.path.dirname(current_dir)
    output_dir = current_dir + '\output'
    whole_name = main_dir.split('\\')[-1]
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    ## Define volumes (unit pL)
    V_myo = 34.4
    V_e = 5.182     # external volume
    V_SR = V_myo*0.035 # SR volume
    V_di = V_myo*0.0539 # diadic volume
    V = dict()
    V['V_myo'] = V_myo
    V['V_SR'] = V_SR
    V['V_di'] = V_di
    V['V_e'] = V_e

    ## Load stoichiometric matrices and kinetic rate constants
    subsystem_names = ['individual'] #['cAMP', 'LRGbinding_B1AR', 'B1AR', 'PKA', 'PLB', 'Inhib1', 'GsProtein']
    num_subsystems = len(subsystem_names)
    sys_struct = {c:{} for c in subsystem_names}
    rxnIDs = []
    Knames = []
    Kname_modules = dict()
    for i_system in range(num_subsystems):
        sys_name = subsystem_names[i_system]
        sys_dir = current_dir + '\\' + sys_name +'\\'
        os.chdir(sys_dir)
        forward_mat_path = 'data\\all_forward_matrix.txt'
        reverse_mat_path = 'data\\all_reverse_matrix.txt'
        N_f = load_matrix(forward_mat_path)
        N_r = load_matrix(reverse_mat_path)
        sys_struct[sys_name]['N_f'] = N_f
        sys_struct[sys_name]['N_r'] = N_r
        
    #     print(subsystem_names[i_system])
        dims = dict()
        dims['num_rows'] = len(N_f)
        dims['num_cols'] = len(N_f[0])
        I = np.identity(dims['num_cols'])
        M = np.append(np.append(I, np.transpose(N_f),1), np.append(I, np.transpose(N_r),1),0)
        sys.path.append(sys_dir)
        globals()['kp_' + sys_name] = importlib.import_module('kinetic_parameters_' + sys_name)
        [k_kinetic, N_cT, K_C, W] = globals()['kp_' + sys_name].kinetic_parameters(dims, V)
        sys_struct[sys_name]['kfkr'] = k_kinetic
        sys_struct[sys_name]['Kc'] = K_C
        sys_struct[sys_name]['N_cT'] = N_cT
    #         sys_struct[i_system].W = W(dims.num_cols+1:end)
        rxnID = read_IDs('data\\rxnID.txt')
        rxnIDs.extend(rxnID)
        sys_struct[sys_name]['rxnID'] = rxnID
        Kname = read_IDs('data\\Kname.txt')
        Knames.extend(Kname)
        sys_struct[sys_name]['Kname'] = Kname

    Kunique = []
    for ik in Knames:
        # if ~any(strcmp(Kunique,ik)):
        if ik not in Kunique:
            Kunique.append(ik)

    os.chdir(current_dir)
    # relations between submodule to whole module

    for name in subsystem_names:
        ids = [Kunique.index(kid) for kid in sys_struct[name]['Kname']]
        sys_struct[name]['I_vec'] = ids
    num_rows = max(sys_struct[subsystem_names[-1]]['I_vec'])+1

    N_f = []
    N_r = []

    for sys_name in subsystem_names:
        # print(sys_name)
        T = calcT(sys_struct[sys_name]['I_vec'],num_rows)
        sys_struct[sys_name]['T'] = T

        new_f = np.matmul(T,sys_struct[sys_name]['N_f'])
        new_r = np.matmul(T,sys_struct[sys_name]['N_r'])

        if not len(N_f):
            N_f = new_f
            N_r = new_r
        else:
            N_f = np.append(N_f, new_f,1)
            N_r = np.append(N_r, new_r,1)

    N_fT = np.transpose(N_f)
    N_rT = np.transpose(N_r)

    N = N_r - N_f
    N_T = N_rT - N_fT

    num_cols = len(N[0])
    I = np.identity(num_cols)

    M = np.append(np.append(I, N_fT,1), np.append(I, N_rT,1),0)
    M_rref = sympy.Matrix(M).rref()

    ## Set up the vectors for kinetic rate constants

    kf = []
    kr = []
    for sys_name in subsystem_names:
        nrx = int(len(sys_struct[sys_name]['kfkr'])/2)
        kf.extend(sys_struct[sys_name]['kfkr'][:nrx])
        kr.extend(sys_struct[sys_name]['kfkr'][nrx:])
    k_kinetic = kf +kr

    W = list(np.append([1]*len(N[0]), [V_myo]*num_rows))

    lambda_expo = np.matmul(np.linalg.pinv(M), [math.log(k) for k in k_kinetic])
    lambdaW = [math.exp(l) for l in lambda_expo]
    lambdak = [lambdaW[i]/W[i] for i in range(len(W))]
    kappa = lambdak[:len(N[0])]
    K = lambdak[len(N[0]):]

    file = open(output_dir + '/all_parameters_out.json', 'w')
    data = { "K": K, "kappa": kappa, "k_kinetic": k_kinetic }
    json.dump(data, file)

    # Checks
    N_rref = sympy.Matrix(N).rref()
    R = nsimplify(Matrix(N), rational=True).nullspace() #rational_nullspace(N, max_denom=len(N[0]))
    if R:
        R = np.transpose(np.array(R).astype(np.float64))[0]
    # Check that there is a detailed balance constraint
    Z = nsimplify(Matrix(M), rational=True).nullspace() #rational_nullspace(M, 2)
    if Z:
        Z = np.transpose(np.array(Z).astype(np.float64))[0]

    k_est = np.matmul(M,[math.log(k) for k in lambdaW])
    k_est = [math.exp(k) for k in k_est]
    diff = [(k_kinetic[i] - k_est[i])/k_kinetic[i] for i in range(len(k_kinetic))]
    error = np.sum([abs(d) for d in diff])

    K_eq = [kf[i]/kr[i] for i in range(len(kr))]
    
    try:
        zero_est = np.matmul(np.transpose(R),K_eq)
        zero_est_log = np.matmul(np.transpose(R),[math.log(k) for k in K_eq])
    except:
        print('undefined R nullspace')


    # ### print outputs ###
    for ik in range(len(kappa)):
        print('var kappa_%s: fmol_per_sec {init: %g, pub: out};' %(rxnIDs[ik],kappa[ik]))
    for ik in range(len(Kunique)):
        print('var K_%s: per_fmol {init: %g, pub: out};' %(Kunique[ik],K[ik]))

    print('error = ', error)

    # initialise struct for storing modules contributing to a given K
    for ik in range(len(Kunique)):
        Kname_modules[Kunique[ik]] = []

    for sys_name in subsystem_names:
        modKname = sys_struct[sys_name]['Kname']
        for ik in range(len(modKname)):
            Kname_modules[modKname[ik]].append(sys_name)

    # write out CellML code
    if True:
        cellmlfilepath = output_dir + '\\TEMP.cellml.txt'
        with open(cellmlfilepath, 'w') as cid:

            cid.write('def model %s as\n def import using "units_and_constants/units_BG.cellml" for\n\
            unit mM using unit mM;\nunit fmol using unit fmol;\nunit per_fmol using unit per_fmol;\n\
            unit J_per_mol using unit J_per_mol;\nunit fmol_per_sec using unit fmol_per_sec;\n\
            unit C_per_mol using unit C_per_mol;\n  unit J_per_C using unit J_per_C;\n\
            unit microm3 using unit microm3;\n  unit fF using unit fF;\n\
            unit fC using unit fC;\n  unit fA using unit fA;\n\
            unit per_second using unit per_second;\n  unit millivolt using unit millivolt;\n\
            unit per_sec using unit per_sec;\n  unit J_per_K_per_mol using unit J_per_K_per_mol;\n\
            unit fmol_per_L using unit fmol_per_L;\n  unit fmol_per_L_per_sec using unit fmol_per_L_per_sec;\n\
            unit per_sec_per_fmol_per_L using unit per_sec_per_fmol_per_L;\n  unit uM using unit uM;\n\
            unit mM_per_sec using unit mM_per_sec;\n  unit uM_per_sec using unit uM_per_sec;\n\
            unit pL using unit pL;\n  unit m_to_u using unit m_to_u;\n enddef;\n' %(whole_name))
            cid.write('def import using "units_and_constants/constants_BG.cellml" for\n\
                comp constants using comp constants;\nenddef;\n\n')
            for module in subsystem_names:
                cid.write('def import using "%s/BG_%s.cellml" for\ncomp %s using comp %s;\nenddef;\n' % (
                module, module, module, module))
            cid.write('\ndef comp BG_parameters as\n')
            for ik in range(len(kappa)):
                cid.write('var kappa_%s: fmol_per_sec {init: %g, pub: out};\n' % (rxnIDs[ik], kappa[ik]))
            for ik in range(len(Kunique)):
                cid.write('var K_%s: per_fmol {init: %g, pub: out};\n' % (Kunique[ik], K[ik]))
            cid.write('enddef;\n')
            cid.write('    def comp environment as\n\
                var time: second {pub: out};\n\
                var vol_myo: pL {init: 34.4, pub: out};\n\
                var freq: dimensionless {init: 500};\n\
                // stimulus\n\
                // ramp UP and ramp DOWN\n\
                var stimSt: second {init: 3.5e-4};\n\
                var stimDur: second {init: 0.25e-4};\n\
                var tRamp: second {init: 1.8e-4};\n\
                var stimMag: fmol {init: 1e1};\n\
                var stimHolding: fmol {init: 1e-5};  \n\
                var m: fmol_per_sec;  \n\
                m = stimMag/tRamp;   \n\
                q_L_B1_init = sel          \n\
                    case (time < stimSt) and (time > stimSt-tRamp):   \n\
                        stimHolding+m*(time-stimSt+tRamp);\n\
                    case (time >= stimSt) and (time < stimSt+stimDur): \n\
                        stimMag+stimHolding;   \n\
                    case (time < stimSt+tRamp+stimDur) and (time >= stimSt+stimDur):  \n\
                        stimHolding+-m*(time-stimSt-tRamp-stimDur);   \n\
                    otherwise:             \n\
                        stimHolding;       \n\
                endsel;\n')
            for j in range(len(K)):
                cid.write('var q_%s_init: fmol {init: 1e-888};\n' % Kunique[j])

            cid.write('\n// mass conservation checks\n')
            cid.write(' var L_B1_T: fmol;\n\
            var R_B1_T: fmol;\n var Gs_T: fmol;\n var adenosine_T: fmol;\n')
            cid.write('        L_B1_T = q_L_B1+q_LR_B1Gs+q_LR_B1+q_LR_B1_aby+q_LR_B1_aby_T;\n\
            R_B1_T = q_R_B1+q_R_B1Gs+q_LR_B1+q_LR_B1Gs+q_R_B1_aby+q_R_B1_aby_T+q_LR_B1_aby+q_LR_B1_aby_T;\n\
            Gs_T = q_Gs+q_R_B1Gs+q_LR_B1Gs+q_a_Gs_GTP+q_a_Gs_GDP+q_R_B1_aby+q_R_B1_aby_T+q_LR_B1_aby+q_LR_B1_aby_T;\n\
            adenosine_T = q_cAMP+q_PDE_cAMP+q_five_AMP+q_ATP+q_AC_ATP+q_a_Gs_GTP_AC_ATP+q_FSK_AC_ATP;\n')
            cid.write('\n// Global value\n')
            for j in range(len(K)):
                cid.write('var q_%s: fmol {pub: out};\n' % Kunique[j])
            for module in subsystem_names:
                modKname = sys_struct[module]['Kname']
                cid.write('\n// %s imports\n' % module)
                for j in modKname:
                    cid.write('var q_%s_m%s: fmol {pub: in};\n' % (j, module))
                cid.write('\n')
            cid.write('\n')
            for kun in Kunique:
                cid.write('q_%s = q_%s_init' % (kun, kun))
                for mod in Kname_modules[kun]:
                    cid.write(' + q_%s_m%s ' % (kun, mod))
                cid.write(';\n')
            cid.write('enddef;\n')

            cid.write('\n')
            for module in subsystem_names:
                modKname = sys_struct[module]['Kname']
                cid.write('def map between environment and %s for\n' % module)
                cid.write('vars time and time;\n')
                for mod in modKname:
                    cid.write('vars q_%s_m%s and q_%s;\n' % (mod, module, mod))
                    cid.write('vars q_%s and q_%s_global;\n' % (mod, mod))
                cid.write('enddef;\n\n')

            for module in subsystem_names:
                modKname = sys_struct[module]['Kname']
                modrxnID = sys_struct[module]['rxnID']
                cid.write('def map between BG_parameters and %s for\n' % (module))
                for ik in modrxnID:
                    cid.write('vars kappa_%s and kappa_%s;\n' % (ik, ik))
                for mod in modKname:
                    cid.write('vars K_%s and K_%s;\n' % (mod, mod))
                cid.write('enddef;\n')
            cid.write('\n')
            for module in subsystem_names:
                cid.write('def map between constants and %s for\n' % (module))
                cid.write('\tvars R and R;\n\tvars T and T;\nenddef;\n')

            cid.write('\nenddef;\n')
        cid.close()
