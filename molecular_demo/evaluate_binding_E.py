####################################################################################
# This script will help add fake CB atoms for Glycines
#
# Written by Xingcheng Lin, 06/10/2020
####################################################################################

import math
import subprocess
import os
import time
import sys


import numpy as np
import matplotlib.pyplot as plt

from common_function import *

################################################


def my_lt_range(start, end, step):
    while start < end:
        yield start
        start += step


def my_le_range(start, end, step):
    while start <= end:
        yield start
        start += step


def RepresentsFloat(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
###########################################


def phi_pairwise_contact_well(res_list_tmonly, res_list_entire, neighbor_list, parameter_list, TCRmodeling=False, TCR_name='IDK'):

    r_min, r_max, kappa, min_seq_sep = parameter_list
    r_min = float(r_min)
    r_max = float(r_max)
    kappa = float(kappa)
    min_seq_sep = int(min_seq_sep)
    phi_pairwise_contact_well = np.zeros((20, 20))
    for res1globalindex, res1 in enumerate(res_list_entire):

        res1index = get_local_index(res1)
        res1chain = get_chain(res1)

        # For TCR modeling, we only need the sequence in the peptide;
        if TCRmodeling:


            if (res1 in res_list_tmonly):
                for res2 in get_neighbors_within_radius(neighbor_list, res1, r_max + 2.0):
                    res2index = get_local_index(res2)
                    res2chain = get_chain(res2)
                    res2globalindex = get_global_index(res_list_entire, res2)
                    # Here, we strictly consider only between the peptide and the TCR (chain C & D):
                    # Res1 through tm_only, is already in chain E, we only need to control the res2 to
                    # be in chian C or D;
                    # The chain ID varies from one TCR to the other, be careful!!!
                    # 2B4: chain C or D; 5CC7: chain D or E; 226: chain C or D;
                    if (TCR_name == '3qib'):                       
                        if (res2chain == 'C' or res2chain == 'D'):
                            res1type = get_res_type(res_list_entire, res1)
                            res2type = get_res_type(res_list_entire, res2)
                            rij = get_interaction_distance(res1, res2)
                            phi_pairwise_contact_well[res1type][res2type] += interaction_well(
                                rij, r_min, r_max, kappa)
                            if not res1type == res2type:
                                phi_pairwise_contact_well[res2type][res1type] += interaction_well(
                                    rij, r_min, r_max, kappa)

            else:
                continue

        else:

            for res2 in get_neighbors_within_radius(neighbor_list, res1, r_max + 2.0):
                res2index = get_local_index(res2)
                res2chain = get_chain(res2)
                res2globalindex = get_global_index(res_list_entire, res2)
                if (res1chain == res2chain and res2index - res1index >= min_seq_sep) or (res1chain != res2chain and res2globalindex > res1globalindex):
                    res1type = get_res_type(res_list_entire, res1)
                    res2type = get_res_type(res_list_entire, res2)
                    rij = get_interaction_distance(res1, res2)
                    phi_pairwise_contact_well[res1type][res2type] += interaction_well(
                        rij, r_min, r_max, kappa)
                    if not res1type == res2type:
                        phi_pairwise_contact_well[res2type][res1type] += interaction_well(
                            rij, r_min, r_max, kappa)

    phis_to_return = []
    for i in range(20):
        for j in range(i, 20):
            phis_to_return.append(phi_pairwise_contact_well[i][j])

    return phis_to_return

def evaluate_phis_for_decoy_protein(protein, phi_list, decoy_method, max_decoys, tm_only=False, TCRmodeling=False, TCR_name='IDK'):
    print(protein)
    structure = parse_pdb(os.path.join(decoy_structures_directory, protein))

    # Two lists of res_list, one for the peptide (selected by the .tm file), one for the entire list
    res_list_tmonly = get_res_list(structure, tm_only=True)
    res_list_entire = get_res_list(structure, tm_only=False)
    # Here, we are going to take every residues close to the pMHC peptide, so there is no restriction (tm_only) on what is going to be taken;
    neighbor_list = get_neighbor_list(structure, tm_only=False)

    sequence = get_sequence_from_structure(structure)

    for phi, parameters in phi_list:
        phi = globals()[phi]
        parameters_string = get_parameters_string(parameters)
        # check to see if the decoys are already generated
#        number_of_lines_in_file = get_number_of_lines_in_file(os.path.join(
#            phis_directory, "%s_%s_decoy_%s" % (phi.__name__, protein, parameters_string)))
#        if not number_of_lines_in_file >= 1:
        output_file = open(os.path.join(phis_directory, "%s_%s_decoy_%s" % (
            phi.__name__, protein, parameters_string)), 'w')
        phis_to_write = phi(res_list_tmonly, res_list_entire,
                            neighbor_list, parameters, TCRmodeling=TCRmodeling, TCR_name=TCR_name)
        output_file.write(str(phis_to_write).strip(
            '[]').replace(',', '') + '\n')
        output_file.close()

def evaluate_phis_over_training_set_for_decoy_structures(decoy_training_set_file, phi_list_file_name, decoy_method, max_decoys, tm_only=False, num_processors=1, TCRmodeling=False, TCR_name='IDK'):
    phi_list = read_phi_list(phi_list_file_name)
    print(phi_list)
    training_set = read_column_from_file(decoy_training_set_file, 1)
    print(training_set)

    function_to_evaluate = functools.partial(
        evaluate_phis_for_decoy_protein, tm_only=tm_only, TCRmodeling=TCRmodeling, TCR_name=TCR_name)
    arguments_lists = [training_set, [phi_list] * len(training_set), [
        decoy_method] * len(training_set), [max_decoys] * len(training_set)]
    call_independent_functions_on_n_processors(
        function_to_evaluate, arguments_lists, num_processors)

def validate_hamiltonian_decoy_structures_provided(hamiltonian, native_training_set_file, decoy_training_set_file, training_decoy_method, native_test_set_file=None, decoy_test_set_file=None, test_decoy_method=None, use_filtered_gammas=False):
    if native_test_set_file == None:
        native_test_set_file = native_training_set_file
    if decoy_test_set_file == None:
        decoy_test_set_file = decoy_training_set_file
    if test_decoy_method == None:
        test_decoy_method = training_decoy_method
    native_test_set = read_column_from_file(native_test_set_file, 1)
    decoy_test_set = read_column_from_file(decoy_test_set_file, 1)
    z_scores = []
    e_natives = []
    e_mgs = []
    e_mg_stds = []
    for protein in native_test_set:
        en = evaluate_hamiltonian_native_structures_provided(
            protein, hamiltonian, native_test_set_file, training_decoy_method, test_decoy_method, use_filtered_gammas)
        e_natives.append(en)
    for protein in decoy_test_set:
        e_decoy = evaluate_hamiltonian_decoy_structures_provided(
            protein, hamiltonian, native_test_set_file, decoy_test_set_file, training_decoy_method, test_decoy_method, use_filtered_gammas)
        e_mgs.append(e_decoy)

    # Calculate the real averaged <E>mg, <E>n, std(Emg) and zScore;
    average_e_mg = np.average(e_mgs)
    std_e_mgs = np.std(e_mgs)
    average_e_native = np.average(e_natives)
    z_score_forall = (average_e_mg - average_e_native) / std_e_mgs
    return z_score_forall, average_e_native, average_e_mg, std_e_mgs, e_mgs, e_natives


def evaluate_hamiltonian_native_structures_provided(protein, hamiltonian, native_training_set_file, training_decoy_method, test_decoy_method, use_filtered_gammas=True):
    phi_list = read_phi_list(hamiltonian)
    native_training_set = read_column_from_file(native_training_set_file, 1)

    # read in Hamiltonian;
    # Find out how many total phi_i there are and get full parameter string
    total_phis, full_parameters_string, num_phis = get_total_phis_and_parameter_string(
        phi_list, native_training_set)

    # read in corresponding gammas
    if use_filtered_gammas:
        gamma_file_name = "%s%s_%s_gamma_filtered" % (
            gammas_directory, native_training_set_file.split('/')[-1].split('.')[0], full_parameters_string)
    else:
        gamma_file_name = "%s%s_%s_gamma" % (gammas_directory, native_training_set_file.split(
            '/')[-1].split('.')[0], full_parameters_string)

    # Need to filter out the complex number if in the "filtered" mode;
    if use_filtered_gammas:
        gamma = np.loadtxt(gamma_file_name, dtype=complex, converters={
                           0: lambda s: complex(s.decode().replace('+-', '-'))})
    else:
        gamma = np.loadtxt(gamma_file_name)

    # Read in corresponding phis (native);
    phi_native = read_native_phi(protein, phi_list, total_phis)
        
    # perform dot products to get native energies
    e_native = np.dot(gamma, phi_native)
    return e_native

def evaluate_hamiltonian_decoy_structures_provided(protein, hamiltonian, native_training_set_file, decoy_training_set_file, training_decoy_method, test_decoy_method, use_filtered_gammas=True):
    phi_list = read_phi_list(hamiltonian)
    decoy_training_set = read_column_from_file(decoy_training_set_file, 1)

    # read in Hamiltonian;
    # Find out how many total phi_i there are and get full parameter string
    total_phis, full_parameters_string, num_phis = get_total_phis_and_parameter_string_decoy_structures_provided(
        phi_list, decoy_training_set)    

    # read in corresponding gammas, note there is only one gamma for both the native and decoy sets;
    if use_filtered_gammas:
        gamma_file_name = "%s%s_%s_gamma_filtered" % (
            gammas_directory, native_training_set_file.split('/')[-1].split('.')[0], full_parameters_string)
    else:
        gamma_file_name = "%s%s_%s_gamma" % (gammas_directory, native_training_set_file.split(
            '/')[-1].split('.')[0], full_parameters_string)

    # Need to filter out the complex number if in the "filtered" mode;
    if use_filtered_gammas:
        gamma = np.loadtxt(gamma_file_name, dtype=complex, converters={
                           0: lambda s: complex(s.decode().replace('+-', '-'))})
    else:
        gamma = np.loadtxt(gamma_file_name)

    # Read in corresponding phis (decoy);
    phi_i_decoy = read_decoy_phi_structures_provided(protein, phi_list, total_phis)
    
    # perform dot products to get decoy energies 
    e_decoy = np.dot(gamma, phi_i_decoy)
    return e_decoy



def draw_plot(data, pos, edge_color, fill_color, labels):
    bp = ax.boxplot(data, positions=pos, patch_artist=True, showfliers=False, widths=0.75, showmeans=False, meanline=False, medianprops=dict(linewidth=2.5), labels=labels)


    for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
        plt.setp(bp[element], color=edge_color)

    for patch in bp['boxes']:
        patch.set(facecolor=fill_color)
    
    first_legend = ax.legend([bp["boxes"][0]], ['Weaker Binders'], loc='best')

###################################################################################################################################################

decoy_structures_directory = "./test_structures_pdbs_with_virtual_cbs"
gammas_directory = "./gammas/randomized_decoy/01022019/direct_contact/"


os.chdir(decoy_structures_directory)
add_virtual_glycines_list("testSetFiles.txt")
os.chdir('../')


evaluate_phis_over_training_set_for_decoy_structures("testSetFiles.txt", "phi1_list.txt", decoy_method='TCR_modeling', max_decoys=1e+5,
                                                     num_processors=1, TCRmodeling=True, TCR_name='3qib')


z_score, ave_e_native, ave_e_mgs, e_mg_std, e_mgs, e_natives = validate_hamiltonian_decoy_structures_provided(hamiltonian='phi1_list.txt', native_training_set_file='proteins_list.txt', decoy_training_set_file='testSetFiles.txt', training_decoy_method='TCR_modeling', decoy_test_set_file='testSetFiles.txt', use_filtered_gammas=True)

np.savetxt('evaluated_binding_E/epitopeE.txt',np.real(e_natives), fmt='%1.3f')
np.savetxt('evaluated_binding_E/non-epitopeE.txt',np.real(e_mgs), fmt='%1.3f')
print(z_score)

E_native_list = [np.real(e_natives)]
E_decoy_list = [np.real(e_mgs)]

# Draw the strong- and weak-binding energies
fig, ax = plt.subplots()

bp = ax.boxplot(E_decoy_list, positions=[0], patch_artist=True, showfliers=False, widths=0.75, showmeans=False, meanline=False, medianprops=dict(linewidth=2.5), labels=["Binding Energies"])


for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
    plt.setp(bp[element], color='red')

for patch in bp['boxes']:
    patch.set(facecolor='tan')
    
first_legend = ax.legend([bp["boxes"][0]], ['Weak Binder'], loc='lower left')

ax.scatter(0, E_native_list[0], s=100, marker='o', color='b', zorder=1000, label='Strong Binders')
second_legend = ax.legend(loc='lower right')
ax.add_artist(first_legend)
plt.ylabel('Binding Energies')

plt.savefig('./Evaluated_bindingE.png')

