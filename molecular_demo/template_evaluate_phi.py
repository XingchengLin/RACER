####################################################################################
# This script will generate decoy sequences based on your given native sequence
#
# Written by Xingcheng Lin, 06/10/2020
####################################################################################

import math
import subprocess
import os
import time
import sys

import numpy as np


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
                    if (TCR_name == 'TCR_NAME'):                       
                        if (res2chain == 'TCR_ALPHACHAIN' or res2chain == 'TCR_BETACHAIN'):
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


def evaluate_phis_over_training_set(training_set_file, phi_list_file_name, decoy_method, max_decoys, tm_only=False, num_processors=1, TCRmodeling=False, TCR_name='IDK'):
    phi_list = read_phi_list(phi_list_file_name)
    print(phi_list)
    training_set = read_column_from_file(training_set_file, 1)
    print(training_set)

    # for protein in training_set:
    evaluate_phis_for_protein(training_set, phi_list, decoy_method, max_decoys, tm_only=tm_only, TCRmodeling=TCRmodeling, TCR_name=TCR_name)


def evaluate_phis_for_protein(training_set, phi_list, decoy_method, max_decoys, tm_only=False, TCRmodeling=False, TCR_name='IDK'):
    # Because there is only one protein in the training set; if there are multiple proteins, the script could be different!
    protein = training_set[0]

    print(native_structures_directory)
    structure = parse_pdb(os.path.join(native_structures_directory, protein))

    # Two lists of res_list, one for the peptide (selected by the .tm file), one for the entire list
    res_list_tmonly = get_res_list(structure, tm_only=True)
    res_list_entire = get_res_list(structure, tm_only=False)
    # Here, we are going to take every residues close to the pMHC peptide, so there is no restriction (tm_only) on what is going to be taken;
    neighbor_list = get_neighbor_list(structure, tm_only=False)

    sequence = get_sequence_from_structure(structure)

    # Before the iteration, we need to store the native res_list_tmonly and res_list_entire; because after mutation for each phi, both of them 
    # are changed afterwards to the decoy sequences; but we still need to evaluate the native phi for different phis;
    res_list_tmonly_native = res_list_tmonly
    res_list_entire_native = res_list_entire

    for phi, parameters in phi_list:

        phi = globals()[phi]
        parameters_string = get_parameters_string(parameters)
        # check to see if the decoys are already generated
#        number_of_lines_in_file = get_number_of_lines_in_file(os.path.join(
#            phis_directory, "%s_%s_native_%s" % (phi.__name__, protein, parameters_string)))
#        if not number_of_lines_in_file >= 1:
        output_file = open(os.path.join(phis_directory, "%s_%s_native_%s" % (
            phi.__name__, protein, parameters_string)), 'w')
        phis_to_write = phi(res_list_tmonly_native, res_list_entire_native,
                            neighbor_list, parameters, TCRmodeling=TCRmodeling, TCR_name=TCR_name)
        output_file.write(str(phis_to_write).strip(
            '[]').replace(',', '') + '\n')
        output_file.close()
#        number_of_lines_in_file = get_number_of_lines_in_file(os.path.join(
#            phis_directory, "%s_%s_decoys_%s_%s" % (phi.__name__, protein, decoy_method, parameters_string)))
#        if not number_of_lines_in_file >= max_decoys:
        output_file = open(os.path.join(phis_directory, "%s_%s_decoys_%s_%s" % (
            phi.__name__, protein, decoy_method, parameters_string)), 'w')
        decoy_sequences = read_decoy_sequences(os.path.join(
            decoys_root_directory, "%s/%s.decoys" % (decoy_method, protein)))
        for i_decoy, decoy_sequence in enumerate(decoy_sequences):
            if i_decoy >= max_decoys:
                break
            mutate_whole_sequence(res_list_entire, decoy_sequence)
            # Note after this mutation, both res_list_entire and res_list_tmonly have been changed accordingly, because this is a change by reference;

            phis_to_write = phi(res_list_tmonly, res_list_entire,
                                neighbor_list, parameters, TCRmodeling=TCRmodeling, TCR_name=TCR_name)
            output_file.write(str(phis_to_write).strip(
                '[]').replace(',', ' ') + '\n')
        output_file.close()


############################################

native_structures_directory = "./native_structures_pdbs_with_virtual_cbs/"
phis_directory = "./phis/"
decoys_root_directory = "./sequences/"

evaluate_phis_over_training_set("proteins_list.txt", "phi1_list.txt", decoy_method='TCR_randomization', 
                                max_decoys=1000, tm_only=False, num_processors=1, TCRmodeling=True, TCR_name='TCR_NAME')
