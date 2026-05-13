#!/usr/bin/env python3

"""
EffectorGeneP: gene annotation in pathogen genomes

Copyright (c) 2026 Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.

This work is freely available for non-commercial scientific research, non-commercial education, 
or non-commercial research projects, under the CSIRO Software License Agreement (https://github.com/JanaSperschneider/EffectorGeneP/blob/main/LICENCE). 
"""
# -----------------------------------------------------------------------------------------------------------
import os
import sys
import re
import random
import functions
import subprocess as sub

import subprocess
import errno
import uuid
import shutil
import tempfile
# -----------------------------------------------------------------------------------------------------------    
def get_model_predictions(WEKA_INPUT_FILE, WEKA_PATH, RESULTS_PATH, MODELS, CLASSIFIERS):

    for model, CLASSIFIER in zip(MODELS, CLASSIFIERS):
        ParamList = ['java', '-cp', WEKA_PATH, CLASSIFIER, '-l', model, '-T', WEKA_INPUT_FILE, '-p', 'first-last']
        #print(ParamList)

        with open(RESULTS_PATH + 'Predictions.txt', 'wb') as out:
            try:
                Process = subprocess.Popen(ParamList, shell=False, stdout=out)
                sts = Process.wait()
                cstdout, cstderr = Process.communicate()

                if Process.returncode:
                    raise Exception("Calling WEKA returned %s"%Process.returncode)
                if cstdout:
                    pass
                elif cstderr:
                    sys.exit()
            except:
                e = sys.exc_info()[1]
                print("Error calling WEKA: %s" % e)
                sys.exit(1)

        file_input = RESULTS_PATH + 'Predictions.txt'

        predictions = parse_weka_output(file_input)

    return predictions
# -----------------------------------------------------------------------------------------------------------        
def parse_weka_output(file_input):
    """ Function: parse_weka_output()
        Purpose:  Given the WEKA output file and the query identifiers and sequences, 
                  parse the predicted class for each protein from the WEKA output. 
              
        Input:    WEKA output file and the query identifiers and sequences.                  
    
        Return:   Predictions. 
    """    
    predictions = []

    with open(file_input) as f:

        content = f.readlines()

        for line in content:
            if line.strip() and ('pos' in line or 'neg' in line):
                position = line.split()[0]
                prediction = line.split()[2]
                prob = float(line.split()[3])        
 
                # WEKA output counts from position 1, our identifiers are counted from zero
                short_ident = 'seq' + str(int(position) - 1)

                if ':neg' in prediction:                                
                    predictions.append((short_ident, 'Negative', prob))
                else:                                                              
                    predictions.append((short_ident, 'Positive', prob))

    return predictions
# -----------------------------------------------------------------------------------------------------------
def reverse_classification(classification):
    
    reverse_classification = []

    for (short_ident, classification, prob) in classification:
        if classification == 'Positive':
            reverse_classification.append( (short_ident, 'Negative', prob) )
        elif classification == 'Negative':
            reverse_classification.append( (short_ident, 'Positive', prob) )

    return reverse_classification
# -----------------------------------------------------------------------------------------------------------
def ensembl_voting(classification, PREDICTIONS):

    yes_prob_average, no_prob_average = 0.0, 0.0
    yes_prob, no_prob = {}, {}

    for (short_ident, prediction, prob) in classification:

        if prediction == 'Positive':
            if short_ident in yes_prob:
                yes_prob[short_ident] = yes_prob[short_ident] + [prob] 
                no_prob[short_ident] = no_prob[short_ident] + [1.0 - prob] 
            else:
                yes_prob[short_ident] = [prob]
                no_prob[short_ident] = [1.0 - prob]

        if prediction == 'Negative':
            if short_ident in yes_prob:
                no_prob[short_ident] = no_prob[short_ident] + [prob] 
                yes_prob[short_ident] = yes_prob[short_ident] + [1.0 - prob]             
            else:
                no_prob[short_ident] = [prob]
                yes_prob[short_ident] = [1.0 - prob]


    for (short_ident, probability_list_yes) in yes_prob.items():
        probability_list_no = no_prob[short_ident]

        # Soft voting: argmax of the sum of predicted probabilities
        yes_prob_average = round(sum(probability_list_yes)/float(len(probability_list_yes)),3)
        no_prob_average = round(sum(probability_list_no)/float(len(probability_list_no)),3)

        if short_ident in PREDICTIONS:
            PREDICTIONS[short_ident] = PREDICTIONS[short_ident] + [yes_prob_average, no_prob_average]
        else:
            PREDICTIONS[short_ident] = [yes_prob_average, no_prob_average]


    return PREDICTIONS    
# -----------------------------------------------------------------------------------------------------------
def run_EffectorGeneP(FASTA_FILE, RESULTS_PATH, PATH_TO_MODELS, WEKA_PATH):

    hexamer_frequencies_intergenic = functions.read_kmer_frequencies_genome(PATH_TO_MODELS + '/intergenic.hexamers.txt', 6)
    hexamer_frequencies_intron = functions.read_kmer_frequencies_genome(PATH_TO_MODELS + '/introns.hexamers.txt', 6)
    hexamer_frequencies_CDS_secreted = functions.read_kmer_frequencies_genome(PATH_TO_MODELS + '/CDS_secreted.hexamers.txt', 6)
    hexamer_frequencies_CDS_notsecreted = functions.read_kmer_frequencies_genome(PATH_TO_MODELS + '/CDS_notsecreted.hexamers.txt', 6)
    hexamer_frequencies_CDS_effectors = functions.read_kmer_frequencies_genome(PATH_TO_MODELS + '/CDS_secreted_effectors.hexamers.txt', 6)
    hexamer_frequencies_UTRs3 = functions.read_kmer_frequencies_genome(PATH_TO_MODELS + '/UTRs3.hexamers.txt', 6)
    hexamer_frequencies_UTRs5 = functions.read_kmer_frequencies_genome(PATH_TO_MODELS + '/UTRs5.hexamers.txt', 6)
    hexamer_frequencies_kozak = functions.read_kmer_frequencies_genome(PATH_TO_MODELS + '/kozak.hexamers.txt', 6)

    codon_frequencies_intergenic = functions.read_kmer_frequencies_genome(PATH_TO_MODELS + '/intergenic.trimers.txt', 3)
    codon_frequencies_intron = functions.read_kmer_frequencies_genome(PATH_TO_MODELS + '/introns.trimers.txt', 3)
    codon_frequencies_CDS_secreted = functions.read_kmer_frequencies_genome(PATH_TO_MODELS + '/CDS_secreted.trimers.txt', 3)
    codon_frequencies_CDS_notsecreted = functions.read_kmer_frequencies_genome(PATH_TO_MODELS + '/CDS_notsecreted.trimers.txt', 3)
    codon_frequencies_CDS_effectors = functions.read_kmer_frequencies_genome(PATH_TO_MODELS + '/CDS_secreted_effectors.trimers.txt', 3)
    codon_frequencies_UTRs3 = functions.read_kmer_frequencies_genome(PATH_TO_MODELS + '/UTRs3.trimers.txt', 3)
    codon_frequencies_UTRs5 = functions.read_kmer_frequencies_genome(PATH_TO_MODELS + '/UTRs5.trimers.txt', 3)
    codon_frequencies_kozak = functions.read_kmer_frequencies_genome(PATH_TO_MODELS + '/kozak.trimers.txt', 3)
    # -----------------------------------------------------------------------------------------------------------
    not_secreted_candidates_CDS = PATH_TO_MODELS +'/trainingset.not_secreted_TMHMMs.homology_reduced.CDS.fa'
    secreted_candidates_CDS = PATH_TO_MODELS +'/trainingset.secreted_noTMHMMs.homology_reduced.CDS.fa'
    secreted_effector_candidates_CDS = PATH_TO_MODELS +'/trainingset.secreted_noTMHMMs.effectors.homology_reduced.CDS.fa'

    CodonsDict_Usage_NotSecreted = functions.read_codon_usage_dic(PATH_TO_MODELS + '/trainingset.not_secreted_CODON_USAGE.txt')
    CodonsDict_Usage_Secreted = functions.read_codon_usage_dic(PATH_TO_MODELS + '/trainingset.secreted_CODON_USAGE.txt')
    CodonsDict_Usage_Effectors = functions.read_codon_usage_dic(PATH_TO_MODELS + '/trainingset.effectors_CODON_USAGE.txt')

    # -----------------------------------------------------------------------------------------------------------
    # -----------------------------------------------------------------------------------------------------------
    # Check if FASTA file exists
    try:
        open(FASTA_FILE, 'r') 
    except OSError as e:
        print("Unable to open FASTA file:", FASTA_FILE)  #Does not exist OR no read permissions
        print("I/O error({0}): {1}".format(e.errno, e.strerror))
        sys.exit(1)
    # -----------------------------------------------------------------------------------------------------------
    # Extract the identifiers and sequences from input FASTA file
    ORIGINAL_IDENTIFIERS, SEQUENCES = [], []
    for identifier, sequence in functions.SimpleFastaParser(open(FASTA_FILE, 'r')):
      ORIGINAL_IDENTIFIERS.append(identifier)
      SEQUENCES.append(sequence.replace('*',''))

    SEQUENCES = [seq.upper() for seq in SEQUENCES]
    # -----------------------------------------------------------------------------------------------------------
    SHORT_IDENTIFIERS = ['seq' + str(index) for index, identifier in enumerate(ORIGINAL_IDENTIFIERS)]
    SHORT_IDENTIFIER_DIC = {}
    for short_ident, long_ident in zip(SHORT_IDENTIFIERS, ORIGINAL_IDENTIFIERS):
        SHORT_IDENTIFIER_DIC[short_ident] = long_ident
    # -----------------------------------------------------------------------------------------------------------
    # Write the WEKA arff file for classification of the input FASTA file
    weka_input = RESULTS_PATH + 'weka.arff'    
    functions.write_weka_input(weka_input, SEQUENCES, hexamer_frequencies_intergenic, hexamer_frequencies_intron, hexamer_frequencies_CDS_secreted, hexamer_frequencies_CDS_notsecreted, hexamer_frequencies_CDS_effectors, hexamer_frequencies_UTRs3, hexamer_frequencies_UTRs5, codon_frequencies_intergenic, codon_frequencies_intron, codon_frequencies_CDS_secreted, codon_frequencies_CDS_notsecreted, codon_frequencies_CDS_effectors, codon_frequencies_UTRs3, codon_frequencies_UTRs5, CodonsDict_Usage_NotSecreted, CodonsDict_Usage_Secreted, CodonsDict_Usage_Effectors)
    # -----------------------------------------------------------------------------------------------------------
    weka_input_kozak = RESULTS_PATH + 'weka_kozak.arff'    
    functions.write_weka_input_kozak(weka_input_kozak, SEQUENCES)
    # -----------------------------------------------------------------------------------------------------------
    ALGORITHM_TREE = 'weka.classifiers.trees.J48'
    ALGORITHM_LOGISTIC = 'weka.classifiers.functions.Logistic'
    # -----------------------------------------------------------------------------------------------------------
    PREFIX_TREE = '/TrainingData_J48_'
    PREFIX_LOGISTIC = '/TrainingData_Logistic_'
    PREDICTIONS = {}   
    # -----------------------------------------------------------------------------------------------------------
    ### Random sequence predictions
    classification_Random_ORFs_TREE = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_TREE + 'Random_ORFs.MODEL'], [ALGORITHM_TREE])
    classification_Random_Intron_TREE = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_TREE + 'Random_Intron.MODEL'], [ALGORITHM_TREE])
    classification_Random_Intergenic_TREE = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_TREE + 'Random_Intergenic.MODEL'], [ALGORITHM_TREE])
    classification_Random_NotSecreted_TREE = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_TREE + 'Random_NotSecreted.MODEL'], [ALGORITHM_TREE])
    classification_Random_Secreted_TREE = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_TREE + 'Random_Secreted.MODEL'], [ALGORITHM_TREE])
    classification_Random_SecretedEffectors_TREE = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_TREE + 'Random_SecretedEffectors.MODEL'], [ALGORITHM_TREE])
    classification_Random_UTRs3_TREE = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_TREE + 'Random_UTRs3.MODEL'], [ALGORITHM_TREE])
    classification_Random_UTRs5_TREE = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_TREE + 'Random_UTRs5.MODEL'], [ALGORITHM_TREE])

    classification_Random_ORFs_LOGISTIC = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_LOGISTIC + 'Random_ORFs.MODEL'], [ALGORITHM_LOGISTIC])
    classification_Random_Intron_LOGISTIC = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_LOGISTIC + 'Random_Intron.MODEL'], [ALGORITHM_LOGISTIC])
    classification_Random_Intergenic_LOGISTIC = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_LOGISTIC + 'Random_Intergenic.MODEL'], [ALGORITHM_LOGISTIC])
    classification_Random_NotSecreted_LOGISTIC = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_LOGISTIC + 'Random_NotSecreted.MODEL'], [ALGORITHM_LOGISTIC])
    classification_Random_Secreted_LOGISTIC = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_LOGISTIC + 'Random_Secreted.MODEL'], [ALGORITHM_LOGISTIC])
    classification_Random_SecretedEffectors_LOGISTIC = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_LOGISTIC + 'Random_SecretedEffectors.MODEL'], [ALGORITHM_LOGISTIC])
    classification_Random_UTRs3_LOGISTIC = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_LOGISTIC + 'Random_UTRs3.MODEL'], [ALGORITHM_LOGISTIC])
    classification_Random_UTRs5_LOGISTIC = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_LOGISTIC + 'Random_UTRs5.MODEL'], [ALGORITHM_LOGISTIC])

    classification_ORFs_Random_TREE = reverse_classification(classification_Random_ORFs_TREE)
    classification_Introns_Random_TREE = reverse_classification(classification_Random_Intron_TREE)
    classification_Intergenic_Random_TREE = reverse_classification(classification_Random_Intergenic_TREE)
    classification_NotSecreted_Random_TREE = reverse_classification(classification_Random_NotSecreted_TREE)
    classification_Secreted_Random_TREE = reverse_classification(classification_Random_Secreted_TREE)
    classification_SecretedEffectorCandidates_Random_TREE = reverse_classification(classification_Random_SecretedEffectors_TREE)
    classification_UTRs3_Random_TREE = reverse_classification(classification_Random_UTRs3_TREE)
    classification_UTRs5_Random_TREE = reverse_classification(classification_Random_UTRs5_TREE)

    classification_ORFs_Random_LOGISTIC = reverse_classification(classification_Random_ORFs_LOGISTIC)
    classification_Introns_Random_LOGISTIC = reverse_classification(classification_Random_Intron_LOGISTIC)
    classification_Intergenic_Random_LOGISTIC = reverse_classification(classification_Random_Intergenic_LOGISTIC)
    classification_NotSecreted_Random_LOGISTIC = reverse_classification(classification_Random_NotSecreted_LOGISTIC)
    classification_Secreted_Random_LOGISTIC = reverse_classification(classification_Random_Secreted_LOGISTIC)
    classification_SecretedEffectorCandidates_Random_LOGISTIC = reverse_classification(classification_Random_SecretedEffectors_LOGISTIC)
    classification_UTRs3_Random_LOGISTIC = reverse_classification(classification_Random_UTRs3_LOGISTIC)
    classification_UTRs5_Random_LOGISTIC = reverse_classification(classification_Random_UTRs5_LOGISTIC)
    # -----------------------------------------------------------------------------------------------------------
    ### Non-secreted CDS predictions
    # -----------------------------------------------------------------------------------------------------------
    classification_NotSecreted_ORFs_TREE = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_TREE + 'NotSecreted_ORFs.MODEL'], [ALGORITHM_TREE])
    classification_NotSecreted_Intergenic_TREE = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_TREE + 'NotSecreted_Intergenic.MODEL'], [ALGORITHM_TREE])
    classification_NotSecreted_Introns_TREE = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_TREE + 'NotSecreted_Introns.MODEL'], [ALGORITHM_TREE])
    classification_NotSecreted_UTRs3_TREE = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_TREE + 'NotSecreted_UTRs3.MODEL'], [ALGORITHM_TREE])
    classification_NotSecreted_UTRs5_TREE = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_TREE + 'NotSecreted_UTRs5.MODEL'], [ALGORITHM_TREE])

    classification_NotSecreted_ORFs_LOGISTIC = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_LOGISTIC + 'NotSecreted_ORFs.MODEL'], [ALGORITHM_LOGISTIC])
    classification_NotSecreted_Intergenic_LOGISTIC = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_LOGISTIC + 'NotSecreted_Intergenic.MODEL'], [ALGORITHM_LOGISTIC])
    classification_NotSecreted_Introns_LOGISTIC = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_LOGISTIC + 'NotSecreted_Introns.MODEL'], [ALGORITHM_LOGISTIC])
    classification_NotSecreted_UTRs3_LOGISTIC = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_LOGISTIC + 'NotSecreted_UTRs3.MODEL'], [ALGORITHM_LOGISTIC])
    classification_NotSecreted_UTRs5_LOGISTIC = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_LOGISTIC + 'NotSecreted_UTRs5.MODEL'], [ALGORITHM_LOGISTIC])

    classification = classification_NotSecreted_ORFs_TREE + classification_NotSecreted_Intergenic_TREE + classification_NotSecreted_Introns_TREE + classification_NotSecreted_UTRs3_TREE + classification_NotSecreted_UTRs5_TREE
    classification += classification_NotSecreted_ORFs_LOGISTIC + classification_NotSecreted_Intergenic_LOGISTIC + classification_NotSecreted_Introns_LOGISTIC + classification_NotSecreted_UTRs3_LOGISTIC + classification_NotSecreted_UTRs5_LOGISTIC
    classification += classification_NotSecreted_Random_TREE + classification_NotSecreted_Random_LOGISTIC

    PREDICTIONS = ensembl_voting(classification, PREDICTIONS)
    # -----------------------------------------------------------------------------------------------------------
    ### Secreted CDS predictions
    # -----------------------------------------------------------------------------------------------------------
    classification_Secreted_ORFs_TREE = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_TREE + 'Secreted_ORFs.MODEL'], [ALGORITHM_TREE])
    classification_Secreted_Intergenic_TREE = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_TREE + 'Secreted_Intergenic.MODEL'], [ALGORITHM_TREE])
    classification_Secreted_Introns_TREE = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_TREE + 'Secreted_Introns.MODEL'], [ALGORITHM_TREE])
    classification_Secreted_UTRs3_TREE = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_TREE + 'Secreted_UTRs3.MODEL'], [ALGORITHM_TREE])
    classification_Secreted_UTRs5_TREE = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_TREE + 'Secreted_UTRs5.MODEL'], [ALGORITHM_TREE])

    classification_Secreted_ORFs_LOGISTIC = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_LOGISTIC + 'Secreted_ORFs.MODEL'], [ALGORITHM_LOGISTIC])
    classification_Secreted_Intergenic_LOGISTIC = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_LOGISTIC + 'Secreted_Intergenic.MODEL'], [ALGORITHM_LOGISTIC])
    classification_Secreted_Introns_LOGISTIC = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_LOGISTIC + 'Secreted_Introns.MODEL'], [ALGORITHM_LOGISTIC])
    classification_Secreted_UTRs3_LOGISTIC = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_LOGISTIC + 'Secreted_UTRs3.MODEL'], [ALGORITHM_LOGISTIC])
    classification_Secreted_UTRs5_LOGISTIC = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_LOGISTIC + 'Secreted_UTRs5.MODEL'], [ALGORITHM_LOGISTIC])

    classification = classification_Secreted_ORFs_TREE + classification_Secreted_Intergenic_TREE + classification_Secreted_Introns_TREE + classification_Secreted_UTRs3_TREE + classification_Secreted_UTRs5_TREE
    classification += classification_Secreted_ORFs_LOGISTIC + classification_Secreted_Intergenic_LOGISTIC + classification_Secreted_Introns_LOGISTIC + classification_Secreted_UTRs3_LOGISTIC + classification_Secreted_UTRs5_LOGISTIC
    classification += classification_Secreted_Random_TREE + classification_Secreted_Random_LOGISTIC

    PREDICTIONS = ensembl_voting(classification, PREDICTIONS)
    # -----------------------------------------------------------------------------------------------------------
    ### Effector CDS predictions
    # -----------------------------------------------------------------------------------------------------------
    classification_SecretedEffectorCandidates_ORFs_TREE = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_TREE + 'SecretedEffectorCandidates_ORFs.MODEL'], [ALGORITHM_TREE])
    classification_SecretedEffectorCandidates_Intergenic_TREE = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_TREE + 'SecretedEffectorCandidates_Intergenic.MODEL'], [ALGORITHM_TREE])
    classification_SecretedEffectorCandidates_Introns_TREE = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_TREE + 'SecretedEffectorCandidates_Introns.MODEL'], [ALGORITHM_TREE])
    classification_SecretedEffectorCandidates_UTRs3_TREE = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_TREE + 'SecretedEffectorCandidates_UTRs3.MODEL'], [ALGORITHM_TREE])
    classification_SecretedEffectorCandidates_UTRs5_TREE = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_TREE + 'SecretedEffectorCandidates_UTRs5.MODEL'], [ALGORITHM_TREE])

    classification_SecretedEffectorCandidates_ORFs_LOGISTIC = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_LOGISTIC + 'SecretedEffectorCandidates_ORFs.MODEL'], [ALGORITHM_LOGISTIC])
    classification_SecretedEffectorCandidates_Intergenic_LOGISTIC = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_LOGISTIC + 'SecretedEffectorCandidates_Intergenic.MODEL'], [ALGORITHM_LOGISTIC])
    classification_SecretedEffectorCandidates_Introns_LOGISTIC = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_LOGISTIC + 'SecretedEffectorCandidates_Introns.MODEL'], [ALGORITHM_LOGISTIC])
    classification_SecretedEffectorCandidates_UTRs3_LOGISTIC = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_LOGISTIC + 'SecretedEffectorCandidates_UTRs3.MODEL'], [ALGORITHM_LOGISTIC])
    classification_SecretedEffectorCandidates_UTRs5_LOGISTIC = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_LOGISTIC + 'SecretedEffectorCandidates_UTRs5.MODEL'], [ALGORITHM_LOGISTIC])

    classification = classification_SecretedEffectorCandidates_ORFs_TREE + classification_SecretedEffectorCandidates_Intergenic_TREE + classification_SecretedEffectorCandidates_Introns_TREE + classification_SecretedEffectorCandidates_UTRs3_TREE + classification_SecretedEffectorCandidates_UTRs5_TREE
    classification += classification_SecretedEffectorCandidates_ORFs_LOGISTIC + classification_SecretedEffectorCandidates_Intergenic_LOGISTIC + classification_SecretedEffectorCandidates_Introns_LOGISTIC + classification_SecretedEffectorCandidates_UTRs3_LOGISTIC + classification_SecretedEffectorCandidates_UTRs5_LOGISTIC
    classification += classification_SecretedEffectorCandidates_Random_TREE + classification_SecretedEffectorCandidates_Random_LOGISTIC

    PREDICTIONS = ensembl_voting(classification, PREDICTIONS)
    # -----------------------------------------------------------------------------------------------------------
    ### Intron predictions
    # -----------------------------------------------------------------------------------------------------------
    classification_Introns_ORFs_TREE = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_TREE + 'Introns_ORFs.MODEL'], [ALGORITHM_TREE])
    classification_Introns_Intergenic_TREE = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_TREE + 'Introns_Intergenic.MODEL'], [ALGORITHM_TREE])
    classification_Introns_NotSecreted_TREE = reverse_classification(classification_NotSecreted_Introns_TREE)
    classification_Introns_Secreted_TREE = reverse_classification(classification_Secreted_Introns_TREE)
    classification_Introns_SecretedEffectorCandidates_TREE = reverse_classification(classification_SecretedEffectorCandidates_Introns_TREE)
    classification_Introns_UTRs3_TREE = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_TREE + 'Introns_UTRs3.MODEL'], [ALGORITHM_TREE])
    classification_Introns_UTRs5_TREE = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_TREE + 'Introns_UTRs5.MODEL'], [ALGORITHM_TREE])

    classification_Introns_ORFs_LOGISTIC = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_LOGISTIC + 'Introns_ORFs.MODEL'], [ALGORITHM_LOGISTIC])
    classification_Introns_Intergenic_LOGISTIC = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_LOGISTIC + 'Introns_Intergenic.MODEL'], [ALGORITHM_LOGISTIC])
    classification_Introns_NotSecreted_LOGISTIC = reverse_classification(classification_NotSecreted_Introns_LOGISTIC)
    classification_Introns_Secreted_LOGISTIC = reverse_classification(classification_Secreted_Introns_LOGISTIC)
    classification_Introns_SecretedEffectorCandidates_LOGISTIC = reverse_classification(classification_SecretedEffectorCandidates_Introns_LOGISTIC)
    classification_Introns_UTRs3_LOGISTIC = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_LOGISTIC + 'Introns_UTRs3.MODEL'], [ALGORITHM_LOGISTIC])
    classification_Introns_UTRs5_LOGISTIC = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_LOGISTIC + 'Introns_UTRs5.MODEL'], [ALGORITHM_LOGISTIC])

    classification = classification_Introns_ORFs_TREE + classification_Introns_Intergenic_TREE + classification_Introns_NotSecreted_TREE 
    classification += classification_Introns_Secreted_TREE + classification_Introns_SecretedEffectorCandidates_TREE + classification_Introns_UTRs3_TREE + classification_Introns_UTRs5_TREE
    classification += classification_Introns_ORFs_LOGISTIC + classification_Introns_Intergenic_LOGISTIC + classification_Introns_NotSecreted_LOGISTIC 
    classification += classification_Introns_Secreted_LOGISTIC + classification_Introns_SecretedEffectorCandidates_LOGISTIC + classification_Introns_UTRs3_LOGISTIC + classification_Introns_UTRs5_LOGISTIC
    classification += classification_Introns_Random_TREE + classification_Introns_Random_LOGISTIC

    PREDICTIONS = ensembl_voting(classification, PREDICTIONS)
    # -----------------------------------------------------------------------------------------------------------
    ### Intergenic predictions
    # -----------------------------------------------------------------------------------------------------------
    classification_Intergenic_Intron_TREE = reverse_classification(classification_Introns_Intergenic_TREE)
    classification_Intergenic_ORFs_TREE = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_TREE + 'Intergenic_ORFs.MODEL'], [ALGORITHM_TREE])
    classification_Intergenic_NotSecreted_TREE = reverse_classification(classification_NotSecreted_Intergenic_TREE)
    classification_Intergenic_Secreted_TREE = reverse_classification(classification_Secreted_Intergenic_TREE)
    classification_Intergenic_SecretedEffectorCandidates_TREE = reverse_classification(classification_SecretedEffectorCandidates_Intergenic_TREE)
    classification_Intergenic_UTRs3_TREE = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_TREE + 'Intergenic_UTRs3.MODEL'], [ALGORITHM_TREE])
    classification_Intergenic_UTRs5_TREE = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_TREE + 'Intergenic_UTRs5.MODEL'], [ALGORITHM_TREE])

    classification_Intergenic_Intron_LOGISTIC = reverse_classification(classification_Introns_Intergenic_LOGISTIC)
    classification_Intergenic_ORFs_LOGISTIC = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_LOGISTIC + 'Intergenic_ORFs.MODEL'], [ALGORITHM_LOGISTIC])
    classification_Intergenic_NotSecreted_LOGISTIC = reverse_classification(classification_NotSecreted_Intergenic_LOGISTIC)
    classification_Intergenic_Secreted_LOGISTIC = reverse_classification(classification_Secreted_Intergenic_LOGISTIC)
    classification_Intergenic_SecretedEffectorCandidates_LOGISTIC = reverse_classification(classification_SecretedEffectorCandidates_Intergenic_LOGISTIC)
    classification_Intergenic_UTRs3_LOGISTIC = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_LOGISTIC + 'Intergenic_UTRs3.MODEL'], [ALGORITHM_LOGISTIC])
    classification_Intergenic_UTRs5_LOGISTIC = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_LOGISTIC + 'Intergenic_UTRs5.MODEL'], [ALGORITHM_LOGISTIC])

    classification = classification_Intergenic_Intron_TREE + classification_Intergenic_ORFs_TREE + classification_Intergenic_NotSecreted_TREE
    classification += classification_Intergenic_Secreted_TREE + classification_Intergenic_SecretedEffectorCandidates_TREE + classification_Intergenic_UTRs3_TREE + classification_Intergenic_UTRs5_TREE
    classification += classification_Intergenic_Intron_LOGISTIC + classification_Intergenic_ORFs_LOGISTIC + classification_Intergenic_NotSecreted_LOGISTIC
    classification += classification_Intergenic_Secreted_LOGISTIC + classification_Intergenic_SecretedEffectorCandidates_LOGISTIC + classification_Intergenic_UTRs3_LOGISTIC + classification_Intergenic_UTRs5_LOGISTIC
    classification += classification_Intergenic_Random_TREE + classification_Intergenic_Random_LOGISTIC

    PREDICTIONS = ensembl_voting(classification, PREDICTIONS)
    # -----------------------------------------------------------------------------------------------------------
    ### 3' UTR predictions
    # -----------------------------------------------------------------------------------------------------------
    classification_UTRs3_Intron_TREE = reverse_classification(classification_Introns_UTRs3_TREE)
    classification_UTRs3_Intergenic_TREE = reverse_classification(classification_Intergenic_UTRs3_TREE)
    classification_UTRs3_NotSecreted_TREE = reverse_classification(classification_NotSecreted_UTRs3_TREE)
    classification_UTRs3_Secreted_TREE = reverse_classification(classification_Secreted_UTRs3_TREE)
    classification_UTRs3_SecretedEffectorCandidates_TREE = reverse_classification(classification_SecretedEffectorCandidates_UTRs3_TREE)
    classification_UTRs3_ORFs_TREE = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_TREE + 'UTRs3_ORFs.MODEL'], [ALGORITHM_TREE])
    classification_UTRs3_UTRs5_TREE = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_TREE + 'UTRs3_UTRs5.MODEL'], [ALGORITHM_TREE])

    classification_UTRs3_Intron_LOGISTIC = reverse_classification(classification_Introns_UTRs3_LOGISTIC)
    classification_UTRs3_Intergenic_LOGISTIC = reverse_classification(classification_Intergenic_UTRs3_LOGISTIC)
    classification_UTRs3_NotSecreted_LOGISTIC = reverse_classification(classification_NotSecreted_UTRs3_LOGISTIC)
    classification_UTRs3_Secreted_LOGISTIC = reverse_classification(classification_Secreted_UTRs3_LOGISTIC)
    classification_UTRs3_SecretedEffectorCandidates_LOGISTIC = reverse_classification(classification_SecretedEffectorCandidates_UTRs3_LOGISTIC)
    classification_UTRs3_ORFs_LOGISTIC = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_LOGISTIC + 'UTRs3_ORFs.MODEL'], [ALGORITHM_LOGISTIC])
    classification_UTRs3_UTRs5_LOGISTIC = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_LOGISTIC + 'UTRs3_UTRs5.MODEL'], [ALGORITHM_LOGISTIC])

    classification = classification_UTRs3_Intron_TREE + classification_UTRs3_Intergenic_TREE + classification_UTRs3_NotSecreted_TREE
    classification += classification_UTRs3_Secreted_TREE + classification_UTRs3_SecretedEffectorCandidates_TREE + classification_UTRs3_ORFs_TREE + classification_UTRs3_UTRs5_TREE
    classification += classification_UTRs3_Intron_LOGISTIC + classification_UTRs3_Intergenic_LOGISTIC + classification_UTRs3_NotSecreted_LOGISTIC
    classification += classification_UTRs3_Secreted_LOGISTIC + classification_UTRs3_SecretedEffectorCandidates_LOGISTIC + classification_UTRs3_ORFs_LOGISTIC + classification_UTRs3_UTRs5_LOGISTIC
    classification += classification_UTRs3_Random_TREE + classification_UTRs3_Random_LOGISTIC

    PREDICTIONS = ensembl_voting(classification, PREDICTIONS)
    # -----------------------------------------------------------------------------------------------------------
    ### 5' UTR predictions
    # -----------------------------------------------------------------------------------------------------------
    classification_UTRs5_Intron_TREE = reverse_classification(classification_Introns_UTRs5_TREE)
    classification_UTRs5_Intergenic_TREE = reverse_classification(classification_Intergenic_UTRs5_TREE)
    classification_UTRs5_NotSecreted_TREE = reverse_classification(classification_NotSecreted_UTRs5_TREE)
    classification_UTRs5_Secreted_TREE = reverse_classification(classification_Secreted_UTRs5_TREE)
    classification_UTRs5_SecretedEffectorCandidates_TREE = reverse_classification(classification_SecretedEffectorCandidates_UTRs5_TREE)
    classification_UTRs5_ORFs_TREE = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_TREE + 'UTRs5_ORFs.MODEL'], [ALGORITHM_TREE])
    classification_UTRs5_UTRs3_TREE = reverse_classification(classification_UTRs3_UTRs5_TREE)

    classification_UTRs5_Intron_LOGISTIC = reverse_classification(classification_Introns_UTRs5_LOGISTIC)
    classification_UTRs5_Intergenic_LOGISTIC = reverse_classification(classification_Intergenic_UTRs5_LOGISTIC)
    classification_UTRs5_NotSecreted_LOGISTIC = reverse_classification(classification_NotSecreted_UTRs5_LOGISTIC)
    classification_UTRs5_Secreted_LOGISTIC = reverse_classification(classification_Secreted_UTRs5_LOGISTIC)
    classification_UTRs5_SecretedEffectorCandidates_LOGISTIC = reverse_classification(classification_SecretedEffectorCandidates_UTRs5_LOGISTIC)
    classification_UTRs5_ORFs_LOGISTIC = get_model_predictions(weka_input, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_LOGISTIC + 'UTRs5_ORFs.MODEL'], [ALGORITHM_LOGISTIC])
    classification_UTRs5_UTRs3_LOGISTIC = reverse_classification(classification_UTRs3_UTRs5_LOGISTIC)

    classification = classification_UTRs5_Intron_TREE + classification_UTRs5_Intergenic_TREE + classification_UTRs5_NotSecreted_TREE
    classification += classification_UTRs5_Secreted_TREE + classification_UTRs5_SecretedEffectorCandidates_TREE + classification_UTRs5_ORFs_TREE + classification_UTRs5_UTRs3_TREE
    classification += classification_UTRs5_Intron_LOGISTIC + classification_UTRs5_Intergenic_LOGISTIC + classification_UTRs5_NotSecreted_LOGISTIC
    classification += classification_UTRs5_Secreted_LOGISTIC + classification_UTRs5_SecretedEffectorCandidates_LOGISTIC + classification_UTRs5_ORFs_LOGISTIC + classification_UTRs5_UTRs3_LOGISTIC
    classification += classification_UTRs5_Random_TREE + classification_UTRs5_Random_LOGISTIC

    PREDICTIONS = ensembl_voting(classification, PREDICTIONS)
    # -----------------------------------------------------------------------------------------------------------
    ### ORF predictions
    # -----------------------------------------------------------------------------------------------------------
    classification_ORFs_Intron_TREE = reverse_classification(classification_Introns_ORFs_TREE)
    classification_ORFs_Intergenic_TREE = reverse_classification(classification_Intergenic_ORFs_TREE)
    classification_ORFs_NotSecreted_TREE = reverse_classification(classification_NotSecreted_ORFs_TREE)
    classification_ORFs_Secreted_TREE = reverse_classification(classification_Secreted_ORFs_TREE)
    classification_ORFs_SecretedEffectorCandidates_TREE = reverse_classification(classification_SecretedEffectorCandidates_ORFs_TREE)
    classification_ORFs_UTRs5_TREE = reverse_classification(classification_UTRs5_ORFs_TREE)
    classification_ORFs_UTRs3_TREE = reverse_classification(classification_UTRs3_ORFs_TREE)

    classification_ORFs_Intron_LOGISTIC = reverse_classification(classification_Introns_ORFs_LOGISTIC)
    classification_ORFs_Intergenic_LOGISTIC = reverse_classification(classification_Intergenic_ORFs_LOGISTIC)
    classification_ORFs_NotSecreted_LOGISTIC = reverse_classification(classification_NotSecreted_ORFs_LOGISTIC)
    classification_ORFs_Secreted_LOGISTIC = reverse_classification(classification_Secreted_ORFs_LOGISTIC)
    classification_ORFs_SecretedEffectorCandidates_LOGISTIC = reverse_classification(classification_SecretedEffectorCandidates_ORFs_LOGISTIC)
    classification_ORFs_UTRs5_LOGISTIC = reverse_classification(classification_UTRs5_ORFs_LOGISTIC)
    classification_ORFs_UTRs3_LOGISTIC = reverse_classification(classification_UTRs3_ORFs_LOGISTIC)

    classification = classification_ORFs_Intron_TREE + classification_ORFs_Intergenic_TREE + classification_ORFs_NotSecreted_TREE
    classification += classification_ORFs_Secreted_TREE + classification_ORFs_SecretedEffectorCandidates_TREE + classification_ORFs_UTRs5_TREE + classification_ORFs_UTRs3_TREE
    classification += classification_ORFs_Intron_LOGISTIC + classification_ORFs_Intergenic_LOGISTIC + classification_ORFs_NotSecreted_LOGISTIC
    classification += classification_ORFs_Secreted_LOGISTIC + classification_ORFs_SecretedEffectorCandidates_LOGISTIC + classification_ORFs_UTRs5_LOGISTIC + classification_ORFs_UTRs3_LOGISTIC
    classification += classification_ORFs_Random_TREE + classification_ORFs_Random_LOGISTIC

    PREDICTIONS = ensembl_voting(classification, PREDICTIONS)
    # -----------------------------------------------------------------------------------------------------------
    MODELS_KOZAK = [
    PATH_TO_MODELS + PREFIX_LOGISTIC + 'Kozak_Random.MODEL',
    PATH_TO_MODELS + PREFIX_TREE + 'Kozak_Random.MODEL'
    ]

    classification_Kozak_Random_TREE = get_model_predictions(weka_input_kozak, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_TREE + 'Kozak_Random.MODEL'], [ALGORITHM_TREE])
    classification_Kozak_Random_LOGISTIC = get_model_predictions(weka_input_kozak, WEKA_PATH, RESULTS_PATH, [PATH_TO_MODELS + PREFIX_LOGISTIC + 'Kozak_Random.MODEL'], [ALGORITHM_LOGISTIC])

    classification = classification_Kozak_Random_TREE + classification_Kozak_Random_LOGISTIC

    PREDICTIONS = ensembl_voting(classification, PREDICTIONS)
    # -----------------------------------------------------------------------------------------------------------
    # Lastly, add the random sequence predictions
    # -----------------------------------------------------------------------------------------------------------
    classification = classification_Random_ORFs_TREE + classification_Random_Intron_TREE + classification_Random_Intergenic_TREE
    classification += classification_Random_NotSecreted_TREE + classification_Random_Secreted_TREE + classification_Random_SecretedEffectors_TREE + classification_Random_UTRs3_TREE + classification_Random_UTRs5_TREE
    classification += classification_Random_ORFs_LOGISTIC + classification_Random_Intron_LOGISTIC + classification_Random_Intergenic_LOGISTIC
    classification += classification_Random_NotSecreted_LOGISTIC + classification_Random_Secreted_LOGISTIC + classification_Random_SecretedEffectors_LOGISTIC + classification_Random_UTRs3_LOGISTIC + classification_Random_UTRs5_LOGISTIC

    PREDICTIONS = ensembl_voting(classification, PREDICTIONS)
    # -----------------------------------------------------------------------------------------------------------
    f = open(RESULTS_PATH + 'EffectorGeneP_Predictions.txt', 'w')
    # -----------------------------------------------------------------------------------------------------------
    output = '#Identifier\tCDS_NotSecreted\tCDS_Secreted\tCDS_Effector\tIntron\tIntergenic\t3UTR\t5UTR\tORF\tKozak\tRandom\n'
    # -----------------------------------------------------------------------------------------------------------
    for short_ident, long_ident in zip(SHORT_IDENTIFIERS, ORIGINAL_IDENTIFIERS):

        CDS_NotSec_yes_prob = PREDICTIONS[short_ident][0]
        CDS_NotSec_no_prob = PREDICTIONS[short_ident][1]

        CDS_Sec_yes_prob = PREDICTIONS[short_ident][2]
        CDS_Sec_no_prob = PREDICTIONS[short_ident][3]

        CDS_Eff_yes_prob = PREDICTIONS[short_ident][4]
        CDS_Eff_no_prob = PREDICTIONS[short_ident][5]

        Intron_yes_prob = PREDICTIONS[short_ident][6]
        Intron_no_prob = PREDICTIONS[short_ident][7]

        Intergenic_yes_prob = PREDICTIONS[short_ident][8]
        Intergenic_no_prob = PREDICTIONS[short_ident][9]

        UTR3_yes_prob = PREDICTIONS[short_ident][10]
        UTR3_no_prob = PREDICTIONS[short_ident][11]

        UTR5_yes_prob = PREDICTIONS[short_ident][12]
        UTR5_no_prob = PREDICTIONS[short_ident][13]    

        ORF_yes_prob = PREDICTIONS[short_ident][14]
        ORF_no_prob = PREDICTIONS[short_ident][15]    

        # These scores should only be used for sequences in the KOZAK length range (15 nts)
        Kozak_yes_prob = PREDICTIONS[short_ident][16]
        Kozak_no_prob = PREDICTIONS[short_ident][17]   

        Random_yes_prob = PREDICTIONS[short_ident][18]
        Random_no_prob = PREDICTIONS[short_ident][19]

        output += long_ident + '\t'
        output += str(CDS_NotSec_yes_prob) + '\t' + str(CDS_Sec_yes_prob) + '\t' + str(CDS_Eff_yes_prob) + '\t'
        output += str(Intron_yes_prob) + '\t' + str(Intergenic_yes_prob) + '\t' + str(UTR3_yes_prob) + '\t' + str(UTR5_yes_prob) + '\t' + str(ORF_yes_prob) + '\t' + str(Kozak_yes_prob) + '\t' + str(Random_yes_prob) + '\n'
        
    f.writelines(output.rstrip('\n'))
    f.close()

    return


