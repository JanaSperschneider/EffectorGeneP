#!/usr/bin/env python3

"""
EffectorGeneP: gene annotation in pathogen genomes

Copyright (c) 2025 Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.

This work is freely available for non-commercial scientific research, non-commercial education, 
or non-commercial research projects, under the CSIRO Non-Commercial License (https://github.com/JanaSperschneider/EffectorGeneP/blob/main/LICENCE). 
"""
# -----------------------------------------------------------------------------------------------------------
import os
import sys
import subprocess

from scipy.stats import norm

import variables
# -----------------------------------------------------------------------------------------------------------
def run_signalP(signalp_input, RESULTS_PATH, SIGNALP_PATH):
    ID_DIC, SIGNALP = {}, {}

    if signalp_input != []:
	    f = open(RESULTS_PATH + 'SignalP_Input.fasta', 'w')
	    for index, ORF in enumerate(list(set(signalp_input))):
	        f.writelines('>' "seq" + str(index) + '\n')
	        f.writelines(ORF[0] + '\n')
	        ID_DIC["seq" + str(index)] = ORF[3]
	    f.close()

	    ParamList = [SIGNALP_PATH + ' -t euk -f short -u 0.34 -U 0.34 ' + RESULTS_PATH + 'SignalP_Input.fasta']

	    with open(RESULTS_PATH + 'SignalP_Output.txt', 'wb') as out:
	        try:
	            Process = subprocess.Popen(ParamList, shell=True, stdout=out)
	            sts = Process.wait()
	            cstdout, cstderr = Process.communicate()

	            if Process.returncode:
	                raise Exception("Calling SignalP 4.1 returned %s"%Process.returncode)
	            if cstdout:
	                pass
	            elif cstderr:
	                sys.exit()
	        except:
	            e = sys.exc_info()[1]
	            print("Error calling SignalP: %s" % e)
	            sys.exit(1)

	    f = open(RESULTS_PATH + 'SignalP_Output.txt', 'r')
	    content = f.readlines()
	    f.close()

	    for line in content:	    	
	        if line.startswith('#'):
	            pass
	        else:
	            predictions = line.split(' ')
	            predictions = [entry for entry in predictions if entry != '']
	            chopped_id = predictions[0]
	            true_id = ID_DIC[chopped_id]

	            if float(predictions[8]) > 0.340:
	            	# In line with SignalP, it has to be > and not >=
	                SIGNALP[true_id] = float(predictions[8])

    return SIGNALP
# -----------------------------------------------------------------------------------------------------------
def run_tmhmm(signalp_input, SIGNALP, RESULTS_PATH, TMHMM_PATH):
	ID_DIC, TMHMM = {}, {}

	if signalp_input != []:
	    f = open(RESULTS_PATH + 'TMHMM_Input.fasta', 'w')
	    for index, ORF in enumerate(list(set(signalp_input))):
	    	if ORF[3] in SIGNALP:
		        f.writelines('>' "seq" + str(index) + '\n')
		        f.writelines(ORF[0] + '\n')
		        ID_DIC["seq" + str(index)] = ORF[3]
	    f.close()

	    ParamList = [TMHMM_PATH + ' -noplot < ' + RESULTS_PATH + 'TMHMM_Input.fasta > ' + RESULTS_PATH + 'TMHMM_Output.txt']

	    with open(RESULTS_PATH + 'TMHMM_Output.txt', 'wb') as out:
	        try:
	            Process = subprocess.Popen(ParamList, shell=True, stdout=out)
	            sts = Process.wait()
	            cstdout, cstderr = Process.communicate()

	            if Process.returncode:
	                raise Exception("Calling TMHMM 2.0 returned %s"%Process.returncode)
	            if cstdout:
	                pass
	            elif cstderr:
	                sys.exit()
	        except:
	            e = sys.exc_info()[1]
	            print("Error calling TMHMM: %s" % e)
	            sys.exit(1)

	    f = open(RESULTS_PATH + 'TMHMM_Output.txt', 'r')
	    content = f.readlines()
	    f.close()

	    for index, line in enumerate(content):
	        # Now check if the start of the protein where the signal peptide is could also be a transmembrane domain
	        if 'Number of predicted TMHs:' in line:
	            ident = line.split('Number of predicted TMHs:')[0].replace('#','').strip()
	            number_of_tmhmms = int(line.split('Number of predicted TMHs:')[1])             

	            if number_of_tmhmms > 0:                
	                cache = content[index:]

	                if '# ' + ident + ' POSSIBLE N-term signal sequence\n' in cache:
	                    number_of_tmhmms -= 1

	            if number_of_tmhmms > 0:
	            	true_id = ID_DIC[ident]
	            	TMHMM[true_id] = number_of_tmhmms

	return TMHMM
# -----------------------------------------------------------------------------------------------------------
def update_scores_signalp(ORFs_per_ISOFORM_DIC, SIGNALP, TMHMM, cds_scores_effector, CONSERVATIVE):

    ORFs_FINAL_SCORES = {}

    for transcript_id, ORFs in ORFs_per_ISOFORM_DIC.items():

        for this_ORF in ORFs:
            full_id = this_ORF[3]

            if full_id in SIGNALP:
                signalp_score = SIGNALP[full_id]
            else:
                signalp_score = 0.0

            if full_id in TMHMM:
                tmhmm_number = TMHMM[full_id]
            else:
                tmhmm_number = 0

            # ORFs with signal peptides get a higher gene score, but only if they have a relatively high Kozak score and are possible effectors
            if CONSERVATIVE == False:
            	MIN_EFFECTOR_CDS_SCORE_SP = variables.MIN_EFFECTOR_CDS_SCORE_SP
            else:
            	MIN_EFFECTOR_CDS_SCORE_SP = variables.CONSERVATIVE_MIN_EFFECTOR_CDS_SCORE_SP

            if signalp_score > 0.340 and tmhmm_number == 0 and cds_scores_effector[full_id][0] > MIN_EFFECTOR_CDS_SCORE_SP:
                # In line with SignalP, it has to be > and not >=
                # Kozak score filtering here
                if this_ORF[10] >= variables.KOZAK_MINIMUM:
                    updated_score = this_ORF[1] + variables.SIGNALP_SCORE
                    updated_ORF = (this_ORF[0], updated_score, this_ORF[2], full_id, this_ORF[4], this_ORF[5], this_ORF[6], this_ORF[7], this_ORF[8], this_ORF[9], this_ORF[10], this_ORF[11])   
                else:
                    updated_ORF = this_ORF
            else:
                updated_ORF = this_ORF

            if transcript_id in ORFs_FINAL_SCORES:
                ORFs_FINAL_SCORES[transcript_id] = ORFs_FINAL_SCORES[transcript_id] + [updated_ORF]
            else:
                ORFs_FINAL_SCORES[transcript_id] = [updated_ORF]

    return ORFs_FINAL_SCORES
# -----------------------------------------------------------------------------------------------------------
def parse_EffectorGeneP_CDS_scores(RESULTS_PATH, mean_not_secreted, std_dev_not_secreted, mean_secreted, std_dev_secreted, mean_effector, std_dev_effector, CONSERVATIVE):

	f = open(RESULTS_PATH + 'EffectorGeneP_Predictions.txt', 'r')
	content = f.readlines()
	f.close()
	# -----------------------------------------------------------------------------------------------------------
	cds_scores, cds_scores_effector = {}, {}

	cds_to_keep = []
	# -----------------------------------------------------------------------------------------------------------
	# Collect the scores and the coordinates of the genomic feature
	ALL_SCORES_FOR_CDS_CANDIDATES = {}
	# -----------------------------------------------------------------------------------------------------------
	for line in content:
	    if line.startswith('#'):
	        pass
	    else:
	        #print(line)
	        # ---------------------------------------------------------------------------------------------------        
	        identifier_field = line.split('\t')[0]
	        # ---------------------------------------------------------------------------------------------------
	        CDS_NotSecreted = float(line.split('\t')[1])
	        CDS_Secreted = float(line.split('\t')[2])
	        CDS_Effector = float(line.split('\t')[3])
	        Intron = float(line.split('\t')[4]) 
	        Intergenic = float(line.split('\t')[5]) 
	        UTR3 = float(line.split('\t')[6])
	        UTR5 = float(line.split('\t')[7])
	        orf_score = float(line.split('\t')[8])	   
	        kozak_score = float(line.split('\t')[9])	 
	        random_score = float(line.split('\t')[10])	 
	        # ---------------------------------------------------------------------------------------------------                
	        if '_CDS_' in identifier_field:

	            candidate_id = line.split('\t')[0].split('_CDS_')[0] + '_' + line.split('\t')[0].split('_')[-1]      
	            transcript_id = line.split('\t')[0].split('_ORF')[0]

	            ### Record the best CDS score ###
	            cds_score = max(CDS_NotSecreted, CDS_Secreted, CDS_Effector)
	            pvalue_cds_score = 1.0

	            if max(CDS_NotSecreted, CDS_Secreted, CDS_Effector) == CDS_NotSecreted:
	            	pvalue_cds_score = pvalue_from_score(CDS_NotSecreted, mean_not_secreted, std_dev_not_secreted)

	            if max(CDS_NotSecreted, CDS_Secreted, CDS_Effector) == CDS_Secreted:
	            	pvalue_cds_score = pvalue_from_score(CDS_Secreted, mean_secreted, std_dev_secreted)

	            if max(CDS_NotSecreted, CDS_Secreted, CDS_Effector) == CDS_Effector:
	            	pvalue_cds_score = pvalue_from_score(CDS_Effector, mean_effector, std_dev_effector)	    

	            # For transcript fusion events we want to keep a record of all the scores to use later
	            ALL_SCORES_FOR_CDS_CANDIDATES[candidate_id] = [cds_score, Intron, Intergenic, UTR3, UTR5, orf_score, random_score]

	            if CONSERVATIVE == False:
	                MIN_CDS_SCORE = variables.MIN_CDS_SCORE
	            else:
	                MIN_CDS_SCORE = variables.CONSERVATIVE_MIN_CDS_SCORE

	            if cds_score > MIN_CDS_SCORE and pvalue_cds_score < variables.PVALUE_THRESHOLD: 

	                # These are 1-based coordinates of the CDS on the genome, including introns
	                start = int(line.split('\t')[0].split('_')[-3])
	                end = int(line.split('\t')[0].split('_')[-2])           
	                cds_scores[candidate_id] = (cds_score, start, end, line.split('_CDS_')[0].split('_ORF')[0])
	                cds_scores_effector[candidate_id] = (CDS_Effector, start, end, line.split('_CDS_')[0].split('_ORF')[0])	    

	                cds_to_keep.append(identifier_field)        

	return cds_scores, cds_scores_effector, ALL_SCORES_FOR_CDS_CANDIDATES, cds_to_keep
# -----------------------------------------------------------------------------------------------------------
def parse_EffectorGeneP_gene_scores(RESULTS_PATH):

	f = open(RESULTS_PATH + 'EffectorGeneP_Predictions.txt', 'r')
	content = f.readlines()
	f.close()
	# -----------------------------------------------------------------------------------------------------------
	intron_scores = {}
	utr3_scores = {}
	utr5_scores = {}
	kozak_scores = {}
	# -----------------------------------------------------------------------------------------------------------
	# Collect the scores and the coordinates of the genomic feature
	ALL_SCORES_FOR_CDS_CANDIDATES = {}
	# -----------------------------------------------------------------------------------------------------------
	for line in content:
	    if line.startswith('#'):
	        pass
	    else:
	        # ---------------------------------------------------------------------------------------------------        
	        identifier_field = line.split('\t')[0] 
	        # ---------------------------------------------------------------------------------------------------                
	        CDS_NotSecreted = float(line.split('\t')[1])
	        CDS_Secreted = float(line.split('\t')[2])
	        CDS_Effector = float(line.split('\t')[3])
	        Intron = float(line.split('\t')[4]) 
	        Intergenic = float(line.split('\t')[5]) 
	        UTR3 = float(line.split('\t')[6])
	        UTR5 = float(line.split('\t')[7])
	        orf_score = float(line.split('\t')[8])	   
	        kozak_score = float(line.split('\t')[9])	 
	        random_score = float(line.split('\t')[10])	 

	        if '_INTRON_' in identifier_field:
	            candidate_id = line.split('_INTRON_')[0] + '_' + line.split('\t')[0].split('_')[-1]
	            intron_score = Intron
	            # These are 1-based coordinates
	            intron_start = int(line.split('\t')[0].split('_')[-3])
	            intron_end = int(line.split('\t')[0].split('_')[-2])
	            intron_strand = line.split('\t')[0].split('_')[-1].strip()	                  
	            	
	            if candidate_id in intron_scores:
	                intron_scores[candidate_id] = intron_scores[candidate_id] + [(intron_score, intron_start, intron_end, intron_strand)]
	            else:
	                intron_scores[candidate_id] = [(intron_score, intron_start, intron_end, intron_strand)]

	        elif '_UTR3_' in identifier_field:
	            candidate_id = line.split('\t')[0].split('_UTR3_')[0] + '_' + line.split('\t')[0].split('_')[-1]
	            utr_score = UTR3

	            # These are 1-based
	            utr_start = int(line.split('\t')[0].split('_')[-3])
	            utr_end = int(line.split('\t')[0].split('_')[-2])
	            utr3_scores[candidate_id] = utr_score

	        elif '_UTR5_' in identifier_field:
	            candidate_id = line.split('\t')[0].split('_UTR5_')[0] + '_' + line.split('\t')[0].split('_')[-1]
	            utr_score = UTR5       

	            # These are 1-based
	            utr_start = int(line.split('\t')[0].split('_')[-3])
	            utr_end = int(line.split('\t')[0].split('_')[-2])
	            utr5_scores[candidate_id] = utr_score                 

	        elif '_KOZAK_' in identifier_field:
	            candidate_id = line.split('\t')[0].split('_KOZAK_')[0] + '_' + line.split('\t')[0].split('_')[-1]
	            kozak_scores[candidate_id] = kozak_score    

	return intron_scores, utr3_scores, utr5_scores, kozak_scores
# -----------------------------------------------------------------------------------------------------------
def collect_score(full_id, cds_score, cds_start, cds_end, intron_id, intron_scores, utr3_score, utr5_score, kozak_score, ORF_blocks, blacklisted_introns, INTRON_MIN_PROB, CONSERVATIVE):

    gene_score = 0.0
    # -----------------------------------------------------------------------------------------------------------
    intron_score_list = []
    ORF_blocks[full_id] = [('CDS', cds_start, cds_end)] 
    # -----------------------------------------------------------------------------------------------------------     
    if intron_id in intron_scores:
        for (score, intron_start, intron_end, intron_strand) in intron_scores[intron_id]:
            #     i.......i
            #  c.................c
            if intron_start > cds_start and intron_start < cds_end and intron_end < cds_end and intron_end > cds_start:
                intron_score_list.append(score)
                ORF_blocks[full_id] = ORF_blocks[full_id] + [('Intron', intron_start, intron_end)]
                if (intron_id, intron_start, intron_end, intron_strand) in blacklisted_introns:
                    intron_score_list = [-1000000.0]

    if intron_score_list != []:
        # This gene has introns, do they have low average probability?
        intron_score = sum(intron_score_list)/len(intron_score_list)  
        if intron_score < INTRON_MIN_PROB:
            gene_score = 0.0
        else:              
            gene_score = round(100.0*(variables.CODING_WEIGHT*cds_score + variables.INTRON_WEIGHT*intron_score + variables.UTR3_WEIGHT*utr3_score + variables.UTR5_WEIGHT*utr5_score + variables.KOZAK_WEIGHT*kozak_score)/(variables.CODING_WEIGHT + variables.INTRON_WEIGHT + variables.UTR3_WEIGHT + variables.UTR5_WEIGHT + variables.KOZAK_WEIGHT), 2)
    else:
        intron_score = 0.0 
        gene_score = round(100.0*(variables.CODING_WEIGHT*cds_score + variables.UTR3_WEIGHT*utr3_score + variables.UTR5_WEIGHT*utr5_score + variables.KOZAK_WEIGHT*kozak_score)/(variables.CODING_WEIGHT + variables.UTR3_WEIGHT + variables.UTR5_WEIGHT + variables.KOZAK_WEIGHT), 2)

    if CONSERVATIVE == False:
        MIN_GENE_SCORE = variables.MIN_GENE_SCORE
    else:
        MIN_GENE_SCORE = variables.CONSERVATIVE_MIN_GENE_SCORE

    if gene_score <= MIN_GENE_SCORE:
        gene_score = 0.0

    return gene_score, ORF_blocks, intron_score, intron_score_list
# -----------------------------------------------------------------------------------------------------------
def weighted_gene_score(cds_length, intron_length, gene_score):

	if gene_score <= variables.MIN_GENE_SCORE:
		gene_score_length = 0.0
	else:
		gene_score_length = gene_score*(cds_length + intron_length)/100.0
		
	return gene_score_length
# -----------------------------------------------------------------------------------------------------------
def pvalue_from_score(cds_score, mean, std_dev):

	if std_dev > 0.0:
		z_score = (cds_score - mean) / std_dev
		p_value = norm.sf(z_score)
	else:
		p_value = 0.0

	return p_value
# -----------------------------------------------------------------------------------------------------------
