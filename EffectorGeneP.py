#!/usr/bin/env python3

"""
EffectorGeneP: gene annotation in pathogen genomes

Copyright (c) 2025 Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.

This work is freely available for non-commercial scientific research, non-commercial education, 
or non-commercial research projects, under the CSIRO Non-Commercial License (https://github.com/JanaSperschneider/EffectorGeneP/blob/main/LICENCE). 
"""
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
import os
import sys
import math
import subprocess
import errno
import uuid
import shutil
import tempfile
# -----------------------------------------------------------------------------------------------------------
import process_gff3
import functions
import scores 
import write_gff3 
import classify
import transcript_fusions
import orfs
import variables
import version
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
SCRIPT_PATH = sys.path[0]
WEKA_PATH = SCRIPT_PATH  + '/target/weka-csiro-agriculture-and-food-0.0.1-SNAPSHOT.jar'
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# Check that the path to the WEKA software exists
path_exists = os.access(WEKA_PATH, os.F_OK)
if not path_exists:
    print()
    print("Path to WEKA software does not exist!")
    print("Check the installation and the given path to the WEKA software %s in EffectorGeneP.py (line 27)." % WEKA_PATH)
    print()
    sys.exit(1)
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
def write_candidate_sequences_for_EffectorGeneP(OUTPUT_FILE, candidate_list, replace_string, mode):

    with open(OUTPUT_FILE, mode) as f:

        for identifier, seq in candidate_list:
            if replace_string != "":
                f.writelines('>' + identifier.replace('CDS', replace_string) + '\n')
                f.writelines(seq + '\n')
            else:
                f.writelines('>' + identifier + '\n')
                f.writelines(seq + '\n')              

    return
# -----------------------------------------------------------------------------------------------------------
def write_candidate_sequences_introns(OUTPUT_FILE, candidate_introns):

    with open(OUTPUT_FILE, 'a') as f:

        for identifier, seq in candidate_introns:
            f.writelines(identifier)
            f.writelines(seq + '\n')

    return
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# Main Program starts here
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
commandline = sys.argv[1:]
# -----------------------------------------------------------------------------------------------------------
if commandline:
    GENOME, TRANSCRIPT_GFF3, PATH_TO_MODELS, out_file, MIN_PROTEIN_LENGTH, STRANDEDNESS, TRANSCRIPT_PADDING, SIGNALP_PATH, TMHMM_PATH, START_SITE_DISTANCE_MAX, INTRON_MIN_PROB, CONSERVATIVE = functions.scan_arguments(commandline)
    # If none of the mandatory arguments were provided
    if not GENOME:
        print()
        print('Please specify a genome FASTA input file using the -g option!')
        functions.usage()
    if not TRANSCRIPT_GFF3:
        print()
        print('Please specify a transcript GFF3 input file (formatted by gffread) using the -t option!')
        functions.usage()
    if not PATH_TO_MODELS:
        print()
        print('Please specify the full path to the EffectorGeneP model files using the -m option!\nDownload species models from here: https://effectorp.csiro.au/effectorgenep.html')
        functions.usage()
    if not out_file:
        print()
        print('Please specify the name to the EffectorGeneP output GFF3 file using the -o option!')
        functions.usage()      

    # Check that the path to the softwares exist
    if SIGNALP_PATH != None:
        path_exists = os.access(SIGNALP_PATH, os.F_OK)
        if not path_exists:
            print()
            print("Given path to SIGNALP 4.1 executable does not exist!")
            print("Check the installation and the given path to the SIGNALP 4.1 software (e.g. /yyy/xxx/signalp-4.1/signalp). Provided path is currently: %s" % SIGNALP_PATH)
            print()
            sys.exit(1)
    if TMHMM_PATH != None:
        path_exists = os.access(TMHMM_PATH, os.F_OK)
        if not path_exists:
            print()
            print("Given path to TMHMM 2.0 executable does not exist!")
            print("Check the installation and the given path to the TMHMM 2.0 software (e.g. /yyy/xxx/tmhmm-2.0c/bin/tmhmm). Provided path is currently: %s" % TMHMM_PATH)
            print()
            sys.exit(1)

else:
    functions.usage() 
# -----------------------------------------------------------------------------------------------------------
print("---------------------------------")
print("EffectorGeneP is running with these options:")
print("---------------------------------")
print("Genome:", GENOME)
print("Transcripts:", TRANSCRIPT_GFF3)
print("Model files:", PATH_TO_MODELS)
print("Output file:", out_file)
print("Minimum protein length:", MIN_PROTEIN_LENGTH, "aas")
print("Strand stringency:", STRANDEDNESS)
print("Transcript padding (bps):", TRANSCRIPT_PADDING)
print("Path to SignalP 4.1:", SIGNALP_PATH)
print("Path to TMHMM 2.0:", TMHMM_PATH)
print("How far upstream to look for another translation start site (aas):", START_SITE_DISTANCE_MAX)
print("Minimum average intron probability for multi-exon genes:", INTRON_MIN_PROB)
if CONSERVATIVE == True:
    print("Genes are annotated in conservative mode.")
print("EffectorGeneP version is:", str(version.VERSION))
print("---------------------------------")

gff3_commandline = "#"
gff3_commandline += " EffectorGeneP.py -t " + TRANSCRIPT_GFF3 + " -g " + GENOME + " -m " + PATH_TO_MODELS  + " -o " + out_file 
gff3_commandline += " -l " + str(MIN_PROTEIN_LENGTH) + " -s " + STRANDEDNESS + " -p " + str(TRANSCRIPT_PADDING) 
gff3_commandline += " -d " + str(START_SITE_DISTANCE_MAX) + " -i " + str(INTRON_MIN_PROB)
if CONSERVATIVE == True:
    gff3_commandline += " -c" 
gff3_commandline += " --SIGNALP4 " + str(SIGNALP_PATH) + " --TMHMM " + str(TMHMM_PATH) 
gff3_commandline += "\n# EffectorGeneP " + "version " + str(version.VERSION)
# -----------------------------------------------------------------------------------------------------------
if TRANSCRIPT_PADDING > 0:
    USE_PADDING = "YES"
else:
    USE_PADDING = "NO"
# -----------------------------------------------------------------------------------------------------------
# Temporary folder name identifier that will be used to store results
RESULTS_PATH = tempfile.mkdtemp() + '/'
# -----------------------------------------------------------------------------------------------------------
# Try to create folder where results will be stored
try:
    os.mkdir(RESULTS_PATH)
except OSError as exception:        
    if exception.errno != errno.EEXIST:
        raise
# -----------------------------------------------------------------------------------------------------------
FORWARD_SEQUENCES = {}
# Read in all the contigs
for identifier, sequence in functions.SimpleFastaParser(open(GENOME, 'r')):
    identifier = identifier.replace('>','').split()[0]
    sequence = sequence.replace('*','').upper()
    FORWARD_SEQUENCES[identifier] = sequence
# -----------------------------------------------------------------------------------------------------------
print("---------------------------------")
print("Derive mean and variance for random CDSs with similar codon content to genome for p-value calculations.")
print("---------------------------------")
# -----------------------------------------------------------------------------------------------------------
# Check if the file exists
if os.path.exists(PATH_TO_MODELS + '/random_ORF_composition.EffectorGeneP_Predictions.txt'):

    print("Now parse the EffectorGeneP classification scores from file:", PATH_TO_MODELS + '/random_ORF_composition.EffectorGeneP_Predictions.txt')
    cds_not_secreted_scores, cds_secreted_scores, cds_effector_scores = [], [], []
    utr3_scores, utr5_scores, intron_scores = [], [], []

    f = open(PATH_TO_MODELS + '/random_ORF_composition.EffectorGeneP_Predictions.txt', 'r')
    content = f.readlines()
    f.close()
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

            cds_not_secreted_scores.append(CDS_NotSecreted)
            cds_secreted_scores.append(CDS_Secreted)
            cds_effector_scores.append(CDS_Effector)

    mean_not_secreted = (sum(cds_not_secreted_scores)/len(cds_not_secreted_scores))
    var = sum((x - mean_not_secreted) ** 2 for x in cds_not_secreted_scores) / len(cds_not_secreted_scores)
    std_dev_not_secreted = (var ** 0.5)

    mean_secreted = (sum(cds_secreted_scores)/len(cds_secreted_scores))
    var = sum((x - mean_secreted) ** 2 for x in cds_secreted_scores) / len(cds_secreted_scores)
    std_dev_secreted = (var ** 0.5)

    mean_effector = (sum(cds_effector_scores)/len(cds_effector_scores))
    var = sum((x - mean_effector) ** 2 for x in cds_effector_scores) / len(cds_effector_scores)
    std_dev_effector = (var ** 0.5)

else:
    e = sys.exc_info()[1]
    print("The training file", PATH_TO_MODELS + '/random_ORF_composition.EffectorGeneP_Predictions.txt', "does not exist. %s" % e)
    sys.exit(1)
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
print("---------------------------------")
# -----------------------------------------------------------------------------------------------------------
print("Parse transcript information from gff3 file (gffread formatted).")
# -----------------------------------------------------------------------------------------------------------
f = open(TRANSCRIPT_GFF3, 'r')
content = f.readlines()
f.close()
# -----------------------------------------------------------------------------------------------------------
# Go through the gff3 file and collect all the information about transcripts (1-based), exons (0-based)
# -----------------------------------------------------------------------------------------------------------
exons, transcript_information, TRANSCRIPT_COORDINATES_PER_CONTIG, TRANSCRIPT_IDs_MAPPING = process_gff3.gff3_read_exons(content)
# -----------------------------------------------------------------------------------------------------------
# Now add the introns (0-based)
introns = process_gff3.gff3_derive_introns_from_exons(exons)
# -----------------------------------------------------------------------------------------------------------
# Now go through each contig and add the padding for transcripts
# -----------------------------------------------------------------------------------------------------------
TRANSCRIPT_COORDINATES_WITH_PADDING, exons_with_padding = process_gff3.add_padding(TRANSCRIPT_COORDINATES_PER_CONTIG, exons, transcript_information, FORWARD_SEQUENCES, TRANSCRIPT_PADDING)
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
print("Done parsing transcript information, now run ORF predictions for this many transcripts:", len(transcript_information))
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
sequences_to_classify = 0.0
transcripts_that_were_padded = []
# -----------------------------------------------------------------------------------------------------------
ORFs_TO_CLASSIFY = {}
# -----------------------------------------------------------------------------------------------------------
all_candidate_cds, all_candidate_utr3s, all_candidate_utr5s, all_candidate_kozaks = [], [], [], []
# -----------------------------------------------------------------------------------------------------------
for transcript_id, (contig, transcript_start, transcript_end, gene_id, transcript_strand) in transcript_information.items():
    # -----------------------------------------------------------------------------------------------------------
    # Concatenate the exons to form the mRNA transcript
    # -----------------------------------------------------------------------------------------------------------
    SEQUENCE = ""
    for (exon_start, exon_end, contig, strand) in sorted(exons[transcript_id]):
        SEQUENCE += FORWARD_SEQUENCES[contig][exon_start:exon_end+1]
    SEQUENCE_REV = functions.rev_comp(SEQUENCE)

    SEQUENCE_PADDED = ""
    for (exon_start, exon_end, contig, strand) in sorted(exons_with_padding[transcript_id]):
        SEQUENCE_PADDED += FORWARD_SEQUENCES[contig][exon_start:exon_end+1]       
    SEQUENCE_REV_PADDED = functions.rev_comp(SEQUENCE_PADDED) 
    # -----------------------------------------------------------------------------------------------------------
    # -----------------------------------------------------------------------------------------------------------
    ### Note that single-exon genes in stringtie have no strand information, so investigate both strands for ORFs
    # -----------------------------------------------------------------------------------------------------------
    # -----------------------------------------------------------------------------------------------------------
    # Here, assume that the strand information is not necessarily correct
    # -----------------------------------------------------------------------------------------------------------     
    if STRANDEDNESS == "NON_STRINGENT":
        # -----------------------------------------------------------------------------------------------------------
        strand = '+'
        if USE_PADDING == "NO":
            ORFs_TO_CLASSIFY, candidate_cds, candidate_utr3s, candidate_utr5s, candidate_kozaks = orfs.ORFS_per_transcript(transcript_id, SEQUENCE, exons, ORFs_TO_CLASSIFY, strand, TRANSCRIPT_PADDING, TRANSCRIPT_COORDINATES_WITH_PADDING, transcript_information, MIN_PROTEIN_LENGTH, START_SITE_DISTANCE_MAX)
        # -----------------------------------------------------------------------------------------------------------                        
        if USE_PADDING == "YES":
            ORFs_TO_CLASSIFY, candidate_cds, candidate_utr3s, candidate_utr5s, candidate_kozaks = orfs.ORFS_per_transcript(transcript_id, SEQUENCE_PADDED, exons_with_padding, ORFs_TO_CLASSIFY, strand, TRANSCRIPT_PADDING, TRANSCRIPT_COORDINATES_WITH_PADDING, transcript_information, MIN_PROTEIN_LENGTH, START_SITE_DISTANCE_MAX)
            transcripts_that_were_padded.append(transcript_id)                   
        # -----------------------------------------------------------------------------------------------------------      
        all_candidate_cds += candidate_cds
        all_candidate_utr3s += candidate_utr3s 
        all_candidate_utr5s += candidate_utr5s   
        all_candidate_kozaks += candidate_kozaks
        # -----------------------------------------------------------------------------------------------------------
        # -----------------------------------------------------------------------------------------------------------
        strand = '-'
        if USE_PADDING == "NO":
            ORFs_TO_CLASSIFY, candidate_cds, candidate_utr3s, candidate_utr5s, candidate_kozaks = orfs.ORFS_per_transcript(transcript_id, SEQUENCE_REV, exons, ORFs_TO_CLASSIFY, strand, TRANSCRIPT_PADDING, TRANSCRIPT_COORDINATES_WITH_PADDING, transcript_information, MIN_PROTEIN_LENGTH, START_SITE_DISTANCE_MAX)   
        # -----------------------------------------------------------------------------------------------------------
        if USE_PADDING == "YES":
            ORFs_TO_CLASSIFY, candidate_cds, candidate_utr3s, candidate_utr5s, candidate_kozaks = orfs.ORFS_per_transcript(transcript_id, SEQUENCE_REV_PADDED, exons_with_padding, ORFs_TO_CLASSIFY, strand, TRANSCRIPT_PADDING, TRANSCRIPT_COORDINATES_WITH_PADDING, transcript_information, MIN_PROTEIN_LENGTH, START_SITE_DISTANCE_MAX) 
            transcripts_that_were_padded.append(transcript_id)
        # -----------------------------------------------------------------------------------------------------------      
        all_candidate_cds += candidate_cds
        all_candidate_utr3s += candidate_utr3s 
        all_candidate_utr5s += candidate_utr5s   
        all_candidate_kozaks += candidate_kozaks        
    # -----------------------------------------------------------------------------------------------------------             
    # Here, assume that the strand information is correct
    # -----------------------------------------------------------------------------------------------------------         
    if STRANDEDNESS == "STRINGENT":
        # -----------------------------------------------------------------------------------------------------------        
        # Now use the reverse sequence if on the - strand
        if transcript_strand == '-':
            if USE_PADDING == "NO":
                ORFs_TO_CLASSIFY, candidate_cds, candidate_utr3s, candidate_utr5s, candidate_kozaks = orfs.ORFS_per_transcript(transcript_id, SEQUENCE_REV, exons, ORFs_TO_CLASSIFY, transcript_strand, TRANSCRIPT_PADDING, TRANSCRIPT_COORDINATES_WITH_PADDING, transcript_information, MIN_PROTEIN_LENGTH, START_SITE_DISTANCE_MAX)
            # -----------------------------------------------------------------------------------------------------------                        
            if USE_PADDING == "YES": 
                ORFs_TO_CLASSIFY, candidate_cds, candidate_utr3s, candidate_utr5s, candidate_kozaks = orfs.ORFS_per_transcript(transcript_id, SEQUENCE_REV_PADDED, exons_with_padding, ORFs_TO_CLASSIFY, transcript_strand, TRANSCRIPT_PADDING, TRANSCRIPT_COORDINATES_WITH_PADDING, transcript_information, MIN_PROTEIN_LENGTH, START_SITE_DISTANCE_MAX)
                transcripts_that_were_padded.append(transcript_id)
            # -----------------------------------------------------------------------------------------------------------      
            all_candidate_cds += candidate_cds
            all_candidate_utr3s += candidate_utr3s 
            all_candidate_utr5s += candidate_utr5s   
            all_candidate_kozaks += candidate_kozaks            
        # -----------------------------------------------------------------------------------------------------------
        if transcript_strand == '+':
            if USE_PADDING == "NO":
                ORFs_TO_CLASSIFY, candidate_cds, candidate_utr3s, candidate_utr5s, candidate_kozaks = orfs.ORFS_per_transcript(transcript_id, SEQUENCE, exons, ORFs_TO_CLASSIFY, transcript_strand, TRANSCRIPT_PADDING, TRANSCRIPT_COORDINATES_WITH_PADDING, transcript_information, MIN_PROTEIN_LENGTH, START_SITE_DISTANCE_MAX)
            # -----------------------------------------------------------------------------------------------------------
            if USE_PADDING == "YES":
                ORFs_TO_CLASSIFY, candidate_cds, candidate_utr3s, candidate_utr5s, candidate_kozaks = orfs.ORFS_per_transcript(transcript_id, SEQUENCE_PADDED, exons_with_padding, ORFs_TO_CLASSIFY, transcript_strand, TRANSCRIPT_PADDING, TRANSCRIPT_COORDINATES_WITH_PADDING, transcript_information, MIN_PROTEIN_LENGTH, START_SITE_DISTANCE_MAX)
                transcripts_that_were_padded.append(transcript_id)
            # -----------------------------------------------------------------------------------------------------------      
            all_candidate_cds += candidate_cds
            all_candidate_utr3s += candidate_utr3s 
            all_candidate_utr5s += candidate_utr5s   
            all_candidate_kozaks += candidate_kozaks              
        # -----------------------------------------------------------------------------------------------------------        
        if transcript_strand == '.':
            # -----------------------------------------------------------------------------------------------------------
            strand = '+'
            if USE_PADDING == "NO":
                ORFs_TO_CLASSIFY, candidate_cds, candidate_utr3s, candidate_utr5s, candidate_kozaks = orfs.ORFS_per_transcript(transcript_id, SEQUENCE, exons, ORFs_TO_CLASSIFY, strand, TRANSCRIPT_PADDING, TRANSCRIPT_COORDINATES_WITH_PADDING, transcript_information, MIN_PROTEIN_LENGTH, START_SITE_DISTANCE_MAX)
            # -----------------------------------------------------------------------------------------------------------                        
            if USE_PADDING == "YES":
                ORFs_TO_CLASSIFY, candidate_cds, candidate_utr3s, candidate_utr5s, candidate_kozaks = orfs.ORFS_per_transcript(transcript_id, SEQUENCE_PADDED, exons_with_padding, ORFs_TO_CLASSIFY, strand, TRANSCRIPT_PADDING, TRANSCRIPT_COORDINATES_WITH_PADDING, transcript_information, MIN_PROTEIN_LENGTH, START_SITE_DISTANCE_MAX)
                transcripts_that_were_padded.append(transcript_id)             
            # -----------------------------------------------------------------------------------------------------------      
            all_candidate_cds += candidate_cds
            all_candidate_utr3s += candidate_utr3s 
            all_candidate_utr5s += candidate_utr5s   
            all_candidate_kozaks += candidate_kozaks              
            # -----------------------------------------------------------------------------------------------------------
            # -----------------------------------------------------------------------------------------------------------
            strand = '-'
            if USE_PADDING == "NO":
                ORFs_TO_CLASSIFY, candidate_cds, candidate_utr3s, candidate_utr5s, candidate_kozaks = orfs.ORFS_per_transcript(transcript_id, SEQUENCE_REV, exons, ORFs_TO_CLASSIFY, strand, TRANSCRIPT_PADDING, TRANSCRIPT_COORDINATES_WITH_PADDING, transcript_information, MIN_PROTEIN_LENGTH, START_SITE_DISTANCE_MAX)   
            # -----------------------------------------------------------------------------------------------------------
            if USE_PADDING == "YES":
                ORFs_TO_CLASSIFY, candidate_cds, candidate_utr3s, candidate_utr5s, candidate_kozaks = orfs.ORFS_per_transcript(transcript_id, SEQUENCE_REV_PADDED, exons_with_padding, ORFs_TO_CLASSIFY, strand, TRANSCRIPT_PADDING, TRANSCRIPT_COORDINATES_WITH_PADDING, transcript_information, MIN_PROTEIN_LENGTH, START_SITE_DISTANCE_MAX) 
                transcripts_that_were_padded.append(transcript_id)
            # -----------------------------------------------------------------------------------------------------------      
            all_candidate_cds += candidate_cds
            all_candidate_utr3s += candidate_utr3s 
            all_candidate_utr5s += candidate_utr5s   
            all_candidate_kozaks += candidate_kozaks              
        # -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
print("Done with ORF collection for all the transcripts.")
# -----------------------------------------------------------------------------------------------------------
# Write the CDS sequences first
write_candidate_sequences_for_EffectorGeneP(RESULTS_PATH + "Candidates.fasta", all_candidate_cds, "", 'w')
# -----------------------------------------------------------------------------------------------------------
print("---------------------------------")
print("Run EffectorGeneP machine learning classification, this might take some time.")
print()
print("Classifying this many CDS sequences with EffectorGeneP:", str(len(all_candidate_cds)))
print("---------------------------------")
# -----------------------------------------------------------------------------------------------------------
if len(all_candidate_cds) != 0.0:
    classify.run_EffectorGeneP(RESULTS_PATH + 'Candidates.fasta', RESULTS_PATH, PATH_TO_MODELS, WEKA_PATH)
else:
    print('Nothing to do, quit EffectorGeneP.')
    sys.exit()
# -----------------------------------------------------------------------------------------------------------
cds_scores, cds_scores_effector, ALL_SCORES_FOR_CDS_CANDIDATES, cds_to_keep = scores.parse_EffectorGeneP_CDS_scores(RESULTS_PATH, mean_not_secreted, std_dev_not_secreted, mean_secreted, std_dev_secreted, mean_effector, std_dev_effector, CONSERVATIVE)
print("Of these CDS sequences,", len(cds_to_keep), "will be considered further as gene candidates.")
print("---------------------------------")
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
if len(cds_scores) == 0.0:
    print('Nothing to do, quit EffectorGeneP.')
    sys.exit()

print("Now add the UTRs, introns and Kozak sequences for classification.")
filtered_candidate_utr3s = [(identifier, utr3) for identifier, utr3 in all_candidate_utr3s if identifier in cds_to_keep]
filtered_candidate_utr5s = [(identifier, utr5) for identifier, utr5 in all_candidate_utr5s if identifier in cds_to_keep]
filtered_candidate_kozaks = [(identifier, kozak) for identifier, kozak in all_candidate_kozaks if identifier in cds_to_keep]

write_candidate_sequences_for_EffectorGeneP(RESULTS_PATH + "Candidates.fasta", filtered_candidate_utr3s, "UTR3", 'w')
write_candidate_sequences_for_EffectorGeneP(RESULTS_PATH + "Candidates.fasta", filtered_candidate_utr5s, "UTR5", 'a')
write_candidate_sequences_for_EffectorGeneP(RESULTS_PATH + "Candidates.fasta", filtered_candidate_kozaks, "KOZAK", 'a')
# -----------------------------------------------------------------------------------------------------------
# Now add the intron sequences to the FASTA file for classification, add all strands here
all_candidate_introns = []
blacklisted_introns = []

for transcript, (contig, transcript_start, transcript_end, gene_id, transcript_strand) in transcript_information.items():
    if transcript in introns:
        for (start, end, contig, strand) in sorted(introns[transcript]):
            SEQUENCE = FORWARD_SEQUENCES[contig][start:end+1]
            SEQUENCE_REV = functions.rev_comp(SEQUENCE)

            if transcript_strand == '.':

                if SEQUENCE != "":
                    # 1-based coordinates
                    strand = '+'
                    identifier= '>' + transcript.replace('"','') + '_INTRON_' + contig + '_' + str(start + 1) + '_' + str(end + 1) + '_' + strand + '\n'
                    all_candidate_introns.append( (identifier, str(SEQUENCE)) )
                if SEQUENCE != "":
                    strand = '-'
                    # 1-based coordinates
                    identifier= '>' + transcript.replace('"','') + '_INTRON_' + contig + '_' + str(start + 1) + '_' + str(end + 1) + '_' + strand + '\n'
                    all_candidate_introns.append( (identifier, str(SEQUENCE_REV)) )

            if transcript_strand == '+':

                if SEQUENCE != "":
                    # 1-based coordinates
                    strand = '+'
                    identifier= '>' + transcript.replace('"','') + '_INTRON_' + contig + '_' + str(start + 1) + '_' + str(end + 1) + '_' + strand + '\n'
                    all_candidate_introns.append( (identifier, str(SEQUENCE)) )

                if SEQUENCE != "":
                    strand = '-'
                    # 1-based coordinates
                    identifier= '>' + transcript.replace('"','') + '_INTRON_' + contig + '_' + str(start + 1) + '_' + str(end + 1) + '_' + strand + '\n'
                    all_candidate_introns.append( (identifier, str(SEQUENCE_REV)) )                    
                    blacklisted_introns.append( (transcript.replace('"','') + "_" + strand, int(start + 1), int(end + 1), strand) )


            if transcript_strand == '-':                    
                if SEQUENCE != "":
                    strand = '-'
                    # 1-based coordinates
                    identifier= '>' + transcript.replace('"','') + '_INTRON_' + contig + '_' + str(start + 1) + '_' + str(end + 1) + '_' + strand + '\n'
                    all_candidate_introns.append( (identifier, str(SEQUENCE_REV)) )

                if SEQUENCE != "":
                    # 1-based coordinates
                    strand = '+'
                    identifier= '>' + transcript.replace('"','') + '_INTRON_' + contig + '_' + str(start + 1) + '_' + str(end + 1) + '_' + strand + '\n'
                    all_candidate_introns.append( (identifier, str(SEQUENCE)) )                    
                    blacklisted_introns.append( (transcript.replace('"','') + "_" + strand, int(start + 1), int(end + 1), strand) )


write_candidate_sequences_introns(RESULTS_PATH + "Candidates.fasta", all_candidate_introns)    
# -----------------------------------------------------------------------------------------------------------
print('''-----------------------------------------------------''')
print("Run EffectorGeneP machine learning classification, this might take some time.")
print()
print("Classifying this many UTR, intron and Kozak sequences with EffectorGeneP:", str(len(filtered_candidate_utr3s) + len(filtered_candidate_utr5s) + len(filtered_candidate_kozaks) + len(all_candidate_introns)))
print('''-----------------------------------------------------''')
if len(filtered_candidate_utr3s) + len(filtered_candidate_utr5s) + len(filtered_candidate_kozaks) + len(all_candidate_introns) == 0.0:
    intron_scores, utr3_scores, utr5_scores, kozak_scores = [], [], [], []
else:       
    classify.run_EffectorGeneP(RESULTS_PATH + 'Candidates.fasta', RESULTS_PATH, PATH_TO_MODELS, WEKA_PATH)
    intron_scores, utr3_scores, utr5_scores, kozak_scores = scores.parse_EffectorGeneP_gene_scores(RESULTS_PATH)
# -----------------------------------------------------------------------------------------------------------
print("Now parse the EffectorGeneP classification scores")
# -----------------------------------------------------------------------------------------------------------
# Now calculate the combined score for each ORF (CDS + introns + UTRs)
# -----------------------------------------------------------------------------------------------------------
ORF_blocks = {}
# -----------------------------------------------------------------------------------------------------------
ORFs_per_ISOFORM_DIC = {}
# -----------------------------------------------------------------------------------------------------------
# 1st step: collect scores for the ORFs for each transcript
# -----------------------------------------------------------------------------------------------------------
for full_id, (cds_score, cds_start, cds_end, intron_id) in cds_scores.items():

    transcript_id = full_id.split('_ORF')[0]
    strand = full_id.split('_ORF')[1].split('_')[1]
    contig, transcript_start, transcript_end, gene_id, transcript_strand = transcript_information[transcript_id]
    # ---------------------------------------------------------------------------------------------------
    if CONSERVATIVE == False:
        MIN_CDS_SCORE = variables.MIN_CDS_SCORE
    else:
        MIN_CDS_SCORE = variables.CONSERVATIVE_MIN_CDS_SCORE

    if cds_score >= MIN_CDS_SCORE:
        # -----------------------------------------------------------------------------------------------------------        
        if full_id in utr3_scores:
            utr3_score = utr3_scores[full_id]
        else:
            utr3_score = 0.95 # In case the UTR was too short, but do not set at 1 because of transcript fusions

        if full_id in utr5_scores:
            utr5_score = utr5_scores[full_id]
        else:
            utr5_score = 0.95 # In case the UTR was too short, but do not set at 1 because of transcript fusions     

        if full_id in kozak_scores:
            kozak_score = kozak_scores[full_id]
        else:
            kozak_score = 1.0 # In case the ORF started at position 1-9 in the transcript    
        # ----------------------------------------------------------------------------------------------------------- 
        intron_scores_strand = {}
        intron_id = transcript_id + '_' + strand

        if intron_id in intron_scores:
            intron_scores_strand_values = intron_scores[intron_id]
            intron_scores_strand[intron_id] = intron_scores_strand_values

        gene_score, ORF_blocks, intron_score, intron_score_list = scores.collect_score(full_id, cds_score, cds_start, cds_end, intron_id, intron_scores_strand, utr3_score, utr5_score, kozak_score, ORF_blocks, blacklisted_introns, INTRON_MIN_PROB, CONSERVATIVE)

        # If this is a strand conflict, penalize the putative gene for overall gene length score
        if strand != transcript_strand and transcript_strand != '.':
            if intron_score < INTRON_MIN_PROB and intron_score_list != []:
                gene_score = 0.0 
            else:
                gene_score = gene_score - variables.STRAND_PENALTY

        encoded_protein = ORFs_TO_CLASSIFY[full_id][2]
        transcript_coverage = ORFs_TO_CLASSIFY[full_id][5]
        # -----------------------------------------------------------------------------------------------------------        
        # Now check if padding was applied
        if transcript_start > cds_start:
            transcript_start = cds_start
        if transcript_end < cds_end:
            transcript_end = cds_end
        # -----------------------------------------------------------------------------------------------------------        
        cds_length = cds_end-cds_start+1
        intron_length = sum([entry[2]-entry[1]+1 for entry in ORF_blocks[full_id] if entry[0] == 'Intron'])
        # -----------------------------------------------------------------------------------------------------------        
        gene_score_length = scores.weighted_gene_score(cds_length, intron_length, gene_score)
        #print(full_id, cds_score, intron_score, utr3_score, utr5_score, "kozak:", kozak_score, "effector cds:", cds_scores_effector[full_id][0])
        #print(gene_score, cds_length, intron_length, gene_score_length)
        #print(encoded_protein)
        #print('---------------')
        # -----------------------------------------------------------------------------------------------------------        
        if gene_score_length > 0.0:
            if round(math.log(gene_score_length), 2) > variables.MIN_GENE_LENGTH_SCORE_LOG:

                if transcript_id in ORFs_per_ISOFORM_DIC:
                    ORFs_per_ISOFORM_DIC[transcript_id] = ORFs_per_ISOFORM_DIC[transcript_id] + [(encoded_protein, gene_score, cds_score, full_id, transcript_id, transcript_start, transcript_end, cds_start, cds_end, strand, kozak_score, transcript_coverage)]
                else:
                    ORFs_per_ISOFORM_DIC[transcript_id] = [(encoded_protein, gene_score, cds_score, full_id, transcript_id, transcript_start, transcript_end, cds_start, cds_end, strand, kozak_score, transcript_coverage)]                                    
    # -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
signalp_input = []
 # -----------------------------------------------------------------------------------------------------------
for transcript_id, ORFs in ORFs_per_ISOFORM_DIC.items():
    for this_ORF in ORFs:
        signalp_input.append(this_ORF) 
# -----------------------------------------------------------------------------------------------------------
if SIGNALP_PATH != None and TMHMM_PATH != None:
    # Also assign a higher score to ORFs with a signal peptide
    # -----------------------------------------------------------------------------------------------------------
    print('Now run signal peptide search')
    SIGNALP = scores.run_signalP(signalp_input, RESULTS_PATH, SIGNALP_PATH)
    print('Now run transmembrane domain search')
    TMHMM = scores.run_tmhmm(signalp_input, SIGNALP, RESULTS_PATH, TMHMM_PATH)
else:
    SIGNALP, TMHMM = {}, {}

ORFs_FINAL_SCORES = scores.update_scores_signalp(ORFs_per_ISOFORM_DIC, SIGNALP, TMHMM, cds_scores_effector, CONSERVATIVE)
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# Now find the best ORF for the transcript
# -----------------------------------------------------------------------------------------------------------
best_ORF_for_transcript = {}
# -----------------------------------------------------------------------------------------------------------
for transcript_id, ORFs in ORFs_FINAL_SCORES.items():

    contig, transcript_start, transcript_end, gene_id, transcript_strand = transcript_information[transcript_id]

    best_ORF_score = 0.0 

    for (encoded_protein, gene_score, cds_score, full_id, short_id, transcript_start, transcript_end, cds_start, cds_end, strand, kozak_score, transcript_coverage) in ORFs: 
        # -----------------------------------------------------------------------------------------------------------        
        cds_length = cds_end-cds_start+1
        intron_length = sum([entry[2]-entry[1]+1 for entry in ORF_blocks[full_id] if entry[0] == 'Intron'])
        # -----------------------------------------------------------------------------------------------------------                       
        gene_score_length = scores.weighted_gene_score(cds_length, intron_length, gene_score)
        # -----------------------------------------------------------------------------------------------------------                 
        if gene_score_length > best_ORF_score:
            best_ORF_score = gene_score_length
            best_ORF_gene_score = gene_score
            best_ORF_gene_strand = strand
            best_ORF_full_id = full_id

            best_ORF_utr3 = ORFs_TO_CLASSIFY[full_id][3]
            best_ORF_utr5 = ORFs_TO_CLASSIFY[full_id][4]

            best_ORF_cds_start = cds_start
            best_ORF_cds_end = cds_end              

    if best_ORF_score != 0.0:
        best_ORF_for_transcript[transcript_id] = best_ORF_score, best_ORF_gene_score, best_ORF_gene_strand, best_ORF_full_id, len(best_ORF_utr5), len(best_ORF_utr3), best_ORF_cds_start, best_ORF_cds_end
# -----------------------------------------------------------------------------------------------------------
# Here, select the ORFs that achieve the highest score for the transcript interval
# This deals with transcript fusion cases where multiple genes might be part of a fused transcript
# -----------------------------------------------------------------------------------------------------------
GENE_ISOFORMS_DIC_CHOSEN = {} 
COORDINATES = {}

for transcript_id, ORFs in ORFs_FINAL_SCORES.items():
    # ------------------------------------------------------------
    gene_id = transcript_information[transcript_id][3]
    # ------------------------------------------------------------
    highest_score = best_ORF_for_transcript[transcript_id][0]
    strand_of_highest_score = best_ORF_for_transcript[transcript_id][2]
    ORF_ID_highest_score = best_ORF_for_transcript[transcript_id][3]
    highest_score_utr5_length = best_ORF_for_transcript[transcript_id][4]
    highest_score_utr3_length = best_ORF_for_transcript[transcript_id][5]
    cds_start = best_ORF_for_transcript[transcript_id][6]
    cds_end = best_ORF_for_transcript[transcript_id][7]
    # ------------------------------------------------------------
    # Transcript boundaries without padding
    contig, transcript_start, transcript_end, gene_id, transcript_strand = transcript_information[transcript_id]

    # If UTR is too short there will be no score in the dictionary
    if ORF_ID_highest_score in utr3_scores:
        utr3_score_of_best_ORF = utr3_scores[ORF_ID_highest_score]
    else:
        utr3_score_of_best_ORF = 0.0 

    if ORF_ID_highest_score in utr5_scores:
        utr5_score_of_best_ORF = utr5_scores[ORF_ID_highest_score]
    else:
        utr5_score_of_best_ORF = 0.0
    # ------------------------------------------------------------
    investigate_transcript_fusions = False
    investigate_transcript_fusions_3utr_plus, investigate_transcript_fusions_5utr_plus = False, False
    investigate_transcript_fusions_3utr_minus, investigate_transcript_fusions_5utr_minus = False, False
    transcript_fusion_regions = []

    if (transcript_end - transcript_start + 1) > variables.TRANSCRIPT_MINLENGTH_FUSION and 100.0*(cds_end - cds_start + 1)/(transcript_end - transcript_start + 1) < variables.TRANSCRIPT_FUSION_BEST_GENE_COVERAGE:
        if strand_of_highest_score == '+':

            if highest_score_utr3_length > variables.MIN_TRANSCRIPT_LENGTH:
                # If it has a very high UTR score, it is unlikely a transcript fusion
                if utr3_score_of_best_ORF > 0.98:
                    investigate_transcript_fusions_3utr_plus = False 
                else:
                    investigate_transcript_fusions_3utr_plus = True
                    transcript_fusion_regions.append([cds_end + variables.CDS_SPACER, transcript_end + TRANSCRIPT_PADDING])                    

            if highest_score_utr5_length > variables.MIN_TRANSCRIPT_LENGTH:
                # If it has a very high UTR score, it is unlikely a transcript fusion
                if utr5_score_of_best_ORF > 0.98:
                    investigate_transcript_fusions_5utr_plus = False 
                else:
                    investigate_transcript_fusions_5utr_plus = True  
                    transcript_fusion_regions.append([transcript_start - TRANSCRIPT_PADDING, cds_start - variables.CDS_SPACER])

            if investigate_transcript_fusions_3utr_plus == True or investigate_transcript_fusions_5utr_plus == True:
                investigate_transcript_fusions = True


        if strand_of_highest_score == '-':

            if highest_score_utr3_length > variables.MIN_TRANSCRIPT_LENGTH:
                # If it has a very high UTR score, it is unlikely a transcript fusion
                if utr3_score_of_best_ORF > 0.98:
                    investigate_transcript_fusions_3utr_minus = False 
                else:
                    investigate_transcript_fusions_3utr_minus = True
                    transcript_fusion_regions.append([transcript_start - TRANSCRIPT_PADDING, cds_start - variables.CDS_SPACER])                    

            if highest_score_utr5_length > variables.MIN_TRANSCRIPT_LENGTH:
                # If it has a very high UTR score, it is unlikely a transcript fusion
                if utr5_score_of_best_ORF > 0.98:
                    investigate_transcript_fusions_5utr_minus = False 
                else:
                    investigate_transcript_fusions_5utr_minus = True  
                    transcript_fusion_regions.append([cds_end + variables.CDS_SPACER, transcript_end + TRANSCRIPT_PADDING])

            if investigate_transcript_fusions_3utr_minus == True or investigate_transcript_fusions_5utr_minus == True:
                investigate_transcript_fusions = True
    # ------------------------------------------------------------
    intervals = []

    if investigate_transcript_fusions == True:

        # Include the ORF with the highest score by default
        for (encoded_protein, gene_score, cds_score, full_id, short_id, transcript_start, transcript_end, cds_start, cds_end, strand, kozak_score, transcript_coverage) in ORFs:

            if full_id == ORF_ID_highest_score:
                intervals.append( transcript_fusions.GeneInterval(cds_start - variables.CDS_SPACER, cds_end + variables.CDS_SPACER, 1000000.0) )
                COORDINATES[(cds_start - variables.CDS_SPACER, cds_end + variables.CDS_SPACER, 1000000.0)] = (encoded_protein, gene_score, cds_score, full_id, short_id, transcript_start, transcript_end, cds_start, cds_end, strand, kozak_score, transcript_coverage)
    
        for (encoded_protein, gene_score, cds_score, full_id, short_id, transcript_start, transcript_end, cds_start, cds_end, strand, kozak_score, transcript_coverage) in ORFs:

            valid_transcript_fusion_candidate = False 
            for (region_start, region_end) in transcript_fusion_regions:
                if cds_start >= region_start and cds_start <= region_end and cds_end >= region_start and cds_end <= region_end:
                    valid_transcript_fusion_candidate = True

            if valid_transcript_fusion_candidate == True:

                # Choose genes that maximize gene scores in the transcript
                # Note that the UTR scores are not ideal and should be re-calculated in the future
                # ------------------------------------------------------------
                ## Make sure to use overall length including introns here
                # ------------------------------------------------------------
                cds_length = cds_end-cds_start+1
                intron_length = sum([entry[2]-entry[1]+1 for entry in ORF_blocks[full_id] if entry[0] == 'Intron'])
                cds_end_with_introns = cds_end + intron_length

                gene_score_length = scores.weighted_gene_score(cds_length, intron_length, gene_score)

                COORDINATES[(cds_start - variables.CDS_SPACER, cds_end + variables.CDS_SPACER, gene_score_length)] = (encoded_protein, gene_score, cds_score, full_id, short_id, transcript_start, transcript_end, cds_start, cds_end, strand, kozak_score, transcript_coverage)

                if gene_score_length > 0.0:
                    if round(math.log(gene_score_length), 2) > variables.MIN_GENE_LENGTH_SCORE_LOG_TRANSCRIPT_FUSION:

                        contig = transcript_information[transcript_id][0]
                        contig_length = len(FORWARD_SEQUENCES[contig])

                        # This are the scores for the CDS to be something other, i.e. intron/intergenic/UTR/ORF
                        cds_score, intron_score, intergenic_score, UTR3_score, UTR5_score, ORF_score, random_score = ALL_SCORES_FOR_CDS_CANDIDATES[full_id]

                        ### This is a stringent filtering process here
                        if cds_score < intron_score or cds_score < UTR3_score or cds_score < UTR5_score:
                            # Rescue if it has a high Kozak score
                            if kozak_score >= 0.5:
                                intervals.append( transcript_fusions.GeneInterval(cds_start - variables.CDS_SPACER, cds_end + variables.CDS_SPACER, gene_score_length) )
                            else:
                                pass 
                        else:                            
                            intervals.append( transcript_fusions.GeneInterval(cds_start - variables.CDS_SPACER, cds_end + variables.CDS_SPACER, gene_score_length) )

        if intervals:
            genes_chosen = transcript_fusions.findMaxGeneScoreJobs(intervals)

            ### First, calculate the coverage of this transcript with genes
            total_transcript_coverage = 0.0
            for (start, end, score) in genes_chosen:          
                isoform = COORDINATES[(start, end, score)]
                total_transcript_coverage += isoform[-1]

            if float(total_transcript_coverage) > len(genes_chosen) * variables.MIN_TRANSCRIPT_COVERAGE:

                ### Now need to add these to the final gene lists!
                for (start, end, score) in genes_chosen:                    
                    isoform = COORDINATES[(start, end, score)]

                    if transcript_id in GENE_ISOFORMS_DIC_CHOSEN:
                        GENE_ISOFORMS_DIC_CHOSEN[transcript_id] = GENE_ISOFORMS_DIC_CHOSEN[transcript_id] + [isoform]
                    else:
                        GENE_ISOFORMS_DIC_CHOSEN[transcript_id] = [isoform]

    else:
        highest_score = 0.0

        for (encoded_protein, gene_score, cds_score, full_id, short_id, transcript_start, transcript_end, cds_start, cds_end, strand, kozak_score, transcript_coverage) in ORFs:

            cds_length = cds_end-cds_start+1
            intron_length = sum([entry[2]-entry[1]+1 for entry in ORF_blocks[full_id] if entry[0] == 'Intron'])

            gene_score_length = scores.weighted_gene_score(cds_length, intron_length, gene_score)    
            
            if highest_score < gene_score_length:
                highest_score = gene_score_length

        for (encoded_protein, gene_score, cds_score, full_id, short_id, transcript_start, transcript_end, cds_start, cds_end, strand, kozak_score, transcript_coverage) in ORFs:

            isoform = (encoded_protein, gene_score, cds_score, full_id, short_id, transcript_start, transcript_end, cds_start, cds_end, strand, kozak_score, transcript_coverage)

            cds_length = cds_end-cds_start+1
            intron_length = sum([entry[2]-entry[1]+1 for entry in ORF_blocks[full_id] if entry[0] == 'Intron'])

            gene_score_length = scores.weighted_gene_score(cds_length, intron_length, gene_score)

            if gene_score_length == highest_score:            

                if float(transcript_coverage) > variables.MIN_TRANSCRIPT_COVERAGE:

                    if transcript_id in GENE_ISOFORMS_DIC_CHOSEN:
                        GENE_ISOFORMS_DIC_CHOSEN[transcript_id] = GENE_ISOFORMS_DIC_CHOSEN[transcript_id] + [isoform]
                    else:
                        GENE_ISOFORMS_DIC_CHOSEN[transcript_id] = [isoform]     

                    break 
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
print("---------------------------------")
print("Now write the gff3 file.")
print("---------------------------------")
# -----------------------------------------------------------------------------------------------------------
out_handle = open(out_file, "w")

transcripts_to_delete = write_gff3.write_gff3_noUTRs_to_file(GENE_ISOFORMS_DIC_CHOSEN, transcript_information, transcripts_that_were_padded, TRANSCRIPT_COORDINATES_WITH_PADDING, ORFs_TO_CLASSIFY, ORF_blocks, cds_scores_effector, out_handle, TRANSCRIPT_IDs_MAPPING, gff3_commandline)
if SIGNALP_PATH != None and TMHMM_PATH != None:
    if out_file.endswith(".gff3"):
        out_handle_secretome = open(out_file.replace(".gff3", ".secretome.gff3"), "w")
    else:
        out_handle_secretome = open(out_file + ".secretome.gff3", "w")

    print("---------------------------------")
    print("Now write the gff3 file for the secretome.")
    print("---------------------------------")
    write_gff3.write_gff3_noUTRs_secretome_to_file(transcripts_to_delete, GENE_ISOFORMS_DIC_CHOSEN, transcript_information, transcripts_that_were_padded, TRANSCRIPT_COORDINATES_WITH_PADDING, ORFs_TO_CLASSIFY, ORF_blocks, cds_scores_effector, out_handle_secretome, SIGNALP, TMHMM, TRANSCRIPT_IDs_MAPPING, gff3_commandline)
# -----------------------------------------------------------------------------------------------------------
print("---------------------------------")
print("All done, bye.")
# Clean up and delete temporary folder that was created
shutil.rmtree(RESULTS_PATH)
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
