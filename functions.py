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
import io
import random
import getopt
from itertools import product
import math 

#from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
# -----------------------------------------------------------------------------------------------------------
def rev_comp(seq):
    
    complement = {'A': 'T', 
                  'C': 'G', 
                  'G': 'C', 
                  'T': 'A',
                  'a': 't', 
                  'c': 'g', 
                  'g': 'c', 
                  't': 'a'}

    # This function returns a reverse complement of a nucleotide sequence
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))

    return reverse_complement
# -----------------------------------------------------------------------------------------------------------
# Copied this BioPython function to ensure consistent results across versions
# -----------------------------------------------------------------------------------------------------------
def lcc_simp(seq):
    """Calculate Local Composition Complexity (LCC) for a sequence.

    seq - an unambiguous DNA sequence (a string or Seq object)

    Returns the Local Composition Complexity (LCC) value for the entire
    sequence (as a float).

    Reference:
    Andrzej K Konopka (2005) Sequence Complexity and Composition
    https://doi.org/10.1038/npg.els.0005260
    """
    wsize = len(seq)
    try:
        # Assume its a string
        upper = seq.upper()
    except AttributeError:
        # Should be a Seq object then
        upper = str(seq).upper()
    l2 = math.log(2)
    if 'A' not in seq:
        term_a = 0
        # Check to avoid calculating the log of 0.
    else:
        term_a = ((upper.count('A')) / float(wsize)) * \
                 ((math.log((upper.count('A')) / float(wsize))) / l2)
    if 'C' not in seq:
        term_c = 0
    else:
        term_c = ((upper.count('C')) / float(wsize)) * \
                 ((math.log((upper.count('C')) / float(wsize))) / l2)
    if 'T' not in seq:
        term_t = 0
    else:
        term_t = ((upper.count('T')) / float(wsize)) * \
                 ((math.log((upper.count('T')) / float(wsize))) / l2)
    if 'G' not in seq:
        term_g = 0
    else:
        term_g = ((upper.count('G')) / float(wsize)) * \
                 ((math.log((upper.count('G')) / float(wsize))) / l2)
    return -(term_a + term_c + term_t + term_g)
# -----------------------------------------------------------------------------------------------------------
def usage():
    """ Function: usage()
        Purpose:  Print helpful information for the user.        
        
        Input:    None.
    
        Return:   Print options for running EffectorGeneP.       
    """
    print('''
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# EffectorGeneP: gene annotation in pathogen genomes
# Copyright (c) 2025 Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.

# This work is freely available for non-commercial scientific research, non-commercial education, 
# or non-commercial research projects, under the CSIRO Non-Commercial License (https://github.com/JanaSperschneider/EffectorGeneP/blob/main/LICENCE). 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ''')
    print("Usage for EffectorGeneP: ")
    print()    
    print("python EffectorGeneP.py [OPTIONS] -g <genome fasta file> -t <transcripts in GFF3 format (gffread formatted)> -m <path to EffectorGeneP model files> -o <EffectorGeneP GFF3 output file>")
    print()
    print("Essential arguments are:")
    print("-g/--genome : genome FASTA file")       
    print("-t/--transcript : assembled transcripts in GFF3 format (have to be gffread formatted)")       
    print("-m/--model_files : path to EffectorGeneP model files, e.g. ./EffectorGeneP_Models/Fusarium_oxysporum_f_sp_lycopersici/")       
    print("-o/--output : EffectorGeneP GFF3 output file")       

    print()
    print("Options are:")
    print("-l/--length <int>: genes have to encode proteins with minimum length of <int> (default: 50 aas)")       
    print("-s/--stringent : assume strand information of transcripts is always correct (not recommended for unstranded RNA-seq data)")   
    print("-p/--padding <int> : add <int> bps to each transcript at both 5' and 3' ends to capture start/stops outside the transcript boundaries. Useful for low coverage data or overlapping UTRs (default: 200)")     
    print("-c/--conservative : annotate genes conservatively")

    print()
    print("to enable sensitive search for genes encoding secreted proteins these two paths need to be provided:")
    print("--SIGNALP4 : path to signalp 4.1 executable")
    print("--TMHMM : path to tmhmm 2.0 executable")
    print()
    print("Advanced options are:")
    print("-d/--distance <int> : how far upstream to look for another translation start site, measured in encoded protein length (default: 70 aas, equates 210 nts)")       
    print("-i/--intron <float> : for a multi-exon gene, average intron probability has to be >= <float> (default: 0.6, range from [0:1])")
    print()
    print("-h/--help : show brief help on usage")

    sys.exit()    
    return
# -----------------------------------------------------------------------------------------------------------
def scan_arguments(commandline):
    """ Function: scan_arguments()
        Purpose:  Scan the input options given to the EffectorGeneP program.        
        
        Input:    Input options given by the user.
    
        Return:   Parsed options.
    """
    try:
        opts, args = getopt.getopt(commandline, "hscl:g:t:m:o:p:d:i:", ["help", "conservative", "length=", "genome=", "transcript=", "model_files=", "output=", "padding=", "stringent", "SIGNALP4=", "TMHMM=", "distance=", "intron="])        
    except getopt.GetoptError as err:
        # print help information and exit:
        print(str(err)) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    GENOME, TRANSCRIPT_GFF3, PATH_TO_MODELS, out_file = None, None, None, None
    SIGNALP4, TMHMM = None, None

    MIN_PROTEIN_LENGTH = 50
    STRANDEDNESS = "NON_STRINGENT"
    TRANSCRIPT_PADDING = 200
    INTRON_MIN_PROB = 0.6
    START_SITE_DISTANCE_MAX = 70 # 70 aas equates 210 nts
    CONSERVATIVE = False

    l_count, g_count, t_count, m_count, o_count, pad_count, s_count = 0, 0, 0, 0, 0, 0, 0
   
    for opt, arg in opts:

        if opt in ['-l', '--length']:
            MIN_PROTEIN_LENGTH = int(arg)
            l_count += 1

        elif opt in ['-g', '--genome']:
            GENOME = arg
            g_count += 1

        elif opt in ['-t', '--transcript']:
            TRANSCRIPT_GFF3 = arg
            t_count += 1

        elif opt in ['-m', '--model_files']:
            PATH_TO_MODELS = arg
            m_count += 1

        elif opt in ['-o', '--output']:
            out_file = arg
            o_count += 1

        elif opt in ['-p', '--padding']:
            TRANSCRIPT_PADDING = int(arg)
            pad_count += 1

        elif opt in ['-s', '--stringent']:
            STRANDEDNESS = "STRINGENT"  
            s_count += 1

        elif opt in ['-c', '--conservative']:
            CONSERVATIVE = True  

        elif opt in ['--SIGNALP4']:
            SIGNALP4 = arg

        elif opt in ['--TMHMM']:
            TMHMM = arg

        elif opt in ['-d', '--distance']:
            START_SITE_DISTANCE_MAX = int(arg)
            if START_SITE_DISTANCE_MAX < 0:
                print()
                print ("Option -d/--distance has to be a number >= 0.")
                usage()

        elif opt in ['-i', '--intron']:
            INTRON_MIN_PROB = float(arg)
            if INTRON_MIN_PROB < 0 or INTRON_MIN_PROB > 1.0:
                print()
                print ("Option -i/--intron has to be a float in the range from [0:1].")
                usage()

        elif opt in ['-h', '--help']:
            usage()

        else:
            print()
            print ("Commandline option was supplied that was not recognized!")
            usage()

    if g_count > 1 or t_count > 1 or m_count > 1 or o_count > 1 or l_count > 1 or pad_count > 1 or s_count > 1:
       usage()

    return GENOME, TRANSCRIPT_GFF3, PATH_TO_MODELS, out_file, MIN_PROTEIN_LENGTH, STRANDEDNESS, TRANSCRIPT_PADDING, SIGNALP4, TMHMM, START_SITE_DISTANCE_MAX, INTRON_MIN_PROB, CONSERVATIVE
# -----------------------------------------------------------------------------------------------------------
ARFF_HEADER = '''@RELATION effectors
@ATTRIBUTE GC NUMERIC
@ATTRIBUTE THIRDA NUMERIC
@ATTRIBUTE THIRDT NUMERIC
@ATTRIBUTE THIRDG NUMERIC
@ATTRIBUTE THIRDC NUMERIC
@ATTRIBUTE HEXAMERSCORE_INTERGENIC NUMERIC
@ATTRIBUTE HEXAMERSCORE_INTRON NUMERIC
@ATTRIBUTE HEXAMERSCORE_CDS_SECRETED NUMERIC
@ATTRIBUTE HEXAMERSCORE_CDS_NOTSECRETED NUMERIC
@ATTRIBUTE HEXAMERSCORE_CDS_EFFECTORS NUMERIC
@ATTRIBUTE HEXAMERSCORE_CDS_UTR3 NUMERIC
@ATTRIBUTE HEXAMERSCORE_CDS_UTR5 NUMERIC
@ATTRIBUTE CODONSCORE_INTERGENIC NUMERIC
@ATTRIBUTE CODONSCORE_INTRON NUMERIC
@ATTRIBUTE CODONSCORE_CDS_SECRETED NUMERIC
@ATTRIBUTE CODONSCORE_CDS_NOTSECRETED NUMERIC
@ATTRIBUTE CODONSCORE_CDS_EFFECTORS NUMERIC
@ATTRIBUTE CODONSCORE_CDS_UTR3 NUMERIC
@ATTRIBUTE CODONSCORE_CDS_UTR5 NUMERIC
@ATTRIBUTE MELTING NUMERIC
@ATTRIBUTE LENGTH NUMERIC
@ATTRIBUTE A NUMERIC
@ATTRIBUTE C NUMERIC
@ATTRIBUTE D NUMERIC
@ATTRIBUTE E NUMERIC
@ATTRIBUTE F NUMERIC
@ATTRIBUTE G NUMERIC
@ATTRIBUTE H NUMERIC
@ATTRIBUTE I NUMERIC
@ATTRIBUTE K NUMERIC
@ATTRIBUTE L NUMERIC
@ATTRIBUTE M NUMERIC
@ATTRIBUTE N NUMERIC
@ATTRIBUTE P NUMERIC
@ATTRIBUTE Q NUMERIC
@ATTRIBUTE R NUMERIC
@ATTRIBUTE S NUMERIC
@ATTRIBUTE T NUMERIC
@ATTRIBUTE V NUMERIC
@ATTRIBUTE W NUMERIC
@ATTRIBUTE Y NUMERIC
@ATTRIBUTE LCC NUMERIC
@ATTRIBUTE CODNOTSEC NUMERIC
@ATTRIBUTE CODSEC NUMERIC
@ATTRIBUTE CODEFF NUMERIC
@ATTRIBUTE class {pos,neg}
@DATA
'''
# -----------------------------------------------------------------------------------------------------------
ARFF_HEADER_KOZAK = '''@RELATION effectors
@ATTRIBUTE GC NUMERIC
@ATTRIBUTE base9 NUMERIC
@ATTRIBUTE base8 NUMERIC
@ATTRIBUTE base7 NUMERIC
@ATTRIBUTE base6 NUMERIC
@ATTRIBUTE base5 NUMERIC
@ATTRIBUTE base4 NUMERIC
@ATTRIBUTE base3 NUMERIC
@ATTRIBUTE base2 NUMERIC
@ATTRIBUTE base1 NUMERIC
@ATTRIBUTE baseafterM NUMERIC
@ATTRIBUTE class {pos,neg}
@DATA
'''
# -----------------------------------------------------------------------------------------------------------
SynonymousCodons = {
    'C': ['TGT', 'TGC'],
    'D': ['GAT', 'GAC'],
    'S': ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT'],
    'Q': ['CAA', 'CAG'],
    'N': ['AAC', 'AAT'],
    'P': ['CCT', 'CCG', 'CCA', 'CCC'],
    'K': ['AAG', 'AAA'],
    'T': ['ACC', 'ACA', 'ACG', 'ACT'],
    'F': ['TTT', 'TTC'],
    'A': ['GCA', 'GCC', 'GCG', 'GCT'],
    'G': ['GGT', 'GGG', 'GGA', 'GGC'],
    'I': ['ATC', 'ATA', 'ATT'],
    'L': ['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA'],
    'H': ['CAT', 'CAC'],
    'R': ['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA'],
    'W': ['TGG'],
    'V': ['GTA', 'GTC', 'GTG', 'GTT'],
    'E': ['GAG', 'GAA'],
    'Y': ['TAT', 'TAC']}
# -----------------------------------------------------------------------------------------------------------
def read_codon_usage_dic(INPUT_FILE):

    CodonsDict_Usage = {}

    f = open(INPUT_FILE)
    content = f.readlines()
    f.close()

    # Format is: CodonsDict_Usage[amino_acid] = freq_of_codons
    for line in content:
        amino_acid = line.split()[0]
        freq = float(line.split()[1])

        if amino_acid in CodonsDict_Usage:
            CodonsDict_Usage[amino_acid] = CodonsDict_Usage[amino_acid] + [freq]
        else:
            CodonsDict_Usage[amino_acid] = [freq]

    return CodonsDict_Usage
# -----------------------------------------------------------------------------------------------------------
def codon_usage_deviation(sequence, CodonsDict_Usage):

    codons_in_sequence = []
    amino_acid_deviation = []

    for i in range(0, len(sequence.upper()), 3):
        codons_in_sequence.append(sequence.upper()[i:i + 3])

    for amino_acid, list_of_codons in SynonymousCodons.items():
        number_of_codons= []
        for codon in list_of_codons:
            number_of_codons.append(codons_in_sequence.count(codon))
        freq_of_codons = [x/sum(number_of_codons) if sum(number_of_codons) > 0 else 0.0 for x in number_of_codons]

        for i in range(len(freq_of_codons)):
            if CodonsDict_Usage[amino_acid][i] > 0.0:
                amino_acid_deviation.append(freq_of_codons[i] / CodonsDict_Usage[amino_acid][i])

    average_amino_acid_deviation = [x for x in amino_acid_deviation if x > 0.0]

    if len(average_amino_acid_deviation) > 0.0:
        deviation = sum(average_amino_acid_deviation)/len(average_amino_acid_deviation)
    else: 
        deviation = 0.0

    return deviation
# -----------------------------------------------------------------------------------------------------------
def kmer_frequencies(sequence, word_size, kmer_frequencies_intergenic, kmer_frequencies_intron, kmer_frequencies_CDS_secreted, kmer_frequencies_CDS_notsecreted, kmer_frequencies_CDS_effectors, kmer_frequencies_UTRs3, kmer_frequencies_UTRs5):

    # All possible hexamers
    all_kmers = [''.join(c) for c in product('ACGT', repeat=word_size)]
    kmer_dic = {}

    for position in range(0, len(sequence), 1):
            kmer = sequence[position:position+word_size]

            if len(kmer) == word_size:
                if kmer in kmer_dic:
                    kmer_dic[kmer] += 1
                else:
                    kmer_dic[kmer] = 1

    # Total number of k-mers
    sum_of_counts = sum([count for kmer, count in kmer_dic.items()])

    kmer_score_intergenic, kmer_score_intron, kmer_score_CDS_secreted, kmer_score_CDS_notsecreted, kmer_score_CDS_effectors = 0.0, 0.0, 0.0, 0.0, 0.0
    kmer_score_UTRs3, kmer_score_UTRs5 = 0.0, 0.0

    for kmer in all_kmers:
        if kmer in kmer_dic:

            # It can happen that the kmer has no counts, thus handle division by zero here!
            if kmer_frequencies_intergenic[kmer] != 0.0:
                kmer_score_intergenic += (kmer_dic[kmer]/sum_of_counts)/kmer_frequencies_intergenic[kmer] 
            if kmer_frequencies_intron[kmer] != 0.0:
                kmer_score_intron += (kmer_dic[kmer]/sum_of_counts)/kmer_frequencies_intron[kmer] 
            if kmer_frequencies_CDS_secreted[kmer] != 0.0:                
                kmer_score_CDS_secreted += (kmer_dic[kmer]/sum_of_counts)/kmer_frequencies_CDS_secreted[kmer] 
            if kmer_frequencies_CDS_notsecreted[kmer] != 0.0:          
                kmer_score_CDS_notsecreted += (kmer_dic[kmer]/sum_of_counts)/kmer_frequencies_CDS_notsecreted[kmer] 
            if kmer_frequencies_CDS_effectors[kmer] != 0.0:         
                kmer_score_CDS_effectors += (kmer_dic[kmer]/sum_of_counts)/kmer_frequencies_CDS_effectors[kmer] 
            if kmer_frequencies_UTRs3[kmer] != 0.0:         
               kmer_score_UTRs3 += (kmer_dic[kmer]/sum_of_counts)/kmer_frequencies_UTRs3[kmer] 
            if kmer_frequencies_UTRs5[kmer] != 0.0:         
                kmer_score_UTRs5 += (kmer_dic[kmer]/sum_of_counts)/kmer_frequencies_UTRs5[kmer] 

    kmer_score_intergenic = float(kmer_score_intergenic/len(all_kmers))
    kmer_score_intron = float(kmer_score_intron/len(all_kmers))
    kmer_score_CDS_secreted = float(kmer_score_CDS_secreted/len(all_kmers))
    kmer_score_CDS_notsecreted = float(kmer_score_CDS_notsecreted/len(all_kmers))
    kmer_score_CDS_effectors = float(kmer_score_CDS_effectors/len(all_kmers))
    kmer_score_UTRs3 = float(kmer_score_UTRs3/len(all_kmers))
    kmer_score_UTRs5 = float(kmer_score_UTRs5/len(all_kmers))

    return kmer_score_intergenic, kmer_score_intron, kmer_score_CDS_secreted, kmer_score_CDS_notsecreted, kmer_score_CDS_effectors, kmer_score_UTRs3, kmer_score_UTRs5
# -----------------------------------------------------------------------------------------------------------
def third_base_in_codon(sequence):

    sequence = sequence.upper()

    third_codon_A, third_codon_T, third_codon_G, third_codon_C = 0.0, 0.0, 0.0, 0.0  

    for position in range(2, len(sequence), 3):
        third_codon = sequence[position]  
        if third_codon == 'A':
            third_codon_A += 1.0
        if third_codon == 'T':
            third_codon_T += 1.0
        if third_codon == 'G':
            third_codon_G += 1.0
        if third_codon == 'C':
            third_codon_C += 1.0     

    if (third_codon_A + third_codon_T + third_codon_G + third_codon_C) == 0.0:
        third_codon_A_freq, third_codon_T_freq, third_codon_G_freq, third_codon_C_freq  = 0.0, 0.0, 0.0, 0.0  
    else:
        third_codon_A_freq = third_codon_A/(third_codon_A + third_codon_T + third_codon_G + third_codon_C)                       
        third_codon_T_freq = third_codon_T/(third_codon_A + third_codon_T + third_codon_G + third_codon_C)                       
        third_codon_G_freq = third_codon_G/(third_codon_A + third_codon_T + third_codon_G + third_codon_C)                       
        third_codon_C_freq = third_codon_C/(third_codon_A + third_codon_T + third_codon_G + third_codon_C)    
    
    return third_codon_A_freq, third_codon_T_freq, third_codon_G_freq, third_codon_C_freq 
# -----------------------------------------------------------------------------------------------------------
def GC_content(sequence):

    sequence = sequence.upper()
    
    if ( sequence.count('A') + sequence.count('T') + sequence.count('G') + sequence.count('C') ) == 0.0:
        gc_content = 0.0
    else:
        gc_content = (sequence.count('G') + sequence.count('C'))/(sequence.count('A') + sequence.count('T') + sequence.count('G') + sequence.count('C'))

    return gc_content
# -----------------------------------------------------------------------------------------------------------
def amino_acid_frequencies(protein_translation):

	# Amino acid frequencies in the sequence
	amino_acid_frequencies = []
	amino_acid_frequencies.append(100.0*protein_translation.count('A')/len(protein_translation))
	amino_acid_frequencies.append(100.0*protein_translation.count('C')/len(protein_translation))
	amino_acid_frequencies.append(100.0*protein_translation.count('D')/len(protein_translation))
	amino_acid_frequencies.append(100.0*protein_translation.count('E')/len(protein_translation))
	amino_acid_frequencies.append(100.0*protein_translation.count('F')/len(protein_translation))
	amino_acid_frequencies.append(100.0*protein_translation.count('G')/len(protein_translation))
	amino_acid_frequencies.append(100.0*protein_translation.count('H')/len(protein_translation))
	amino_acid_frequencies.append(100.0*protein_translation.count('I')/len(protein_translation))
	amino_acid_frequencies.append(100.0*protein_translation.count('K')/len(protein_translation))
	amino_acid_frequencies.append(100.0*protein_translation.count('L')/len(protein_translation))
	amino_acid_frequencies.append(100.0*protein_translation.count('M')/len(protein_translation))
	amino_acid_frequencies.append(100.0*protein_translation.count('N')/len(protein_translation))
	amino_acid_frequencies.append(100.0*protein_translation.count('P')/len(protein_translation))
	amino_acid_frequencies.append(100.0*protein_translation.count('Q')/len(protein_translation))
	amino_acid_frequencies.append(100.0*protein_translation.count('R')/len(protein_translation))
	amino_acid_frequencies.append(100.0*protein_translation.count('S')/len(protein_translation))
	amino_acid_frequencies.append(100.0*protein_translation.count('T')/len(protein_translation))
	amino_acid_frequencies.append(100.0*protein_translation.count('V')/len(protein_translation))
	amino_acid_frequencies.append(100.0*protein_translation.count('W')/len(protein_translation))
	amino_acid_frequencies.append(100.0*protein_translation.count('Y')/len(protein_translation))

	return amino_acid_frequencies
# -----------------------------------------------------------------------------------------------------------
def read_kmer_frequencies_genome(input_file, word_size):

    f = open(input_file, 'r')
    content = f.readlines()
    f.close()

    hexamer_frequencies = {}

    all_kmers = [''.join(c) for c in product('ACGT', repeat=word_size)]

    for line in content:

        if line.startswith('#'):
            pass
        else:
            try:
                if line.split()[0] in all_kmers:
                    hexamer_frequencies[line.split()[0]] = float(line.split()[2])
            except:
                pass

    return hexamer_frequencies
# -----------------------------------------------------------------------------------------------------------
def SimpleFastaParser(handle):
    for line in handle:
        if line[0] == ">":
            title = line[1:].rstrip()
            break

    lines = []
    for line in handle:
        if line[0] == ">":
            yield title, "".join(lines).replace(" ", "").replace("\r", "")
            lines = []
            title = line[1:].rstrip()
            continue
        lines.append(line.rstrip())

    yield title, "".join(lines).replace(" ", "").replace("\r", "")
# -----------------------------------------------------------------------------------------------------------   
def translate(seq):
	
	table = {
		'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
		'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
		'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
		'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',				
		'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
		'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
		'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
		'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
		'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
		'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
		'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
		'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
		'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
		'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
		'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
		'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
	}
	protein =""
	if len(seq)%3 == 0:
		for i in range(0, len(seq), 3):
			codon = seq[i:i + 3]
			try:
				protein += table[codon]
			except:
				protein += '_'
	else:
		for i in range(0, len(seq), 3):
			codon = seq[i:i + 3]
			try:
				protein += table[codon]
			except:
				protein += '_'

	return protein
# -----------------------------------------------------------------------------------------------------------   
def write_weka_input(weka_input, SEQUENCES, hexamer_frequencies_intergenic, hexamer_frequencies_intron, hexamer_frequencies_CDS_secreted, hexamer_frequencies_CDS_notsecreted, hexamer_frequencies_CDS_effectors, hexamer_frequencies_UTRs3, hexamer_frequencies_UTRs5, codon_frequencies_intergenic, codon_frequencies_intron, codon_frequencies_CDS_secreted, codon_frequencies_CDS_notsecreted, codon_frequencies_CDS_effectors, codon_frequencies_UTRs3, codon_frequencies_UTRs5, CodonsDict_Usage_NotSecreted, CodonsDict_Usage_Secreted, CodonsDict_Usage_Effectors):
    with open(weka_input, 'w') as f:
        # Create a list of features for each protein
        X = [[] for __ in range(len(SEQUENCES))]

        for protein_position, sequence in enumerate(SEQUENCES):

            # If the sequence contains 'N', replace with random selection of ACGT
            letters = 'ACGT'
            sequence = ''.join(random.choice(letters) if c == 'N' else c for c in sequence)

            # Translate sequence into mature protein
            coding_dna = Seq(sequence)
            protein_translation = translate(sequence)
            # -----------------------------------------------------------------------------------------------------------
            # Amino acid frequencies in the sequence
            amino_acid_frequencies_list = amino_acid_frequencies(protein_translation)
            # -----------------------------------------------------------------------------------------------------------
            # GC content
            gc_content = GC_content(sequence)
            # -----------------------------------------------------------------------------------------------------------
            # Melting temperature
            melting_temperature = mt.Tm_NN(sequence)
            # -----------------------------------------------------------------------------------------------------------
            # Local composition complexity
            lcc = lcc_simp(sequence)
            # -----------------------------------------------------------------------------------------------------------
            # Frequency of third base in codons
            third_codon_A_freq, third_codon_T_freq, third_codon_G_freq, third_codon_C_freq = third_base_in_codon(sequence)  
            # Hexamer frequencies in the sequence
            kmer6_intergenic, kmer6_intron, kmer6_CDS_secreted, kmer6_CDS_notsecreted, kmer6_CDS_effectors, kmer6_score_UTRs3, kmer6_score_UTRs5 = kmer_frequencies(sequence, 6, hexamer_frequencies_intergenic, hexamer_frequencies_intron, hexamer_frequencies_CDS_secreted, hexamer_frequencies_CDS_notsecreted, hexamer_frequencies_CDS_effectors, hexamer_frequencies_UTRs3, hexamer_frequencies_UTRs5)
            # -----------------------------------------------------------------------------------------------------------
            # Codon frequencies in the sequence
            kmer3_intergenic, kmer3_intron, kmer3_CDS_secreted, kmer3_CDS_notsecreted, kmer3_CDS_effectors, kmer3_score_UTRs3, kmer3_score_UTRs5 = kmer_frequencies(sequence, 3, codon_frequencies_intergenic, codon_frequencies_intron, codon_frequencies_CDS_secreted, codon_frequencies_CDS_notsecreted, codon_frequencies_CDS_effectors, codon_frequencies_UTRs3, codon_frequencies_UTRs5)
            # -----------------------------------------------------------------------------------------------------------
            # Codon usage
            codon_deviation_not_secreted = codon_usage_deviation(sequence, CodonsDict_Usage_NotSecreted)
            codon_deviation_secreted = codon_usage_deviation(sequence, CodonsDict_Usage_Secreted)
            codon_deviation_effector = codon_usage_deviation(sequence, CodonsDict_Usage_Effectors)            
            # -----------------------------------------------------------------------------------------------------------
            feature_vector = [gc_content, third_codon_A_freq, third_codon_T_freq, third_codon_G_freq, third_codon_C_freq]
            feature_vector += [kmer6_intergenic, kmer6_intron, kmer6_CDS_secreted, kmer6_CDS_notsecreted, kmer6_CDS_effectors, kmer6_score_UTRs3, kmer6_score_UTRs5]     
            feature_vector += [kmer3_intergenic, kmer3_intron, kmer3_CDS_secreted, kmer3_CDS_notsecreted, kmer3_CDS_effectors, kmer3_score_UTRs3, kmer3_score_UTRs5]  
            feature_vector += [melting_temperature, len(sequence.strip())]
            feature_vector += amino_acid_frequencies_list + [lcc] + [codon_deviation_not_secreted, codon_deviation_secreted, codon_deviation_effector]
             
            X[protein_position] = feature_vector

        # Write protein feature data to WEKA arff file 
        f.writelines(ARFF_HEADER)
        for index, vector in enumerate(X):
            for feature in vector:
                f.writelines(str(feature) + ',')
            f.writelines('?\n')

    return
# -----------------------------------------------------------------------------------------------------------
def base_to_number(base):

    if base.upper() == 'A':
        number = 1
    elif base.upper() == 'C':
        number = 2
    elif base.upper() == 'G':
        number = 3
    elif base.upper() == 'T':
        number = 4
    else:
        number = 0 
    
    return number
# -----------------------------------------------------------------------------------------------------------   
def write_weka_input_kozak(weka_input, SEQUENCES):
    with open(weka_input, 'w') as f:
        # Create a list of features for each protein
        X = [[] for __ in range(len(SEQUENCES))]

        for protein_position, sequence in enumerate(SEQUENCES):

            # If the sequence contains 'N', replace with random selection of ACGT
            letters = 'ACGT'
            sequence = ''.join(random.choice(letters) if c == 'N' else c for c in sequence)
            # -----------------------------------------------------------------------------------------------------------
            # GC content
            gc_content = GC_content(sequence)

            if len(sequence) >= 13.0:
                base9 = base_to_number(sequence[0])
                base8 = base_to_number(sequence[1])
                base7 = base_to_number(sequence[2])
                base6 = base_to_number(sequence[3]) 
                base5 = base_to_number(sequence[4]) 
                base4 = base_to_number(sequence[5]) 
                base3 = base_to_number(sequence[6]) 
                base2 = base_to_number(sequence[7]) 
                base1 = base_to_number(sequence[8]) 
                baseafterM = base_to_number(sequence[12])

                feature_vector = [gc_content, base9, base8, base7, base6, base5, base4, base3, base2, base1, baseafterM]    
            else:
                feature_vector = [gc_content, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]    


            X[protein_position] = feature_vector

        # Write protein feature data to WEKA arff file 
        f.writelines(ARFF_HEADER_KOZAK)
        for index, vector in enumerate(X):
            for feature in vector:
                f.writelines(str(feature) + ',')
            f.writelines('?\n')

    return
# -----------------------------------------------------------------------------------------------------------    