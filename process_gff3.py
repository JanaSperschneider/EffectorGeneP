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
# -----------------------------------------------------------------------------------------------------------
def gff3_read_exons(gff3_file_content):
	# -------------------------------------------------------------------------------------------------------	
	# Read the exon information and other transcript-related information from the gff3 file (gffread input format)
	# -------------------------------------------------------------------------------------------------------	
	exons = {}
	transcript_information = {}
	TRANSCRIPT_IDs_MAPPING = {}
	# -------------------------------------------------------------------------------------------------------
	transcript_counter = 0
	# -------------------------------------------------------------------------------------------------------	
	for line in gff3_file_content:
	    if line.startswith('#') or not line.strip():
	        pass
	    else:
	        contig = line.split('\t')[0]
	        ### Note that single-exon genes in stringtie have no strand information!
	        strand = line.split('\t')[6]    
	        identifier_line = line.split('\t')[8]

	        if line.split('\t')[2] == 'transcript' or line.split('\t')[2] == 'mRNA':

	            transcript_counter += 1
	        	                
	            transcript_id = identifier_line.split('ID=')[1].split(';')[0].strip()
	            if 'geneID=' in identifier_line:
	            	gene_id = identifier_line.split('geneID=')[1].split(';')[0].strip()
	            else:
	            	gene_id = transcript_id

	            # Watch out in case transcripts have identical names, should not happen but might	   
	            TRANSCRIPT_IDs_MAPPING[transcript_id + '_UniqueID_' + str(transcript_counter)] = transcript_id
	            transcript_id = transcript_id + '_UniqueID_' + str(transcript_counter)	            

	            # These are 1-based coordinates from gff3
	            transcript_start = int(line.split('\t')[3])
	            transcript_end = int(line.split('\t')[4])
	            transcript_information[transcript_id] = (contig, transcript_start, transcript_end, gene_id, strand)

	        elif line.split('\t')[2] == 'exon' or line.split('\t')[2] == 'CDS':
	            transcript_id = identifier_line.split('Parent=')[1].split(';')[0].strip()
	            # Note: here we convert to 0-based 
	            exon_start = int(line.split()[3]) - 1 
	            exon_end = int(line.split()[4]) - 1 

	            # Watch out in case transcripts have identical names, should not happen but might	
	            TRANSCRIPT_IDs_MAPPING[transcript_id + '_UniqueID_' + str(transcript_counter)] = transcript_id   
	            transcript_id = transcript_id + '_UniqueID_' + str(transcript_counter)

	            if transcript_id in exons:
	                exons[transcript_id] = exons[transcript_id] + [(exon_start, exon_end, contig, strand)]
	            else:
	                exons[transcript_id] = [(exon_start, exon_end, contig, strand)]   
	# -------------------------------------------------------------------------------------------------------
	TRANSCRIPT_COORDINATES_PER_CONTIG = {}
	# -------------------------------------------------------------------------------------------------------		
	for transcript_id, (contig, transcript_start, transcript_end, gene_id, strand) in transcript_information.items():

	    if contig in TRANSCRIPT_COORDINATES_PER_CONTIG:
	        TRANSCRIPT_COORDINATES_PER_CONTIG[contig] = TRANSCRIPT_COORDINATES_PER_CONTIG[contig] + [ [transcript_start, transcript_end, transcript_id] ]
	    else:
	        TRANSCRIPT_COORDINATES_PER_CONTIG[contig] = [ [transcript_start, transcript_end, transcript_id] ]
	# -------------------------------------------------------------------------------------------------------	
	return exons, transcript_information, TRANSCRIPT_COORDINATES_PER_CONTIG, TRANSCRIPT_IDs_MAPPING
# -----------------------------------------------------------------------------------------------------------
def gff3_derive_introns_from_exons(exons):

	introns = {}

	for transcript_id, list_of_exons in exons.items():
	    for index in range(len(list_of_exons) - 1):
	        this_exon, next_exon = list_of_exons[index], list_of_exons[index + 1]
	        contig, strand = this_exon[2], this_exon[3]
	        intron_start, intron_end = this_exon[1] + 1, next_exon[0] - 1

	        if transcript_id in introns:
	            introns[transcript_id] = introns[transcript_id] + [(intron_start, intron_end, contig, strand)]
	        else:   
	            introns[transcript_id] = [(intron_start, intron_end, contig, strand)]

	return introns
# -----------------------------------------------------------------------------------------------------------
def add_padding(TRANSCRIPT_COORDINATES_PER_CONTIG, exons, transcript_information, FORWARD_SEQUENCES, TRANSCRIPT_PADDING):

	TRANSCRIPT_COORDINATES_WITH_PADDING = {}
	exons_with_padding = {}

	for contig, transcript_coordinates in TRANSCRIPT_COORDINATES_PER_CONTIG.items():

	    contig_length = len(FORWARD_SEQUENCES[contig])

	    for transcript in transcript_coordinates:
	    	TRANSCRIPT_COORDINATES_WITH_PADDING[transcript[2]] = max(1, transcript[0] - TRANSCRIPT_PADDING), min(transcript[1] + TRANSCRIPT_PADDING, contig_length)	
	# -----------------------------------------------------------------------------------------------------------
	# Now add the padding to the first exon, and the last exon
	for transcript_id, exon_list in exons.items():

	    gene_id = transcript_id#transcript_information[transcript_id][3]
	    exon_list_with_padding = []

	    # These are the 1-based coordinates with padding added if it was possible
	    if transcript_id in transcript_information:
	        # These are the boundaries in which padding can be applied
	        transcript_start, transcript_end = TRANSCRIPT_COORDINATES_WITH_PADDING[transcript_id][0], TRANSCRIPT_COORDINATES_WITH_PADDING[transcript_id][1] 

	        # Note that the exons are 0-based
	        if len(exon_list) >= 2.0:

	            first_exon = exon_list[0]
	            middle_exons = exon_list[1:-1]
	            last_exon = exon_list[-1]

	            # Get the padded coordinates from the transcript
	            exon_list_with_padding.append( (transcript_start - 1, exon_list[0][1], exon_list[0][2], exon_list[0][3]) )
	            if middle_exons != []:
	                exon_list_with_padding += middle_exons
	            exon_list_with_padding.append( (last_exon[0], transcript_end - 1, exon_list[-1][2], exon_list[-1][3]) )

	        else:
	            exon_list_with_padding.append( (transcript_start - 1, transcript_end - 1,  exon_list[0][2], exon_list[0][3]) )

	        exons_with_padding[transcript_id] = exon_list_with_padding

	return TRANSCRIPT_COORDINATES_WITH_PADDING, exons_with_padding
# -----------------------------------------------------------------------------------------------------------
