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
import clean_gff3
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
def write_gff3_noUTRs_to_file(GENE_ISOFORMS_DIC_CHOSEN, transcript_information, transcripts_that_were_padded, TRANSCRIPT_COORDINATES_WITH_PADDING, ORF_CONTIG, ORF_blocks, cds_scores_effector, out_handle, TRANSCRIPT_IDs_MAPPING, gff3_commandline):
	# -------------------------------------------------------------------------------------------------------	
	# Write the EffectorGeneP GFF3 file
	gff3_output_lines = ["##gff-version 3"]	
	gff3_output_lines.append(gff3_commandline)
	# -------------------------------------------------------------------------------------------------------	
	# Now write to gff3 file
	for gene_id, isoforms_to_keep in GENE_ISOFORMS_DIC_CHOSEN.items():

	    unique_proteins_encoded, unique_isoforms = [], []

	    # Find the gene boundaries first to write to the gff3
	    lowest_transcript_start, highest_transcript_end = isoforms_to_keep[0][5], isoforms_to_keep[0][6]

	    for (encoded_protein, gene_score, cds_score, full_id, transcript_id, transcript_start, transcript_end, cds_start, cds_end, strand, kozak_score, transcript_coverage) in isoforms_to_keep:
	        # This includes the padding, which was added to the transcript boundaries & 1st/last exons
	        transcript_start = transcript_information[transcript_id][1] 
	        transcript_end = transcript_information[transcript_id][2]
	        
	        if transcript_start < lowest_transcript_start:
	            lowest_transcript_start = transcript_start

	        if transcript_end > highest_transcript_end:
	            highest_transcript_end = transcript_end

	        # Now check if padding was applied
	        if transcript_id in transcripts_that_were_padded:
	            transcript_start, transcript_end = TRANSCRIPT_COORDINATES_WITH_PADDING[gene_id][0], TRANSCRIPT_COORDINATES_WITH_PADDING[gene_id][1] 
	            lowest_transcript_start, highest_transcript_end = TRANSCRIPT_COORDINATES_WITH_PADDING[gene_id][0], TRANSCRIPT_COORDINATES_WITH_PADDING[gene_id][1] 

	        if encoded_protein not in unique_proteins_encoded:
	            unique_isoforms.append((encoded_protein, gene_score, cds_score, full_id, transcript_id, transcript_start, transcript_end, cds_start, cds_end, strand))
	            unique_proteins_encoded.append(encoded_protein)

	    for (encoded_protein, gene_score, cds_score, full_id, transcript_id, transcript_start, transcript_end, cds_start, cds_end, strand) in unique_isoforms:
	        # -----------------------------------------------------------------------------------------------------------
	        # -----------------------------------------------------------------------------------------------------------
	        # Now write GFF3 output, write the gene/mRNA entry line
	        transcript_gff3 = ""        
	        # -----------------------------------------------------------------------------------------------------------
	        # -----------------------------------------------------------------------------------------------------------
	        # -----------------------------------------------------------------------------------------------------------            
	        cds_score = cds_scores_effector[full_id][0]

	        transcript_id_to_write = TRANSCRIPT_IDs_MAPPING[gene_id]

	        mRNA_line = ORF_CONTIG[full_id][0] + '\t' + "EffectorGeneP" + '\t' + 'mRNA' + '\t' + str(cds_start) + '\t' + str(cds_end) + '\t' + str(cds_score) + '\t' + ORF_CONTIG[full_id][1] + '\t' + '.' + '\t' 
	        mRNA_line += 'ID=' + full_id.replace(gene_id, transcript_id_to_write) + ';geneID=' + transcript_id_to_write #+ ';Note=EffectorGeneP score ' + str(cds_score)

	        gff3_output_lines.append(mRNA_line)
	        # -----------------------------------------------------------------------------------------------------------     
	        # -----------------------------------------------------------------------------------------------------------
	        # Now write the CDS blocks
	        # -----------------------------------------------------------------------------------------------------------
	        # This is for a gene without introns, the phase is 0
	        # -----------------------------------------------------------------------------------------------------------
	        if len(ORF_blocks[full_id]) == 1:

	            index = 1
	            block = ORF_blocks[full_id][0][0]

	            if block == 'CDS':
	                cds_start = cds_start                
	                cds_end = cds_end
	                length = cds_end - cds_start + 1

	                cds_line = ORF_CONTIG[full_id][0] + '\t' + "EffectorGeneP" + '\t' + 'CDS' + '\t' + str(cds_start) + '\t' + str(cds_end) + '\t' + '.' + '\t' + ORF_CONTIG[full_id][1] + '\t' + '0' + '\t' 
	                cds_line += 'ID=cds.' + full_id.replace(gene_id, transcript_id_to_write) + ';Parent=' + full_id.replace(gene_id, transcript_id_to_write)

	                gff3_output_lines.append(cds_line)
	        # -----------------------------------------------------------------------------------------------------------
	        # This is for a gene with introns, the phase needs to be determined
	        # -----------------------------------------------------------------------------------------------------------  
	        # The CDS block can start in any of the exons, at any position. 
	        # After the start, concatenate the remaining exons until the stop is reached.
	        # -----------------------------------------------------------------------------------------------------------           
	        if len(ORF_blocks[full_id]) > 1:
	        	# -----------------------------------------------------------------------------------------------------------           
	            if ORF_CONTIG[full_id][1] == '+':
	        	# -----------------------------------------------------------------------------------------------------------           
	                start = cds_start

	                this_phase, carryover = 0, 0

	                for index, (block, intron_start, intron_end) in enumerate(ORF_blocks[full_id][1:]):

	                    if block == 'Intron':
	                        cds_start = start       
	                        cds_end_temp = intron_start - 1

	                        cds_line = ORF_CONTIG[full_id][0] + '\t' + "EffectorGeneP" + '\t' + 'CDS' + '\t' + str(cds_start) + '\t' + str(cds_end_temp) + '\t' + '.' + '\t' + ORF_CONTIG[full_id][1] + '\t' + str(this_phase) + '\t' 
	                        cds_line += 'ID=cds.' + full_id.replace(gene_id, transcript_id_to_write) + ';Parent=' + full_id.replace(gene_id, transcript_id_to_write)

	                        gff3_output_lines.append(cds_line)
	                        start = intron_end + 1 

	                        carryover = ((cds_end_temp-cds_start+1)%3)
	                        this_phase += (3-carryover)%3
	                        this_phase = this_phase%3

	                cds_line = ORF_CONTIG[full_id][0] + '\t' + "EffectorGeneP" + '\t' + 'CDS' + '\t' + str(intron_end + 1) + '\t' + str(cds_end) + '\t' + '.' + '\t' + ORF_CONTIG[full_id][1] + '\t' + str(this_phase) + '\t' 
	                cds_line += 'ID=cds.' + full_id.replace(gene_id, transcript_id_to_write) + ';Parent=' + full_id.replace(gene_id, transcript_id_to_write)
	                gff3_output_lines.append(cds_line)
	        	# -----------------------------------------------------------------------------------------------------------           
	            # On the minus strand the phase needs to be calculated differently
	            if ORF_CONTIG[full_id][1] == '-':
	        	# -----------------------------------------------------------------------------------------------------------           
	                cds_blocks = []

	                start = cds_start

	                for index, (block, intron_start, intron_end) in enumerate(ORF_blocks[full_id][1:]):

	                    if block == 'Intron':
	                        cds_start = start       
	                        cds_end_temp = intron_start - 1

	                        cds_blocks.append( (ORF_CONTIG[full_id][0], cds_start, cds_end_temp, full_id.replace(gene_id, transcript_id_to_write)) )

	                        start = intron_end + 1 

	                cds_blocks.append( (ORF_CONTIG[full_id][0], intron_end + 1, cds_end, full_id.replace(gene_id, transcript_id_to_write)) )
	            
	                previous_phase = 0
	                this_phase, carryover = 0, 0
	                output_lines = []

	                for (cds_block_contig, cds_block_start, cds_block_end, cds_block_id ) in reversed(cds_blocks):                

	                    cds_line = cds_block_contig + '\t' + "EffectorGeneP" + '\t' + 'CDS' + '\t' + str(cds_block_start) + '\t' + str(cds_block_end) + '\t' + '.' + '\t' + '-' + '\t' + str(this_phase) + '\t' 
	                    cds_line += 'ID=cds.' + cds_block_id + ';Parent=' + cds_block_id
	                    output_lines.append(cds_line)

	                    carryover = ((cds_block_end-cds_block_start+1)%3)
	                    this_phase += (3-carryover)%3
	                    this_phase = this_phase%3

	                for cds_line in reversed(output_lines):
	                    gff3_output_lines.append(cds_line)     

	# -----------------------------------------------------------------------------------------------------------
	# Here, invoke the cleaning of the gff3 function
	transcripts_to_delete = clean_gff3.clean_gff3(gff3_output_lines)
	gff3_output_lines_clean = []
	# -----------------------------------------------------------------------------------------------------------
	for line in gff3_output_lines:
		to_be_deleted = False
		for black_listed_transcript in transcripts_to_delete:
			if black_listed_transcript in line:
				to_be_deleted = True
				break
		if to_be_deleted == False:
			gff3_output_lines_clean.append(line)

	print("EffectorGeneP will now write", len(gff3_output_lines_clean), "lines in the output gff3 file.")
	
	for line in gff3_output_lines_clean:
		out_handle.writelines(line + '\n')

	return transcripts_to_delete
# -----------------------------------------------------------------------------------------------------------
def write_gff3_noUTRs_secretome_to_file(transcripts_to_delete, GENE_ISOFORMS_DIC_CHOSEN, transcript_information, transcripts_that_were_padded, TRANSCRIPT_COORDINATES_WITH_PADDING, ORF_CONTIG, ORF_blocks, cds_scores_effector, out_handle, SIGNALP, TMHMM, TRANSCRIPT_IDs_MAPPING, gff3_commandline):
	# -------------------------------------------------------------------------------------------------------	
	# Write the EffectorGeneP GFF3 file
	gff3_output_lines = ["##gff-version 3"]	
	gff3_output_lines.append(gff3_commandline)
	# -------------------------------------------------------------------------------------------------------	
	# Now write to gff3 file
	for gene_id, isoforms_to_keep in GENE_ISOFORMS_DIC_CHOSEN.items():

	    unique_proteins_encoded, unique_isoforms = [], []

	    for (encoded_protein, gene_score, cds_score, full_id, transcript_id, transcript_start, transcript_end, cds_start, cds_end, strand, kozak_score, transcript_coverage) in isoforms_to_keep:

	        # This includes the padding, which was added to the transcript boundaries & 1st/last exons
	        transcript_start = transcript_information[transcript_id][1] 
	        transcript_end = transcript_information[transcript_id][2]
	        
	        # Now check if padding was applied
	        if transcript_id in transcripts_that_were_padded:
	            transcript_start, transcript_end = TRANSCRIPT_COORDINATES_WITH_PADDING[gene_id][0], TRANSCRIPT_COORDINATES_WITH_PADDING[gene_id][1] 

	        if encoded_protein not in unique_proteins_encoded:
	            unique_isoforms.append((encoded_protein, gene_score, cds_score, full_id, transcript_id, transcript_start, transcript_end, cds_start, cds_end, strand))
	            unique_proteins_encoded.append(encoded_protein)

	    # Is this a secreted protein?
	    secreted = []
	    for (encoded_protein, gene_score, cds_score, full_id, transcript_id, transcript_start, transcript_end, cds_start, cds_end, strand) in unique_isoforms:
	        signalp_score, tmhmm_number = 0.0, 0.0
	        if full_id in SIGNALP:
	            signalp_score = SIGNALP[full_id]

	        if full_id in TMHMM:
	            tmhmm_number = TMHMM[full_id]

	        if signalp_score > 0.340 and tmhmm_number == 0:
	            secreted.append((encoded_protein, gene_score, cds_score, full_id, transcript_id, transcript_start, transcript_end, cds_start, cds_end, strand))

	    for (encoded_protein, gene_score, cds_score, full_id, transcript_id, transcript_start, transcript_end, cds_start, cds_end, strand) in secreted:

	        # -----------------------------------------------------------------------------------------------------------
	        # -----------------------------------------------------------------------------------------------------------
	        # Now write GFF3 output, write the gene/mRNA entry line
	        transcript_gff3 = ""        
	        # -----------------------------------------------------------------------------------------------------------
	        # -----------------------------------------------------------------------------------------------------------
	        # -----------------------------------------------------------------------------------------------------------            
	        cds_score = cds_scores_effector[full_id][0]

	        transcript_id_to_write = TRANSCRIPT_IDs_MAPPING[gene_id]
	        
	        mRNA_line = ORF_CONTIG[full_id][0] + '\t' + "EffectorGeneP" + '\t' + 'mRNA' + '\t' + str(cds_start) + '\t' + str(cds_end) + '\t' + str(cds_score) + '\t' + ORF_CONTIG[full_id][1] + '\t' + '.' + '\t' 
	        mRNA_line += 'ID=' + full_id.replace(gene_id, transcript_id_to_write) + ';geneId=' + transcript_id_to_write #+ ';Note=EffectorGeneP score ' + str(cds_score)

	        gff3_output_lines.append(mRNA_line)
	        # -----------------------------------------------------------------------------------------------------------     
	        # -----------------------------------------------------------------------------------------------------------
	        # Now write the CDS blocks
	        # -----------------------------------------------------------------------------------------------------------
	        # This is for a gene without introns, the phase is 0
	        # -----------------------------------------------------------------------------------------------------------
	        if len(ORF_blocks[full_id]) == 1:

	            index = 1
	            block = ORF_blocks[full_id][0][0]

	            if block == 'CDS':
	                cds_start = cds_start                
	                cds_end = cds_end
	                length = cds_end - cds_start + 1

	                cds_line = ORF_CONTIG[full_id][0] + '\t' + "EffectorGeneP" + '\t' + 'CDS' + '\t' + str(cds_start) + '\t' + str(cds_end) + '\t' + '.' + '\t' + ORF_CONTIG[full_id][1] + '\t' + '0' + '\t' 
	                cds_line += 'ID=cds.' + full_id.replace(gene_id, transcript_id_to_write) + ';Parent=' + full_id.replace(gene_id, transcript_id_to_write)

	                gff3_output_lines.append(cds_line)
	        # -----------------------------------------------------------------------------------------------------------
	        # This is for a gene with introns, the phase needs to be determined
	        # -----------------------------------------------------------------------------------------------------------  
	        # The CDS block can start in any of the exons, at any position. 
	        # After the start, concatenate the remaining exons until the stop is reached.
	        # -----------------------------------------------------------------------------------------------------------           
	        if len(ORF_blocks[full_id]) > 1:
	        	# -----------------------------------------------------------------------------------------------------------           
	            if ORF_CONTIG[full_id][1] == '+':
	        	# -----------------------------------------------------------------------------------------------------------           
	                start = cds_start

	                this_phase, carryover = 0, 0

	                for index, (block, intron_start, intron_end) in enumerate(ORF_blocks[full_id][1:]):

	                    if block == 'Intron':
	                        cds_start = start       
	                        cds_end_temp = intron_start - 1

	                        cds_line = ORF_CONTIG[full_id][0] + '\t' + "EffectorGeneP" + '\t' + 'CDS' + '\t' + str(cds_start) + '\t' + str(cds_end_temp) + '\t' + '.' + '\t' + ORF_CONTIG[full_id][1] + '\t' + str(this_phase) + '\t' 
	                        cds_line += 'ID=cds.' + full_id.replace(gene_id, transcript_id_to_write) + ';Parent=' + full_id.replace(gene_id, transcript_id_to_write)

	                        gff3_output_lines.append(cds_line)
	                        start = intron_end + 1 

	                        carryover = ((cds_end_temp-cds_start+1)%3)
	                        this_phase += (3-carryover)%3
	                        this_phase = this_phase%3

	                cds_line = ORF_CONTIG[full_id][0] + '\t' + "EffectorGeneP" + '\t' + 'CDS' + '\t' + str(intron_end + 1) + '\t' + str(cds_end) + '\t' + '.' + '\t' + ORF_CONTIG[full_id][1] + '\t' + str(this_phase) + '\t' 
	                cds_line += 'ID=cds.' + full_id.replace(gene_id, transcript_id_to_write) + ';Parent=' + full_id.replace(gene_id, transcript_id_to_write)
	                gff3_output_lines.append(cds_line)
	        	# -----------------------------------------------------------------------------------------------------------           
	            # On the minus strand the phase needs to be calculated differently
	            if ORF_CONTIG[full_id][1] == '-':
	        	# -----------------------------------------------------------------------------------------------------------           
	                cds_blocks = []

	                start = cds_start

	                for index, (block, intron_start, intron_end) in enumerate(ORF_blocks[full_id][1:]):

	                    if block == 'Intron':
	                        cds_start = start       
	                        cds_end_temp = intron_start - 1

	                        cds_blocks.append( (ORF_CONTIG[full_id][0], cds_start, cds_end_temp, full_id.replace(gene_id, transcript_id_to_write)) )

	                        start = intron_end + 1 

	                cds_blocks.append( (ORF_CONTIG[full_id][0], intron_end + 1, cds_end, full_id.replace(gene_id, transcript_id_to_write)) )
	            
	                previous_phase = 0
	                this_phase, carryover = 0, 0
	                output_lines = []

	                for (cds_block_contig, cds_block_start, cds_block_end, cds_block_id ) in reversed(cds_blocks):                

	                    cds_line = cds_block_contig + '\t' + "EffectorGeneP" + '\t' + 'CDS' + '\t' + str(cds_block_start) + '\t' + str(cds_block_end) + '\t' + '.' + '\t' + '-' + '\t' + str(this_phase) + '\t' 
	                    cds_line += 'ID=cds.' + cds_block_id + ';Parent=' + cds_block_id
	                    output_lines.append(cds_line)

	                    carryover = ((cds_block_end-cds_block_start+1)%3)
	                    this_phase += (3-carryover)%3
	                    this_phase = this_phase%3

	                for cds_line in reversed(output_lines):
	                    gff3_output_lines.append(cds_line)     

	    # -----------------------------------------------------------------------------------------------------------
	# -----------------------------------------------------------------------------------------------------------
	# -----------------------------------------------------------------------------------------------------------
	# Here, invoke the cleaning of the gff3 function with the transcripts flagged for deletion from whole annotation
	gff3_output_lines_clean = []
	# -----------------------------------------------------------------------------------------------------------
	for line in gff3_output_lines:
		to_be_deleted = False
		for black_listed_transcript in transcripts_to_delete:
			if black_listed_transcript in line:
				to_be_deleted = True
				break
		if to_be_deleted == False:
			gff3_output_lines_clean.append(line)

	print("EffectorGeneP will now write", len(gff3_output_lines_clean), "lines in the output gff3 file.")
	
	for line in gff3_output_lines_clean:
		out_handle.writelines(line + '\n')
		
	return
