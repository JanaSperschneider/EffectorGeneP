#!/usr/bin/env python3

"""
EffectorGeneP: gene annotation in pathogen genomes

Copyright (c) 2026 Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.

This work is freely available for non-commercial scientific research, non-commercial education, 
or non-commercial research projects, under the CSIRO Software License Agreement (https://github.com/JanaSperschneider/EffectorGeneP/blob/main/LICENCE). 
"""
# -----------------------------------------------------------------------------------------------------------
from Bio.Seq import Seq
import variables
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
def ORFS_per_transcript(transcript_id, SEQUENCE, exons, ORFs_TO_CLASSIFY, strand, TRANSCRIPT_PADDING, TRANSCRIPT_COORDINATES_WITH_PADDING, transcript_information, MIN_PROTEIN_LENGTH, START_SITE_DISTANCE_MAX):
    # -----------------------------------------------------------------------------------------------------------
    candidate_cds, candidate_utr3s, candidate_utr5s, candidate_kozaks = [], [], [], []
    # -----------------------------------------------------------------------------------------------------------
    # Get the exon coordinates relative to the genome sequence, i.e. with introns present
    # -----------------------------------------------------------------------------------------------------------    
    exon_lengths = []
    # -----------------------------------------------------------------------------------------------------------    
    for (genome_start, genome_end, contig, strand_exon) in exons[transcript_id]:
        exon_lengths.append(genome_end-genome_start+1)
    # -----------------------------------------------------------------------------------------------------------
    exon_positions = [x for x in range(sum(exon_lengths))]
    # -----------------------------------------------------------------------------------------------------------
    # Now produce the corresponding genomic positions
    ##   1.................................................n
    ##   ===========>......==========>..........===========>   exons
    ##          ---->......---------->..........---->          CDS
    # ----------------------------------------------------------------------------------------------------------- 
    ##   n.................................................1       
    ##   <===========......<==========..........<===========   exons
    ##          <----......<----------..........<----          CDS    
    # -----------------------------------------------------------------------------------------------------------        
    genomic_positions, genomic_positions_minus_strand = [], []

    for (genome_start, genome_end, contig, strand_exon) in exons[transcript_id]:
        exon_length = genome_end-genome_start+1
        for x in range(exon_length):
            genomic_positions.append(genome_start+x)

    for (start, end, contig, strand_exon) in reversed(exons[transcript_id]):
        exon_length = end-start+1
        for x in range(exon_length):
            genomic_positions_minus_strand.append(end-x)        
    # -----------------------------------------------------------------------------------------------------------            
    exon_position_dic, exon_position_dic_minus_strand = {}, {}

    for x, y in zip(exon_positions,genomic_positions):
        exon_position_dic[x] = y
    for x, y in zip(exon_positions,genomic_positions_minus_strand):
        exon_position_dic_minus_strand[x] = y
    # -----------------------------------------------------------------------------------------------------------
    possible_ORFs, intervals = ORFs_in_transcript(SEQUENCE, TRANSCRIPT_PADDING, transcript_id, TRANSCRIPT_COORDINATES_WITH_PADDING, transcript_information, MIN_PROTEIN_LENGTH)
    possible_ORFs_filtered = filter_ORFs(possible_ORFs, intervals, START_SITE_DISTANCE_MAX)
    # -----------------------------------------------------------------------------------------------------------
    ORFs_to_be_scored, ORFs_transcript_fusion_to_be_scored = [], []
    # -----------------------------------------------------------------------------------------------------------
    if possible_ORFs_filtered != []:
        # -------------------------------------------------------------------        
        # Append the longest ORF first as the best ORF for this locus
        longest_ORF = sorted(possible_ORFs_filtered)[-1]
        longest_ORF_protein_sequence = longest_ORF[1]
        longest_ORF_start, longest_ORF_stop = longest_ORF[3], longest_ORF[4]
        longest_ORF_transcript_coverage = longest_ORF[8]
        ORFs_to_be_scored.append(longest_ORF)   
        # -------------------------------------------------------------------
        # Also add ORFs that are very similar in length to the best one, best score wins later
        # -------------------------------------------------------------------
        additional_long_ORFs = []

        for entry in possible_ORFs_filtered:
            if entry != longest_ORF:
                protein_length_possible_ORF = entry[0]
                entry_start, entry_end = entry[3], entry[4]

                if 100.0*protein_length_possible_ORF/len(longest_ORF_protein_sequence) > 50.0:
                    additional_long_ORFs.append(entry)   
        # -------------------------------------------------------------------
        ORFs_to_be_scored += additional_long_ORFs
        # -------------------------------------------------------------------
        # Now investigate if there could be a transcript fusion
        # -------------------------------------------------------------------  
        LOOK_LEFT, LOOK_RIGHT = False, False
        # -------------------------------------------------------------------  
        free_bases_left_of_ORF = longest_ORF_start - variables.CDS_SPACER
        free_bases_right_of_ORF = len(SEQUENCE) - longest_ORF_stop + 1 + variables.CDS_SPACER
        # -------------------------------------------------------------------
        # For long transcripts, investigate transcript fusions
        # If the isoform occupies < x% of the transcript, try and add more ORFs
        # -------------------------------------------------------------------   
        if len(SEQUENCE) - 2.0*TRANSCRIPT_PADDING > variables.TRANSCRIPT_MINLENGTH_FUSION and free_bases_left_of_ORF > variables.MIN_TRANSCRIPT_LENGTH:
            LOOK_LEFT = True

        if len(SEQUENCE) - 2.0*TRANSCRIPT_PADDING > variables.TRANSCRIPT_MINLENGTH_FUSION and free_bases_right_of_ORF > variables.MIN_TRANSCRIPT_LENGTH:
            LOOK_RIGHT = True 
        # -------------------------------------------------------------------       
        if LOOK_LEFT == True:
            # Go through the ORFs sorted by coordinate in the transcript, non-overlapping ORFs that can be added will be added
            for entry in possible_ORFs_filtered:
                if entry not in ORFs_to_be_scored:
                    start, stop = entry[3], entry[4]           

                    # Add the ORF if it is located in the left region of the selected ORF
                    # xxxxxxxxxxxxxxxxxxx.......
                    # .....................SooooooooooooooE
                    if start < longest_ORF_start and stop < longest_ORF_start:                    
                        ORFs_transcript_fusion_to_be_scored.append(entry)                            
        # -------------------------------------------------------------------       
        if LOOK_RIGHT == True:
            # Go through the ORFs sorted by coordinate in the transcript, non-overlapping ORFs that can be added will be added
            for entry in possible_ORFs_filtered:                
                if entry not in ORFs_to_be_scored:
                    start, stop = entry[3], entry[4]

                    # Add the ORF if it is located in the right region of the selected ORF
                    # SooooooooooooooE.......
                    # .....................xxxxxxxxxxxxxxxxxxx
                    if start > longest_ORF_stop and stop > longest_ORF_stop:
                        ORFs_transcript_fusion_to_be_scored.append(entry)   
        # -------------------------------------------------------------------   
        # Note that the longest ORF is at position 0, but the others are unordered by size      
        # -------------------------------------------------------------------    
        for index, (protein_length, protein, seq, start, stop, utr3, utr5, kozak, transcript_coverage) in enumerate(ORFs_to_be_scored + ORFs_transcript_fusion_to_be_scored):
            # Now translate the CDS coordinates in the transcript to genomic coordinates
            # Note: these are 1-based coordinates

            if strand == '-':
                # On the minus strand, the start of the ORF is calculated relative to the last exon
                genome_end = exon_position_dic_minus_strand[start] + 1
                genome_start = exon_position_dic_minus_strand[stop] + 1

                identifier = transcript_id.replace('"','') + '_ORF' + str(index) + '_CDS_' + contig + '_' + str(genome_start) + '_' + str(genome_end) + '_' + strand

                candidate_cds.append((identifier, str(seq)))  
                if len(utr3) >= 6.0:
                    candidate_utr3s.append((identifier, utr3))
                if len(utr5) >= 6.0:
                    candidate_utr5s.append((identifier, utr5))
                if len(kozak) >= 6.0:
                    candidate_kozaks.append((identifier, kozak))

                ORFs_TO_CLASSIFY[transcript_id.replace('"','') + '_ORF' + str(index) + '_' + strand] = (contig, strand, protein, utr3, utr5, transcript_coverage)

            if strand == '+':    
                genome_start = exon_position_dic[start] + 1
                genome_end = exon_position_dic[stop] + 1

                identifier = transcript_id.replace('"','') + '_ORF' + str(index) + '_CDS_' + contig + '_' + str(genome_start) + '_' + str(genome_end) + '_' + strand

                candidate_cds.append((identifier, str(seq)))  
                if len(utr3) >= 6.0:
                    candidate_utr3s.append((identifier, utr3))
                if len(utr5) >= 6.0:
                    candidate_utr5s.append((identifier, utr5))
                if len(kozak) >= 6.0:
                    candidate_kozaks.append((identifier, kozak))

                ORFs_TO_CLASSIFY[transcript_id.replace('"','') + '_ORF' + str(index) + '_' + strand] = (contig, strand, protein, utr3, utr5, transcript_coverage)
        # -------------------------------------------------------------------    
        # -------------------------------------------------------------------                         
    else:
        pass

    return ORFs_TO_CLASSIFY, candidate_cds, candidate_utr3s, candidate_utr5s, candidate_kozaks
# -----------------------------------------------------------------------------------------------------------
def ORFs_in_transcript(SEQUENCE, TRANSCRIPT_PADDING, transcript_id, TRANSCRIPT_COORDINATES_WITH_PADDING, transcript_information, MIN_PROTEIN_LENGTH):
    # -----------------------------------------------------------------------------------------------------------
    # Find all potential starts and stops in the transcript, note that these are 0-based
    start_sites = [pos for pos, base in enumerate(SEQUENCE) if SEQUENCE[pos:pos+3] == 'ATG']
    stop_sites = [pos for pos, base in enumerate(SEQUENCE) if SEQUENCE[pos:pos+3] == 'TGA' or SEQUENCE[pos:pos+3] == 'TAA' or SEQUENCE[pos:pos+3] == 'TAG']
    kozak_sequences = []
    # -----------------------------------------------------------------------------------------------------------
    possible_ORFs = []
    intervals = []

    for start in start_sites:
        # Get the ORF until the next stop
        for stop in stop_sites:
            if stop > start and (stop-start)%3 == 0:

                # Translate this CDS sequence
                seq = Seq(SEQUENCE[start:stop+3])
                protein = seq.translate()
                internal_stops = protein[:-1].count('*')   

                padded_length = TRANSCRIPT_COORDINATES_WITH_PADDING[transcript_id][1] - TRANSCRIPT_COORDINATES_WITH_PADDING[transcript_id][0] + 1
                original_length = transcript_information[transcript_id][2] - transcript_information[transcript_id][1] + 1
                TRANSCRIPT_PADDING_APPLIED = padded_length - original_length
                transcript_coverage = 100.0*len(seq)/(len(SEQUENCE) - TRANSCRIPT_PADDING_APPLIED)      

                if internal_stops == 0:
                    if len(protein) >= MIN_PROTEIN_LENGTH:
                        # Extract the 3' UTR for this CDS - do not include the padding             
                        if stop+3 < len(SEQUENCE) - TRANSCRIPT_PADDING:
                            utr3 = SEQUENCE[stop+3:len(SEQUENCE) - TRANSCRIPT_PADDING]
                        else:
                            utr3 = "".join((SEQUENCE[stop+3:]))

                        # Extract the 5' UTR for this CDS - do not include the padding
                        if start > TRANSCRIPT_PADDING:
                            utr5 = SEQUENCE[TRANSCRIPT_PADDING:start]
                        else:
                            utr5 = "".join((SEQUENCE[:start]))

                        # Extract the Kozak sequence for this CDS
                        if start >= 9:
                            kozak = SEQUENCE[start-9:start+4]
                        else:
                            kozak = ""

                        possible_ORFs.append((len(protein), str(protein), seq, start, stop+2, utr3, utr5, kozak, transcript_coverage))
                        intervals.append((start, stop+2))

                else:
                    # If there is a stop in-frame, no need to look for another stop further along the sequence
                    break

    return possible_ORFs, intervals
# -----------------------------------------------------------------------------------------------------------
def filter_ORFs(possible_ORFs, intervals, START_SITE_DISTANCE_MAX):

    # Now delete those smaller ORFs that have the same stop, but a different start site
    # M..........M..........M.......M...............STOP

    possible_ORFs_filtered = []

    all_stops = list(set([stop for start, stop in intervals]))
    longest_intervals = []

    for all_stop in all_stops:
        # Keep the start with the lowest position
        for start, stop in sorted(intervals):
            if all_stop == stop:
                longest_intervals.append((start, all_stop))                
                break

    # Append the longest ORF with that stop first
    for ORF in possible_ORFs:
        if (ORF[3], ORF[4]) in longest_intervals:
            possible_ORFs_filtered.append(ORF)

    # Here keep those ORFs that have a start site within a certain distance
    for ORF in possible_ORFs:
        for lowest_start, highest_stop in longest_intervals:
            if ORF[4] == highest_stop and ORF not in possible_ORFs_filtered:
                # START_SITE_DISTANCE_MAX * 3 as this is now nts
                if ORF[3] - lowest_start > 3.0 and ORF[3] - lowest_start < START_SITE_DISTANCE_MAX*3:
                    possible_ORFs_filtered.append(ORF)

    return possible_ORFs_filtered
# -----------------------------------------------------------------------------------------------------------
