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
import re
import random
import functions
import subprocess as sub

from Bio import SeqIO
from itertools import product

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from Bio.SeqUtils import ProtParam
from itertools import product

import subprocess
import errno
import uuid
import shutil
import tempfile
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
def clean_gff3(content):
    # Input is the EffectorGeneP output GFF3 file
    # Now read in the mRNAs, CDS blocks
    CDS_DIC = {}
    transcripts_per_contig = {}
    transcript_count = 0
    contig_DIC = {}
    # -----------------------------------------------------------------------------------------------------------
    for line in content:
        if line.startswith('#'):
            pass
        else:
            contig = line.split('\t')[0]
            feature = line.split('\t')[2]
            start, end = int(line.split('\t')[3]), int(line.split('\t')[4])
            strand = line.split('\t')[6]
            phase = line.split('\t')[7]
            identifier_line = line.split('\t')[8]

            if feature == 'mRNA':
                transcript_count += 1
                transcript_id = identifier_line.split('ID=')[1].split(';')[0].strip()

                if contig not in transcripts_per_contig:
                    transcripts_per_contig[contig] = [transcript_id]
                else:
                    transcripts_per_contig[contig] = transcripts_per_contig[contig] + [transcript_id]

                if contig not in contig_DIC:
                    contig_DIC[contig] = [transcript_id]
                else:
                    contig_DIC[contig] = contig_DIC[contig] + [transcript_id]   

            if feature == 'CDS':
                transcript_id = identifier_line.split('Parent=')[1].split(';')[0].strip()

                if transcript_id in CDS_DIC:
                    CDS_DIC[transcript_id] = CDS_DIC[transcript_id] + [(contig, start, end, strand)]
                else:
                    CDS_DIC[transcript_id] = [(contig, start, end, strand)]                
    # -----------------------------------------------------------------------------------------------------------
    if transcript_count == 1:
        print("There is", transcript_count, 'gene.')   
    else:  
        print("There are", transcript_count, 'genes.')             
    # -----------------------------------------------------------------------------------------------------------
    transcripts_to_delete = []
    # -----------------------------------------------------------------------------------------------------------
    ### 1st case: find contained CDS - single-exon genes only
    for contig in contig_DIC:

        transcripts = contig_DIC[contig]
        
        for transcript_id in transcripts:
            cds_blocks_of_this_transcript = CDS_DIC[transcript_id]

            if len(cds_blocks_of_this_transcript) == 1.0:

                for other_transcript_id in transcripts:
                    cds_blocks_of_other_transcript = CDS_DIC[other_transcript_id]

                    contained = False

                    for cds_block_of_this_transcript in cds_blocks_of_this_transcript:
                        for cds_block_of_other_transcript in cds_blocks_of_other_transcript:
                            start_this_transcript = cds_block_of_this_transcript[1] 
                            end_this_transcript = cds_block_of_this_transcript[2] 
                            start_other_transcript = cds_block_of_other_transcript[1] 
                            end_other_transcript = cds_block_of_other_transcript[2] 

                            # The CDS of the first transcript is fully contained in another one 
                            if start_this_transcript >= start_other_transcript and end_this_transcript <= end_other_transcript:
                                if start_this_transcript != start_other_transcript or end_this_transcript != end_other_transcript:
                                    contained = True
                                    break

                            # The CDS of the first transcript is partially contained in another one but it is shorter
                            if start_this_transcript > start_other_transcript and start_this_transcript <= end_other_transcript:
                                # ...................ttttttttttttttttttttt...
                                # ....ooooooooooooooooooooooooooooooooo.....
                                if end_this_transcript > end_other_transcript:
                                    if end_this_transcript - start_this_transcript + 1 < end_other_transcript - start_other_transcript + 1:
                                        contained = True
                                        break

                            if start_this_transcript <= start_other_transcript and end_this_transcript < end_other_transcript:
                                # ......ttttttttttttttttttttt...
                                # ........ooooooooooooooooooooooooooooooooo.....
                                if end_this_transcript > start_other_transcript:
                                    if end_this_transcript - start_this_transcript + 1 < end_other_transcript - start_other_transcript + 1:
                                        contained = True
                                        break

                    if contained == True:
                        if transcript_id in transcripts_to_delete:
                            pass
                        else:
                            transcripts_to_delete.append(transcript_id)
    # -----------------------------------------------------------------------------------------------------------
    ### 2nd case: find exact duplicates
    for contig, transcript_list in transcripts_per_contig.items():
        for transcript_id in transcript_list:
            cds_blocks_of_this_transcript = CDS_DIC[transcript_id]

            for other_transcript_id in transcript_list:
                if transcript_id != other_transcript_id:
                    cds_blocks_of_other_transcript = CDS_DIC[other_transcript_id]

                    # 1st case: the two sets are identical
                    if cds_blocks_of_this_transcript == cds_blocks_of_other_transcript:


                        if transcript_id in transcripts_to_delete or other_transcript_id in transcripts_to_delete:
                            pass                        
                        else:
                            transcripts_to_delete.append(transcript_id)


    print("This many genes will be deleted:", len(transcripts_to_delete))

    if transcript_count - len(transcripts_to_delete) > 1:
        print("In total, EffectorGeneP returned", transcript_count - len(transcripts_to_delete), "genes.")
    else:
        print("In total, EffectorGeneP returned", transcript_count - len(transcripts_to_delete), "gene.")

    return transcripts_to_delete
