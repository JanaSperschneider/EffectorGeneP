#!/usr/bin/env python3

"""
EffectorGeneP: gene annotation in pathogen genomes

Copyright (c) 2025 Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.

This work is freely available for non-commercial scientific research, non-commercial education, 
or non-commercial research projects, under the CSIRO Non-Commercial License (https://github.com/JanaSperschneider/EffectorGeneP/blob/main/LICENCE). 
"""
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# GLOBAL VARIABLES
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# The minimum number of free bases required on the left/right of the best gene for investigating transcript fusion events
MIN_TRANSCRIPT_LENGTH = 450
# The minimum total transcript length that will be investigated for transcript fusion events
TRANSCRIPT_MINLENGTH_FUSION = 1000
# -----------------------------------------------------------------------------------------------------------
# For transcript fusions, there needs to be at least this many bps between the adjacent CDS sequences
CDS_SPACER = 20
# -----------------------------------------------------------------------------------------------------------
MIN_GENE_LENGTH_SCORE_LOG = 5.0
MIN_GENE_LENGTH_SCORE_LOG_TRANSCRIPT_FUSION = 5.7
# -----------------------------------------------------------------------------------------------------------
STRAND_PENALTY = 20.0
MIN_TRANSCRIPT_COVERAGE = 5.0
TRANSCRIPT_FUSION_BEST_GENE_COVERAGE = 75.0
# -----------------------------------------------------------------------------------------------------------
CODING_WEIGHT = 3.0
INTRON_WEIGHT = 6.0
UTR3_WEIGHT = 3.0
UTR5_WEIGHT = 3.0
KOZAK_WEIGHT = 3.0
# -----------------------------------------------------------------------------------------------------------
PVALUE_THRESHOLD = 0.05
KOZAK_MINIMUM = 0.1
SIGNALP_SCORE = 40.0
# -----------------------------------------------------------------------------------------------------------
MIN_EFFECTOR_CDS_SCORE_SP = 0.5
MIN_CDS_SCORE = 0.50
MIN_GENE_SCORE = 50.0
# -----------------------------------------------------------------------------------------------------------
# In conservative mode, these are used
CONSERVATIVE_MIN_EFFECTOR_CDS_SCORE_SP = 0.7
CONSERVATIVE_MIN_CDS_SCORE = 0.60
CONSERVATIVE_MIN_GENE_SCORE = 60.0
# -----------------------------------------------------------------------------------------------------------
