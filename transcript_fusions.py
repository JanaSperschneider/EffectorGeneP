#!/usr/bin/env python3

"""
EffectorGeneP: gene annotation in pathogen genomes

Copyright (c) 2025 Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.

This work is freely available for non-commercial scientific research, non-commercial education, 
or non-commercial research projects, under the CSIRO Non-Commercial License (https://github.com/JanaSperschneider/EffectorGeneP/blob/main/LICENCE). 
"""
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# A class to store a gene interval
class GeneInterval:
    def __init__(self, start, end, gene_score):
        self.start = start
        self.end = end
        self.gene_score = gene_score
# -----------------------------------------------------------------------------------------------------------
# Function to perform a binary search on the given genes, which are sorted
# by end coordinate. The function returns the index of the last gene, which
# doesn't conflict with the given gene, i.e., whose end coordinate is
# less than or equal to the given gene's start coordinate.
def findLastNonConflictingGene(intervals, n):
 
    # search space
    (low, high) = (0, n)
 
    # iterate till the search space is exhausted
    while low <= high:
        mid = (low + high) // 2
        if intervals[mid].end <= intervals[n].start:
            if intervals[mid + 1].end <= intervals[n].start:
                low = mid + 1
            else:
                return mid
        else:
            high = mid - 1
 
    # return the negative index if no non-conflicting gene is found
    return -1

# -----------------------------------------------------------------------------------------------------------
# Function to print the non-overlapping genes involved in score maximization using dynamic programming
def findMaxGeneScoreJobs(intervals):
 
    # base case
    if not intervals:
        return 0
 
    # sort intervals in increasing order of their end coordinates
    intervals.sort(key=lambda x: x.end)
 
    # get the number of intervals
    n = len(intervals)
 
    # maxGeneScore[i] stores the maximum gene scores possible for the first i intervals, and
    # tasks[i] stores the index of genes involved in the maximum gene scores
    maxGeneScore = [None] * n
    tasks = [[] for _ in range(n)]
 
    # initialize maxGeneScore and tasks
    maxGeneScore[0] = intervals[0].gene_score
    tasks[0].append(0)
 
    for i in range(1, n):
 
        # find the index of the last non-conflicting gene with the current gene
        index = findLastNonConflictingGene(intervals, i)
 
        # include the current gene with its non-conflicting intervals
        currentGeneScore = intervals[i].gene_score
        if index != -1:
            currentGeneScore += maxGeneScore[index]
 
        # if including the current gene leads to the maximum gene scores so far
        if maxGeneScore[i - 1] < currentGeneScore:
            maxGeneScore[i] = currentGeneScore
 
            if index != -1:
                tasks[i] = tasks[index][:]
            tasks[i].append(i)
 
        # else if excluding the current gene leads to the maximum gene scores so far
        else:
            tasks[i] = tasks[i - 1][:]
            maxGeneScore[i] = maxGeneScore[i - 1]
 
    genes_chosen = []
    for i in tasks[n - 1]:
        genes_chosen.append((intervals[i].start, intervals[i].end, intervals[i].gene_score))

    return genes_chosen
# -----------------------------------------------------------------------------------------------------------
