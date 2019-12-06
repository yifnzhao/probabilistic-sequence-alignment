#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 20:46:14 2019

@author: yifan
"""


import numpy as np
from re import finditer
from math import log
from queryGenerator import queryGenerator
import time


def loadData(genome_path, prob_path):
    allNu = ['A','C','T','G']
    matrix = dict()
    # read genome data
    f_genome = open(genome_path)
    genome = f_genome.read()
    f_genome.close()
    # read probability data
    f_prob = open(prob_path)
    prob = f_prob.read()
    f_prob.close()
    prob = prob.split()
    genome = list(genome)
    for nu in allNu:
        matrix[nu] = []
    highest_prob_genome = []
    for i in range(len(genome)):
        # update matrix and highestProbGenome
        #  for the nucleotide with highest probability
        maxProb = -log(float(prob[i]))
        matrix[genome[i]].append(maxProb)
        highest_prob_genome.append(genome[i])
        
        # update matrix for the rest
        rest_nu = [n for n in allNu if n!=genome[i]]
        rest_prob = -log((1-maxProb) / float(3)) 
        for n in rest_nu:
            matrix[n].append(rest_prob)
    
    highest_prob_genome_str =''.join(highest_prob_genome)
    return matrix, highest_prob_genome_str

def seedS(querySeq, highest_prob_genome_str, k = 8):
    '''
    k = word size to be searched
    reutrns a list of seed positions [(x,y)，seed]
    '''
    start = time.time()
    seedPosList = []
    for i in range(len(querySeq)-k):
        query = querySeq[i:i+k]
        for match in finditer(query, highest_prob_genome_str):            
            seedPosList.append([match.span(), match.group(), (i, i+k)])
            # i is the start position of the seed in query
    end = time.time()
    timeElapsed = end - start
    return timeElapsed, seedPosList


def seedM_fix_n(query, matrix, k = 8, n = 100):
    '''
    k: word size
    n: number of seeds to look for (fixed)
    returns seedPosList which is a list of seed positions (x,y)
    '''
    start = time.time()
    seedPosList = [0 for i in range(n)]
    seedProbList = [100 for i in range(n)]
    #loop over all k-words in querySeq
    for i in range(len(query)-k): 
        k_word = query[i:i+k]
        #loop over all positions in matrix
        for pos in range(len(matrix['A'])-k): 
            #find probability : thisQuery_prob
            thisQuery_prob = 0    
            for j in range(k): 
                query_n = k_word[j]
                thisQuery_prob += matrix[query_n][pos+j]
                # TODO: argmax is probably the most time consuming step
                # find index of lowest probability 
                index_min = np.argmax(seedProbList) 
                # find the lowest seed probability seen so far
                minSeedProb = seedProbList[index_min] 
            if thisQuery_prob < minSeedProb: # comparing neg log prob
                seedPosList[index_min] = [(pos, pos+k), k_word,  (i, i+k)]
                # i is the start position of the seed in query
                seedProbList[index_min] = thisQuery_prob
    end = time.time()
    timeElapsed = end - start
    return timeElapsed, seedPosList,seedProbList

def seedM_fix_min_prob(query, matrix, k = 8, min_prob = 0.1):
    '''
    k: word size
    n: number of seeds to look for 
    returns seedPosList which is a list of seed positions (x,y)
    '''
    start = time.time()
    seedPosList = []
    seedProbList = []
    for i in range(len(query)-k): #loop over all k-words in querySeq
        k_word = query[i:i+k]
        for pos in range(len(matrix['A'])-k): #loop over all positions in Matrix
            #find probability : thisQuery_prob
            thisQuery_prob = 0    
            for j in range(k): 
                query_n = k_word[j]
                thisQuery_prob += matrix[query_n][pos+j]
            if thisQuery_prob < -log(min_prob): # comparing neg log prob
                seedPosList.append([(pos, pos+k), k_word,  (i, i+k)])
                # i is the start position of the seed in query
                seedProbList.append(thisQuery_prob)
    end = time.time()
    timeElapsed = end - start
    return timeElapsed, seedPosList,seedProbList      


def seed(gpath, ppath, matrix, highest_prob_genome_str, kword = 8):
    '''
    output[0]: querySeeds, a dictionary querySeeds["generationMethod"] = [list of seeds for each query]
    output[1]: the same output as the queryGenerator, contains query and its true prob 
    '''
    query = queryGenerator(gpath, ppath)
    print("Start querying...")
    querySeeds = {}
    for k, v in query.items(): # iterate through each method of query generation
        if k[0] == "q": #check if value is a query list (could be true prob list)
            #k: a str indicates how this list of query sequences was generated, refer to queryGenerator.py
            querySeeds[k] = []
            for query_seq in v:
                seed_dict = {}
                
                t, seedPosListS = seedS(query_seq, highest_prob_genome_str, k=kword)
                seed_dict["S"] = [t, seedPosListS]
                
                t, seedPosListM_fix_n, seedProbListM_fix_n = seedM_fix_n(query_seq, matrix , k=kword)
                seed_dict["M_fix_min_100"] = [t, seedPosListM_fix_n, seedProbListM_fix_n]
                
                t, seedPosListM_fix_min_prob, seedProbListM_fix_min_prob \
                    = seedM_fix_min_prob(query_seq, matrix, k=kword)
                seed_dict["M_fix_min_prob_10"] = [t, seedPosListM_fix_min_prob, seedProbListM_fix_min_prob]

                # -- vary number of seeds (if fix n) --- 
                t, seedPosListM_fix_n, seedProbListM_fix_n = seedM_fix_n(query_seq, matrix, n = 75, k=kword)
                seed_dict["M_fix_min_75"] = [t, seedPosListM_fix_n, seedProbListM_fix_n]
                t, seedPosListM_fix_n, seedProbListM_fix_n = seedM_fix_n(query_seq, matrix, n = 50, k=kword)
                seed_dict["M_fix_min_50"] = [t, seedPosListM_fix_n, seedProbListM_fix_n]
                t, seedPosListM_fix_n, seedProbListM_fix_n = seedM_fix_n(query_seq, matrix, n = 25, k=kword)
                seed_dict["M_fix_min_25"] = [t, seedPosListM_fix_n, seedProbListM_fix_n]
                
                # --- vary prob threshold --- 
                t, seedPosListM_fix_min_prob, seedProbListM_fix_min_prob \
                    = seedM_fix_min_prob(query_seq, matrix, min_prob = 0.05, k=kword)
                seed_dict["M_fix_min_prob_05"] = [t, seedPosListM_fix_min_prob, seedProbListM_fix_min_prob]
                t, seedPosListM_fix_min_prob, seedProbListM_fix_min_prob \
                    = seedM_fix_min_prob(query_seq, matrix, min_prob = 0.03, k=kword)
                seed_dict["M_fix_min_prob_03"] = [t, seedPosListM_fix_min_prob, seedProbListM_fix_min_prob]
                t, seedPosListM_fix_min_prob, seedProbListM_fix_min_prob \
                    = seedM_fix_min_prob(query_seq, matrix, min_prob = 0.01, k=kword)
                seed_dict["M_fix_min_prob_01"] = [t, seedPosListM_fix_min_prob, seedProbListM_fix_min_prob]
            

                querySeeds[k].append([query_seq, seed_dict])
    
    
    return querySeeds, query


if __name__ == "__main__":
    # example usage
    gpath = "../data/chr22.maf.ancestors.42000000.complete.boreo.fa.txt"
    ppath = "../data/chr22.maf.ancestors.42000000.complete.boreo.conf.txt"
    matrix, hpg_str = loadData(gpath, ppath)
    querySeeds, query = seed(gpath, ppath, matrix, hpg_str)
    
#    query_seq = 'AGGAATTCCGAAC' #tester, to be commented
    
    

 


