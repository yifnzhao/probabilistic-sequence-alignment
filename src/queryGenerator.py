#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 15:41:46 2019

@author: yifan
"""

import random
from datetime import datetime
from math import log
import pickle

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

def findProb(matrix, seedStart, seed):
    '''
    seedStart: starting position of the seed wrt matrix
    output is a probability of the seed aligning to the matrix gaplessly
    '''
    prob = 0
    for index in range(len(seed)):
        n = seed[index]
        prob += matrix[n][seedStart]
        seedStart+=1
    return prob

def generateHighestScoreQuery(gpath, ppath, n_query=5, length = -1):
    random.seed(datetime.now())
    matrix, highest_prob_genome_str = loadData(gpath, ppath)
    len_genome = len(highest_prob_genome_str)
    queryList = []
    trueProbList = []
    for i in range(n_query):
        start_pos = random.randint(0, len_genome-1)
        maxLen = len_genome - start_pos
        # -- find query length ---
        if length == -1: # unspecified
            query_len = random.randint(6, 200) #TODO: 200 is chosen arbitrarily, what would be a more statistically relevant value?
        else:
            query_len = length
        if query_len > maxLen:
            continue
        query = highest_prob_genome_str[start_pos: start_pos+query_len]
        # --- find true prob ---
        prob = findProb(matrix, start_pos, query)
        queryList.append(query)
        trueProbList.append(prob)
    return queryList, trueProbList
    

def generateQueryGapless(gpath, ppath, weight = 0.75, n_query=5, length=-1):
    # 75% probability select highest scoring nucleotide by default 
    random.seed(datetime.now())
    allNu = ['A','C','T','G']
    matrix, highest_prob_genome_str = loadData(gpath, ppath)
    len_genome = len(highest_prob_genome_str)
    queryList = []
    trueProbList = []
    for i in range(n_query):        
        start_pos = random.randint(0, len_genome-1)
        maxLen = len_genome - start_pos
        # -- find query length ---
        if length == -1: # unspecified
            query_len = random.randint(6, 200) #TODO: 200 is chosen arbitrarily, what would be a more statistically relevant value?
        else:
            query_len = length
        if query_len > maxLen:
            continue        
        query = ''
        for j in range(query_len):
            sample = random.uniform(0,1)
            highestScoreNu = highest_prob_genome_str[start_pos + j]
            otherNu = [n for n in allNu if n!=highestScoreNu]
            if sample > weight: # do random nucleotide selection
                query+=random.choice(otherNu)
            else:
                query+=highestScoreNu
        # --- find true prob ---
        prob = findProb(matrix, start_pos, query)
        queryList.append(query)
        trueProbList.append(prob)
    return queryList, trueProbList


def generateQueryGapless_shortQueryHighScore(gpath, ppath, n_query=5, length=-1):
    random.seed(datetime.now())
    allNu = ['A','C','T','G']
    matrix, highest_prob_genome_str = loadData(gpath, ppath)
    len_genome = len(highest_prob_genome_str)
    queryList = []
    trueProbList = []
    for i in range(n_query):        
        start_pos = random.randint(0, len_genome-1)
        maxLen = len_genome - start_pos
        # -- find query length ---
        if length == -1: # unspecified
            query_len = random.randint(6, 200) #TODO: 200 is chosen arbitrarily, what would be a more statistically relevant value?
        else:
            query_len = length
        if query_len > maxLen:
            continue        
        # -- find weight of choosing most probable nu ---
        weight = 1/float(query_len) #weight is inversely proportional to query len
        query = ''
        for j in range(query_len):
            sample = random.uniform(0,1)
            highestScoreNu = highest_prob_genome_str[start_pos + j]
            otherNu = [n for n in allNu if n!=highestScoreNu]
            if sample > weight: # do random nucleotide selection
                query+=random.choice(otherNu)
            else:
                query+=highestScoreNu
        # --- find true prob ---
        prob = findProb(matrix, start_pos, query)
        queryList.append(query)
        trueProbList.append(prob)
    return queryList, trueProbList

def queryGenerator(gpath, ppath):
    '''
    q_100_rand: "q" for query list, 
                "100" specifies weight (i.e., 1.0 possibility of choosing most probable nu)
                "rand" for random, specifies length (i.e. random length)
    p_100_rand: "p" for true probability list
    
    returns a dict of all queries and true probabilities
    
    '''
    query = {}
    
    # -- generate 5 highest scoring query sequences of variable lengths b/t 6 and 200 ---

    q_100_rand, p_100_rand = generateHighestScoreQuery(gpath, ppath)
    query["q_100_rand"] = q_100_rand
    query["p_100_rand"] = p_100_rand
    # -- generate 5 highest scoring query sequences of fixed lengths b/t 6 and 200 --- 
    q_100_25, p_100_25 = generateHighestScoreQuery(gpath, ppath, length = 25)
    q_100_50, p_100_50 = generateHighestScoreQuery(gpath, ppath, length = 50)
    q_100_200, p_100_200 = generateHighestScoreQuery(gpath, ppath, length = 200)
    query["q_100_25"] = q_100_25
    query["q_100_50d"] = q_100_50
    query["q_100_200"] = q_100_200
    query["p_100_25"] = p_100_25
    query["p_100_50"] = p_100_50
    query["p_100_200"] = p_100_200
    
    # -- generate 5 query sequences of variable lengths b/t 6 and 200, with specified (fixed) weights for choosing high prob nucleotide--- 
    q_75_rand, p_75_rand= generateQueryGapless(gpath, ppath, weight = 0.75)
    q_50_rand, p_50_rand = generateQueryGapless(gpath, ppath, weight = 0.50)
    q_25_rand, p_25_rand = generateQueryGapless(gpath, ppath, weight = 0.25)
    q_rand_rand, p_rand_rand = generateQueryGapless(gpath, ppath, weight = 0)
    
    query["q_75_rand"] = q_75_rand
    query["q_50_rand"] = q_50_rand   
    query["q_25_rand"] = q_25_rand 
    query["q_rand_rand"] =  q_rand_rand   
    query["p_75_rand"] = p_75_rand    
    query["p_50_rand"] = p_50_rand    
    query["p_25_rand"] = p_25_rand  
    query["p_rand_rand"] = p_rand_rand

    # -- generate 25 query sequences of specified (fixed) lengths, with specified (fixed) weights for choosing high prob nucleotide--- 
    # ----- weight = 0.75 ------
    q_75_25, p_75_25= generateQueryGapless(gpath, ppath, weight = 0.75, length = 25)
    q_75_50, p_75_50= generateQueryGapless(gpath, ppath, weight = 0.75, length = 50)
    q_75_200, p_75_200= generateQueryGapless(gpath, ppath, weight = 0.75, length = 200)
    
    query["q_75_25"] = q_75_25
    query["q_75_50"] = q_75_50
    query["q_75_200"] = q_75_200
    query["p_75_25"] = p_75_25
    query["p_75_50"] = p_75_50
    query["p_75_200"] = p_75_200
    
    # ----- weight = 0.50 ------
    q_50_25, p_50_25= generateQueryGapless(gpath, ppath, weight = 0.50, length = 25)
    q_50_50, p_50_50= generateQueryGapless(gpath, ppath, weight = 0.50, length = 50)
    q_50_200, p_50_200= generateQueryGapless(gpath, ppath, weight = 0.50, length = 200) 
    
    query["q_50_25"] = q_50_25
    query["q_50_50"] = q_50_50
    query["q_50_200"] = q_50_200
    query["p_50_25"] = p_50_25
    query["p_50_50"] = p_50_50
    query["p_50_200"] = p_50_200
    
    # ----- weight = 0.25 ------
    q_25_25, p_25_25= generateQueryGapless(gpath, ppath, weight = 0.25, length = 25)
    q_25_50, p_25_50= generateQueryGapless(gpath, ppath, weight = 0.25, length = 50)
    q_25_200, p_25_200= generateQueryGapless(gpath, ppath, weight = 0.25, length = 200)
    
    query["q_25_25"] = q_25_25
    query["q_25_50"] = q_25_50
    query["q_25_200"] = q_25_200
    query["p_25_25"] = p_25_25
    query["p_25_50"] = p_25_50
    query["p_25_200"] = p_25_200
    
    # ----- weight = 0 ---------      
    q_rand_25, p_rand_25= generateQueryGapless(gpath, ppath, weight = 0, length = 25)
    q_rand_50, p_rand_50= generateQueryGapless(gpath, ppath, weight = 0, length = 50)
    q_rand_200, p_rand_200= generateQueryGapless(gpath, ppath, weight = 0, length = 200) 
    
    query["q_rand_25"] = q_rand_25
    query["q_rand_50"] = q_rand_50
    query["q_rand_200"] = q_rand_200
    query["p_rand_25"] = p_rand_25
    query["p_rand_50"] = p_rand_50
    query["p_rand_200"] = p_rand_200

    #--- now, let the weight be inversely proportional to query length --- 
    # --- variable length (randomly generated) --- 
    q_prop_rand, p_prop_rand = generateQueryGapless_shortQueryHighScore(gpath, ppath)
    query["q_prop_rand"] = q_prop_rand
    query["p_prop_rand"] = p_prop_rand
    
    # --- specifies length, weight always proportional to length, i.e. "prop" --- 
    q_prop_25, p_prop_25 = generateQueryGapless_shortQueryHighScore(gpath, ppath, length = 25)
    q_prop_50, p_prop_50 = generateQueryGapless_shortQueryHighScore(gpath, ppath, length = 50)
    q_prop_100, p_prop_100 = generateQueryGapless_shortQueryHighScore(gpath, ppath, length = 100)
    q_prop_200, p_prop_200 = generateQueryGapless_shortQueryHighScore(gpath, ppath, length = 200)

    query["q_prop_25"] = q_prop_25
    query["q_prop_50"] = q_prop_50
    query["q_prop_100"] = q_prop_100
    query["q_prop_200"] = q_prop_200
    query["p_prop_25"] = p_prop_25
    query["p_prop_50"] = p_prop_50
    query["p_prop_100"] = p_prop_100
    query["p_prop_100"] = q_prop_200
            
    
    with open('../data/query.txt', 'wb') as handle:
        pickle.dump(query, handle)
    return query      
                
            
    

if __name__ == "__main__":
    #example usage
    
    gpath = "../data/chr22.maf.ancestors.42000000.complete.boreo.fa.txt"
    ppath = "../data/chr22.maf.ancestors.42000000.complete.boreo.conf.txt"
    query = queryGenerator(gpath, ppath)
    
  
    