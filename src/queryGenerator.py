#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 15:41:46 2019

@author: yifan
"""

from seed import loadData
import random
from datetime import datetime




def generateHighestScoreQuery(g_path, p_path, n_query=100):
    random.seed(datetime.now())

    matrix, highest_prob_genome_str = loadData(genome_path, prob_path)
    len_genome = len(highest_prob_genome_str)
    queryList = []
    for i in range(n_query):
        start_pos = random.randint(0, len_genome-1)
        maxLen = len_genome - start_pos
        query_len = random.randint(0, 200) #TODO: 200 is chosen arbitrarily, what would be a more statistically relevant value?
        if query_len > maxLen:
            continue
        query = highest_prob_genome_str[start_pos: start_pos+query_len]
        queryList.append(query)
    return queryList
    

def generateQueryGapless(g_path, p_path, weight = 0.75, n_query=100):
    # 75% probability select highest scoring nucleotide by default 
    random.seed(datetime.now())
    allNu = ['A','C','T','G']
    matrix, highest_prob_genome_str = loadData(genome_path, prob_path)
    len_genome = len(highest_prob_genome_str)
    queryList = []
    for i in range(n_query):        
        start_pos = random.randint(0, len_genome-1)
        maxLen = len_genome - start_pos
        query_len = random.randint(0, 200) #TODO: 200 is chosen arbitrarily, what would be a more statistically relevant value?
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
        queryList.append(query)
    return queryList


if __name__ == "__main__":
    #example usage
    genome_path = "../data/chr22.maf.ancestors.42000000.complete.boreo.fa.txt"
    prob_path = "../data/chr22.maf.ancestors.42000000.complete.boreo.conf.txt"
    highest_score_query_list= generateHighestScoreQuery(genome_path, prob_path)
    query_list_75_gapless = generateQueryGapless(genome_path, prob_path, weight = 75)
    query_list_50_gapless = generateQueryGapless(genome_path, prob_path, weight = 50)
    query_list_25_gapless = generateQueryGapless(genome_path, prob_path, weight = 25)
    query_list_random_gapless = generateQueryGapless(genome_path, prob_path, weight = 0)
    




'''
One potential approach to evaluate your algorithm is to use the probabilistic
 genome to generate a short query sequence by picking a starting position 
 randomly in the genome, generating the query sequence by selecting nucleotides
 randomly according to the probabilities at each position, and possibly adding
 a few additional random substitutions and indels. You can then feed this query 
 sequence to your algorithm and see if it succeeds at identifying the correct 
 portion of the probabilistic genome.
'''
