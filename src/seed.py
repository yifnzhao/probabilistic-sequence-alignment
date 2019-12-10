#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 20:46:14 2019

@author: yifan
"""


import numpy as np
from re import finditer
from math import log
import queryGenerator
import time
import pickle
import seaborn as sns
from statistics import mean, stdev
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import PercentFormatter

def seedS(querySeq, highest_prob_genome_str, k = 8):
    '''
    k = word size to be searched
    reutrns a list of seed positions [(x,y)，seed]
    '''
    print("seedS in progress...")
    start = time.time()
    seedPosList = []
    for i in range(len(querySeq)-k):
        query = querySeq[i:i+k]
        for match in finditer(query, highest_prob_genome_str):            
            seedPosList.append([match.span(), match.group(), (i, i+k)])
            # i is the start position of the seed in query
    end = time.time()
    timeElapsed = end - start
    print("seedS: " + str(timeElapsed))
    return timeElapsed, seedPosList


def seedM_fix_n(querySeq, matrix, k = 8, n = 100):
    '''
    k: word size
    n: number of seeds to look for (fixed)
    returns seedPosList which is a list of seed positions (x,y)
    '''
    print("seedM_fix_n in progress...")
    start = time.time()
    seedPosList = [0 for i in range(n)]
    seedProbList = [100 for i in range(n)]
    #loop over all k-words in querySeq
    for i in range(len(querySeq)-k): 

        k_word = querySeq[i:i+k]
#        print (k_word)
        #loop over all positions in matrix
        for pos in range(len(matrix['A'])-k): 
            #find probability : thisQuery_prob
            thisQuery_prob = 0    
            for j in range(k): 
                query_n = k_word[j]
                thisQuery_prob += matrix[query_n][pos+j]
                # TODO: argmax is probably the most time consuming step
                # find index of lowest probability 
            
            # find min prob
            index_min = np.argmax(seedProbList) 
            minSeedProb = seedProbList[index_min] 
            
            # compare this query prob with min prob
            if thisQuery_prob < minSeedProb: # comparing neg log prob
                seedPosList[index_min] = [(pos, pos+k), k_word,  (i, i+k)]
                # i is the start position of the seed in query
                seedProbList[index_min] = thisQuery_prob
    
    end = time.time()
    timeElapsed = end - start
    print("seedM_fix_n: " + str(timeElapsed))
    return timeElapsed, seedPosList,seedProbList

def seedM_fix_min_prob(querySeq, matrix, k = 8, min_prob = 0.1):
    '''
    k: word size
    n: number of seeds to look for 
    returns seedPosList which is a list of seed positions (x,y)
    '''
    print("seedM_fix_min_prob in progress...")
    start = time.time()
    seedPosList = []
    seedProbList = []
    for i in range(len(querySeq)-k): #loop over all k-words in querySeq
        k_word = querySeq[i:i+k]
        for pos in range(len(matrix['A'])-k): #loop over all positions in Matrix
            # -- find probability : thisQuery_prob --
            thisQuery_prob = 0    
            for j in range(k): 
                query_n = k_word[j]
                thisQuery_prob += matrix[query_n][pos+j]
            # -- compare this query prob with the min prob --
            if thisQuery_prob < -log(min_prob): # comparing neg log prob
                seedPosList.append([(pos, pos+k), k_word,  (i, i+k)])
                # i is the start position of the seed in query
                seedProbList.append(thisQuery_prob)
    end = time.time()
    timeElapsed = end - start
    print("seedM_fix_min_prob: " + str(timeElapsed))
    return timeElapsed, seedPosList,seedProbList      


def seedGenerator(query, gpath, ppath, matrix, highest_prob_genome_str, kword = 8):
    '''
    output: querySeeds, a dictionary querySeeds["generationMethod"] = [list of seeds for each query]
    '''
    
    print("Start seeding...")
    querySeeds = {}
    for k, v in query.items(): # iterate through each method of query generation
        
#        if k in ['q_25_rand','q_50_rand','q_75_25','q_100_50d',
#            'q_75_50', 'q_75_rand', 'q_100_25', 'q_100_50',
#            'q_100_200', 'q_100_rand', 'q_rand_rand']:
#            continue
        print("Seeding " + k)  
        
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
                seed_dict["M_fix_maxSeeds_75"] = [t, seedPosListM_fix_n, seedProbListM_fix_n]
                t, seedPosListM_fix_n, seedProbListM_fix_n = seedM_fix_n(query_seq, matrix, n = 50, k=kword)
                seed_dict["M_fix_maxSeeds_50"] = [t, seedPosListM_fix_n, seedProbListM_fix_n]
                t, seedPosListM_fix_n, seedProbListM_fix_n = seedM_fix_n(query_seq, matrix, n = 25, k=kword)
                seed_dict["M_fix_maxSeeds_25"] = [t, seedPosListM_fix_n, seedProbListM_fix_n]
                
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
            eval(k)
            
            with open( '../data/' + k + '.txt', 'wb') as handle:
                pickle.dump(querySeeds[k], handle)
        
        
    return querySeeds

def eval(queryGenerationMethod):
    
    print("Reading seeds and query...")
    with open('../data/' +queryGenerationMethod+ '.txt','rb') as handle:
        seeds = pickle.loads(handle.read())
    print(queryGenerationMethod) # method of query generation

    for k, v in seeds[0][1].items(): 
        if k == 'S':
            S_p = []
            for s in v[1]:
                seedStart = s[0][0]
                seedSeq = s[1]
                seed_p = queryGenerator.findProb(matrix, seedStart, seedSeq)
                S_p.append(seed_p)
            v.append(S_p)
        
        print(k) # method of seed generation
        print(v[0]) # time
        
        # read probability info, find distribution
        
        x = v[2]
        
        plt.hist(x, bins=30)
        plt.xlabel('- log(p) distribution of '+ queryGenerationMethod + ';' + k);
        plt.tight_layout()

        plt.savefig('../data/' + queryGenerationMethod + '_' + k + '.png',dpi=500)
        plt.show()
         # calculate statistics : mean and std dev
        avg_p = mean(v[2])
        stdev_p = stdev(v[2])
        print("Mean: " +str(avg_p))
        print("Std Deviation: " +str(stdev_p))
        v.append(avg_p)
        v.append(stdev_p)
        
#    with open( '../data/' + queryGenerationMethod + '.txt', 'wb') as handle:
#                pickle.dump(seeds, handle)
#        
#    
#    

if __name__ == "__main__":
    # example usage
    gpath = "../data/chr22.maf.ancestors.42000000.complete.boreo.fa.txt"
    ppath = "../data/chr22.maf.ancestors.42000000.complete.boreo.conf.txt"
    matrix, hpg_str = queryGenerator.loadData(gpath, ppath)
    
    print("Reading query...")
#    with open('../data/query.txt','rb') as handle:
#        query = pickle.loads(handle.read())
#    querySeeds = seedGenerator(query, gpath, ppath, matrix, hpg_str)
#    
    
   #  --- evaluate seeds quality  --- 
   # example usage

#    eval('q_100_50')
   
   
        
   
    with open('../data/hspList_1210_runtime.txt','rb') as handle:
        runtime = pickle.loads(handle.read())
    
    x = [pos[2] for pos in runtime]
    y = [pos[0] for pos in runtime]
    plt.xlabel('Seed start position relative to query');
    plt.ylabel('Ungapped viterbi runtime (s)')
    plt.plot(x, y, '.');
    
    plt.savefig('../diagrams/' + 'fig2.png',dpi=500)
    plt.show()
    
    plt.hist(x, bins=20)
    plt.xlabel('Seed start position relative to query');
    # total: 83 seeds
    
    plt.savefig('../diagrams/' + 'fig2(2).png',dpi=500)
