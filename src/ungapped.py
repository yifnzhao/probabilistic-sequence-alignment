#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 11:03:26 2019

@author: yifan
"""
from seed import seedGenerator, seedS, seedM_fix_n
from queryGenerator import loadData, findProb
import random
from datetime import datetime
import time
import pickle

def findBaselineProb(matrix, hspStart, hsp):
    allNu = ['A','C','T','G']
    prob = 0
    for index in range(len(hsp)):
        n = random.choice(allNu)
        prob += matrix[n][hspStart]
        hspStart+=1
    return prob
    

def ungapped_viterbi(obs, seed, matrix):
    
    '''
    seed: an element from seedPosList
    obs: query sequence
    matrix: probability matrix of the genome
    V[0]: current probability (neg log)
    V[1]: current position wrt matrix, tuple (start, end)
    V[2]: current position wrt query sequence, tuple (start, end)
    V[3]: current state (seed, Mi or Mj), str
    V[4]: previous position wrt matrix
    V[5]: previous position wrt query
    V[6]: previous state (seed, Mi or Mj)
    '''
    random.seed(datetime.now())
    start_p = {'seed': 1,
              'Mi': 0,
              'Mj' : 0}
    
    # TODO: naive implementation, what would be more statistically relevant?
    trans_p = {'seed': {'seed':0, 'Mi':0.5, 'Mj': 0.5},
               'Mi': {'seed':0, 'Mi':0.5, 'Mj': 0.5},
               'Mj': {'seed':0, 'Mi':0.5, 'Mj': 0.5}}
    
    states = ('Mi', 'Mj') # no need to check for the seed after initialization
    
    
    
    seedPos_m, seedStr, seedPos_q = seed
    seedStart_m,seedEnd_m = seedPos_m
    seedStart_q,seedEnd_q = seedPos_q
    length_hsp = len(seedStr)
    
    #initialization
    V = [[]] * (len(obs)+1)       
    currState = 'seed'
    currPos_m = seedPos_m
    currPos_q = seedPos_q
    current_start_prob = start_p[currState]
    current_emit_prob = findProb(matrix, seedStart_m, seedStr)
    info_thisPos_thisState = [current_start_prob + current_emit_prob,
                  currPos_m, currPos_q, currState, None, None, None]
    V[length_hsp] = [info_thisPos_thisState]
    
    
    while length_hsp < len(obs)-1:
        V[length_hsp] = sorted(V[length_hsp]) 
        # --- checking baseline -----
        hspPos_m = V[length_hsp][0][1]
        hspPos_q = V[length_hsp][0][2]
        hspStart = hspPos_m[0]
        hsp = obs[hspPos_q[0]: hspPos_q[1]]
        bestProb = V[length_hsp][0][0]
        baselineProb = findBaselineProb(matrix, hspStart, hsp)
        if bestProb > baselineProb:
            break

        for info_thisPos_thisState in V[length_hsp]:  
            prevProb = info_thisPos_thisState[0]
            prevPos_m = info_thisPos_thisState[1]
            prevPos_q = info_thisPos_thisState[2]
            prevState = info_thisPos_thisState[3]
            # ---- loop over all possible current states ----
            for possible_currState in states:
                 if possible_currState == 'Mi':
                     if prevPos_q[0]==0:
                         continue
                     nu_pos_m = prevPos_m[0]-1
                     currPos_m = (nu_pos_m, prevPos_m[1])
                     nu_pos_q = prevPos_q[0]-1
                     currPos_q = (nu_pos_q, prevPos_q[1])
                     nu = obs[nu_pos_q]
                     length_hsp = currPos_m[1] - currPos_m[0]
                 elif possible_currState == 'Mj':
                     if prevPos_q[1]==len(obs)-1:
                         continue
                     nu_pos_m = prevPos_m[1]+1
                     currPos_m = (prevPos_m[0], nu_pos_m)
                     nu_pos_q = prevPos_q[1]+1
                     currPos_q = (prevPos_q[0], nu_pos_q)
                     nu = obs[nu_pos_q]
                     length_hsp = currPos_m[1] - currPos_m[0]
                 emission_p = matrix[nu][nu_pos_m]
                 transition_p = trans_p[prevState][possible_currState]
                 prob = prevProb + transition_p + emission_p 
        
                    
                 # --- update V ----
                 info = [prob, currPos_m, currPos_q, possible_currState, 
                         prevPos_m, prevPos_q, prevState]
                 if V[length_hsp] == []:
                    V[length_hsp] = [info]
                 else:
                    V[length_hsp].append(info)
    
    V[length_hsp] = sorted(V[length_hsp])
    return V[length_hsp][0][0:3] 


if __name__ == '__main__':
    
    # Calculate the most probable state path using the Viterbi algorithm. 
    gpath = "../data/chr22.maf.ancestors.42000000.complete.boreo.fa.txt"
    ppath = "../data/chr22.maf.ancestors.42000000.complete.boreo.conf.txt"
    m, hpg_str = loadData(gpath, ppath)
#    # --- toy example --- 
#    querySeqTester = 'AAGGGGTTTTACGGAATTCCGAAC'
#    seedTester = "AAG" 
#    hsp = ungapped_viterbi(querySeqTester, seedTester, m)
#    
    # in hsp:
    #  0: prob; 1:pos wrt matrix; 2: pos wrt query

    # --- generate seeds with various k-word length ---
#    start = time.time()
#    
#    end = time.time()
#    timeElapsed = end - start #TODO: record this somewhere
#

    with open('../data/query.txt', 'rb') as handle:
        query = pickle.loads(handle.read())
#    
    
    querySeqList = query['q_100_50'] #  a list of 25 queries 
    queryAnswer = query['p_100_50']
    
    with open('../data/query_seeds_1209_test.txt', 'rb') as handle:
        querySeeds = pickle.loads(handle.read())
    counter = 0
    hspList = []   
#    for querySeq in querySeeds:
    # reading the first query
    querySeq = querySeeds[0]
    seq = querySeq[0]
    seedDict = querySeq[1]
    hspS = []
    hspM = []
    runtime = []
    print(seq)
    
#    print("seedM")
#    for seed in seedDict["M_fix_maxSeeds_25"][1]:
#        counter+=1
#        print(counter)
#        hsp = ungapped_viterbi(seq, seed, m)
#        hspM.append(hsp)
#        
        
    print("seedS")
    for seed in seedDict["S"][1]:
        start_g = seed[0][0]
        start_q = seed[2][0]
        counter+=1
        print(counter)
        start = time.time()
        hsp = ungapped_viterbi(seq, seed, m)
        end = time.time()
        timeElapsed = end - start 
        hspS.append(hsp)
        runtime.append([timeElapsed, start_g, start_q])
    
    
    hspList.append((querySeq, hspS, hspM))


        
    with open('../data/hspList_1210_runtime.txt','wb') as handle:
        pickle.dump(runtime, handle)
        
        
        
    
#    querySeeds = []
#    for querySeq in querySeqList: # iterate through each query sequence (25 in total)
#
#        print("Seeding a new sequence...")  
#        kword = 8
#        seed_dict = {}
#        t, seedPosListS = seedS(querySeq, hpg_str, k=kword) 
#        seed_dict["S"] = [t, seedPosListS]
#        t, seedPosListM_fix_n, seedProbListM_fix_n = seedM_fix_n(querySeq, m, n = 25, k=kword)
#        seed_dict["M_fix_maxSeeds_25"] = [t, seedPosListM_fix_n, seedProbListM_fix_n]
#        querySeeds.append([querySeq, seed_dict])
#            
   
    
    

#
#
#
#    with open('../data/q_100_50.txt','rb') as handle:
#        seeds = pickle.loads(handle.read())
#    
#    