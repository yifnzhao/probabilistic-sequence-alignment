#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 11:03:26 2019

@author: yifan
"""
from seed import seedS, loadData
import random
from datetime import datetime





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
    genome_path = "../data/chr22.maf.ancestors.42000000.complete.boreo.fa.txt"
    prob_path = "../data/chr22.maf.ancestors.42000000.complete.boreo.conf.txt"
    matrix, highest_prob_genome_str = loadData(genome_path, prob_path)
    querySeqTester = 'AAGGGGTTTTACGGAATTCCGAAC' #tester
    seedPosList = seedS(querySeqTester, highest_prob_genome_str, k=7)
    # toy example
    obs = querySeqTester # query sequence
    seed = seedPosList[100]
    hsp = ungapped_viterbi(obs, seed, matrix)
    # in hsp:
    #  0: prob; 1:pos wrt matrix; 2: pos wrt query


        
    