#!/usr/bin/python
# -*- coding: UTF-8 -*-
from __future__ import division

import time
from random import choice
import operator
import random
import math
from functools import reduce

class Solver(object):
    def __init__(self, ligands, receptors,configs):

        random.seed(int(time.time()))

        self.bandit = ligands
        self.receptors = receptors
        self.configs = configs
        self.t = 0
        self.counts = {}
        for item in self.bandit.items:
            self.counts[item] = 0  
        self.actions = []  # A list of machine ids, 0 to bandit.n-1.

    @property
    def estimated_probas(self):
        raise NotImplementedError

    def run_one_step(self):
        """Return the machine index to take action on."""
        raise NotImplementedError

    def run(self, num_steps):
        assert self.bandit is not None
        for x in range(num_steps):
            item = self.run_one_step()
            if x%10==0:
                print("%s dock %s times"%(self.name,x))
            self.t += 1
            self.counts[item] += 1
            self.actions.append(item)

def stddev(lst):
    mean = float(sum(lst)) / len(lst)
    return math.sqrt(float(reduce(lambda x, y: x + y, map(lambda x: (x - mean) ** 2, lst))) / len(lst))


class rUCB(Solver):
    def __init__(self, ligands, receptors, configs, percent, c=1, init_proba=0.0):
        super().__init__(ligands, receptors, configs)
        from main import virtual
        self.name = "rUCB"
        self.c = c
        if percent<1:
            self.winner = math.ceil(percent * self.bandit.length )  # if percent < 1 change it to int according the bandit length 
        else:
            self.winner = int(percent)
        self.estimates = {}
        self.scores = {}
        self.explore = {}   
        self.exploit = {}  
        self.logKas = {}
        self.sigmas = {}
        self.virtual = virtual
        self.c = math.exp(1.96-1.9*((self.winner/self.bandit.length)**0.16))   # if you have a better factor you can change it
        self.process = []
        for item in self.bandit.items:
            self.estimates[item] = init_proba
            self.scores[item] = init_proba
            self.explore[item] = 0
            self.exploit[item] = 0 
            self.logKas[item] = []
            self.sigmas[item] = 0.35
    @property
    def estimated_probas(self):
        return self.scores
    
    @property
    def exploit_counts(self):
        return sorted(self.exploit.items(),key = lambda one_item : one_item[1],reverse = True)
    
    @property
    def explore_counts(self):
        return sorted(self.explore.items(),key = lambda one_item : one_item[1],reverse = True)

    def update_all_item_sorted(self):
        # Docking all the ligands for the first time and sorting according to the scores after docking
        self.remain_items = sorted([(one_item[0],(math.log10(one_item[1])+ self.c*self.sigmas[one_item[0]]*math.sqrt(1 / (self.counts[one_item[0]]))))for one_item in self.estimates.items()],key=lambda one_item:one_item[1],reverse=True)
        self.winner_items = self.remain_items[:self.winner]
        self.sorted_items = sorted([( one_item[0],(math.log10(self.estimates[one_item[0]]) - self.c*self.sigmas[one_item[0]]*math.sqrt(1 / (self.counts[one_item[0]]))) ) for one_item in self.winner_items ],key=lambda one_item:one_item[1], reverse=True)
    def update_one_item_sorted(self):
        # after docking every ligands for once,then every dock, use this function to sort.
        new_score_add = math.log10(self.estimates[self.item])+ self.c*self.sigmas[self.item]*math.sqrt(1 / (self.counts[self.item]))
        flag = -1
        for i in range(self.bandit.length):
            if(self.remain_items[i][0]==self.item):
                self.remain_items.pop(i)
                break
        for i in range(self.bandit.length-1):
            if new_score_add > self.remain_items[i][1]:
                self.remain_items.insert(i,(self.item,new_score_add))
                flag = i
                break        
        if  flag== -1:
            self.remain_items.append((self.item,new_score_add))
        for i in range(self.winner):
            if(self.sorted_items[i][0]==self.item):
                self.sorted_items.pop(i)
                break
        if flag>=self.winner:
            self.item = self.remain_items[self.winner-1][0]
        flag = -1
        new_score_minus = math.log10(self.estimates[self.item])- self.c*self.sigmas[self.item]*math.sqrt(1 / (self.counts[self.item])) 
        for i in range(self.winner-1):
            if(new_score_minus>self.sorted_items[i][1]):
                self.sorted_items.insert(i,(self.item,new_score_minus))
                flag=i
                break 
        if flag==-1:
            self.sorted_items.append((self.item,new_score_minus))
        assert len(self.remain_items) == self.bandit.length
        assert len(self.sorted_items) == self.winner

    def run_one_step(self):
        '''
            this function was called every dock
        '''
        if self.t < self.bandit.length:
            item = self.bandit.items[self.t]
            self.explore[item] += 1
        else:
            if self.t == self.bandit.length:
                self.update_all_item_sorted()
            else:
                self.update_one_item_sorted()
            self.item = self.sorted_items[-1][0]
            item = self.item
            self.exploit[item] += 1
        if self.virtual:
            ligand = self.bandit.choose_item(item)
            r = self.bandit.generate_reward(ligand)
        else:
            ligand = self.bandit.choose_item(item)
            receptor = self.receptors.choose_item(item,self.bandit,self.estimates,self.t,self.bandit.length)
            config = self.configs.choose_item(receptor)
            r = self.bandit.generate_reward(ligand,receptor,config)  # call vina to judge this pair of ligand and receptor
        self.estimates[item] += 1. / (self.counts[item] + 1) * (r - self.estimates[item])
        self.scores[item] = self.estimates[item]
        self.logKas[item].append(math.log10(r))
        if self.counts[item] >=1:
            sigma = stddev(self.logKas[item])
            if sigma == 0:
                sigma = 0.35
            self.sigmas[item]=0.35
        if receptor:
            self.process.append((self.t+1,item,receptor.split("/")[-1],r))
        if self.receptors.delete:       # if every receptor was docked with this ligand , then the ligand no longer participate in the next docks
            del(self.estimates[item])
        return item
