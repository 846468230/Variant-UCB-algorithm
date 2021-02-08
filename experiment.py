#!/usr/bin/python
# -*- coding: UTF-8 -*-
import shutil,os
import subprocess
import time
from random import choice
import math
import random

standard_score_initialed = False
standard_scores = {}

class Files(object):
    def __init__(self, path):
        """
        Files the ligands receptors or configs.
        """
        assert os.path.isdir(path)
        self.path = path
        random.seed(int(time.time()))
        
        self.items = gci(self.path)
        self.length = len(self.items)
    
    def choose_item(self):
        raise NotImplementedError
        

class Proteins(Files):
    def __init__(self,path):
        super(Proteins, self).__init__(path)

    def choose_item(self,item):
        return os.path.join(self.path,item)

    def generate_reward(self,ligand,receptor,config):
        from main import mode
        if mode == 'offline' :
            return judgeOffline(ligand,receptor,config)
        return judge(ligand,receptor,config)

class Virtual_ligands(Proteins):
    def __init__(self,n_ligands):
        random.seed(int(time.time()))
        self.items = [ i for i in range(n_ligands)]
        self.length = len(self.items)
        self.avedG = {}
        self.sigmadG = {}
        self.aveKa = {}
        self.init_ligand_aveKa()

    def init_ligand_aveKa(self,sigma_avedG=0.65,ave_avedG=-5.1,sigma_sigmadG=0.08,ave_sigmadG=0.44,mu=0,sigma=1,saved=True):
        self.mu = mu
        self.sigma =sigma
        for item in self.items:
            self.avedG[item]= random.gauss(self.mu,self.sigma)*sigma_avedG+ave_avedG
            self.sigmadG[item]= random.gauss(self.mu,self.sigma)*sigma_sigmadG+ave_sigmadG
            if(self.sigmadG[item]<0.0):
                self.sigmadG[item]=0.0
            self.aveKa[item]= math.exp(-self.avedG[item]*4.18e3/(8.31*300)+0.5*math.pow(self.sigmadG[item]*4.18e3/(8.31*300),2.0))
        if saved:
            import json
            with open('Virtual_ligand_standard.txt','w') as f:
                f.write(json.dumps(self.aveKa))
    
    def choose_item(self, item):
        return item

    def generate_reward(self,ligand):
        dG = random.gauss(self.mu,self.sigma)*self.sigmadG[ligand] + self.avedG[ligand] 
        Ka = math.exp(-dG*4.18e3/(8.31*300))
        return Ka

class Receptors(Files):
    def __init__(self, path):
        super().__init__(path)
        self.remain_items = {}
        self.delete = False
        self.base_item = choice(self.items)
    
    def choose_item(self,item,ligands,estimates,t,ligand_length):       #i is the i'th receptor every ligand dock with one receptor no more once
        self.delete = False
        if ligand_length < 100 * self.length and self.length > 500:  # The settings can be changed manually, so that a large number of receptors will not be repeatedly docked because of the difference
            if t<ligand_length:
                temp_item = self.base_item 
            else:
                temp_item = choice(self.items) 
            return os.path.join(self.path,temp_item)
        else:
            try:
                remain_items = self.remain_items[item]       #In the first round, every item is a new slot machine and must be created to record
                temp_item = choice(remain_items)
            except KeyError:
                self.remain_items[item] = self.items[:]      # The i-th ligand has not been docked with the protein
                temp_item = self.base_item                   # In the first round, each one is docked with the selected basic protein              
            finally:
                self.remain_items[item].remove(temp_item)
                if len(self.remain_items[item]) == 0:
                    ligands.items.remove(item)
                    self.delete = True
                return os.path.join(self.path,temp_item)
        

class Configs(Files):
    def __init__(self,path):
        super(Configs, self).__init__(path) 
    
    def choose_item(self,names):
        names = names.split("/")[-1][:-6]+'.txt'
        return os.path.join(self.path,names)

def gci(filepath):  
    files = os.listdir(filepath)
    names = []  
    for fi in files:    
        fi_d = os.path.join(filepath,fi)    
        if os.path.isdir(fi_d):
            gci(fi_d)    
        else:
            if fi.endswith(".pdbqt"):
                names.append(fi)
    return names

def judgeOffline(ligand,receptor,conf,out="./out/",standard_score_path="standard_scores"):
    import re
    global standard_score_initialed
    global standard_scores
    ligand_name = int(re.findall(r'\d+',ligand)[0])
    receptor_name = int(re.findall(r'\d+',receptor)[0]) 
    score = 0
    if  standard_score_initialed:
        score = standard_scores[ligand_name].get(receptor_name,-4)
    else:
        files = os.listdir(standard_score_path)
        for fi in files:
            temp_ligand_name=int(re.findall(r'\d+',fi)[0])
            standard_scores[temp_ligand_name]={}
            with open(os.path.join(standard_score_path,fi),'r') as f:
                text=f.readlines()
                for item in text:
                    key = int(item.split(":")[0])
                    value = float(item.split(":")[1])
                    standard_scores[temp_ligand_name][key]=value
        standard_score_initialed = True
        score = standard_scores[ligand_name].get(receptor_name,-4)
    return math.exp(score * -1.6782342233982843)

def judge(ligand,receptor,conf,out="./out/"):
    """
        call vina to score this pair
    """
    out = os.path.join(out,"out"+ligand.split('/')[-1])
    p = subprocess.Popen('./vina --ligand %s --receptor %s --config %s --out %s --cpu 20 --seed 5'%(ligand,receptor,conf,out), shell=True,stdout=subprocess.PIPE)
    p = subprocess.Popen(["grep"," 1 "],stdin=p.stdout, stdout=subprocess.PIPE)
    score = p.communicate()
    try:
        temp=bytes.decode(score[0])
        score = float(temp.split()[1])
    except IndexError as ex:
        print(ligand.split('/')[-1])
        print(receptor.split('/')[-1])
        print(temp)
        print()
        score = -3.0
    return math.exp(score * -1.6782342233982843)


class Receptors_for_random(Receptors):
    def __init__(self, path):
        super().__init__(path)
        self.time = 0
        
    def choose_item(self,item,ligands,estimates):       #i is the i'th receptor every ligand dock with one receptor no more once
        self.delete = False
        if self.time % ligands.length == 0 :
            self.base_item = choice(self.items)
            self.items.remove(self.base_item)
        self.time += 1
        return os.path.join(self.path,self.base_item) 
    
    

if __name__=="__main__":
    try:
        t = judge('ligand_all/Specs_6858.pdbqt','receptor/3.pdbqt','conf/3.txt','out/out.pdbqt') 
    except :
        t = -3
    print(t)
    