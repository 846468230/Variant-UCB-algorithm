#!/usr/bin/python
# -*- coding: UTF-8 -*-

import sys, getopt
from experiment import Proteins,Configs,Receptors,Receptors_for_random,Virtual_ligands
from solvers import UCB_Variation
from config import search_bind
import os
import operator
import json
import csv
mode = "online" # "offline" 
virtual = False # True
def experiment(N,ligands_path,receptors_path,configs_path,percent = 0.15,ucb_coefficient = 1,n_virtual_ligands=1000):
    """
    Run a small experiment on solving a molecular docking problem,
    each with a reward judged by vina.

    Args:
        N (int): number of time steps to try.
        ligands_path (str) : path of the ligand files stored . (suffix of the ligands file is pdbqt)
        receptors_path (str) :path of the receptor files store . (suffix of the receptor file is pdbqt) 
        configs_path (str) : path of the configs files store . (suffix of the configs file is txt,every receptor's config file's name should the same with receptor's) 
        percent (int or float) : percentage or numbers of proteins to be selected through dock.
        ucb_coefficient (int) : The coefficient of the formula
    """
    global virtual
    if not os.path.exists("./out"):
        os.system("mkdir out")
    receptors = Receptors(receptors_path)
    configs = Configs(configs_path)
    if virtual:
        ucb_variation_ligands = Virtual_ligands(n_virtual_ligands)   # numbers of virtual ligands to be created
    else:
        ucb_variation_ligands = Proteins(ligands_path)

    test_solvers = [
        UCB_Variation(ucb_variation_ligands,receptors,configs,percent,c=ucb_coefficient)
    ]

    for s in test_solvers:
        s.run(N)
        item_list = sorted(s.estimated_probas.items(),key=operator.itemgetter(1),reverse=True)
        #print('\n%sligands_score:'%(s.name),item_list)
        with open('%s_order.txt'%(s.name),'w') as f:
            f.write(json.dumps(s.actions))
        with open('%s_exploit.txt'%(s.name),'w') as f:
            f.write(json.dumps(s.exploit_counts))
        with open('%s_explore.txt'%(s.name),'w') as f:
            f.write(json.dumps(s.explore_counts))
        with open('./%s_result.txt'%(s.name), 'w') as f:
            f.write("filename\t score\t \n")
            for i in range(len(item_list)):
                f.write("%s\t %s\t \n"%(item_list[i][0],item_list[i][1]))
        if hasattr(s,"process"):
            with open('%s_process.csv'%(s.name),"w") as f: 
                writer = csv.writer(f)
                writer.writerow(["id","ligand","receptor","PKa"])
                writer.writerows(s.process)
def main(argv):
   global mode
   counts = None
   percent = None
   ucb_coefficient = None
   config_path = None
   solve = None
   n_virtual_ligands = None
   try:
        opts, args = getopt.getopt(argv,"hc:p:u:f:o:m",["counts=","percent=","ucb=","config=","solve","mode="])
   except getopt.GetoptError:
        print ('python main.py -c <dock counts> -p <percentage> -u <ucb_coefficient> -f <config_path> -o <don\'t solve>')
        sys.exit(2)
   for opt, arg in opts:
        if opt == '-h':
            print('python main.py -c <dock counts> -p <percentage> -u <ucb_coefficient> -f <config_path> -o <don\'t solve>')
            sys.exit()
        elif opt in ("-c", "--counts"):
            counts = int(arg)
        elif opt in ("-p", "--percent"):
            percent = float(arg)
        elif opt in ("-u", "--ucb"):
            ucb_coefficient = float(arg)
        elif opt in ("-f", "--config"):
            config_path = arg
        elif opt in ("-o", "--solve"):# dont solve the question
            solve = arg
        elif opt in ("-m","--mode"):
            mode = arg
        elif opt in ("-n","--n_virtual_ligands"):
            n_virtual_ligands = arg
        
   return counts,percent,ucb_coefficient,config_path,solve,n_virtual_ligands

if __name__ == '__main__':
    counts,percent,ucb_coefficient,config_path,solve,n_virtual_ligands= main(sys.argv[1:])
    if percent == None:
        percent = 3
    if counts == None:
        counts = 100
    if ucb_coefficient == None:
        ucb_coefficient = 1.0
    if n_virtual_ligands == None:
        n_virtual_ligands = 1000
    if config_path != None:
        search_bind(config_path)
        print("config is completed searched")
    if solve == None:
        print("solving····")
        experiment(counts,"./ligand","./receptor","./conf",percent=percent,ucb_coefficient=ucb_coefficient,n_virtual_ligands=n_virtual_ligands)
