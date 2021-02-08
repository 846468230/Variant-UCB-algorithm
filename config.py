#!/usr/bin/python
# -*- coding: UTF-8 -*-
import os,sys
import shutil,subprocess
def traversal(filepath): 
    filenames = []
    files = os.listdir(filepath)  
    for fi in files:    
        fi_d = os.path.join(filepath,fi)    
        if os.path.isdir(fi_d):
            print(os.path.join(filepath, fi_d))
            traversal(fi_d)    
        else:
            if fi.endswith(".pdb"):
                filenames.append(fi)
    return filenames

def traversal_cavity(filepath,config_path):
    filenames = []
    counts = {}
    files = os.listdir(filepath)  
    for fi in files:    
        fi_d = os.path.join(filepath,fi)    
        if os.path.isdir(fi_d):
            print(os.path.join(filepath, fi_d))
            traversal_cavity(fi_d)    
        else:
            if fi.endswith(".pdb") and "_surface_" in fi:
                name = fi.split("_")[0]
                if name not in filenames:
                    filenames.append(name)
                    counts[name]= 1
                else:
                    counts[name]+= 1
    for item in filenames:
        result = None
        choosed_name = None
        num = 0
        x = 0
        y = 0 
        z = 0
        for i in range(counts[item]):
            texts = []
            with open(os.path.join(filepath,item+"_surface_"+str(i+1)+".pdb")) as f:
                texts = f.readlines() 
            for text in texts:
                if "DrugScore" in text:
                    score = float(text.split()[-1])
                    if result == None:
                        result = score
                        choosed_name = i+1
                    else:
                        if result < score:
                            result = score
                            choosed_name = i+1
                    break
        with open(os.path.join(filepath,item+"_cavity_"+str(choosed_name)+".pdb"),"r") as f:
            datas = f.readlines()
            for data in datas:
                if data.startswith("ATOM"):
                    num+=1
                    x += float(data[30:38])
                    y += float(data[38:46])
                    z += float(data[46:54])
        x /= num
        y /= num
        z /= num
        with open(os.path.join(config_path,item+".txt"),"w") as f:
            f.write('receptor = 1.pdbqt\nligand = 1.pdbqt\n\n\n')
            f.write('center_x = %.3f\n'%(x,))
            f.write('center_y = %.3f\n'%(y,))
            f.write('center_z = %.3f\n'%(z,))
            f.write('\n\n\n')
            f.write('size_x = 22\n')
            f.write('size_y = 22\n')
            f.write('size_z = 22\n')
            f.write('\n\n\n')
            f.write('energy_range = 3\n') 
            f.write('exhaustiveness = 8\n')
            f.write('num_modes = 8')
    

def change_parameter(path,name):
    input = []
    result = []
    with open(path,'r') as f:
        input = f.readlines()
    for item in input:
        if item.startswith("RECEPTOR_FILE"):
            result.append("RECEPTOR_FILE "+name+"\n")
        else:
            result.append(item)
    with open(path,'w+') as f:
        f.writelines(result)


def search_bind(filepath,config_path="./cavity",inputfile_name = "cavity-AA.input",conf_path="./config"):
    '''
        if you don't have the config files about where is the pocket of the receptor.
        this function may helpful for you.
    '''
    '''
    filenames = []
    filenames = traversal(filepath)
    if not os.path.exists(config_path):
        os.system("mkdir %s"%(config_path,))
    if not os.path.exists(conf_path):
        os.system("mkdir %s"%(conf_path,))
    for item in filenames:
        shutil.copy(os.path.join(filepath,item),config_path)
        change_parameter(inputfile_name,os.path.join(config_path,item))
        out_bytes = subprocess.check_output(['./cavity64',inputfile_name])
        #out_bytes = e.output       # Output generated before error
        #code      = e.returncode   # Return code
        #print(code)
        os.remove(os.path.join(config_path,item))
        print("finish1")
    '''
    print("please waiting for complete.")
    traversal_cavity(config_path,conf_path)
    
    
