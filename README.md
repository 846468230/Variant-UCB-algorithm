#  The interpretation of this code
1. The source codes are main.py experiment.py config.py solvers.py
2. The file called vina can only execute in linux. if you want to execute on macos, you should download the corresponding version and repace it
3. To run main.py, you should prepare ligands and receptor files and the same name config files as receptor files. Place them in the path of ./conf,./ligand,./receptor,the symbol './' represents the current directory.  And delete the existing file which were as examples before .
4. My running example
   1. python main.py -p 3 -c 100 (select the first three ligands and the total docking time is 100)
5. After excecuting the codes, five results file will be created.
   1. UCB_Variation_process.csv records the docking process and scores.
   2. UCB_Variation_explore.txt records every ligand was explored times.
   3. UCB_Variation_exploit.txt records every ligand was exploited times.
   4. UCB_Variation_order.txt records ligands docking sequence.
   5. UCB_Variation_result.txt records every ligand's score.
6. Reference: Bin Chong, Yingguang Yang, Zile Wang, and Zhirong Liu,Reinforcement Learning to Boost Molecular Docking upon Protein Conformational Ensemble.
