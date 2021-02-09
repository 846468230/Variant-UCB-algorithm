#  The interpretation of this code
1. The source codes are main.py experiment.py config.py solvers.py
   1. The experiment.py : Three classes of Receptors, Proteins, and Configs are saved in the experiment.py file. In order to save the file name and path for the receptor and ligand directory, so that when calling vina to score, provide the file path parameter.
   2. The config.py : if you don't have the config files about where is the pocket of the receptor.You need to use cavity to compute the pocket files,this function may helpful for you.
   3. The solvers.py : The class rUCB is a implementation of the algorithm which we proposed.
2. The file called vina can only execute in linux. if you want to execute on macos, you should download the corresponding version and repace it. Also you maybe need gave the vina execute permission use ```chmod u+x vina```
3. To run main.py, you should prepare ligands and receptor files and the same name config files as receptor files. Place them in the path of ./conf,./ligand,./receptor,the symbol './' represents the current directory.  And delete the existing file which were as examples before.
4. My running example
   1. ```python main.py -m 20 -T 56958``` (select the first 20 ligands and the total docking time is 56958)
   2. if you want to change the coefficient of the formula rUCB you can delete the line 66 in solvers.py and ```python main.py -m 20 -T 56958 -u <float>```
5. After excecuting the codes, five results file will be created.
   1. rUCB_process.csv records the docking process and predicted Ka (calculated from scores).
   2. rUCB_initialization.txt records every ligand was docked once.
   3. rUCB_loop.txt records every ligand was docked times during loop.
   4. rUCB_order.txt records ligands docking sequence.
   5. rUCB_result.txt records every ligand's \<Ka\>.
6. Reference: Bin Chong, Yingguang Yang, Zi-Le Wang, Han Xing and Zhirong Liu,Reinforcement Learning to Boost Molecular Docking upon Protein Conformational Ensemble.
