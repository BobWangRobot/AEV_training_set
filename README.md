# Introduction
This project can annotate helix by AEV method. We downloaded 10 high resolution PDBs and stored them in **PDB/standard_set**. Then we checked the 10 PDBs, got the best structure manually, and stored the residues' name and sequeen in the **standard_helices.xlsx** file. We use the script **read_xls.py** to read the standard.xls and modified records of PDBs, then store the new structures(PDBs) in the **mannual_set**. We also want to get a training set automacally, so we use **train.py** modified helix records of PDBs from **standard_set** and stored the training structrue in the **train_set**. Finally, we compared the training structures **(train_set)** with the structures modified manually **(manual_set)** by **evaluate.py** to evaluate the training process **(train.py)**  
## script
### read_xls.py
read standard_helices.xlsx and output structures to the manual_set directory.


### train.py
input PDBs in the manual_set and output PDBs to the training_set directory.  


### evaluate.py
compare manual_set with training_set and evaluate the result. The result include 5 parameters: helix_ID, Status, length,  name and sequeen of start residue, name and sequeen of end residue. 
For example:

|ID|Status|length|start_res|end_res|
|--|------|------|---------|-------|
|1 | get  |identical|GLU  59 |ALA  59|
|2 | get  |longer   |ALA  22(ASP  20) |HIS  36(SER  35)|


It means both traing set and manaul set have helix 1 and helix 2, but helix 1 from manual set and training set are identical, helix 2 from manual set and training set are simialr which start residues and end residues are different.

There are 3 parameters of  **Status**: **get, more, less**. **get**, the training strcuture get the helix which is in the manual structure. **more**, the training structure has the 'wrong' helix which is not in the manual structure. **less**, the training structure lost the helix which is in the manual structure. 


There are also 5 parameters of **length**: **same, identical, shoter, longer, None**. **identical**: tow helices are identical. **same**: the length of two helices are identical but start residues are differnt. **shorter**: the helix from structures of training set are shorter than manual set. **longer**: the helix from structures of training set are longer than manual set. The max diffrence of longer and shorter less than 4 reisidues. **None**: can't compare(such as loss).

# Run script
**read_xls.py:** `phenix.python read_xls.py standard_helices.xlsx`  
**train.py:** `phenix.python train.py PDB/standard_set`  
**evaluate.py:** `phenix.python evaluate.py PDB/training_set PDB/manual_set`
