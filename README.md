# Introduction
This project can build a data set of helices. 

We saved ten high-resolution PDBs in **PDB/standard_set** after downloading them from the PDB website. Then we manually checked these 10 PDBs, extracted the good helices (we hope), and saved the residues' name and sequence in **standard_helices.xlsx**  so that we could edit them easily. We used the script **read_xls.py** to read records from the standard.xlsfile,  then adjusted the helix records in the PDB file based on the records read from the xls file we prepared previously,  and then saved the new modified PDBs in the **manual_set**. 

We also want to get a training set automacally, so we use **train.py** modified helices records in PDBs from **standard_set** and saved the modified PDBs in the **train_set**. 

Finally, we compared the training PDBs **(train_set)** with the PDBs modified manually **(manual_set)** by **evaluate.py** to evaluate the training process **(train.py)**  if it is useful.
# script
## read_xls.py
Reading standard_helices.xlsx, changing the helices in the standard_set PDBs, and saving the updated PDBs to the manual_set directory.

## train.py
PDBs in the manual_set are inputed and modified PDBs are outputted to the training_set directory.  To make helices for the train set, we follow these rules:

1. Merging the ksdssp, cablam and dsssp annotating results.
2. Using phi-psi angles to reduce redundancy residues.
3. Checking each residue's hydrogen bond length and cutting those that don't match the theoretical bond length.

## evaluate.py
comparing PDBs in the manual_set with training_set and evaluatting the result. The result include 5 parameters: helix_ID, Status, length,  name and sequeen of start residue, name and sequeen of end residue. 
For example:

|ID|Status|length|start_res|end_res|
|--|------|------|---------|-------|
|1 | get  |identical|GLU  59 |ALA  59|
|2 | get  |longer   |ALA  22(ASP  20) |HIS  36(SER  35)|


It means that both PDBs have helix 1 and helix 2, but helix 1 from the manual set and the training set are identical, whereas helix 2 from the manual set and the training set are similar but have distinct start and end residues.

There are 3 parameters of  **Status**: **get, more, less**. 

|Parameters |Describetion|
|-----------|------------|
|**get**    |the training PDB get the helix which is also in the manual PDB.
|**more**   |the training PDB has the 'wrong' helix which is not in the manual PDB.
|**less**   |the training PDB lost the helix which is in the manual PDB. 

There are also 5 parameters of **length**: **same, identical, shoter, longer, None**.

|Parameters |Describetion|
|-----------|------------|
|**identical**|tow helices are identical.
|**same**   |the length of two helices are identical but start residues are differnt.
|**shorter**|the helix from PDBs in training set is shorter than manual set.
|**longer** |the helix from structures of training set is longer than manual set. 
|**None**   | can't compare(such as loss).

*The max diffrence of longer and shorter less than 4 reisidues.
# Run script
**read_xls.py:** `phenix.python read_xls.py standard_helices.xlsx`  
**train.py:** `phenix.python train.py PDB/standard_set`  
**evaluate.py:** `phenix.python evaluate.py PDB/training_set PDB/manual_set`
