# RNA folding problem
### Creation of an objective function for the RNA folding problem

For a given ribonucleotide chain, the RNA folding problem consists in finding the native fold
among the astronomically large number of possible conformations. The native fold being the
one with the lowest Gibbs free energy, the objective function should be an estimator of this
energy.

**This deposit was made in the context of a [university project](https://github.com/Anthony96p/RNA_folding_problem/blob/master/TP_RNA.pdf) in a RNA bioinformatics course 
taught by <cite>_Mr. Guillaume Postic_</cite>.**


## Installation:
````shell
# 1. Clone Git repository
git clone https://github.com/Anthony96p/RNA_folding_problem

# 2. cd into it !
cd RNA_folding_problem/

# 3. Add permissions (optional)
chmod +rwx TP_RNA_main.py TP_RNAScript1.py TP_RNAScript2.py TP_RNAScript3.py

# 4. Install the necessary dependencies using pip:
pip install -r requirements.txt
````

## Quick Start - using `TP_RNA_main.py`
````shell
./TP_RNA_main.py -input PDB -output ScoringValues -savepng True -test 2jyj.pdb
````
_**n.b** : Make sure you are using ``python 3``._

## 1<sup>st</sup> script : Train function
### Explanatory notes :
In first, the RNA folding problems tool will train the objective function using the 
interatomic distance distribution from a 3D RNA structure dataset of train (`TP_RNAScript1.py`).

This script takes as input the name of the directory where all the pdb files for the training 
are located (e.g. : `-input PDB`).

The second required option will be the name of the directory (existing or not) in which 
the scoring value files should be saved (e.g. : ``-output ScoringValues``).

[//]: # (The third option is not mandatory, it allows to display the training time for each pdb file )

[//]: # (&#40;e.g. : ``-details False``&#41;. Default ``False``.)

### Run only Train function :
```bash
./TP_RNAScript1.py -input PDB -output ScoringValues -details False
```

## 2<sup>nd</sup> script : Plot function
### Explanatory notes :
In the second part, the scoring profile of each pair of residues will be plotted as a function of
the interatomic distance (`TP_RNAScript2.py`).

The only required option will be the name of the directory in which the scoring values files have 
been saved by the script 1  (e.g. : `-input ScoringValues`).

The script 2 has two non required options: :
- The first one allows to save or not the plots and asks as input a boolean (e.g. : ``-savepng False``). 
Default ``False``. Save in the input folder.
- The second one allows to display or not the plots and asks for a boolean as input (e.g. : ``-printpng True``).
Default ``True``.

### Run only Plot function :
```bash
./TP_RNAScript2.py -input ScoringValues -savepng False -printpng True
```

### Expected output :
<img src="https://github.com/Anthony96p/RNA_folding_problem/blob/master/ScoringValues/AA.png" alt=""/>

## 3<sup>rd</sup> script : Test function
### Explanatory notes :
In the third part, we use the objective function to evaluate the predicted structures from another 
3D RNA structure dataset of test (`TP_RNAScript3.py`).

The first required option will be the name of the directory in which the scoring values files have 
been saved by the script 1  (e.g. : `-train ScoringValues`).

The second required option will be the name of the pdb file that contains the 3D structure 
of the test dataset (e.g. : `-test 2jyj.pdb`)

### Run only Test function :
```bash
./TP_RNAScript3.py -train ScoringValues -test 2jyj.pdb
```

## Main script : All functions
### Explanatory notes :
This is the main script that will run the 3 functions (train, plot and test) with a 
single command (`TP_RNA_main.py`).

This script takes as input the name of the directory where all the pdb files for the training 
are located (e.g. : `-input PDB`).

The second required option will be the name of the directory (existing or not) in which 
the scoring value files and plot should be saved (e.g. : ``-output ScoringValues``).

The third required option will be the name of the pdb file that contains the 3D structure 
of the test dataset (e.g. : `-test 2jyj.pdb`)

The script 2 has two non required options:

[//]: # (- The first option allows to display the training time for each pdb file &#40;e.g. : ``-details False``&#41;. )

[//]: # (Default ``False``.)
- The first one allows to save or not the plots and asks as input a boolean (e.g. : ``-savepng False``). 
Default ``False``. Save in the input folder.
- The second one allows to display or not the plots and asks for a boolean as input (e.g. : ``-printpng True``).
Default ``True``.

### Run all functions :
```bash
./TP_RNA_main.py -input PDB -output ScoringValues -details False -savepng False -printpng True -test 4gxy.pdb
```

### Expected output :
Returns a file (``Results_RNA_folding_problem.txt``) which will contain the final gibbs free energy and all 
the process information.

## Others :
### Fonctions script :
Grouping of redundant functions between scripts :

The distance calculation function (`calcul_dist()`) shared by script 1 and 3 and the peer formalization 
function (`pair_res_format()`) shared by script 1 and 3.

**<font size="4"> ⚠ </font> Cannot be executed alone** 

### Constraints

* Only C3' atoms are taken into account
* Only "intrachain" distances are considered
* Only consider residues separated by at least 3 positions on the sequence

### Observed Frequency

The probability (i.e. frequency) of observing two residues *i* and *j* separated by distance bin *r* is calculated as follows:

$$ f_{i,j} ^{OBS}(r) = { N_{i,j}(r) \over N_{i,j} } $$

where $N_{i,j(r)}$ is the count of i and j within the distance bin r, and Ni,j is the count of i and j for
all distance bins. Only distance intervals of 0 to 20 Å are taken into account.

### Reference Frequency

The reference frequency is the same formula except that the different residue types (A, U, C, G) are indistinct ("X"):

$$ f_{X,X} ^{REF}(r) = { N_{X,X}(r) \over N_{X,X}} $$

### Pseudo-Energy

The score for each residue pair is computed as followed:

$$ u_{i,j}(r) = { -log \left( f _{i,j} ^{OBS}(r) \over f_{i,j} ^{REF}(r) \right) } $$


### To install :
- ``numpy``
- ``tqdm``
- ``matplotlib.pyplot``
- ``colorhash``

### Requirements :
- ``colorhash==1.2.1``
- ``cycler==0.11.0``          
- ``fonttools==4.38.0``       
- ``kiwisolver==1.4.4``       
- ``matplotlib==3.5.3``       
- ``numpy==1.21.6``           
- ``packaging==21.3``         
- ``Pillow==9.3.0``           
- ``pyparsing==3.0.9``        
- ``python-dateutil==2.8.2``
- ``tqdm==4.64.1``
- ``six==1.16.0  ``           
- ``typing_extensions==4.4.0``

# Author
PRAGASSAM Anthony

