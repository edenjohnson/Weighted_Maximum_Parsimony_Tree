# Constructing Weighted Maximum Parsimony Tree of an MSA using the Sankoff Algorithm

This program uses BioPython to construct the weighted maximum parsimonious tree. 
This program also returns the closest OTU(s) for a given OTU.


### Prerequisites/Assumptions

Dependencies:
Python 3+ must be installed on the system running the program. 

BioPython must be installed on a user's system before proceeding with this pipeline. 
(MacOS command: pip3 install biopython)

MatPlotLib must be installed on a user's system before proceeding with this pipeline.
(MacOS command: pip3 install matplotlib)

All .py files must be in the same directory in order for this program to execute correctly. The 
clustal omega binary file MUST be present in the same folder as MSA.py in order to perform
an MSA locally. 

Note:
This program performs a BLASTn search over the internet and an MSA locally. Due to these processes, it should be noted that these programs will take a while to execute. (Ex. BLASTn takes ~ 15 mins,MSA takes ~ 30 mins) 


## Author

**Eden Johnson** 