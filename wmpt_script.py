import os

os.system("python3 BLAST.py")
os.system("python3 parse_sequences.py")
os.system("python3 MSA.py")
os.system("python3 tree_building.py")
os.system("python3 closest_OTU.py")