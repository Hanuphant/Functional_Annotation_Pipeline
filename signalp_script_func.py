#!/usr/bin/env python3

import os
from os import listdir 

# os.mkdir(path) # creating the output directory
# mypath = "/home/groupa/final_gene_prediction"
"The tool needs to be added to the path while running the script"
"The command does not have an output flag, hence the scripts moves the required outputs to signalp5_outputs directory"
"For merging, only the output .gff3 files need to be considered"

#get the names of the directories in a list


def signalP(mypath, path):
	dirs = listdir(mypath)
	try:
		os.mkdir(path)
	except:
		print("Folder already exists")
	for i in range(len(dirs)):
			os.system(f" signalp -fasta '{mypath}/{dirs[i]}/consensus_and_fasta/contigs_uniq_consensus.faa' -org gram- -format short -gff3 -prefix sig_'{dirs[i]}' ")
	os.system(f"mv sig_* '{path}'")

if __name__=="__main__":
	directory_name = "signalp5_outputs"
	parent_dir = "/home/groupa/functional_annotation/tools/tool_outputs"
	path = os.path.join(parent_dir,directory_name)
	signalP("/home/groupa/final_gene_prediction", path)

