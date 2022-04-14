# !/usr/bin/env python

import os

def tmhmm(input_directory_path, output_directory_path):
    for file in os.listdir(input_directory_path):
        f = os.path.join(input_directory_path, file)
        # print(f"\n\n\n\n\n\nRunning tmhmm on {f}")
        # for isolate in os.listdir(f):
        # print(f"tmhmm {f}/consensus_and_fasta/contigs_uniq_consensus.faa {output_directory_path}/{file}_tmhmm.out")
        os.system(f"tmhmm '{f}/consensus_and_fasta/contigs_uniq_consensus.faa' > '{output_directory_path}/{file}.out' ")



def converttogff(path):
    for output in os.listdir(path):
        tmhmm_out = open(path + "/" + output, "r")
        tmhmm_gff = []
        for line in tmhmm_out:
            row = line.strip().split()
            tmhmm_gff.append(row)
        source = "TMHMM2.0"
        typ = "TMhelix"
        score = "."
        phase = "."
        attributes = "."
        strand = "."


        for i in range(len(tmhmm_gff)):
            if tmhmm_gff[i][0].startswith("#"):
                continue
            if tmhmm_gff[i][2] == "TMhelix":
                seqid = tmhmm_gff[i][0]
                start = int(tmhmm_gff[i][3])
                end = start + int(tmhmm_gff[int(i)][4])


                with open(output[:-4] + ".gff", 'a') as f:
                    f.write(seqid + "\t" + source + "\t" + typ + "\t" + str(start) + "\t" + str(end) + "\t" + score + "\t" + strand + "\t" + phase + "\t" + attributes + "\n")
                f.close()
if __name__ == "__main__":
    input_directory_path = "/home/groupa/final_gene_prediction/"
    output_directory_path = "/home/groupa/functional_annotation/tools/tool_outputs/tmhmm_results"
    tmhmm(input_directory_path)
    converttogff(output_directory_path)
