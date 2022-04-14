#!/usr/bin/env python3

import os
import sys
import subprocess


def vfdb_module(input_dir, output_dir):
    # takes the names from the directory and looks for the contigs_uniq_consensus.fna file for blastn
    cwd = str(os.getcwd())
    sample_list = os.listdir(input_dir)
    for sample in sample_list:
        fullpath = input_dir
        # constructs the full filepath to be used in the subprocess calls
        for root, dirs, files in os.walk(str(os.path.join(input_dir, str(sample)))):
            if 'contigs_uniq_consensus.fna' in files:
                fullpath = os.path.join(root, 'contigs_uniq_consensus.fna')
        # The generated VFDB database files are located within the ncbi folder
        subprocess.call(
            ['blastn', '-db', '/home/groupa/functional_annotation/ncbi-blast-2.13.0+/bin/VFDB/VFDB', '-query',
             str(fullpath), '-num_threads', '4', '-evalue', '1e-10', '-outfmt', '6', '-max_target_seqs', '1', '-out',
                'VFDB_{}.txt'.format(sample)])
        # generates an empty gff file in the working script's directory
        subprocess.call(['touch', 'VFDB_{}.gff'.format(sample)])
        # moves the blast2gff.py file from the tools directory to the working script's directory
        subprocess.call(['mv', 'blast2gff.py', cwd], cwd='/home/groupa/functional_annotation/tools')
        # generates a gff3 file from the blast txt file, with the stdout getting written into the file
        with open('VFDB_{}.gff'.format(sample), "w") as out_file:
            subprocess.call(['./blast2gff.py', '-b', 'VFDB_{}.txt'.format(sample)], stdout=out_file)
        # removes the now unnecessary txt file
        subprocess.call(['rm', 'VFDB_{}.txt'.format(sample)])
        # moves the blast2gff.py script back to the tools directory
        subprocess.call(['mv', 'blast2gff.py', '/home/groupa/functional_annotation/tools'])
        # adds ##gff-version 3 at the start of each gff file
        with open('VFDB_{}.gff'.format(sample), 'r+') as file:
            content = file.read()
            file.seek(0, 0)
            file.write('##gff-version 3' '\n' + content)
        # moves gff file into the desired directory; make sure the desired directory is present!
        subprocess.call(['mv', 'VFDB_{}.gff'.format(sample), str(output_dir)])


if __name__ == '__main__':
    vfdb_module(sys.argv[1], sys.argv[2])
