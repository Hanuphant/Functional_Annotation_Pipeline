#!/usr/bin/env python3

import os
import sys
import subprocess


def piler_cr_module(input_dir, output_dir):
    # takes the names from the directory and looks for the contigs.fasta file for piler_cr
    cwd = str(os.getcwd())
    sample_list = os.listdir(input_dir)
    for sample in sample_list:
        fullpath = input_dir
        # constructs the full filepath to be used in the subprocess calls
        for root, dirs, files in os.walk(str(os.path.join(input_dir, str(sample)))):
            if 'contigs.fasta' in files:
                fullpath = os.path.join(root, 'contigs.fasta')
        subprocess.call(
            ['pilercr', '-in', str(fullpath), '-out', 'PilerCR_{}.txt'.format(sample)])
        # removes unnecessary lines to allow the CRISPRFileToGFF_1 script to work properly
        subprocess.call(['sed', '-i', '1,75d', 'PilerCR_{}.txt'.format(sample)])
        # moves the CRISPRFileToGFF_1 script to the working directory
        subprocess.call(['mv', 'CRISPRFileToGFF_1.pl', str(cwd)], cwd='/home/groupa/functional_annotation/tools/')
        # generates an empty gff file
        subprocess.call(['touch', 'PilerCR_{}.gff'.format(sample)])
        # generates a gff3 file from the txt file
        subprocess.call(['perl', './CRISPRFileToGFF_1.pl', '-in', 'PilerCR_{}.txt'.format(sample), '-out',
                        'PilerCR_{}.gff'.format(sample)])
        # deletes now unnecessary txt file
        subprocess.call(['rm', 'PilerCR_{}.txt'.format(sample)])
        # adds ##gff-version 3 at the start of each gff file
        with open('PilerCR_{}.gff'.format(sample), 'r+') as file:
            content = file.read()
            file.seek(0, 0)
            file.write('##gff-version 3' '\n' + content)
        # moves gff file into the desired directory; make sure the desired directory is present!
        subprocess.call(['mv', 'PilerCR_{}.gff'.format(sample), str(output_dir)])


if __name__ == '__main__':
    piler_cr_module(sys.argv[1], sys.argv[2])
