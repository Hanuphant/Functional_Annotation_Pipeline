#!/usr/bin/python
import os
import re
import subprocess


def GFF_Merge(indir, outdir, isotaker):
    gffs = []
    for root, dirs, files in os.walk(indir,topdown=False):
        for name in files:
            if re.search(pattern=".gff", string=name):
                gffs.append(os.path.join(root,name))
    reqdirs = []
    for dir in dirs:
        if re.search(r"CARD|Piler|signal|tmhmm|VFDB", dir) is not None:
            reqdirs.append(dir)
    isos = []
    for root2, dirs2, files2 in os.walk(isotaker, topdown=False):
        isos = dirs2
    passer = {}
    for iso, gff in zip(isos * 5,gffs):
        passer[iso] = passer.get(iso, []) + [gff]
    for iso, ps in passer.items():
        p = "--gff " + " --gff ".join(ps)
        command = f"agat_sp_merge_annotations.pl "+p+" --out "+iso+"merged.gff"
        subprocess.check_output(command.split())

if __name__ == "__main__":
    GFF_Merge("..", ".", "/home/groupa/final_gene_prediction")
