import os
import subprocess





def getFile(cwd, out):

    flist = []
    nlist = []
    for dir,sub,file in os.walk(cwd):
        if "CGT" and "consensus" in dir:
            flist.append(dir+"/"+"contigs_uniq_consensus.faa")
        if "CGT" in dir and "consensus" not in dir and "sorted" not in dir and "gff" not in dir and "tmp" not in dir:
            nlist.append(os.path.basename(dir))

    for f,val in enumerate(range(len(flist))):
            eggCommand = 'emapper.py -i {} -o {} --cpu 6 --decorate_gff yes'.format(flist[f],os.path.join(out, nlist[f]))
            subprocess.run(eggCommand, shell=True)
            # print(flist[f],val)
            # print(nlist[f],val)
    
if __name__=="__main__":
    cwd = os.path.abspath('.')
    files = os.listdir(cwd)
    getFile(cwd=cwd,out="")