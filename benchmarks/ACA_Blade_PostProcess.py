import os
import subprocess
import shlex
import string
import distutils.core
import shutil


# copy results from source dir      
currDir =  os.getcwd()                                                                                                                                                             
dst = currDir+'/ACA_Blade';
idx = currDir.find('benchmarks')
src = currDir[:idx]+'build/'+currDir[idx:]+'/ACA_Blade';
distutils.dir_util.copy_tree(src,dst)

data = {}
sizeSet = set([]);
tolSet  = set([]);

for root,dirs,fileNames in os.walk(currDir + '/ACA_Blade'):
        for f in fileNames:
            currFilePath = os.path.join(root,f)
            if 'Tolerance' in f:
                sizeIdx = string.find(f,'_',6)
                size = int(f[6:sizeIdx-1])
                sizeSet.add(size);
                currFile = open(currFilePath)
                for line in currFile.readlines():
                    lineSplit = string.split(line)
                    tol = float(lineSplit[0])
                    tolSet.add(tol)
                    key = (size,tol)
                    if key in data:
                        val = data[key]
                        if 'error' in f:
                            val[0] = float(lineSplit[1])
                        if 'speed' in f:
                            val[1] = float(lineSplit[1])
                    else:
                        val=[0,0]
                        if 'error' in f:
                            val[0] = float(lineSplit[1])
                        if 'speed' in f:
                            val[1] = float(lineSplit[1])
                    data[key] = val;
                currFile.close();

actualSize = {9:9239,12:12936,16:16632}
sizeSet = sorted(sizeSet)
tolSet  = sorted(tolSet)
# create SpeedvsSizePlots
for tol in tolSet:
    fileName = 'blade_speedVsSize_'+str(tol)
    filePath = os.path.join(root,fileName)
    f = open(filePath,'w')
    for size in sizeSet:
        key = (size,tol)
        speed = data[key][1]
        f.write(str(actualSize[size])+'\t'+str(speed)+'\n')
    f.close()
