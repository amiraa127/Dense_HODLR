import os
import subprocess
import shlex
import string
import distutils.core
import shutil


def sortFile(filePath):
	f = open(filePath)
	lines = f.readlines()
	map = {}
	for line in lines:
		lineSplit = string.split(line)
		map[int(lineSplit[0])] = lineSplit[1]
	f.close()
	f = open(filePath,'w')
	for key in sorted(map.keys()):
		f.write(str(key)+'\t'+map[key]+'\n')
	f.write('10\t \n') 
	f.close()

# copy results from source dir      
currDir =  os.getcwd()                                                                                                                                                             
dst = currDir+'/ACA';
idx = currDir.find('benchmarks')
src = currDir[:idx]+'build/'+currDir[idx:]+'/ACA';
distutils.dir_util.copy_tree(src,dst)


# cleanup
for root,dirs,fileNames in os.walk(currDir + '/ACA'):
        for f in fileNames:
            currFilePath = os.path.join(root,f)
            if '.svn' in currFilePath or 'iter' in currFilePath:
                    os.remove(currFilePath)
# work
mySet = set([])
for root,dirs,fileNames in os.walk(currDir + '/ACA'):
        for f in fileNames:
            currFilePath = os.path.join(root,f)
	    if 'fixedSize' in currFilePath and '.' not in currFilePath and 'Detail' not in currFilePath:
		    diagValue = (int(f[string.find(f,'_') + 1]))
		    path = currFilePath[:string.rfind(currFilePath,'/')+1]
		    mySet.add((path,diagValue))

myMap = {'quadratic':0,'multiQuadratic':1,'inverseQuadratic':2,'inverseMultiQuadratic':3,'gaussian':4,'exponential':5,'logarithmic':6,'oneOverR':7,'oneOverSqR':8,'logR':9}
for entry in mySet:
	factor_Path = entry[0]+'Factorization_'+str(entry[1]);
	LR_Path     = entry[0]+'LR-Approximation_'+str(entry[1]);
	solve_Path  = entry[0]+'Solve_'+str(entry[1])
	total_Path  = entry[0]+'Total_'+str(entry[1])
	factor_f = open(factor_Path,'w')
	LR_f     = open(LR_Path,'w')
	solve_f  = open(solve_Path,'w')
	total_f  = open(total_Path,'w')

	for key in myMap.keys():
		currFilePath = entry[0]+key+'_'+str(entry[1])
		f = open(currFilePath)
		currMap = {}
		for line in f.readlines():
			lineSplit = string.split(line)
			currMap[lineSplit[0]] = float(lineSplit[1])
		factor_f.write(str(myMap[key])+'\t'+str(currMap['Factorization'])+'\n')
		LR_f.write(str(myMap[key])+'\t'+str(currMap['LR-Approximation'])+'\n')
		solve_f.write(str(myMap[key])+'\t'+str(currMap['Solve'])+'\n')
		total_f.write(str(myMap[key])+'\t'+str(currMap['Total'])+'\n')
		
	#factor_f.write('10\t0\n')
	#LR_f.write('10\t0\n')
	#solve_f.write('10\t0\n')
	#total_f.write('10\t0\n')
       
	factor_f.close()
	LR_f.close()
	solve_f.close()
	total_f.close()
	
	sortFile(factor_Path)
	sortFile(LR_Path)
	sortFile(solve_Path)
	sortFile(total_Path)
       
