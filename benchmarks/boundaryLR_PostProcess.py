import os
from subprocess import call
import subprocess
import shlex
import string 
import distutils.core
import shutil

def uploadToDropBox(currDir):
    for root,dirs,fileNames in os.walk(currDir):
        for f in fileNames:
            currFilePath = os.path.join(root,f)
            if (('tex' in currFilePath or 'pdf' in currFilePath or 'result' in currFilePath) and 'ACA' not in currFilePath and './' not in currFilePath and 'log' not in currFilePath):
                src = currFilePath;
                if 'boundaryLR' in currFilePath:
                    dst = '/Users/Amir/Dropbox/Amir_Siva/'+currFilePath[string.find(currFilePath,'benchmarks/')+11:]
                else:
                    dst = '/Users/Amir/Dropbox/Amir_Siva/boundaryLR/'+currFilePath[string.find(currFilePath,'benchmarks/')+11:]
                shutil.copyfile(src,dst)

def getTime(stream,token,startIdx):
    timeStartIdx = string.find(stream,token,startIdx);
    timeEndIdx   = string.find(stream,'seconds',timeStartIdx);
    totalTime = string.split(stream[timeStartIdx:timeEndIdx - 1])
    totalTime = totalTime[len(totalTime) - 1];
    return float(totalTime);

def getInfo(stream,token,startIdx):
    infoStartIdx = string.find(stream,token,startIdx);
    infoEndIdx   = string.find(stream,'\n',infoStartIdx)
    info         = stream[infoStartIdx:infoEndIdx]
    info         = string.split(info)
    info         = info[len(info) - 1]
    return info

def insertNumPoints(filePath,numPointsVec):
    f    = open(filePath)
    data = f.readlines()
    f.close();
    f = open(filePath,'w')
    for i in range(len(numPointsVec)):
        f.write(string.split(numPointsVec[i])[1] +'\t'+ string.split(data[i])[1]+'\n')
    f.close();

def parseResultFile(path,nameList):
    f = open(path)
    data = f.read()
    f.close();
    
    size = {}
    directAccuracy = {}
    directTime     = {}
    finalAccuracy  = {} 
    totalTime  = {}
    LUTime     = {}
    numIter    = {}
    for matrixPath in nameList:
        matrix = matrixPath[matrixPath.find('input/') + 6:]
        startIdx = string.find(data,matrix)

        size[matrixPath]          = int(getInfo(data,'Matrix Size',startIdx))
        numIter[matrixPath]       = int(getInfo(data,'Number of Iterations',startIdx))
        finalAccuracy[matrixPath] = float(getInfo(data,'Residual l2 Relative Error',startIdx))
        
        directAccuracyStartIdx = startIdx + len(matrix) + 3
        directAccuracyEndIdx   = string.find(data,'\n',directAccuracyStartIdx)
        directAccuracy[matrixPath] = float(data[directAccuracyStartIdx:directAccuracyEndIdx])
  
        totalTime[matrixPath] = getTime(data,'Total Solve Time',startIdx)
        directTime[matrixPath] = getTime(data,'HODLR Direct Solver Total Time',startIdx)
        LUTime[matrixPath] = getTime(data,'LU Solve Time',startIdx)
  
    return size,directAccuracy,directTime,finalAccuracy,totalTime,LUTime,numIter

def correctTimingCaption(inputFilePath,size):
    matrixName       = inputFilePath[inputFilePath.find('input/') + 6:]
    timingFileName   = matrixName + '_Timing_Plot.tex';
    root = inputFilePath[0:inputFilePath.find('input/')]
    timingFolderPath = root + 'timing_Plots/'
    timingFilePath   = timingFolderPath + timingFileName;
    f = open(timingFilePath)
    data = f.read()
    f.close()
    startIdx = string.find(data,'size of');
    endIdx   = string.find(data,'.',startIdx);
    correctSize = size[inputFilePath];
    correctData = data[:startIdx+8]+str(correctSize)+data[endIdx:]
    f = open(timingFilePath,'w')
    f.write(correctData)
    f.close()
    os.chdir(timingFolderPath)
    proc = subprocess.Popen(shlex.split('pdflatex ' + timingFileName))
    proc.communicate()
    os.unlink(timingFileName[:-4] + '.aux')
    os.unlink(timingFileName[:-4] + '.log')
    
def getMatrixInfo(inputFilePath):
    matrixName = inputFilePath[inputFilePath.find('input/') + 6:]
    nameSplit  = string.split(matrixName,'_');
    if 'blade' in matrixName:
        matrixLevel = 'top'
    else:
        level = int(nameSplit[5]);
        if level == 0:
            matrixLevel = 'top'
        elif level == 1:
            matrixLevel = 'first'
        elif level == 2:
            matrixLevel = 'second'
        else:
            matrixLevel = str(level)+'th'
        
    if 'blade' in matrixName:
        if '9' in matrixName:
            sparseSize = '642K'
        if '12' in matrixName:
            sparseSize = '1.26M'
        if '16' in matrixName:
            sparseSize = '2.14M'
    else:
        sparseSize = str(int(nameSplit[0]))+'K'   
            
    if 'FETI' in inputFilePath:
        matrixType = 'FETI local sparse'
    elif 'stiffness' in inputFilePath:
        matrixType = 'stiffness';

    if '/structured/' in inputFilePath:
        meshType = 'a structured'
    elif '/unstructured' in inputFilePath:
        meshType = 'an unstructured'
        
    if '2D' in inputFilePath:
        meshDim = '2D'
    else:
        meshDim = '3D'

    return matrixLevel,sparseSize,matrixType,meshType,meshDim

def generateLRAnalysisFigure(inputFilePath):
    root = inputFilePath[0:inputFilePath.find('input/')]
    LRFolderPath = root + 'LR_Analysis_Plots/'
    matrixName   = inputFilePath[inputFilePath.find('input/') + 6:]
    LRFileName   = matrixName + '_LR_Plot.tex';
    SVDFileName  = matrixName + '_SingularValueDecay';
    distanceFileName_Col       = matrixName + '_distanceFromBoundaryVsPivotSize_Col'
    distanceFileName_Row       = matrixName + '_distanceFromBoundaryVsPivotSize_Row'
    SVDErrorFileName           = matrixName + '_SVDErrorVsBoundaryDistance'
    boundaryErrorFileName      = matrixName + '_boundaryErrorVsBoundaryDistance'
    boundaryErrorCompFileName1 = matrixName + '_boundaryErrorCompVsBoundaryDistance1'
    boundaryErrorCompFileName3 = matrixName + '_boundaryErrorCompVsBoundaryDistance3'
    boundaryErrorCompFileName5 = matrixName + '_boundaryErrorCompVsBoundaryDistance5'
    numPointsFileName          = matrixName + '_numPointsVsBoundaryDistance'
    SVDFilePath                = root + 'results_LR/' + SVDFileName;
    distanceFilePath_Col       = root + 'results_LR/' + distanceFileName_Col;
    distanceFilePath_Row       = root + 'results_LR/' + distanceFileName_Row;
    SVDErrorFilePath           = root + 'results_LR/' + SVDErrorFileName;
    boundaryErrorFilePath      = root + 'results_LR/' + boundaryErrorFileName;
    boundaryErrorCompFilePath1 = root + 'results_LR/' + boundaryErrorCompFileName1;
    boundaryErrorCompFilePath3 = root + 'results_LR/' + boundaryErrorCompFileName3;
    boundaryErrorCompFilePath5 = root + 'results_LR/' + boundaryErrorCompFileName5;
    LRFilePath                 = LRFolderPath + LRFileName;  
    numPointsFilePath          = root + 'results_LR/' + numPointsFileName;
    width = str(8)
    
    numPoints_f = open(numPointsFilePath)
    numPoints   = numPoints_f.readlines();
    numPoints_f.close();
    
    insertNumPoints(SVDErrorFilePath,numPoints)
    insertNumPoints(boundaryErrorFilePath,numPoints)
    insertNumPoints(boundaryErrorCompFilePath1,numPoints)
    insertNumPoints(boundaryErrorCompFilePath3,numPoints)
    insertNumPoints(boundaryErrorCompFilePath5,numPoints)
    #SVDError_f = open(SVDErrorFilePath)
    #SVDError   = SVDError_f.readlines()
    #SVDError_f.close();
    
    #SVDError_f = open(SVDErrorFilePath,'w')
    #for i in range(len(numPoints)):
    #    SVDError_f.write(string.split(numPoints[i])[1] +'\t'+ string.split(SVDError[i])[1]+'\n')
    #SVDError_f.close();

    #boundaryError_f = open(boundaryErrorFilePath)
    #boundaryError   = boundaryError_f.readlines()
    #boundaryError_f.close();

    #boundaryError_f = open(boundaryErrorFilePath,'w')
    #for i in range(len(numPoints)):
    #    boundaryError_f.write(string.split(numPoints[i])[1] +'\t'+ string.split(boundaryError[i])[1]+'\n')
    #boundaryError_f.close();

    #boundaryErrorComp_f = open(boundaryErrorCompFilePath1)
    #boundaryErrorComp   = boundaryErrorComp_f.readlines()
    #boundaryErrorComp_f.close();
    
    #boundaryErrorComp_f = open(boundaryErrorCompFilePath1,'w')
    #for i in range(len(numPoints)):
    #    boundaryErrorComp_f.write(string.split(numPoints[i])[1] +'\t'+ string.split(boundaryErrorComp[i])[1]+'\n')
    #boundaryErrorComp_f.close();
    
    SVD_f = open(SVDFilePath)
    denseSize = str(len(SVD_f.readlines()));
    SVD_f.close();

    f = open(LRFilePath,'w');
    f.write('\\documentclass[12pt,letterpaper]{article}\n')
    f.write('\n')
    f.write('\\usepackage{tikz}\n')
    f.write('\\usepackage{graphicx,subfigure}\n')
    f.write('\\usepackage{pgfplots}\n')
    f.write('\\usepackage[margin=0.5in]{geometry}\n')
    f.write('\\pgfplotsset{compat=1.3}\n')
    f.write('\n')
    f.write('\\begin{document}\n')
    f.write('\\pagenumbering{gobble}\n')
    f.write('\\begin{figure}\n')
    f.write('\\centering\n')
    f.write('\\subfigure[Row distance from boundary vs pivot size]{\n');
    f.write('\t\\begin{tikzpicture}\n');
    f.write('\t\t\\pgfplotsset{every axis legend/.append style={ at={(0.5,1.03)},anchor=south}}\n')
    f.write('\t\t\\begin{axis}[width='+width+'cm,height=7cm,grid=major,xlabel=Pivot size (largest to smallest),ylabel=Row distance from boundary]\n');
    f.write('\t\t\t\\addplot file{"'+distanceFilePath_Row+'"};\n');
    f.write('\t\t\\end{axis}\n')
    f.write('\t\\end{tikzpicture}\n')
    f.write('}\n')
    f.write('\\subfigure[Col distance from boundary vs pivot size]{\n');
    f.write('\t\\begin{tikzpicture}\n');
    f.write('\t\t\\pgfplotsset{every axis legend/.append style={ at={(0.5,1.03)},anchor=south}}\n')
    f.write('\t\t\\begin{axis}[width='+width+'cm,height=7cm,grid=major,xlabel=Pivot size (largest to smallest),ylabel=Col distance from boundary]\n');
    f.write('\t\t\t\\addplot file{"'+distanceFilePath_Col+'"};\n');
    f.write('\t\t\\end{axis}\n')
    f.write('\t\\end{tikzpicture}\n')
    f.write('}\n')
    f.write('\\subfigure[Singular value decay]{\n');
    f.write('\t\\begin{tikzpicture}\n');
    f.write('\t\t\\pgfplotsset{every axis legend/.append style={ at={(0.5,1.03)},anchor=south}}\n')
    f.write('\t\t\\begin{semilogyaxis}[width='+width+'cm,height=7cm,grid=major,xlabel=Singular value index, ylabel=Singular values]\n');
    f.write('\t\t\t\\addplot file{"'+SVDFilePath+'"};\n')
    f.write('\t\t\\end{semilogyaxis}\n');
    f.write('\t\\end{tikzpicture}\n');
    f.write('}\n')
    f.write('\\subfigure[Error vs NumPoints]{\n');
    f.write('\t\\begin{tikzpicture}\n');
    f.write('\t\t\\begin{semilogyaxis}[width='+width+'cm,height=7cm,grid=major,xlabel=Number of Added Points, ylabel=Relative Error,legend style={font=\\tiny}]\n');
    f.write('\t\t\t\\addplot file{"'+SVDErrorFilePath+'"};\n')
    f.write('\t\t\t\\addplot file{"'+boundaryErrorFilePath+'"};\n')
    f.write('\t\t\t\\addplot file{"'+boundaryErrorCompFilePath1+'"};\n')
    f.write('\t\t\t\\addplot file{"'+boundaryErrorCompFilePath3+'"};\n')
    f.write('\t\t\t\\addplot file{"'+boundaryErrorCompFilePath5+'"};\n')
    f.write('\t\t\t\\legend{SVD,RRLU\\_Boundary\\_1e-15,RRLU\\_Boundary\\_1e-1,RRLU\\_Boundary\\_1e-3,RRLU\\_Boundary\\_1e-5}\n')
    f.write('\t\t\\end{semilogyaxis}\n');
    f.write('\t\\end{tikzpicture}\n');
    f.write('}\n')
    matrixLevel,sparseSize,matrixType,meshType,meshDim = getMatrixInfo(inputFilePath)
    f.write('\\caption{Figures are for the top off-diagonal block of the '+matrixLevel+' level frontal matrix of a '+sparseSize+' '+matrixType+' matrix corresponding to '+meshType+' '+meshDim+' mesh. The off diagonal matrix has a size of '+denseSize+'.}\n')
    f.write('\\end{figure}\n')
    f.write('\\end{document}\n')
    f.close();
    os.chdir(LRFolderPath)
    proc = subprocess.Popen(shlex.split('pdflatex ' + LRFileName))
    proc.communicate()
    os.unlink(LRFileName[:-4] + '.aux')
    os.unlink(LRFileName[:-4] + '.log')

def generateTimeAnalysisFigure(inputFilePath):
    root = inputFilePath[0:inputFilePath.find('input/')]
    timingFolderPath = root + 'timing_Plots/'
    matrixName       = inputFilePath[inputFilePath.find('input/') + 6:]
    timingFileName   = matrixName + '_Timing_Plot.tex';
    timingFilePath   = timingFolderPath + timingFileName;
    recLUFactorFileName  = matrixName + '_recLU_FactorTiming'
    recLUSolveFileName   = matrixName + '_recLU_SolveTiming'
    LRTimingFileName     = matrixName + '_LR_Timing'
    levelRankFileName    = matrixName + '_levelRankAverage'
    iterTimingFileName   = matrixName + '_iter_IterTiming'
    iterAccuracyFileName = matrixName + '_iter_Accuracy' 
    recLUFactorFilePath  = root + 'results_timing/' + recLUFactorFileName;
    recLUSolveFilePath   = root + 'results_timing/' + recLUSolveFileName;
    LRTimingFilePath     = root + 'results_timing/' + LRTimingFileName;
    levelRankFilePath    = root + 'results_timing/' + levelRankFileName;
    iterTimingFilePath   = root + 'results_timing/' + iterTimingFileName;
    iterAccuracyFilePath = root + 'results_timing/' + iterAccuracyFileName;
    width     = str(16)
    iterWidth = str(16)

    recLUFactor_f = open(recLUFactorFilePath)
    lines         = recLUFactor_f.readlines();
    numLevels     = int(string.split(lines[len(lines) - 1])[0])
    val = float(string.split(lines[len(lines) - 2])[1])
    recLUFactor_f.close();
    
    if val > 0:
        recLUFactor_f = open(recLUFactorFilePath,'a')
        recLUFactor_f.write(str(numLevels + 1) +' ' +'0\n')
        recLUFactor_f.close();

    recLUSolve_f = open(recLUSolveFilePath)
    lines        = recLUSolve_f.readlines();
    numLevels    = int(string.split(lines[len(lines) - 1])[0])
    val = float(string.split(lines[len(lines) - 2])[1])
    recLUSolve_f.close();

    if val > 0:
        recLUSolve_f = open(recLUSolveFilePath,'a')
        recLUSolve_f.write(str(numLevels + 1) +' ' +'0\n')
        recLUSolve_f.close();

    LRTiming_f = open(LRTimingFilePath)
    lines      = LRTiming_f.readlines();
    numLevels  = int(string.split(lines[len(lines) - 1])[0])
    val = float(string.split(lines[len(lines) - 2])[1])
    LRTiming_f.close();

    if val > 0:
        LRTiming_f = open(LRTimingFilePath,'a')
        LRTiming_f.write(str(numLevels + 1) +' ' +'0\n')
        LRTiming_f.close();
    
    iter_f  = open(iterAccuracyFilePath)
    lines   = iter_f.readlines();
    numIter = str(len(lines) - 1)
    iter_f.close();

    f = open(timingFilePath,'w');
    f.write('\\documentclass[12pt,letterpaper]{article}\n')
    f.write('\n')
    f.write('\\usepackage{tikz}\n')
    f.write('\\usepackage{graphicx,subfigure}\n')
    f.write('\\usepackage{pgfplots}\n')
    f.write('\\usepackage[margin=0.5in]{geometry}\n')
    f.write('\\pgfplotsset{compat=1.3}\n')
    f.write('\n')
    f.write('\\begin{document}\n')
    f.write('\\pagenumbering{gobble}\n')
    f.write('\\begin{figure}\n')
    f.write('\\subfigure[Detail timings and rank information for direct HODLR solver]{\n')
    f.write('\t\\begin{tikzpicture}\n')
    f.write('\t\\pgfplotsset{every axis legend/.append style={at={(0.5,1.03)},anchor=south},y axis style/.style={yticklabel style=#1,ylabel style=#1,y axis line style=#1,ytick style=#1}}\n')
    f.write('\t\t\\begin{axis}[scale only axis, axis y line*=left,ybar,ybar interval = 0.7,enlargelimits=0.05,width = '+width+'cm, height =7cm,legend style={at={(0.5,1.03)},anchor=south,legend columns=3},ylabel={Time (s)},xlabel=HODLR Level]\n')
    f.write('\t\t\t\\addplot table [x index=0,y index=1,header=false] {' +recLUFactorFilePath + '};\n')
    f.write('\t\t\t\\addplot table [x index=0,y index=1,header=false] {' +recLUSolveFilePath + '};\n')
    f.write('\t\t\t\\addplot table [x index=0,y index=1,header=false] {' +LRTimingFilePath + '};\n')
    f.write('\t\t\t\\legend{Factorization, Solve, Low-Rank Approximation}\n')
    f.write('\t\t\\end{axis}\n')
    f.write('\t\t\\begin{axis}[scale only axis, axis y line*=right,y axis style=blue!75!black,width='+width+'cm, height =7cm,ylabel=Average Rank,xtick=\empty]\n')
    f.write('\t\t\t\\addplot [line legend, thick, sharp plot, stack plots=false,color=blue,mark=*] table[y index =1]{'+levelRankFilePath+'};\n');
    f.write('\t\t\\end{axis}\n')
    f.write('\t\\end{tikzpicture}\n')
    f.write('}\n')
    f.write('\\subfigure[Iteration timing and accuracy ]{\n')
    f.write('\t\\begin{tikzpicture}\n')
    f.write('\t\\pgfplotsset{y axis style/.style={yticklabel style=#1,ylabel style=#1,y axis line style=#1,ytick style=#1}}\n')
    f.write('\t\t\\begin{axis}[scale only axis, axis y line*=left,y axis style=blue!75!black,width='+iterWidth+'cm,height=7cm,grid=major,xlabel=Iteration,ylabel=Time(s),xtick=data,legend columns=4,xmin=0,xmax='+numIter+']\n')
    f.write('\t\t\t\\addplot file{"'+iterTimingFilePath +'"};\n');
    f.write('\t\t\\end{axis}\n')
    f.write('\t\t\\begin{semilogyaxis}[scale only axis, axis y line*=right,y axis style=red!75!black,width='+iterWidth+'cm,height=7cm,ylabel=Relative Error,xtick=data,xmin=0,xmax='+numIter+']\n')
    f.write('\t\t\t\\addplot[color = red, mark = *] file{"'+iterAccuracyFilePath +'"};\n')
    f.write('\t\t\\end{semilogyaxis}\n')
    f.write('\t\\end{tikzpicture}\n')
    f.write('}\n')
    matrixLevel,sparseSize,matrixType,meshType,meshDim = getMatrixInfo(inputFilePath)
    f.write('\t\\caption{Figures are for the '+ matrixLevel +' level frontal matrix of a '+sparseSize+' '+matrixType+' matrix corresponding to '+meshType+' '+meshDim+' mesh. The dense frontal matrix has a size of 7K. Iteration 0 on the iteration plot corresponds to the direct HODLR solver accuracy. The pivot threshold for the direct solver is 0.1. That is, in obtaining the low-rank approximations, we only keep rows and columns that correspond to pivots which have a magnitude greater than 0.1 times the magnitude of the largest pivot. In computing the low-rank approximations, we only considered the points with distance at most one from the boundary.}\n')
    f.write('\end{figure}\n')
    f.write('\end{document}\n')
    f.close()
    os.chdir(timingFolderPath)
    proc = subprocess.Popen(shlex.split('pdflatex ' + timingFileName))
    proc.communicate()
    os.unlink(timingFileName[:-4] + '.aux')
    os.unlink(timingFileName[:-4] + '.log')
    

currDir =  os.getcwd()
# copy results from source dir
dst = currDir+'/boundaryLR';
idx = currDir.find('benchmarks')
src = currDir[:idx]+'build/'+currDir[idx:]+'/boundaryLR';
#distutils.dir_util.copy_tree(src,dst)
nameList = []

# loop through all the files
for root,dirs,fileNames in os.walk(currDir):
    for f in fileNames:
        currFilePath = os.path.join(root,f)
        if 'input/' in currFilePath:
            currResultPath = currFilePath[0:currFilePath.find('input/')]
            timingFolderPath = currResultPath + 'timing_Plots';
            LRFolderPath     = currResultPath + 'LR_Analysis_Plots';
            if not os.path.exists(LRFolderPath):
                os.makedirs(LRFolderPath)
            if not os.path.exists(timingFolderPath):
                os.makedirs(timingFolderPath)
            if ('Graph' not in f):
                if ('level' in f):
                    nameList.append(currFilePath);
                    generateLRAnalysisFigure(currFilePath)
                    generateTimeAnalysisFigure(currFilePath);


# copy output from source dir
src = currDir[:idx]+'build/'+currDir[idx:]+'/results'
dst = currDir+'/results'
#shutil.copyfile(src,dst);
size,directAccuracy,directTime,finalAccuracy,totalTime,LUTime,numIter = parseResultFile(dst,nameList);
# Open summary table file
sumFile = open(currDir+'/summary.tex','w')
sumFile.write('\\documentclass[12pt,letterpaper]{article}\n')
sumFile.write('\n')
sumFile.write('\\usepackage{tikz}\n')
sumFile.write('\\usepackage{graphicx,subfigure}\n')
sumFile.write('\\usepackage{pgfplots}\n')
sumFile.write('\\usepackage[margin=0.5in]{geometry}\n')
sumFile.write('\\pgfplotsset{compat=1.3}\n')
sumFile.write('\\usepackage{multirow}\n')
sumFile.write('\n')
sumFile.write('\\begin{document}\n')
sumFile.write('\\begin{table}\n') 
sumFile.write('\t\\scalebox{0.75}{\n')
sumFile.write('\t\\begin{tabular}{|c|c|c|c||c|c|c|c|c|c|c|c|c|}\n')
sumFile.write('\t\\hline\n')
sumFile.write('\t Matrix & Mesh & Mesh & Elimination & \\multicolumn{2}{c}{Matrix Size}\\vline & \\multicolumn{2}{c}{Direct Solver}\\vline & \\multicolumn{3}{c}{Iterative Solver} \\vline&  Partial Pivoting & Percent\\\  \\cline{5-12}\n')
sumFile.write('\t Type & Type & Dim & Tree Level & Sparse & Dense& Time & Accuracy & Time & Accuracy & num Iterations &LU Time & Speedup\\\ \\hline\n')
# Loop over all files and fill the table
for root,dirs,fileNames in os.walk(currDir):
    for f in fileNames:
        currFilePath = os.path.join(root,f)
        if 'input/' in currFilePath:
            if ('Graph' not in f):
                if ('level' in f):
                    # Correct caption
                    correctTimingCaption(currFilePath,size)
                    
                    # Write to the table file
                    matrixLevel,sparseSize,matrixType,meshType,meshDim = getMatrixInfo(currFilePath)
                    if 'FETI' in matrixType:
                        matrixType = 'FETI local'
                    if 'an' in meshType:
                        meshType = meshType[3:]
                    else:
                        meshType = meshType[2:]
            
                    #speedup = 100.*(LUTime[currFilePath])/totalTime[currFilePath]
                    speedup = (LUTime[currFilePath])/totalTime[currFilePath]
                    sumFile.write('\t'+matrixType+'&'+meshType+'&'+meshDim+'&'+matrixLevel+'&'+sparseSize+'&'+str(size[currFilePath])+'&'+'%.2e' %directTime[currFilePath]+'&'+'%.2e' %directAccuracy[currFilePath]+'&'+'%.2e' %totalTime[currFilePath]+'&'+'%.2e' %finalAccuracy[currFilePath]+'&'+str(numIter[currFilePath])+'&'+'%.2e' %LUTime[currFilePath]+'&'+ '%.2f' %speedup+'\\\ \n')

sumFile.write('\t\\hline\n')
sumFile.write('\t\\end{tabular}}\n')
sumFile.write('\t\\end{table}\n')
sumFile.write('\\end{document}\n')
sumFile.close()
os.chdir(currDir)
proc = subprocess.Popen(shlex.split('pdflatex ' + 'summary.tex'))
proc.communicate()
os.unlink('summary.aux')
os.unlink('summary.log')

#uploadToDropBox(currDir)
