import os.path
import numpy as np
from collections import OrderedDict

#define directories' name

dictList=list(map(lambda x:"NC_0"+x+".2",list(map(str,list(range(19483,19484))))))

#define the base directory

dirName="/home/maulik/data/Shared/Maulik/HetDepthWindow"

#threshold of the standard deviations of mean to be considered to identify windows for beast analysis

stdTh=0.80

#output of the final directory (the final iteration)

rule all:
    input:
        expand("{dirName}/{dictL}/AllSpecies_{dictL}_windowStats.txt",dictL=dictList[-1],dirName=dirName),
        expand("{dirName}/{dictL}/AllSpecies_{dictL}_sorted.txt",dictL=dictList[-1],dirName=dirName)


#run awk command to make the input file to be used in the following python script
rule run_awk:
    output:
        expand("{dirName}/{dictL}/Chrm_Species.txt",dictL=dictList[-1],dirName=dirName)
    run:
        for chrom in dictList:
            newInput=os.path.join(dirName,chrom)
            newOutput=os.path.join(dirName,chrom,"Chrm_Species.txt")
            shell("ls {newInput} |awk '{{split($1,a,\"_\");print $1,a[1]}}' >{newOutput}")


def ExtractWindows(inFile,dirN,chrmN,stdT,outFile1,outFile2):
    dest=open(outFile1,"w")
    dest1=open(outFile2,"w")
    windowDict=OrderedDict()
    sampleHetDict=OrderedDict()
    sampleDepthDict=OrderedDict()
    stdHetDict=OrderedDict()
    stdDepthDict=OrderedDict()
    countWindowDict=OrderedDict()
    sampleList=[]
    with open(inFile) as source:
        for line in source:
            line=line.split()
            species=line[1]
            newInput1=os.path.join(dirN,chrmN,line[0])
            with open(newInput1) as source1:
                for lin in source1:
                    lin=lin.split()
                    count=0
                    window=lin[0]+"_"+lin[1]
                    if window not in windowDict:
                        windowDict[window]=OrderedDict()
                    sampleValuesList=[lin[i:i+3] for i in range(2,len(lin),3)]
                    for sample in sampleValuesList:
                        count+=1
                        sampleName=species+"_"+str(count)
                        if sampleName not in windowDict[window]:
                            windowDict[window][sampleName]=[]
                        if sampleName not in sampleHetDict:
                            sampleHetDict[sampleName]=[]
                            sampleDepthDict[sampleName]=[]
                        try:
                            windowDict[window][sampleName].append(round(int(sample[0])/int(sample[1]),5))
                            windowDict[window][sampleName].append(round(float(sample[2]),3))
                            sampleHetDict[sampleName].append(round(int(sample[0])/int(sample[1]),5))
                            sampleDepthDict[sampleName].append(round(float(sample[2]),3))
                        except:
                            pass
    for sample in sampleDepthDict:
        depthValues=sampleDepthDict[sample]
        hetValues=sampleHetDict[sample]
        stdDepthDict[sample]=[np.mean(depthValues)-np.std(depthValues)*float(stdT),np.std(depthValues)*float(stdT)+np.mean(depthValues)]
        stdHetDict[sample]=[np.mean(hetValues)-np.std(hetValues)*float(stdT),np.std(hetValues)*float(stdT)+np.mean(hetValues)]
    for window in windowDict:
        countWindowDict[window]=[0,0]
        for sample in windowDict[window]:
            sampleList.append(sample)
            try:
                if windowDict[window][sample][1]>=stdDepthDict[sample][0] and windowDict[window][sample][1]<=stdDepthDict[sample][1]:
                    countWindowDict[window][0]+=1
                if windowDict[window][sample][0]>=stdHetDict[sample][0] and windowDict[window][sample][0]<=stdHetDict[sample][1]:
                    countWindowDict[window][1]+=1
            except:pass
    sortedWindowDict=dict(sorted(countWindowDict.items(), key=lambda item: item[1]))
    windowArrangedList=list(sortedWindowDict.keys())[::-1]
    dest.write("start"+"\t"+"end"+"\t"+"\t".join(sampleList))
    for window in windowArrangedList:
        dest1.write(window+"\t"+str(countWindowDict[window][0])+"\t"+str(countWindowDict[window][1])+"\n")
        dest.write("\n")
        dest.write(window.split("_")[0]+"\t"+window.split("_")[1])
        for sample in windowDict[window]:
            try:
                dest.write("\t")
                dest.write(str(windowDict[window][sample][0])+"\t"+str(windowDict[window][sample][1]))
            except:pass
    dest.close()
    dest1.close()

##extract genomic windows with the threshold defined earlier
rule extract_windows:
    input:
        expand("{dirName}/{dictL}/Chrm_Species.txt",dictL=dictList[-1],dirName=dirName)
    output:
        expand("{dirName}/{dictL}/AllSpecies_{dictL}_windowStats.txt",dictL=dictList[-1],dirName=dirName),
        expand("{dirName}/{dictL}/AllSpecies_{dictL}_sorted.txt",dictL=dictList[-1],dirName=dirName)
    run:
        for chrom in dictList:
            newInput=os.path.join(dirName,chrom,"Chrm_Species.txt")
            newOutput1=os.path.join(dirName,chrom,"AllSpecies_"+chrom+"_windowStats.txt")
            newOutput2=os.path.join(dirName,chrom,"AllSpecies_"+chrom+"_sorted.txt")
            ExtractWindows(newInput,dirName,chrom,stdTh,newOutput1,newOutput2)
