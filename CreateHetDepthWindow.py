import sys
from collections import OrderedDict
from pysam import VariantFile

def IdentifyMissingSites(bcfIn,sampleDepthFile,chrmName,windowSize,stepSize,outPrefix):
    dest=open(outPrefix+".txt","w")
    bcf_in=VariantFile(bcfIn)
    sampleList=list((bcf_in.header.samples))
    genoDict={(0,1):"Het",(None,None):"Mis",(0,0):"Hom",(1,1):"Hom"}
    posDepthDict=OrderedDict()
    sampleDepth=OrderedDict()
    chrmLengthDict=OrderedDict()
    indelList=[-9]
    posGenoDict=OrderedDict()
    snpPos=[]
    with open(sampleDepthFile) as infile:
        for line in infile:
            line=line.split()
            sampleDepth[line[0]]=float(line[1])
    for x in bcf_in.header.records:
        if x.key=="contig":
            chrmLengthDict[x.values()[0]]=int(x.values()[1])
    chrmLength=chrmLengthDict[chrmName]+1
    chrmWindowList=[]
    for i in range(1,chrmLength,int(stepSize)):
        if i+int(windowSize)>=chrmLength:
            chrmWindowList.append([i,chrmLength])
        else:
            chrmWindowList.append([i,i+int(windowSize)])
    for rec in bcf_in.fetch():
        posGenoDict[rec.pos]=[]
        posDepthDict[rec.pos]=[]
        if rec.info["MQ"]>=40:
            if len(rec.ref)>1:
                for i,k in enumerate(sampleList):
                    posDepthDict[rec.pos].append(-9)
                    posGenoDict[rec.pos].append("Mis")
            else:
                if rec.alts==None:
                    for i,k in enumerate(sampleList):
                        posDepthDict[rec.pos].append(rec.samples[k]["DP"])
                        posGenoDict[rec.pos].append("Hom")
                else:
                    if len(rec.alts)>1 or len(rec.alts[0])>1:
                        for i,k in enumerate(sampleList):
                            posDepthDict[rec.pos].append(-9)
                            posGenoDict[rec.pos].append("Hom")
                    elif len(rec.alts)==1 and len(rec.alts[0])==1:
                        for sample in sampleList:
                            posDepthDict[rec.pos].append(rec.samples[sample]["DP"])
                            if (rec.samples[sample]["DP"]>float(sampleDepth[sample])*3) or \
                            (rec.samples[sample]["DP"]<float(sampleDepth[sample]/3)):
                                posGenoDict[rec.pos].append("Hom")
                            else:
                                posGenoDict[rec.pos].append(genoDict[rec.samples[sample]["GT"]])
        else:
            for i,k in enumerate(sampleList):
                posDepthDict[rec.pos].append(-9)
                posGenoDict[rec.pos].append("Mis")
    for window in chrmWindowList:
        tmpList=[]
        tmpDepthList=[]
        for sample in sampleList:
            tmpList.append([])
            tmpDepthList.append([])
        for i in range(window[0],window[1]):
            if i in posGenoDict:
                genoList=posGenoDict[i]
                depthList=posDepthDict[i]
                for i,v in enumerate(genoList):
                    tmpList[i].append(v)
                    if depthList[i]!=-9:
                        tmpDepthList[i].append(depthList[i])
        dest.write(str(window[0])+"\t"+str(window[1]))
        for i,v in enumerate(tmpList):
            sampleDepthList=tmpDepthList[i]
            dest.write("\t")
            try:
                dest.write(str(v.count("Het"))+"\t"+str(v.count("Hom"))+"\t"+str(round(sum(sampleDepthList)/len(sampleDepthList),3)))
            except:
                dest.write("0"+"\t"+"0"+"\t"+"0")
        dest.write("\n")
    dest.close()
if __name__=="__main__":
    IdentifyMissingSites(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6])
