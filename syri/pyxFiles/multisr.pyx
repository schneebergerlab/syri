from igraph import *
#from scipy.stats import *
from datetime import datetime
import pandas as pd

getAllLongestPaths, getConnectivityGraph

cimport numpy as np
cimport cython

np.random.seed(1)

def getBlocksData(filePaths, fileTypes):
    genomeData = pd.DataFrame()
    synData = pd.read_table(filePaths[0],header = None)
    synData = synData.loc[synData[0] != "#"]
    synData = synData[[0,1,2,3]]
    synData.columns = ["aStart", "aEnd","bStart","bEnd"]
    synData = synData.astype("int64")
    synData["state"] = fileTypes.pop(0)
    genomeData = genomeData.append(synData)

    for i in filePaths[1:]:
        fileData = pd.read_table(i, header = None)
        fileData = fileData.loc[fileData[0] == "#"]
        fileData = fileData[[1,2,4,5]]
        fileData.columns = ["aStart", "aEnd","bStart","bEnd"]
        fileData = fileData.astype("int64")
        fileData["state"] = fileTypes.pop(0)
        genomeData = genomeData.append(fileData)
    genomeData.sort_values(by = ["aStart", "aEnd","bStart","bEnd"], inplace = True)
    genomeData.index = range(len(genomeData))
    return(genomeData)


def getConservedRegions(dataList, isSyn = True):
    if not isinstance(isSyn, bool):
        raise TypeError("isSyn must be a bool")

    bGenomes = np.unique(dataList.bGenome)
    genomeCount = len(bGenomes)
    if isSyn:
        allCoords = dataList[["start","end","bGenome"]]
    else:
        allCoords = pd.DataFrame()
        for i in bGenomes:
            gData = dataList.loc[dataList.bGenome == i]
            start = list(gData.start)
            end = list(gData.end)
            s = start[0]
            e = end[0]
            endStack = [float('inf'),e]
            region = []
            for j in range(1,len(start)):
                s1 = start[j]
                e1 = end[j]
                if s1 < e:
                    region.append([s,s1])
                    s = s1
                    endStack.append(e1)
                    e = min(endStack)
                elif e <= s1:
                    region.append([s,e])
                    while True:
                        s = e
                        endStack.remove(e)
                        e = min(endStack)
                        if e <= s1:
                            region.append([s,e])
                        else:
                            break
                    if len(endStack) > 1:
                        region.append([s,s1])
                        s = s1
                        endStack.append(e1)
                        e = min(endStack)
                    else:
                        s = s1
                        endStack.append(e1)
                        e = min(endStack)
            while len(endStack) > 1:
                region.append([s,e])
                s = e
                endStack.remove(e)
                e = min(endStack)
            region = [a for a in region if a[0] < a[1]]
            region = pd.DataFrame(region)
            region.columns = ["start","end"]
            region["bGenome"] = i
            allCoords = allCoords.append(region)
        allCoords.sort_values(["start","end"],inplace = True)

    terminalData = pd.DataFrame(data = np.zeros([2,genomeCount], dtype = "int"), index = ["start","end"], columns = bGenomes)

    inGenomeCount = 0
    regions = []
    count = 0
    for row in allCoords.itertuples(index=False):
        count+=1
        if row.end <= terminalData[row.bGenome].end:
            print("values must be sorted. Invalid Entry: ",row, terminalData[row.bGenome])
            sys.exit()
        if row.start <= terminalData[row.bGenome].end:
            terminalData[row.bGenome].start = terminalData[row.bGenome].end + 1
            terminalData[row.bGenome].end = row.end
        else:
            terminalData[row.bGenome].start = row.start
            terminalData[row.bGenome].end = row.end
        if max(terminalData.loc["start"]) < min(terminalData.loc["end"]):
            regions.append((max(terminalData.loc["start"]),min(terminalData.loc["end"])))
    regions = pd.DataFrame(regions, columns = ["start","end"])
    return regions


def getDataList(dataTables, genomeID, identity):
    start = []
    end = []
    bGenome = []
    bStart = []
    bEnd = []
    state = []

    if len(dataTables) != len(genomeID):
        print("need 1 identifier for each table")
        sys.exit()
    else:
        for i in range(len(genomeID)):
            if identity[i] == "a":
                start.extend(dataTables[i].aStart.tolist())
                end.extend(dataTables[i].aEnd.tolist())
                bStart.extend(dataTables[i].bStart.tolist())
                bEnd.extend(dataTables[i].bEnd.tolist())
            elif identity[i] == "b":
                start.extend(dataTables[i].bStart.tolist())
                end.extend(dataTables[i].bEnd.tolist())
                bStart.extend(dataTables[i].aStart.tolist())
                bEnd.extend(dataTables[i].aEnd.tolist())
            state.extend(dataTables[i].state.tolist())
            bGenome.extend([genomeID[i]]*len(dataTables[i]))
        outData = pd.DataFrame({"start":start,"end":end,"bStart":bStart,"bEnd":bEnd,"state":state,"bGenome":bGenome})
        outData.sort_values(["start","end","bStart","bEnd"], inplace = True)
        outData = outData[["start","end","bStart","bEnd","state","bGenome"]]
        outData.index = range(len(outData))
        return(outData)


def getCLQ(adjM,partSize):
    """
    Mirghorbani, M., & Krokhmal, P. (2013). On finding k-cliques in k-partite graphs. Optim Lett, 7, 1155–1165. https://doi.org/10.1007/s11590-012-0536-y
    """
    class startCLQ:
        def __init__(self, adjM, partSize):
            self.clq = []
            self.t = 0
            self.sub = []
            self.BsOut = list(range(len(partSize)))
            self.Bs = []
            self.Z = np.ones([len(partSize),sum(partSize)], dtype = "bool")
            self.Z0 = np.ones([len(partSize), sum(partSize)])
            self.adjM = adjM
            self.partSize = partSize
            self.partIndex = [sum(partSize[:i]) for i in range(len(partSize))]
            self.clqSize = len(partSize)
            self.S = []
            self.Q = []
            self.bitCLQ(self.t)

        def getBits(self, t, b):
            return self.Z[t, self.partIndex[b]:self.partIndex[b]+self.partSize[b]]

        def getBt(self, partSize, t):
            nodeCount = [sum(self.getBits(t,i)) for i in self.BsOut]
            return self.BsOut[nodeCount.index(min(nodeCount))]

        def bitCLQ(self, t):
            bt = self.getBt(self.partSize, t)
            sigBits = np.where(np.array(self.getBits(t,bt)) == True)[0]
            sigBitsLen = len(sigBits)
            count = 0
            for i in sigBits:
                count+=1
                if t == 0:
                    print(count, sigBitsLen, datetime.datetime.now())
                nt = self.partIndex[bt]+i
                self.Z[t,nt] = 0
                self.S.append(nt)
                if len(self.S) == self.clqSize:
                    self.Q.append(self.S.copy())
                    self.S.remove(nt)
                else:
                    self.Z[t+1] = self.Z[t] & self.adjM[nt]
                    self.Bs.append(bt)
                    self.BsOut.remove(bt)
                    P = sum([1 for i in self.BsOut if sum(self.getBits(t,i)) > 0])
                    if len(self.S) + P == self.clqSize:
                        self.bitCLQ(t+1)
                        self.S.remove(nt)
                        self.Bs.remove(bt)
                        self.BsOut.append(bt)
                    else:
                        self.S.remove(nt)
                        self.Bs.remove(bt)
                        self.BsOut.append(bt)

    def filterCLQ(clqList, partIndex):
        clqList = [sorted(i) for i in clqList]
        clqList = [[i[j] - partIndex[j] for j in range(len(i))] for i in clqList]
        return(clqList)
    CLQData = startCLQ(adjM, partSize)
    return(filterCLQ(CLQData.Q, CLQData.partIndex))
