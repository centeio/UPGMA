from matrices import *
from Bio import pairwise2, Phylo, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from io import StringIO
from util.io import *
import sys
import itertools
from heapq import heappush, heappop

DEBUG = True
DM = {}

class Seq(object):
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence
    
    def __repr__(self):
        return self.name

class Tree(object):
    def __init__(self, left = None, right = None, data = None, height = None, clusters = []):
        self.left = left
        self.right = right
        self.data = data
        self.height = height
        self.clusters = clusters

    def __lt__(self, other):
        return self.height * 2 < other.height * 2

    def __eq__(self, other):
        if self.height == 0 or other.height ==0:
            return self.data.name == other.data.name

        return self.left == other.left and self.right == self.right 
             
    
    def __repr__(self):
        #return "data: %s, weight: %s, height: %s" % (self.data, self.weight, self.height)
        return "%s" % (self.data)


def toNewick(tree):
    if(tree is None):
        return ''
    
    if(tree.left is None and tree.right is None):
        return str(tree.data)  

    currentMinusLeft = str(tree.height - tree.left.height)
    currentMinusRight = str(tree.height - tree.right.height)
    
    return "("+toNewick(tree.left)+":"+currentMinusLeft+")"+"("+toNewick(tree.right)+":"+currentMinusRight+")"
    
def d(ci, cj):
    calc = DistanceCalculator(model='blosum62')
    
    alignment = pairwise2.align.globalxx(ci.data.sequence,cj.data.sequence)[0]

    seq1 = alignment[0]
    seq2 = alignment[1]
    
    #An improvement for preventing re-calculating
    try:
        dist = aux[(seq1,seq2)]
    except KeyError as ex:
        aux[(seq1,seq2)] = calc._pairwise(seq1, seq2)
        dist = aux[(seq1,seq2)]

    return dist


def getDistanceMatrix(msaPath):
    aln = AlignIO.read(open(msaPath), 'phylip')
    calculator = DistanceCalculator('blosum62')
    return calculator.get_distance(aln)


def distance(t1,t2):
    ci = t1.clusters
    cj = t2.clusters
    comb = [DM["sp|{}|".format(p),"sp|{}|".format(q)] for p in ci for q in cj]

    length = len(comb)
    
    return sum(comb) / length


def buildtreeheap(sequences):
    heap = []
    bitmap = {}

    clusters = []

    #Dar a cada sequencia o seu proprio cluster
    #definir uma folha para cada sequencia e colocar na altura 0

    for i in range(0, len(sequences)):
        clusters.append(Tree(data = sequences[i], height = 0, clusters=[sequences[i]]))
    if(DEBUG):
        print(clusters)
    #enquanto ha mais que 2 clusters

    while (len(clusters) > 2):
        #fazer tuplos de todos os clusters possíveis com distancia associada -> (dij, Tree(di,dj))
        comb = list(itertools.combinations(clusters, 2))
        trees = [
                Tree(left=x, right=y,
                    data=Seq(x.data.name+"&"+y.data.name, #name
                        x.data.sequence+y.data.sequence #sequence
                        ),
                    height=distance(x,y)/2,
                    clusters= x.clusters + y.clusters                   
                    ) for (x,y) in comb
                ]

        
        #incluir tuplos na heap
        #fazer, paralelamente, um dicionario com correspondencia de validade do tuplo; 'abcd' -> 'on'
        for tree in trees:
            heappush(heap,tree)
            t = tree
            bitmap[t.data.name] = "on"

        newcluster = heappop(heap)

        while(bitmap[newcluster.data.name] != "on"):
            newcluster = heappop(heap)
        
        clusters.append(newcluster)
        clusters.remove(newcluster.left)
        clusters.remove(newcluster.right)

        
        #invalidar ('off') todos os elementos do dicionario cuja chave contenha um ou mais elementos do tuplo de cima
        appears_on = [key for key in bitmap.keys() if bitmap[key] == "on"]
        
        for key in appears_on:
            if newcluster.left.data.name in key or newcluster.right.data.name in key:
                bitmap[key] = "off"
    
    newSeq = Seq(   clusters[0].data.name+"&"+clusters[1].data.name, #name
                    clusters[0].data.sequence+clusters[1].data.sequence #sequence
                )

    final = Tree(clusters[0],clusters[1], newSeq,distance(clusters[0],clusters[1])/2,clusters = clusters[0].clusters + clusters[1].clusters)
    return final
    
if __name__ == '__main__':

    sequencias = [Seq(filename[0:-6],sequence) for (filename,sequence) in read_directory("inputs")]
    DM = getDistanceMatrix('./aligned.phy')
    result = buildtreeheap(sequencias)
    string = toNewick(result)
    tree = Phylo.read(StringIO(string+';'),"newick")
    Phylo.draw(tree)
    
