from matrices import *
from Bio import pairwise2, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from io import StringIO
from util.io import *
import sys
import itertools
from heapq import heappush, heappop

DEBUG = False

class Seq(object):
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence
    
    def __repr__(self):
        return self.name

class Tree(object):
    def __init__(self, left = None, right = None, data = None, height = None):
        self.left = left
        self.right = right
        self.data = data
        self.height = height

    def __lt__(self, other):
        return self.height * 2 < other.height * 2

    def __eq__(self, other):
        if len(self.data) == 1 or len(other.data) == 1:
            return self.data == other.data

        return self.left == other.left and self.right == self.right 
             
    
    def __repr__(self):
        #return "data: %s, weight: %s, height: %s" % (self.data, self.weight, self.height)
        return "%s" % (self.data)


def toNewick(tree):
    
    if(tree is None):
        return ''
    
    if(tree.left is None and tree.right is None):
        return str(tree.data)
        
    return "("+toNewick(tree.left)+")("+toNewick(tree.right)+")"
    
aux = {}

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

def buildtree(sequences):
    clusters = []

    for i in range(0, len(sequences)):
        clusters.append(Tree(data = sequences[i], height = 0))
    if(DEBUG):
        print(clusters)
    
    min_dij = sys.maxsize
    ci = Tree()
    cj = Tree()

    while (len(clusters) > 2):
        comb = list(itertools.combinations(clusters, 2))
        l,r = min(comb,key=lambda e: d(e[0],e[1]))
        clusters.remove(l)
        clusters.remove(r)
        newSeq = Seq(l.data.name+"&"+r.data.name, l.data.sequence+r.data.sequence)
        clusters.append(Tree(l, r, newSeq ,d(l,r)/2))
        if(DEBUG):
            print(clusters)

    
    newSeq = Seq(   clusters[0].data.name+"&"+clusters[1].data.name, #name
                    clusters[0].data.sequence+clusters[1].data.sequence #sequence
                )

    final = Tree(clusters[0],clusters[1], newSeq,d(clusters[0],clusters[1])/2)
    return final

def buildtreeheap(sequences):
    heap = []
    bitmap = {}

    clusters = []

    #Dar a cada sequencia o seu proprio cluster
    #definir uma folha para cada sequencia e colocar na altura 0

    for i in range(0, len(sequences)):
        clusters.append(Tree(data = sequences[i], height = 0))
    if(DEBUG):
        print(clusters)
    #enquanto ha mais que 2 clusters

    

    while (len(clusters) > 2):
        #fazer tuplos de todos os clusters possíveis com distancia associada -> (dij, Tree(di,dj))
        listadetuplos = [Tree(x,y,x.data+y.data,d(x,y),d(x,y)/2) for (x,y) in list(itertools.combinations(clusters, 2))]
        #incluir tuplos na heap
        #fazer, paralelamente, um dicionario com correspondencia de validade do tuplo; 'abcd' -> 'on'
        for tree in listadetuplos:
            heappush(heap,tree)
            t = tree
            bitmap[t.data] = "on"

        newcluster = heappop(heap)

        #fazer pop da heap ate sair um cujo data corresponda a 'on'

        while(bitmap[newcluster.data] != "on"):
            newcluster = heappop(heap)

        
        clusters.append(newcluster)
        clusters.remove(newcluster.left)
        clusters.remove(newcluster.right)

        #invalidar ('off') todos os elementos do dicionario cuja chave contenha um ou mais elementos do tuplo de cima
        appears_on = [key for key in bitmap.keys() if bitmap[key] == "on"]

        for l in newcluster.data:
            for key in appears_on:
                if l in key:
                    bitmap[key] = "off"

    final = Tree(clusters[0],clusters[1], clusters[0].data + clusters[1].data, d(clusters[0],clusters[1])/2)
    
    return final

if __name__ == '__main__':

    sequencias = [Seq(filename[0:-6],sequence) for (filename,sequence) in read_directory("inputs")]
    result = buildtree(sequencias)
    string = toNewick(result)
    tree = Phylo.read(StringIO(string+';'),"newick")
    Phylo.draw(tree)
    
