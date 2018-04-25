from matrices import *
from Bio import pairwise2
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from io import StringIO
from util.io import *
import sys
import itertools

DEBUG = True

class Tree(object):
    def __init__(self, left = None, right = None, data = None, height = None):
        self.left = left
        self.right = right
        self.data = data
        self.height = height
    
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
    
    alignment = pairwise2.align.globalxx(ci.data,cj.data)[0]

    seq1 = alignment[0]
    seq2 = alignment[1]
    
    #An improvement for preventing re-calculating
    try:
        dist = aux[(seq1,seq2)]
    except KeyError as ex:
        aux[(seq1,seq2)] = calc._pairwise(seq1, seq2)
        dist = aux[(seq1,seq2)]

    return dist

def buildtree(seq1):
    clusters = []

    #Dar a cada sequencia o seu proprio cluster
    #definir uma folha para cada sequencia e colocar na altura 0

    for i in range(0, len(seq1)):
        clusters.append(Tree(data = seq1[i], height = 0))
    if(DEBUG):
        print(clusters)
    #enquanto ha mais que 2 clusters
    min_dij = sys.maxsize
    ci = Tree()
    cj = Tree()

    while (len(clusters) >= 2):
        comb = list(itertools.combinations(clusters, 2))
        l,r = min(comb,key=lambda e: d(e[0],e[1]))
        clusters.remove(l)
        clusters.remove(r)
        clusters.append(Tree(l, r, l.data + r.data,d(l,r)/2))
        if(DEBUG):
            print(clusters)
    
    return clusters[0]

if __name__ == '__main__':
    #seq1 = read_fasta(read_file('inputs/Q61743.fasta'))
    seq1 = "MLSRKGIIPEEYVLTRLAEDPAEPRY"
    result = buildtree(seq1)
    string = toNewick(result)
    tree = Phylo.read(StringIO(string+';'),"newick")
    Phylo.draw(tree)
    
