from matrices import *
import sys
import itertools
from heapq import heappush, heappop

matrix_dict =   { 'BLOSUM62': BLOSUM62, 'DNAFULL': DNAfull, 'PAM250': PAM250}
matrix = matrix_dict['BLOSUM62']

class Tree(object):
    def __init__(self, left = None, right = None, data = None, weight = None, height = None):
        self.left = left
        self.right = right
        self.data = data
        self.weight = weight
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

def d(ci, cj):
#    if len(ci.data) == 1:
#        if len(cj.data) == 1: #both have one letter
#            return matrix[ci.data][cj.data]
#        return (d(ci, cj.right) + d(ci, cj.left)) / 2 #only ci has one letter
#    elif len(cj.data) == 1
#        return (d(ci.right, cj) + d(ci.right, cj)) / 2 # only cj has one letter

#    return d(ci.right, cj.) + d(ci.right, cj)) / 2 
    if len(ci.data) == len(cj.data) == 1: #both have one letter
        return matrix[ci.data][cj.data]
    return (ci.weight + cj.weight) / 2


def buildtree(seq1):
    clusters = []

    #Dar a cada sequencia o seu proprio cluster
    #definir uma folha para cada sequencia e colocar na altura 0

    for i in range(0, len(seq1)):
        clusters.append(Tree(data = seq1[i], height = 0, weight = 0))

    #enquanto ha mais que 2 clusters

    print(clusters)

    while (len(clusters) > 2):
        comb = list(itertools.combinations(clusters, 2))
        l,r = min(comb,key=lambda e: d(e[0],e[1]))
        clusters.remove(l)
        clusters.remove(r)
        clusters.append(Tree(l, r, l.data + r.data, 1,d(l,r)/2))
        print(clusters)

    #determinar dois clusters i e j, com a menor dij
    #definir novo cluster Ck com Clusters i e j e criar nos na altura dij/2
    #substituir clusters i e j por Ck

    #dado um novo cluster ck podemos calcular a sua distancia a todos os outros

    final = Tree(clusters[0],clusters[1], clusters[0].data + clusters[1].data, d(clusters[0],clusters[1])/2)
    print(final)

def buildtreeheap(seq1):
    heap = []
    bitmap = {}

    clusters = []

    #Dar a cada sequencia o seu proprio cluster
    #definir uma folha para cada sequencia e colocar na altura 0

    for i in range(0, len(seq1)):
        clusters.append(Tree(data = seq1[i], height = 0, weight = 0))

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

        print(newcluster)
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
    print(final)




    

if __name__ == '__main__':
    seq1 = "ACGCG"

    buildtreeheap(seq1)



