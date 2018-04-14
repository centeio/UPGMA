from matrices import *
import sys
import itertools

matrix_dict =   { 'BLOSUM62': BLOSUM62, 'DNAFULL': DNAfull, 'PAM250': PAM250}
matrix = matrix_dict['BLOSUM62']

class Tree(object):
    def __init__(self, left = None, right = None, data = None, weight = None, height = None):
        self.left = left
        self.right = right
        self.data = data
        self.weight = weight
        self.height = height
    
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
    min_dij = sys.maxsize
    ci = Tree()
    cj = Tree()
    print(clusters)

    while (len(clusters) >= 2):
        comb = list(itertools.combinations(clusters, 2))
        l,r = min(comb,key=lambda e: d(e[0],e[1]))
        clusters.remove(l)
        clusters.remove(r)
        clusters.append(Tree(l, r, l.data + r.data, 1,d(l,r)/2))
        print(clusters)

       # print("current: %s, min: %s" % (current.data, x.data))
#      clusters.append(Tree(current, m, current.data + x.data, d(m,current)))
 #       for i in range(0, len(clusters)):
  #          for j in range(i+1, len(clusters)):
   #             if (d(clusters[i], clusters[j]) < min_dij):
    #                ci = clusters[i]
     #               cj = clusters[j]
      #              min_dij = d(ci, cj) 



    #determinar dois clusters i e j, com a menor dij
    #definir novo cluster Ck com Clusters i e j e criar nos na altura dij/2
    #substituir clusters i e j por Ck

    #dado um novo cluster ck podemos calcular a sua distancia a todos os outros

def print_tree:

    

if __name__ == '__main__':
    seq1 = "AAPGE"

    buildtree(seq1)



