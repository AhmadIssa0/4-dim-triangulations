
from triangulation import *
import itertools

def calcFreeFaces(pents):
    free = []
    for i in range(len(pents)):
        p = pents[i]
        for j in p:
            tet = [m for m in p if m != j]
            is_free = True
            for k in range(len(pents)):
                if k != i and set(tet).issubset(pents[k]):
                    is_free = False
                    break
            if is_free:
                free.append(tet+[16])
    return free
    
def calcInternalFaces(pents):
    internal = []
    for i in range(len(pents)):
        p = pents[i]
        for j in p:
            tet = [m for m in p if m != j]
            for k in range(i+1,len(pents)):
                if set(tet).issubset(pents[k]):
                    internal.append([i, tet])
                    internal.append([k, tet])
                    break
    return internal
    
# simplicial complex of S^2 x D^2 with an ideal vertex.
# First we obtained a simplicial complex of S^2 x D^2 with free faces,
# obtained from simpcomp (GAP package). Then we coned off the boundary.
simplicial_complex = [[2,3,6,9,16], [1,3,6,9,16], [2,3,6,12,16], [1,3,6,12,16], [2,3,9,12,16], [1,3,9,12,16], [2,5,6,9,16], [2,5,6,12,16], [2,5,8,9,16], [1,2,5,8,16], [2,5,11,12,16], [1,2,5,11,16], [2,8,9,12,16], [2,8,11,12,16], [1,2,8,11,16], [1,4,6,9,16], [1,4,6,12,16], [1,4,5,8,16], [1,4,5,11,16], [1,4,7,9,16], [1,4,7,8,16], [1,4,10,12,16], [1,4,10,11,16], [1,7,9,12,16], [1,7,8,11,16], [1,7,10,12,16], [1,7,10,11,16], [5,6,9,15,16], [4,6,9,15,16], [5,6,12,15,16], [4,6,12,15,16], [5,8,9,15,16], [5,8,14,15,16], [4,5,8,14,16], [5,11,12,15,16], [5,11,14,15,16], [4,5,11,14,16], [4,7,9,15,16], [4,7,8,14,16], [4,7,13,15,16], [4,7,13,14,16], [4,10,12,15,16], [4,10,11,14,16], [4,10,13,15,16], [4,10,13,14,16], [8,9,12,15,16], [7,9,12,15,16], [8,11,12,15,16], [8,11,14,15,16], [7,8,11,14,16], [7,10,12,15,16], [7,10,11,14,16], [7,10,13,15,16], [7,10,13,14,16], [1,2,3,6,9], [1,2,3,6,12], [1,2,3,9,12], [1,2,5,6,9], [1,2,5,6,12], [1,2,5,8,9], [1,2,5,11,12], [1,2,8,9,12], [1,2,8,11,12], [1,4,5,6,9], [1,4,5,6,12], [1,4,5,8,9], [1,4,5,11,12], [1,4,7,8,9], [1,4,10,11,12], [1,7,8,9,12], [1,7,8,11,12], [1,7,10,11,12], [4,5,6,9,15], [4,5,6,12,15], [4,5,8,9,15], [4,5,8,14,15], [4,5,11,12,15], [4,5,11,14,15], [4,7,8,9,15], [4,7,8,14,15], [4,7,13,14,15], [4,10,11,12,15], [4,10,11,14,15], [4,10,13,14,15], [7,8,9,12,15], [7,8,11,12,15], [7,8,11,14,15], [7,10,11,12,15], [7,10,11,14,15], [7,10,13,14,15]]
            
gluings = calcInternalFaces(simplicial_complex)
tri = FourTriangulation(len(simplicial_complex))
for i in range(len(gluings)/2):
    n1, f1 = gluings[2*i]
    n2, f2 = gluings[2*i+1]
    f1 = [simplicial_complex[n1].index(j) for j in f1]
    f2 = [simplicial_complex[n2].index(j) for j in f2]
    g = [0]*5
    v1 = FourTriangulation.oppositeVertex(f1)
    g[v1] = FourTriangulation.oppositeVertex(f2)
    for j in range(4):
        g[f1[j]] = f2[j]
    p1 = tri.getPentachoronByNum(n1)
    p2 = tri.getPentachoronByNum(n2)
    p1.joinTo(v1, p2, g)
    
#tri.link().writeSnapPea('C:/pr2011/4dim_triangulation/s2xd2_link.tri')
def random33Move(tri):           
    # returns tri if there are no possible 3-3 moves
    lst = []
    cycles = []
    for p_no in range(tri.numPentachora()):
        for face in map(list, list(itertools.combinations(range(5),3))):
            if [p_no, face] not in lst:
                cycles.append(tri.cycle(p_no, face))
                lst += cycles[-1]
    cycles = filter(lambda x: len(x) == 3, cycles)
    valid_cycles = []
    for c in cycles:
        p_nos = [n for n, f in c]
        if len(set(p_nos)) == 3:
            valid_cycles.append(c)
    if len(valid_cycles) == 0:
        return tri
    else:
        import random
        i = random.randint(0, len(valid_cycles)-1)
        c = valid_cycles[i]
        p_nos = [n for n, f in c]
        faces = [f for n, f in c]
        return tri.pachner(p_nos, faces)
        
def fastSimplify(tri):
    # try any 5-1 moves first, then try 4-2 moves
    lst = []
    for p_no in range(tri.numPentachora()):
        for face in map(list, list(itertools.combinations(range(5),1))):
            if [p_no, face] not in lst:
                c = tri.cycle(p_no, face)
                if len(c) == 5: # found possible 5-1 move
                    p_nos = [n for n, f in c]
                    if len(set(p_nos)) == 5: # 5 distinct pentachora
                        faces = [f for n, f in c]
                        ntri = tri.pachner(p_nos, faces)
                        if ntri != tri:
                            ntri = fastSimplify(ntri)
                        return ntri
                lst += c
                
    lst = []
    for p_no in range(tri.numPentachora()):
        for face in map(list, list(itertools.combinations(range(5),2))):
            if [p_no, face] not in lst:
                c = tri.cycle(p_no, face)
                if len(c) == 4: # found possible 4-2 move
                    p_nos = [n for n, f in c]
                    if len(set(p_nos)) == 4: # 4 distinct pentachora
                        faces = [f for n, f in c]
                        ntri = tri.pachner(p_nos, faces)
                        if ntri != tri:
                            ntri = fastSimplify(ntri)
                        return ntri
                lst += c
    return tri
    
def random15Move(tri):
    lst = []
    for p_no in range(tri.numPentachora()):
        for face in map(list, list(itertools.combinations(range(5),5))):
            if [p_no, face] not in lst:
                c = tri.cycle(p_no, face)
                if len(c) == 1: # found possible 1-5 move
                    p_nos = [n for n, f in c]
                    if len(set(p_nos)) == 1: # 1 distinct pentachora
                        faces = [f for n, f in c]
                        return tri.pachner(p_nos, faces)
                lst += c
    return tri

def random24Move(tri):
    lst = []
    for p_no in range(tri.numPentachora()):
        for face in map(list, list(itertools.combinations(range(5),4))):
            if [p_no, face] not in lst:
                c = tri.cycle(p_no, face)
                if len(c) == 2: # found possible 2-4 move
                    p_nos = [n for n, f in c]
                    if len(set(p_nos)) == 2: # 2 distinct pentachora
                        faces = [f for n, f in c]
                        return tri.pachner(p_nos, faces)
                lst += c
    return tri
        
def foo():
    ntri = random33Move(tri)
    for i in range(500):
        ntri = random33Move(ntri)
        if i % max(5, 2*(i/100)) == 0:
            ntri = fastSimplify(ntri)
            print ntri.numPentachora()
    ntri.link().writeSnapPea('C:/pr2011/4dim_triangulation/s2xd2_pachner_link.tri')

#foo()

tri = FourTriangulation(1)
tri = tri.idealToFinite()
tri.boundaryTriangulation().writeSnapPea('C:/pr2011/4dim_triangulation/test.tri')