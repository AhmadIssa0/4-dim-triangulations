
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
# obtained from simpcomp (GAP package).
simplicial_complex = [[1,2,3,6,9], [1,2,3,6,12], [1,2,3,9,12], [1,2,5,6,9], [1,2,5,6,12], [1,2,5,8,9], [1,2,5,11,12], [1,2,8,9,12], [1,2,8,11,12], [1,4,5,6,9], [1,4,5,6,12], [1,4,5,8,9], [1,4,5,11,12], [1,4,7,8,9], [1,4,10,11,12], [1,7,8,9,12], [1,7,8,11,12], [1,7,10,11,12], [4,5,6,9,15], [4,5,6,12,15], [4,5,8,9,15], [4,5,8,14,15], [4,5,11,12,15], [4,5,11,14,15], [4,7,8,9,15], [4,7,8,14,15], [4,7,13,14,15], [4,10,11,12,15], [4,10,11,14,15], [4,10,13,14,15], [7,8,9,12,15], [7,8,11,12,15], [7,8,11,14,15], [7,10,11,12,15], [7,10,11,14,15], [7,10,13,14,15]]

def triangFromSC(simplicial_complex):
    # return a FourTriangulation of a simplicial complex
    # given in "simpcomp" format
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
    return tri
    
#tri = triangFromSC(simplicial_complex)
#tri.writeToFile('C:/pr2011/4dim_triangulation/s2xd2.tri')

def random23(tri):
    lst = []
    bdryTet = tri.boundaryTetrahedra()
    from random import shuffle
    shuffle(bdryTet)
    for t_no in bdryTet:
        p_no, opp_vert = t_no
        if t_no in lst:
            continue
        verts = [x for x in range(5) if x != opp_vert]
        for face in map(list, list(itertools.combinations(verts,3))):
            if [t_no, face] not in lst:
                c = tri.boundaryCycle(t_no, face)
                if len(c) == 2:
                    t_nos = [n for n, f in c]
                    if len(set(map(tuple, t_nos))) == 2:
                        faces = [f for n, f in c]
                        ntri = tri.boundaryPachner(t_nos, faces)
                        if ntri != tri:
                            return ntri
                lst += c
    return tri
    
def randomBoundaryMove(tri, a, b):
    lst = []
    bdryTet = tri.boundaryTetrahedra()
    from random import shuffle
    shuffle(bdryTet)
    for t_no in bdryTet:
        p_no, opp_vert = t_no
        if t_no in lst:
            continue
        verts = [x for x in range(5) if x != opp_vert]
        for face in map(list, list(itertools.combinations(verts,b))):
            if [t_no, face] not in lst:
                c = tri.boundaryCycle(t_no, face)
                if len(c) == a:
                    t_nos = [n for n, f in c]
                    if len(set(map(tuple, t_nos))) == a:
                        faces = [f for n, f in c]
                        ntri = tri.boundaryPachner(t_nos, faces)
                        if ntri != tri:
                            return ntri
                lst += c
    return tri
    
def validMoves(tri,a,b):
    # a-b pachner move
    lst = []
    moves = []
    for t_no in tri.boundaryTetrahedra():
        p_no, opp_vert = t_no
        if t_no in lst:
            continue
        verts = [x for x in range(5) if x != opp_vert]
        for face in map(list, list(itertools.combinations(verts,b))):
            if [t_no, face] not in lst:
                c = tri.boundaryCycle(t_no, face)
                if len(c) == a:
                    t_nos = [n for n, f in c]
                    if len(set(map(tuple, t_nos))) == a:
                        faces = [f for n, f in c]
                        ntri = tri.boundaryPachner(t_nos, faces)
                        if ntri != tri:
                            moves.append([t_nos, faces])
                lst += c
    return moves
    
def random32(tri):
    lst = []
    for t_no in tri.boundaryTetrahedra():
        p_no, opp_vert = t_no
        if t_no in lst:
            continue
        verts = [x for x in range(5) if x != opp_vert]
        for face in map(list, list(itertools.combinations(verts,2))):
            if [t_no, face] not in lst:
                c = tri.boundaryCycle(t_no, face)
                if len(c) == 3:
                    t_nos = [n for n, f in c]
                    if len(set(map(tuple, t_nos))) == 3:
                        faces = [f for n, f in c]
                        #print '3 distinct', t_nos
                        ntri = tri.boundaryPachner(t_nos, faces)
                        if ntri != tri:
                            return ntri
                lst += c
    return tri
    
def random41(tri, debug=False):
    lst = []
    for t_no in tri.boundaryTetrahedra():
        p_no, opp_vert = t_no
        if t_no in lst:
            continue
        verts = [x for x in range(5) if x != opp_vert]
        for face in map(list, list(itertools.combinations(verts,1))):
            if [t_no, face] not in lst:
                c = tri.boundaryCycle(t_no, face)
                if len(c) == 4:
                    t_nos = [n for n, f in c]
                    if len(set(map(tuple, t_nos))) == 4:
                        faces = [f for n, f in c]
                        ntri = tri.boundaryPachner(t_nos, faces, debug)
                        if ntri != tri:
                            return ntri
                lst += c
    return tri
    
def random14(tri, debug=False):
    lst = []
    for t_no in tri.boundaryTetrahedra():
        p_no, opp_vert = t_no
        if t_no in lst:
            continue
        verts = [x for x in range(5) if x != opp_vert]
        for face in map(list, list(itertools.combinations(verts,4))):
            if [t_no, face] not in lst:
                c = tri.boundaryCycle(t_no, face)
                if len(c) == 1:
                    t_nos = [n for n, f in c]
                    if len(set(map(tuple, t_nos))) == 1:
                        faces = [f for n, f in c]
                        ntri = tri.boundaryPachner(t_nos, faces, debug)
                        if ntri != tri:
                            return ntri
                lst += c
    return tri
    
def fastSimplify(tri):
    ntri = random41(tri)
    ntri = random32(ntri)
    if ntri != tri:
        return fastSimplify(ntri)
    else:
        return tri
        
def reduceWithDepth(tri, depth=2, size=None):
    if depth == 0:
        return fastSimplify(tri)
    if size is None:
        size = tri.numBoundaryTet()
    moves = validMoves(tri,2,3)

    for move in moves:
        ntri = tri.boundaryPachner(move[0], move[1])
        stri = reduceWithDepth(ntri, depth-1, size)
        if stri.numBoundaryTet() < size:
            return stri
    return tri
        
def deepReduce(tri):
    moves1 = validMoves(tri,2,3)
    size = len(tri.boundaryTetrahedra())
    for move in moves1:
        ntri = tri.boundaryPachner(move[0], move[1])
        moves2 = validMoves(ntri,2,3)
        for move2 in moves2:
            nntri = fastSimplify(ntri.boundaryPachner(move2[0], move2[1]))
            if len(nntri.boundaryTetrahedra()) < size:
                return nntri
    return tri
    
def reduce(tri):
    moves = validMoves(tri,2,3)
    size = len(tri.boundaryTetrahedra())
    best_tri = tri
    best_size = size
    for move in moves:
        ntri = fastSimplify(tri.boundaryPachner(move[0], move[1]))
        if len(ntri.boundaryTetrahedra()) < best_size:
            best_size = len(ntri.boundaryTetrahedra())
            best_tri = ntri
            return ntri
    return best_tri
    
def reduceVerts(tri):
    moves = validMoves(tri,2,3)
    size = len(tri.boundaryTetrahedra())
    best_tri = tri
    best_size = size
    for move in moves:
        ntri = random41(tri.boundaryPachner(move[0], move[1]), False)
        if len(ntri.boundaryTetrahedra()) < best_size:
            best_size = len(ntri.boundaryTetrahedra())
            best_tri = ntri
            return ntri
    return best_tri
        
def boundaryRandomize(tri):
    for j in range(5):
        tri = randomBoundaryMove(tri,1,4)
        for i in range(5):
            tri = randomBoundaryMove(tri,2,3)
        for i in range(20):
            tri = randomBoundaryMove(tri,3,2)
        tri = tri.fastSimplify()
    for i in range(20):
        tri = randomBoundaryMove(tri,2,3)
        tri = randomBoundaryMove(tri,3,2)
    return tri
    
def glueAlongBoundaries(tri, tri1, output_fname, debug=False):
    if debug:
        print tri1.numBoundaryTet()
        print tri.numBoundaryTet()
        tri.boundaryTriangulation().writeSnapPea('C:/pr2011/4dim_triangulation/s2xd2_boundary.tri')
        tri.boundaryTriangulation().writeSnapPea('C:/pr2011/4dim_triangulation/s2xd2_boundary2.tri')
    
    all_isos = tri.boundaryIsomorphisms(tri1, debug=debug)
    if debug:
        print 'number of isos', len(all_isos)
        print 'all_isos', all_isos
    
    if len(all_isos) == 0:
        return False # failed
    else:
        iso = all_isos[0]
        c_tri = tri.copy()
        #n_tri = c_tri.glueAlongBoundary(tri, [[[138,2],[138,2],[0,1,2,3,4]], 
        #                            [[139,0],[139,0],[0,1,2,3,4]]])
        n_tri = c_tri.glueAlongBoundary(tri1, iso)
        n_tri = n_tri.fastSimplify()
        n_tri.writeToFile(output_fname)
    
    
# glue two S2 x D2's together, see if you get S2 x S2 or S2 (twisted) S2

# starting_tri = filename of S2 x D2 we begin with
# folder_path = all filenames are searched from folder path
def attemptGluckRetriang(
        starting_tri='s2xd2_gluck_piece1.tri',
        final_tri='s2xd2_rand1.tri',
        bundle_tri='s2xs2.tri',
        folder_path='C:/Users/Ahmad/Dropbox/4dim_triangulation'):
    
    tri = FourTriangulation.loadFromFile(folder_path + '/' + starting_tri)
    tri = boundaryRandomize(tri)
    print 'num of pent', tri.numPentachora()
    print 'boundary tet:', tri.numBoundaryTet()
    
    
    # this definitely reduces so boundary is 2 tetrahedra
    print 'num of tet', tri.numPentachora()
    
    tri = boundaryRandomize(tri)
    print 'num of pent', tri.numPentachora()
    print 'boundary tet:', tri.numBoundaryTet()
    
    tri = tri.fastSimplify()
    for i in range(3):
        tri = reduce(tri)
    for i in range(9):
        print 'deepReduce', i, tri.numBoundaryTet()
        tri = deepReduce(tri)
    for i in range(2):
        print 'depth 3:', i
        tri = reduceWithDepth(tri, 3)
    for i in range(20):
        for j in range(10):
            tri = random23(tri)
        for j in range(10):
            tri = reduce(tri)
        if tri.numBoundaryTet() == 2:
            break
    tri = tri.fastSimplify()
    print 'num of pent', tri.numPentachora()
    print 'boundary tet:', tri.numBoundaryTet()
    
    #tri.writeToFile('C:/pr2011/4dim_triangulation/s2xd2_rand1.tri')
    
    if tri.numBoundaryTet() != 2:
        print 'boundary tet:', tri.numBoundaryTet()
        exit(1)
    else:
        tri.writeToFile(folder_path + '/' + final_tri)
        print 'succeeded'
        tri1 = FourTriangulation.loadFromFile(folder_path + '/' + starting_tri)
        glueAlongBoundaries(tri1, tri, folder_path + '/' + bundle_tri, debug=False)
        tri = FourTriangulation.loadFromFile(folder_path + '/' + bundle_tri)
        print tri.intersectionForm()
    #tri = deepReduce(tri)
    print 'END'

#attemptGluckRetriang()


############ this part glues s2xd2 onto budney triangulation ########
"""
bTet = tri.boundaryTetrahedra()
print bTet
for t_no in bTet:
    p_no, v = t_no
    for j in [x for x in range(5) if x != v]:
        print t_no, 'vert', j, 'glues to', tri.boundaryGluing(t_no, j)
        
ctri = FourTriangulation(2)
P = ctri.pentachora
P[0].joinTo(3, P[0], [1,2,3,0,4])
P[1].joinTo(4, P[1], [0,2,3,4,1])
P[1].joinTo(2, P[0], [0,1,4,2,3])
P[1].joinTo(3, P[0], [0,2,3,1,4])
P[1].joinTo(0, P[0], [2,0,1,3,4])

ctri = ctri.idealToFinite()
for i in range(3):
    ctri = ctri.boundaryReduce()
    
ntri = ctri.glueAlongBoundary(tri, [ [[300,0],[138,2],[2,1,0,3,4]], [[301,0],[139,0],[0,2,4,3,1]] ])
"""
########################################################################

#import cProfile
#cProfile.run('ntri = ntri.fastSimplify()')
#ntri = ntri.simplify()
#print 'size:', ntri.numPentachora()

#ntri.writeToFile('C:/pr2011/4dim_triangulation/homotopy_4sphere.tri')
