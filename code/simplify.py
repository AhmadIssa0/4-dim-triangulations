
from triangulation import *
import sys
import random

#import cProfile
#cProfile.run('tri = tri.fastSimplify()')
#cProfile.run('tri = tri.simplify()')

#print 'valid 3-3 moves:', len(tri.validMoves(3,3))

#tri = tri.simplify()

def randomizeFast(tri):
    # n is how many times to try to randomize
    lst = []
    valid_moves = []
    c_tri = tri.copy()

    from random import shuffle
    pent_nos = tri.getPentachoronList()
    shuffle(pent_nos)

    for p in pent_nos:
        p_no = p.getPentNo()
        for face in map(list, list(itertools.combinations(range(5),3))):
            if [p_no, face] not in lst:
                c = tri.cycle(p_no, face)
                if len(c) == 3: # found possible 3-3 move
                    p_nos = [n for n, f in c]
                    if len(set(p_nos)) == 3: # 3 distinct pentachora
                        faces = [f for n, f in c]
                        ntri = c_tri.pachner(p_nos, faces)
                        if ntri is not False and ntri != tri:
                            return ntri
                            #c_tri = tri.copy()
                            #valid_moves.append([p_nos, faces])
                lst += c
    return tri

def randomize(tri):
    moves = tri.validMoves(3,3)
    if len(moves) == 0:
        return tri
    c = random.randint(0,len(moves)-1)
    moves = moves[c:] + moves[:c]
    used = []
    for m in moves:
        if m[0][0] not in used and m[0][1] not in used and m[0][2] not in used:
            if random.randint(0,1) == 0:
                tri.pachner(m[0], m[1])
                used.extend(m[0])

def smallReduce(tri, allow51=False):
    for j in range(10):
        for i in range(5):
            #randomize(tri)
            tri = randomizeFast(tri)
        tri = tri.fastSimplify(allow51=allow51)
        tri.renumberTriangulation()
        #if j % 5 == 0:
        #    print 'size:', tri.numPentachora()
        #    sys.stdout.flush()
    
    tri = tri.simplify(allow51)
    #print 'size:', tri.numPentachora()
    
    #tri = tri.simplifyWithDepth(2)
    #print 'simplified size:', tri.numPentachora()
    #if tri.numPentachora() <= 36:
    #    tri = tri.simplifyWithDepth(4)
    #    print 'depth 3 size:', tri.numPentachora()
    #print 'valid 2-4 moves:', len(tri.validMoves(2,4))
    #fVector(tri)
    return tri
    
def vertInfo(tri):
    v_cycs = tri.vertexCycles()
    print 'vertices:', len(v_cycs), 'pent:', tri.numPentachora()
    print sorted([len(c) for c in v_cycs])

def fVector(tri):
    print 'fVector:', tri.fVector(), 'complexity:', tri.complexity()

def heighten(tri, m=6, allow51=True, allow15=False):
    n = tri.numPentachora()
    for i in range(m):
        if i == 0 and allow15:
            tri = tri.anyMove(1,5)
        else:
            tri = tri.anyMove(2,4)
        tri = randomizeFast(tri)
        tri = randomizeFast(tri)
        tri = randomizeFast(tri)
        tri = tri.fastSimplify(allow51)
        #tri = tri.simplify()
        #print 'size:', tri.numPentachora()
    
    #v_cycs = tri.vertexCycles()
    #print 'vertices:', len(v_cycs), 'pent:', tri.numPentachora()
    #print sorted([len(c) for c in v_cycs])
    return tri

def simplifyUntilS4(tri, time_in_millis=60000):
    """
    Attempts to simplify the triangulation for a maximum time period of `time_in_millis`
    milliseconds, or until it's the standard two 4-simplex triangulation of S^4.
    
    To achieve this we try to make the triangulation have exactly 5 vertices.

    Returns a pair of a new triangulation (leaves current triangulation unmodified),
    and a boolean of whether returned triangulation is the standard S^4 triangulation.
    """
    tri = tri.copy()
    size = tri.numPentachora()
    print 'size:', size
    import time
    init_time = int(round(time.time() * 1000))
    curr_time = init_time
    loop = 0
    min_tri = tri.copy()
    global_min = tri.copy()
    last_loop = 0 # last local reduction
    heat = 1
    allow51 = False
    firstTimeReached5Verts = False
    reached_min = 0
    def spike(triang, iterations=3, heat=3, allow51=False):
        for i in range(iterations):
            triang = heighten(triang, heat, allow51=allow51, allow15=False)
            triang = randomizeFast(triang)
        for i in range(100*iterations):
            triang = randomizeFast(triang)
        return triang

    while curr_time - init_time < time_in_millis:
        # try to keep 5 vertices
        if global_min.numVertices() < 5:
            while global_min.numVertices() < 5:
                global_min = global_min.anyMove(1,5)

        has5Vertices = (global_min.numVertices() == 5)
        if not firstTimeReached5Verts and has5Vertices:
            firstTimeReached5Verts = True
            min_tri = global_min.copy()
            tri = min_tri.copy()
        allow51 = not has5Vertices
        
        loop += 1
        size = tri.numPentachora()
        tri = smallReduce(tri, allow51)

        # Check if standard triangulation of S^4 (the double of a 4-simplex)
        if tri.numPentachora() == 2 and tri.isoSig() == "cPkbbbbaaaaaaaa":
            return (tri, True)
        
        #vertInfo(tri)
        if tri.numPentachora() >= min_tri.numPentachora():
            if tri.numPentachora() == min_tri.numPentachora():
                min_tri = tri.copy()
            if loop - last_loop > 30:
                if tri.numPentachora() == global_min.numPentachora():
                    reached_min += 1
                if min_tri.numPentachora() < 20:
                    min_tri = spike(min_tri, iterations=3, heat=3, allow51=allow51)
                else:
                    min_tri = spike(min_tri, iterations=3, heat=3, allow51=allow51)
                print 'permanently heightened:', min_tri.numPentachora()
                last_loop = loop
            
            if reached_min % 30 == 29:
                if global_min.numPentachora() < 20:
                    min_tri = spike(min_tri, iterations=10, heat=4, allow51=allow51)
                elif reached_min % 30 == 29:
                    min_tri = spike(min_tri, iterations=60, heat=60, allow51=allow51)
                    reached_min = 0
                else:
                    min_tri = spike(min_tri, iterations=60, heat=60, allow51=allow51)
                    reached_min = 0
                print 'spike heightened:', min_tri.numPentachora()
                reached_min += 1
                
            tri = min_tri.copy()
        elif tri.numPentachora() < min_tri.numPentachora():
            min_tri = tri.copy()
            if tri.numPentachora() < global_min.numPentachora():
                global_min = tri.copy()
                reached_min = 0
            last_loop = loop
        
        if tri.numPentachora() == min_tri.numPentachora():
            n = 2
            print 'heightening, current size:', tri.numPentachora(), 'verts:', tri.numVertices(), 'reached_min', reached_min
            tri = heighten(tri, n, allow51=allow51, allow15=False)
            print 'heightened size:', tri.numPentachora(), 'global min:', global_min.numPentachora(), 'verts in min:', global_min.numVertices()
        curr_time = int(round(time.time() * 1000))
    return (global_min, False)

def simplifyForTime(tri, time_in_millis=60000, allow51=True):
    """
    Attempts to simplify the triangulation for a time period of `time_in_millis`
    milliseconds. Returns a new triangulation (leaves current triangulation unmodified).
    """
    tri = tri.copy()
    size = tri.numPentachora()
    print 'size:', size
    import time
    init_time = int(round(time.time() * 1000))
    curr_time = init_time
    loop = 0
    min_tri = tri.copy()
    global_min = tri.copy()
    last_loop = 0
    heat = 1
    while curr_time - init_time < time_in_millis and tri.numPentachora() > 2:
        loop += 1
        size = tri.numPentachora()
        tri = smallReduce(tri, allow51)
        if tri.numVertices() <= 5:
            allow51 = False
        else:
            allow51 = True
        #vertInfo(tri)
        if tri.numPentachora() >= min_tri.numPentachora():
            if loop - last_loop > heat*50:
                for i in range(10):
                    min_tri = heighten(min_tri, heat, allow51)
                    min_tri = randomizeFast(min_tri)
                for i in range(20):
                    min_tri = randomizeFast(min_tri)
                print 'permanently heightened:', min_tri.numPentachora()
                heat += 1
            tri = min_tri.copy()
            print 'heat:', heat
        elif tri.numPentachora() < min_tri.numPentachora():
            min_tri = tri.copy()
            if tri.numPentachora() < global_min.numPentachora():
                global_min = tri.copy()
                heat = 1
            last_loop = loop
            
        if tri.numPentachora() == min_tri.numPentachora():
            print 'heightening, current size:', tri.numPentachora(), 'global min:', global_min.numPentachora()
            n = 2
            if loop % 20 == 5:
                n = 15
            tri = heighten(tri, n, allow51)
            print 'heightened size:', tri.numPentachora()
        curr_time = int(round(time.time() * 1000))
    return global_min

def simplifyUntil(input_fname, output_fname, pent_size=2):
    round = 0
    tri = FourTriangulation.loadFromFile(input_fname)
    allow51 = True
    
    """
    for i in range(2):
        tri = tri.anyMove(1,5)
    allow51 = False
    print 'cyc len', len(tri.cycleByDim(0)), tri.numPentachora()
    """
    
    vertInfo(tri)
    fVector(tri)
    size = tri.numPentachora()
    print 'size:', size
    
    min_tri = tri.copy()
    loop = 0
    
    while size > pent_size:
        loop += 1
        size = tri.numPentachora()
        tri = smallReduce(tri, allow51) 
        #vertInfo(tri)
        if tri.numPentachora() > min_tri.numPentachora():
            tri = min_tri.copy()
        elif tri.numPentachora() <= min_tri.numPentachora():
            if tri.numPentachora() < min_tri.numPentachora():
                round += 1
                tri.writeToFile(output_fname)
            min_tri = tri.copy()
            
        if tri.numPentachora() == size:
            print 'heightening, current size:', tri.numPentachora()
            n = 2
            if loop % 20 == 5:
                n = 15
            tri = heighten(tri, n, allow51=allow51)
            print 'heightened size:', tri.numPentachora()

    tri.writeToFile(output_fname)
    
#simplifyUntil('C:/pr2011/4dim_triangulation/homotopy_4sphere2.tri', 'C:/pr2011/4dim_triangulation/homotopy_4sphere2_simp.tri')
    
