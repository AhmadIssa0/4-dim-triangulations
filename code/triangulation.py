
import sys
import copy
import os
import random
from snappea_triangulation import *
import itertools
from matrix import *


class FourTriangulation:
    pachnerCount = 0
    debugFile = open("debug.txt", "w")
    # pentachora = list of Pentachora
    # cusps = number of cusps
    
    def __init__(self, n=0, tri_title=None):
        # tri_title is the title of the triangulation as written in the triangulation file
        self.pentachora = []
        self.cusps = 0
        self.pent_by_num = {}
        for i in range(n):
            self.addPentachoron(i)
            #pent.cusp_no = range(5*i+1, 5*i+6)
            
        self.tri_title = tri_title
        
    def addPentachoron(self, num=None):
        if num is None:
            num = self.nextPentIndex()
        pent = Pentachoron(num)
        self.pent_by_num[num] = pent
        self.pentachora.append(pent)
        return pent

    def removePentachoron(self, num):
        # Removes Pentachora with getPentNo() == num
        del self.pent_by_num[num]
        for i in range(self.numPentachora()):
            if self.pentachora[i].getPentNo() == num:
                del self.pentachora[i]
                break
                
    def link(self):
        lnk = Triangulation(5*self.numPentachora())
        #print 'START'
        for p_no in range(self.numPentachora()):
            pent = self.pentachora[p_no]
            for tet_no in range(5):
                tet = oppositeTetrahedron(tet_no)
                for tri_no in tet:
                    # looking at tetrahedron opposite vertex tet_no of pent
                    # want to glue face of tet opposite tri_no
                    # replace tri_no with tet_no in tet
                    # then look at what the defined tet is glued to
                    # use this gluing to see where the triangle should be glued
                    tri = filter(lambda x: x != tri_no, tet)
                    tet_g = tri + [tet_no] # tetrahedral face which will determine gluing of tri
                                           # opposite tri_no
                    pent2 = pent.neighbour[tri_no]
                    gluing = pent.getGluing(tri_no)
                    #gluing[tri_no]
                    #print 'pent', p_no, 'tet opposite:', tet_no, 'face:', tri, 'glued to pent ', \
                    #      pent2.n, 'tet opposite:', gluing[tet_no], 'tri', map(lambda x: gluing[x], tri), 'gluing:', gluing
                    
                    glued_tet = [gluing[i] for i in tet]
                    for i in range(len(glued_tet)):
                        if glued_tet[i] > gluing[tet_no]:
                            glued_tet[i] -= 1
                    #print 'glued_tet', glued_tet
                    t1 = lnk.getTetrahedron(p_no*5 + tet_no)
                    t2 = lnk.getTetrahedron(pent2.n*5 + gluing[tet_no])
                    if tri_no > tet_no:
                        tri_no -= 1
                    t1.joinTo(tri_no, t2, glued_tet)
        #print 'DONE'
        return lnk
        
    def cycleByDim(self, d):
        p_nos = [p.getPentNo() for p in self.getPentachoronList()]
        lst = []
        cyc = []
        for p_no in p_nos:
            for face in map(list, list(itertools.combinations(range(5),1+d))):
                if [p_no, set(face)] not in lst:
                    c = self.cycle(p_no, face)
                    cyc.append(c)
                    lst += [[x[0], set(x[1])] for x in c]
        return cyc

    def edgeDegrees(self):
        return sorted([len(x) for x in self.cycleByDim(1)])
    
    def numVertices(self):
        return len(self.cycleByDim(0))
        
    def fVector(self):
        f = []
        for i in range(5):
            f.append(len(self.cycleByDim(i)))
        return f
        
    def complexity(self):
        return sum(self.fVector())
        
    
    def vertexCycles(self):
        p_nos = [p.getPentNo() for p in self.getPentachoronList()]
        lst = []
        cyc = []
        for p_no in p_nos:
            for i in range(5):
                if [p_no, [i]] not in lst:
                    c = self.cycle(p_no, [i])
                    cyc.append(c)
                    lst += c
        return cyc
        
    def cycle(self, p_no, lst):
        # returns a list of [pentachoron number, [cycle]]
        # which get identified to lst of pentachoron p_no
        cyc = [[p_no, lst[:]]]
        prev_size = 0
        size = 1
        
        while size > prev_size:
            for i in range(prev_size,size):
                n, face = cyc[i]
                # vertex opposite faces of tetrahedra containing face
                opposite_verts = filter(lambda x: x not in face, range(5))
                pent = self.getPentachoronByNum(n)
                for v in opposite_verts:
                    g = pent.getGluing(v)
                    if pent.getNeighbour(v) is None:
                        continue
                    n_nbr = pent.getNeighbour(v).getPentNo()
                    c = [n_nbr, [g[j] for j in face]]
                    if c not in cyc:
                        cyc.append(c)
            prev_size, size = size, len(cyc)
        return cyc
        
    def boundaryCycle(self, t_no, lst):
        # t_no = [p_no, opp_vert]
        # returns a list of [[p_no, opp_vert], cycle]
        # which get identified to lst of t_no, along the boundary of the mfd
        cyc = [[t_no, lst[:]]]
        prev_size = 0
        size = 1
        
        while size > prev_size:
            for i in range(prev_size, size):
                [n, opp_vert], face = cyc[i]
                tet_verts = [x for x in range(5) if x != opp_vert]
                opposite_verts = filter(lambda x: x not in face, tet_verts)
                
                for v in opposite_verts:
                    [nt, nv], tet_opp, g = self.boundaryGluing([n, opp_vert], v)
                    c = [[nt, nv], [g[j] for j in face]]
                    if c not in cyc:
                        cyc.append(c)
            prev_size, size = size, len(cyc)
        return cyc    
                
    def triFaceCycle(self, lst):
        # list will contain [pent_no, vert opposite tet, vert opposite tri, tri] 
        pent_no, tet_no, tri_no, tri = lst[-1]
        
        neighbour = self.pentachora[pent_no].neighbour[tet_no]
        gluing = self.pentachora[pent_no].getGluing(tet_no)
        n_tet_no = gluing[tet_no]
        n_tri_no = gluing[tri_no]
        n_tri = map(lambda x: gluing[x], tri)
        
        elt = [neighbour.n, n_tet_no, n_tri_no, n_tri]
        if elt not in lst:
            lst.append(elt)
        else:
            return lst
        
        elt2 = [neighbour.n, n_tri_no, n_tet_no, n_tri]
        if elt2 not in lst:
            lst.append(elt2)
        else:
            return lst
        
        return self.triFaceCycle(lst)
        
    def boundaryGluing(self, t_no, v):
        # returns [[pent number, opposite vertex], opposite vertex of tet, gluing map]
        p_no, opp_vert = t_no
        pent = self.getPentachoronByNum(p_no)
        
        p_verts = [x for x in range(5) if x != opp_vert] # verts of pent corresponding to tet's [0,1,2,3]
        
        nbr = pent.getNeighbour(v)
        if nbr is None:
            g = range(5)
            g[opp_vert] = v
            g[v] = opp_vert
            return [[p_no, v], opp_vert, g]
            #return [[p_no, v], opp_vert, range(5)]
        else:
            gluing = pent.getGluing(v)
            nbr_v = gluing[v]
            nbr_opp_vert = gluing[opp_vert] # opposite vertex in nbr
            
            #pivot around edge until we hit a free face
            while nbr.getNeighbour(nbr_opp_vert) is not None:
                g = nbr.getGluing(nbr_opp_vert)
                for i in range(5):
                    gluing[i] = g[gluing[i]]
                nbr = nbr.getNeighbour(nbr_opp_vert)
                nbr_opp_vert, nbr_v = g[nbr_v], g[nbr_opp_vert]
                
            gluing[opp_vert] = nbr_opp_vert
            gluing[v] = nbr_v
            t2_no = [nbr.getPentNo(), nbr_opp_vert]
            return [t2_no, nbr_v, gluing]
            
    def boundaryTetrahedra(self):
        pents = self.getPentachoronList()
        
        tets = [] # each tet is [p_no, opposite vertex in pent]
        for p in pents:
            for i in range(5):
                if p.getNeighbour(i) is None:
                    tets.append([p.getPentNo(), i])
        return tets
        
    def boundaryIsomorphisms(self, tri2, debug=False):
        all_isos = self.boundaryIsos(tri2, debug=debug)
        
        # change the format of the result
        # a pairing from boundaryIsos looks like: 
        # [[138, [0, 1, 3, 4]], [138, [0, 1, 3, 4]]]
        # Want this to look like:
        # [[138,2],[138,2], [0,1,2,3,4]]
        # meaning tet [138,2] glues to tet [138,2] via [0,1,2,3,4]
        def complement(lst):
            for i in range(5):
                if i not in lst:
                    return i
                    
        def changeFormat(iso):
            n_iso = []
            for pairing in iso:
                t1, t2 = pairing
                v1 = complement(t1[1])
                v2 = complement(t2[1])
                g = range(5)
                g[v1] = v2
                for i in range(5):
                    for j in range(4):
                        if t1[1][j] == i:
                            g[i] = t2[1][j]
                            break
                n_iso.append([[t1[0], v1], [t2[0], v2], g])
            return n_iso
            
        return map(changeFormat, all_isos)
        
        
    def boundaryIsos(self, tri2, bdryTet1=None, bdryTet2=None, ind=0, iso=None, debug=False):
        # return a list of all possible ways to extend iso (partial isomorphism)
        # ind is the index in bdryTet1 that we're currently trying to extend over
        
        # check consistency of last pairing
        if ind > 0:
            if debug:
                print 'checking consistency of:', iso
            def facePairing(triang, tet, vert):
                p_no = tet[0]
                opp_vert = [x for x in range(5) if x not in tet[1]][0]
                face = [x for x in tet[1] if x != vert]
                bCycle = triang.boundaryCycle([tet[0], opp_vert], face)
                origTet, orig_cyc = bCycle[0]
                gluedTet, cyc = bCycle[1]
                orig_cyc.append(origTet[1])
                cyc.append(gluedTet[1])
                
                i = [x for x in range(5) if x not in orig_cyc][0]
                j = [x for x in range(5) if x not in cyc][0]
                orig_cyc.append(i)
                cyc.append(j)
                
                n = gluedTet[0]
                g = []
                for i in tet[1]:
                    for j in range(5):
                        if orig_cyc[j] == i:
                            g.append(cyc[j])
                return [n, g]
                
            def mappedTo(tet):
                def complement(t):
                    return [x for x in range(5) if x not in t]
                for pairing in iso:
                    t1, t2 = pairing
                    if t1[0] == tet[0] and complement(t1[1]) == complement(tet[1]):
                        g = []
                        for v in tet[1]:
                            for i in range(4):
                                if t1[1][i] == v:
                                    g.append(t2[1][i])
                                    break
                        return [t2[0], g]
                return None
               
            t1, t2 = iso[-1]
            for i in range(len(t1[1])):
                fp1 = facePairing(self, t1, t1[1][i])
                fp2 = facePairing(tri2, t2, t2[1][i])
                im_tet = mappedTo(fp1)
                if im_tet is not None and fp2 != im_tet:
                    if debug:
                        print 'failed consistency', iso
                        print 'bdry1', bdryTet1
                        print 'bdry2', bdryTet2
                        print 'im_tet', im_tet
                        print 'fp2', fp2
                    return []
                    
                    
            #print 'passed consistency', iso
            # if we reached here it's consistent
                
        def complement(t):
            return [t[0], [x for x in range(5) if x != t[1]]]
            
        if iso is None:
            iso = []
        if bdryTet1 is None:
            bdryTet1 = self.boundaryTetrahedra()
            bdryTet2 = tri2.boundaryTetrahedra()
            bdryTet1 = map(complement, bdryTet1)
            bdryTet2 = map(complement, bdryTet2)
        if ind >= len(bdryTet1):
            return [iso]
            
        all_isos = []
        tet = bdryTet1[ind]
        for tet2 in bdryTet2:
            paired = False
            for x in iso:
                def complement(t):
                    return [x for x in range(5) if x not in t]
                if tet2[0] == x[1][0] and complement(tet2[1]) == complement(x[1][1]):
                    paired = True
                    break
            if not paired: # haven't paired with tet2 in iso
                for perm in itertools.permutations(tet2[1]):
                    perm = list(perm)
                    all_isos += self.boundaryIsos(tri2, bdryTet1, bdryTet2, ind+1, iso+[[tet,[tet2[0],perm]]], debug=debug)
        return all_isos
                    
        
    def boundaryTriangulation(self):
        pents = self.getPentachoronList()
        
        tets = self.boundaryTetrahedra() # each tet is [p_no, opposite vertex in pent]
                    
        def complement(lst, n=5):
            return [i for i in range(n) if i not in lst]
            
        tri = Triangulation(len(tets))
        triTets = tri.getTetrahedraList()
        
        def getTet(p_no, opp_vert):
            return triTets[tets.index([p_no, opp_vert])]
        
        for i in range(len(tets)):
            p_no, opp_vert = tets[i]
            pent = self.getPentachoronByNum(p_no)
            p_verts = complement([opp_vert]) # vertices of pent corresponding to tet's [0,1,2,3]
            for v in p_verts:
                nbr = pent.getNeighbour(v)
                if nbr is None:
                    np_verts = complement([v])
                    t1 = getTet(p_no, opp_vert)
                    t2 = getTet(p_no, v)
                    g = [0]*4
                    for j in range(4):
                        if p_verts[j] in np_verts:
                            g[j] = np_verts.index(p_verts[j])
                    g[p_verts.index(v)] = np_verts.index(opp_vert)
                    t1.joinTo(p_verts.index(v), t2, g)
                else:
                    gluing = pent.getGluing(v)
                    nbr_v = gluing[v]
                    nbr_opp_vert = gluing[opp_vert] # opposite vertex in nbr
                    
                    # pivot around edge until we hit an empty face
                    while nbr.getNeighbour(nbr_opp_vert) is not None:
                        g = nbr.getGluing(nbr_opp_vert)
                        for i in range(5):
                            gluing[i] = g[gluing[i]]
                        nbr = nbr.getNeighbour(nbr_opp_vert)
                        nbr_opp_vert, nbr_v = g[nbr_v], g[nbr_opp_vert]
                        
                    gluing[opp_vert] = nbr_opp_vert
                    gluing[v] = nbr_v
                    
                    t1 = getTet(p_no, opp_vert)
                    t2 = getTet(nbr.getPentNo(), nbr_opp_vert)
                    g = [0]*4
                    np_verts = complement([nbr_opp_vert])
                    for j in range(4):
                        if gluing[p_verts[j]] in np_verts:
                            g[j] = np_verts.index(gluing[p_verts[j]])
                    t1.joinTo(p_verts.index(v), t2, g)
        return tri
                    
    def boundaryPachner(self, tet_nos, faces, debug=False):
        # a tet number is of the form [p_no, opposite-vertex]
        nTri = self.copy()

        pent = []
        p_nos = []
        opp_verts = []
        for p_no, opp_vert in tet_nos:
            p_nos.append(p_no)
            opp_verts.append(opp_vert)
            pent.append(nTri.getPentachoronByNum(p_no))
        
        def complement(lst, n=4):
            return [i for i in range(n) if i not in lst]
            
        # vertices of abstract
        init_tets = [[i for i in range(5) if i != j] for j in range(5-len(tet_nos),5)]
        # Explanation for 1 dimension higher:
        # pairing vertices of abstract pentachora to pent[]'s vertices
        # vertex pairing: init_pent_pairing[i][j] <-> init_pents[i][j]
        # We really want an isomorphism between the local configuration
        # of pentachora and a subcomplex of the boundary of a 5-simplex.
        init_tet_pairing = [faces[0] + [x for x in range(5) if x not in faces[0] + [opp_verts[0]]]]

        for i in range(1,len(faces)):
            gluings = []
            if debug:
                print 'searching for', tet_nos[i]
            for j in range(i):
                gluings.append([-1]*5)
            for j in range(i):
                for opp_vert in [x for x in range(5) if x not in faces[j] + [opp_verts[j]]]:
                    [n_p_no, ov], tov, g = nTri.boundaryGluing(tet_nos[j], opp_vert)
                    if debug:
                        print tet_nos[j], 'glued to', [n_p_no, ov], 'via vertex', opp_vert
                        print tov, 'is opposite vertex of', [n_p_no, ov]
                    nbr = nTri.getPentachoronByNum(n_p_no)
                    #if nbr == pent[i]:
                    if [n_p_no, ov] == tet_nos[i]:
                        for k in range(4):
                            if init_tet_pairing[j][k] != opp_vert:
                                abs_v = init_tets[j][k]
                                gluings[j][abs_v] = g[init_tet_pairing[j][k]]
                        v = [x for x in range(5) if x not in init_tets[j]][0]
                        gluings[j][v] = g[opp_vert]
                        
                        if -1 in gluings[j]:
                            gluings[j][gluings[j].index(-1)] = [x for x in range(5) if x not in gluings[j]][0]
                        
                        break
            for j in range(i):
                if gluings[j] != gluings[0]:
                    # pachner move not valid
                    #print gluings
                    if debug:
                        print 'pachner not valid', gluings
                    return self
            init_tet_pairing.append([gluings[0][init_tets[i][j]] for j in range(4)])
            
        for i in range(len(init_tet_pairing)):
            pairing = init_tet_pairing[i]
            if sorted(pairing) != [x for x in range(5) if x != opp_verts[i]]:
                if debug:
                    print 'invalid pachner move', init_tet_pairing
                return self # pachner move not valid
        
        npent = nTri.addPentachoron(nTri.nextPentIndex())
        for i in range(len(tet_nos)):
            p_no, opp_vert = tet_nos[i]
            v = complement(init_tets[i],5)[0]
            g = [-1]*5
            v2 = complement(init_tet_pairing[i],5)[0]
            g[v2] = v
            
            for j in range(len(init_tets[i])):
                g[init_tet_pairing[i][j]] = init_tets[i][j]
                
            pent[i].joinTo(opp_vert, npent, g)
        
        return nTri
        
    def pachner(self, pent_nos, faces, countPachner=True):
        # pent_nos is a list of pentachora numbers
        # face is a list of face vertices of tetrahedron number pent_nos[0]
        #   which is common to all pentachora

        nTri = self.copy()
        #init_fVec = self.fVector()
        pent = []
        for p_no in pent_nos:
            pent.append(nTri.getPentachoronByNum(p_no))
            
        # TO DO: set cusp numbers
        
        # 1,2,3,4,5,6 are vertices of a 5-dimensional simplex
        init_pents = [[i for i in range(1,7) if i != j] for j in range(7-len(pent_nos),7)]
        final_pents = [[i for i in range(1,7) if i != j] for j in range(1,7-len(pent_nos))]
        
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
                        free.append([i, tet])
            return free
                        
        init_pents_free = calcFreeFaces(init_pents)
        init_pents_free.sort(key=lambda x: x[1])
        final_pents_free = calcFreeFaces(final_pents)
        final_pents_free.sort(key=lambda x: x[1])
        
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
            
        final_pents_internal = calcInternalFaces(final_pents)
        
        def complement(lst):
            return [i for i in range(5) if i not in lst]
            
        # pairing vertices of abstract pentachora to pent[]'s vertices
        # vertex pairing: init_pent_pairing[i][j] <-> init_pents[i][j]
        # We really want an isomorphism between the local configuration
        # of pentachora and a subcomplex of the boundary of a 5-simplex.
        
        init_pent_pairing = [faces[0] + complement(faces[0])]
        for i in range(1,len(faces)):
            gluings = []
            for j in range(i):
                gluings.append([-1]*6)
            for j in range(i):
                for opp_vert in complement(faces[j]):
                    nbr = pent[j].getNeighbour(opp_vert)
                    if nbr == pent[i]:
                        g = pent[j].getGluing(opp_vert)
                        for k in range(5):
                            if init_pent_pairing[j][k] != opp_vert:
                                abs_v = init_pents[j][k]
                                gluings[j][abs_v-1] = g[init_pent_pairing[j][k]]
                        v = [x for x in range(1,7) if x not in init_pents[j]][0]-1
                        gluings[j][v] = g[opp_vert]
            for j in range(i):
                if gluings[j] != gluings[0]:
                    # pachner move not valid
                    #print gluings
                    return False
            init_pent_pairing.append([gluings[0][init_pents[i][j]-1] for j in range(5)])
            
        for pairing in init_pent_pairing:
            if sorted(pairing) != range(5):
                print 'invalid pachner move'
                return False # pachner move not valid
                
        # only add pentachora once we know the pachner move is valid
        npent = [] # new pentachora
        for i in range(6-len(pent_nos)):
            n = nTri.nextPentIndex()
            npent.append(nTri.addPentachoron(n))
                             
        # map from free faces of abstract pentachora to faces of actual initial pentachora
        free_abs_to_actual_init = [[p_no, [init_pent_pairing[p_no][init_pents[p_no].index(v)] for v in face]] for p_no, face in init_pents_free]
        
        final_pent_pairing = [range(5) for i in range(len(final_pents))]
        free_abs_to_actual_final = [[t_no, [final_pent_pairing[t_no][final_pents[t_no].index(v)] for v in face]] for t_no, face in final_pents_free]
        internal_abs_to_actual_final = [[t_no, [final_pent_pairing[t_no][final_pents[t_no].index(v)] for v in face]] for t_no, face in final_pents_internal]
        
        def glue_faces(p1, p2, f1, f2):
            # t1 tetrahedron face f1 gets glued to f2 of t2
            v1 = FourTriangulation.oppositeVertex(f1)
            v2 = FourTriangulation.oppositeVertex(f2)
            g = {}
            g[v1] = v2
            for i in range(len(f1)):
                g[f1[i]] = f2[i]
            p1.joinTo(v1, p2, [g[0],g[1],g[2],g[3],g[4]])
            
        def init_pent_to_final_pent(n_init_pent, face):
            for i in range(len(free_abs_to_actual_init)):
                t_no, f = free_abs_to_actual_init[i]
                if t_no == n_init_pent and set(face).issubset(f):
                    new_t_no, new_f = free_abs_to_actual_final[i]
                    dic = {}
                    for j in range(4):
                        dic[f[j]] = new_f[j]
                    return new_t_no, [dic[face[j]] for j in range(4)]
            
        # do the internal gluings of new pentachora
        for i in range(len(internal_abs_to_actual_final)/2):
            nt1, f1 = internal_abs_to_actual_final[2*i]
            nt2, f2 = internal_abs_to_actual_final[2*i+1]
            glue_faces(npent[nt1], npent[nt2], f1, f2)
           
        gluing_arguments = []
        for i in range(len(free_abs_to_actual_init)):
            nt1, f1 = free_abs_to_actual_init[i]
            v = FourTriangulation.oppositeVertex(f1)
            nbr = pent[nt1].getNeighbour(v)
            g = pent[nt1].getGluing(v)
            if nbr is None:
                continue
            outer_face = [g[j] for j in f1]
            # deal with case where external face is of an internal pentachoron
            for j in range(len(pent)):
                if nbr == pent[j]:
                    nt, nf = init_pent_to_final_pent(j, outer_face)
                    nbr = npent[nt]
                    g = nbr.getGluing(FourTriangulation.oppositeVertex(nf))
                    outer_face = nf
                    break
            nt2, f2 = free_abs_to_actual_final[i]
            gluing_arguments.append([nbr, npent[nt2], outer_face, f2])
            #glue_faces(nbr, npent[nt2], outer_face, f2)
            
        for args in gluing_arguments:
            apply(glue_faces, args)
            
        for p_no in pent_nos:
            nTri.removePentachoron(p_no)
        nTri.renumberTriangulation()

        if countPachner:
            FourTriangulation.pachnerCount += 1
            FourTriangulation.debugFile.write(str(FourTriangulation.pachnerCount) + ' ' +
                                              str(nTri.numPentachora()) + '\r\n')
        return nTri
        
    def two_to_four(self, pent_nos, faces):
        # pent_nos is a list of pentachora numbers
        # face is a list of face vertices of tetrahedron number pent_nos[0]
        #   which is common to both pentachora
        
        nTri = self.copy()
        
        pent = []
        for p_no in pent_nos:
            pent.append(nTri.getPentachoronByNum(p_no))
            
        npent = [] # new pentachora
        for i in range(4):
            n = nTri.nextPentIndex()
            npent.append(nTri.addPentachoron(n))
            
        # TO DO: set cusp numbers
        
        # 1,2,3,4,5,6 are vertices of a 5-dimensional simplex
        init_pents = [[i for i in range(1,7) if i != j] for j in range(5,7)]
        final_pents = [[i for i in range(1,7) if i != j] for j in range(1,5)]
        
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
                        free.append([i, tet])
            return free
                        
        init_pents_free = calcFreeFaces(init_pents)
        init_pents_free.sort(key=lambda x: x[1])
        final_pents_free = calcFreeFaces(final_pents)
        final_pents_free.sort(key=lambda x: x[1])
        
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
            
        final_pents_internal = calcInternalFaces(final_pents)
        # pairing vertices of abstract pentachora to pent[]'s vertices
        # currently specific to two_to_four, need to generalise
        init_pent_pairing = [faces[0] + [FourTriangulation.oppositeVertex(faces[0])],
                             faces[1] + [FourTriangulation.oppositeVertex(faces[1])]]
                             
        # map from free faces of abstract pentachora to faces of actual initial pentachora
        free_abs_to_actual_init = [[t_no, [init_pent_pairing[t_no][init_pents[t_no].index(v)] for v in face]] for t_no, face in init_pents_free]
        
        final_pent_pairing = [range(5) for i in range(len(final_pents))]
        free_abs_to_actual_final = [[t_no, [final_pent_pairing[t_no][final_pents[t_no].index(v)] for v in face]] for t_no, face in final_pents_free]
        internal_abs_to_actual_final = [[t_no, [final_pent_pairing[t_no][final_pents[t_no].index(v)] for v in face]] for t_no, face in final_pents_internal]
        
        def glue_faces(p1, p2, f1, f2):
            # t1 tetrahedron face f1 gets glued to f2 of t2
            v1 = FourTriangulation.oppositeVertex(f1)
            v2 = FourTriangulation.oppositeVertex(f2)
            g = {}
            g[v1] = v2
            for i in range(len(f1)):
                g[f1[i]] = f2[i]
            p1.joinTo(v1, p2, [g[0],g[1],g[2],g[3],g[4]])
            
        def init_pent_to_final_pent(n_init_pent, face):
            for i in range(len(free_abs_to_actual_init)):
                t_no, f = free_abs_to_actual_init[i]
                if t_no == n_init_pent and set(face).issubset(f):
                    new_t_no, new_f = free_abs_to_actual_final[i]
                    dic = {}
                    for j in range(4):
                        dic[f[j]] = new_f[j]
                    return new_t_no, [dic[face[j]] for j in range(4)]
            
        # do the internal gluings of new pentachora
        for i in range(len(internal_abs_to_actual_final)/2):
            nt1, f1 = internal_abs_to_actual_final[2*i]
            nt2, f2 = internal_abs_to_actual_final[2*i+1]
            glue_faces(npent[nt1], npent[nt2], f1, f2)
            
        for i in range(len(free_abs_to_actual_init)):
            nt1, f1 = free_abs_to_actual_init[i]
            v = FourTriangulation.oppositeVertex(f1)
            nbr = pent[nt1].getNeighbour(v)
            g = pent[nt1].getGluing(v)
            outer_face = [g[j] for j in f1]
            # deal with case where external face is of an internal pentachoron
            for j in range(len(pent)):
                if nbr == pent[j]:
                    nt, nf = init_pent_to_final_pent(j, outer_face)
                    nbr = npent[nt]
                    g = nbr.getGluing(FourTriangulation.oppositeVertex(nf))
                    outer_face = nf
                    break
            nt2, f2 = free_abs_to_actual_final[i]
            glue_faces(nbr, npent[nt2], outer_face, f2)
            
        for p_no in pent_nos:
            nTri.removePentachoron(p_no)
        nTri.renumberTriangulation()
        return nTri
        
    def idealToFinite(self):
        # Generalisation of the Regina source code (which deals with 3 dimensional case)
        # WARNING: Calls self.renumberTriangulation()
        
        def perm5(i,j):
            p = range(5)
            p[i], p[j] = j, i
            return p
            
        self.renumberTriangulation()
        numOldPent = self.numPentachora()
        oldPents = []
        for i in range(numOldPent):
            oldPents.append(self.getPentachoronByNum(i))
            
        tri = FourTriangulation()
        
        # 4 different types of pentachora in subdivision
        pentInterior = []
        interior = [] # face interior to tetrahedron
        vertex = []
        edge = []

        for k in range(numOldPent):
            pentInterior.append([])
            interior.append([])
            edge.append([])
            vertex.append([])
            
            for i in range(5):
                interior[k].append([None]*5)
                edge[k].append([])
                vertex[k].append([])
                    
                for j in range(5):
                    edge[k][i].append([None]*5)
                    vertex[k][i].append([None]*5)
        
            for i in range(5):
                pentInterior[k].append(tri.addPentachoron())
                
                for j in range(5):
                    if j != i:
                        interior[k][i][j] = tri.addPentachoron()
                    for m in range(5):
                        if i != j and i != m and j != m:
                            edge[k][i][j][m] = tri.addPentachoron()
                            vertex[k][i][j][m] = tri.addPentachoron()
                    
        # glue pentachora inside the same old pentachoron together
        for i in range(numOldPent):      
            for n in range(5):
                for m in range(5):
                    for k in range(5):
                        if m != n and m != k and n != k:
                            edge[i][n][m][k].joinTo(m, edge[i][m][n][k], perm5(m,n)) # ok
                            vertex[i][n][m][k].joinTo(m, vertex[i][m][n][k], perm5(m,n)) # ok
                            pentInterior[i][m].joinTo(n, interior[i][n][m], perm5(m,n)) # ok
                            
                            interior[i][n][m].joinTo(k, vertex[i][n][k][m], perm5(m,k)) # ok
                            edge[i][n][m][k].joinTo(k, edge[i][n][k][m], perm5(k,m)) # ok
                            
                            for j in [x for x in range(5) if x not in [n,m,k]]:
                                q = [x for x in range(5) if x not in [n,m,k,j]][0]
                                p = range(5)
                                p[j], p[q], p[k] = k, j, q
                                edge[i][n][m][k].joinTo(j, vertex[i][n][m][q], p) # ok
                                
        # global gluings
        for i in range(numOldPent):
            oldPent = oldPents[i]
            for j in range(5):
                if oldPent.getNeighbour(j) is not None:
                    oppPent = oldPent.getNeighbour(j).getPentNo()
                    p = oldPent.getGluing(j)
                    
                    for k in range(5):
                        if j != k:
                            interior[i][j][k].joinTo(j, interior[oppPent][p[j]][p[k]], p)
                        for m in range(5):
                            if j != k and k != m and m != j:
                                edge[i][j][k][m].joinTo(j, edge[oppPent][p[j]][p[k]][p[m]], p)
                                vertex[i][j][k][m].joinTo(j, vertex[oppPent][p[j]][p[k]][p[m]], p)
        return tri

    def renumberTriangulation(self):
        self.pent_by_num = {}
        for i in range(self.numPentachora()):
            self.pentachora[i].n = i
            self.pent_by_num[i] = self.pentachora[i]
            
    def numBoundaryTet(self):
        return len(self.boundaryTetrahedra())
            
    def nextPentIndex(self):
        pents = self.getPentachoronList()
        if len(pents) == 0:
            return 0
        return max(map(lambda tet: tet.getPentNo(), pents))+1
        
    def copy(self):
        # CONDITION: Requires the tetrahedra in self to be numbered 0..(n-1)
        nTri = FourTriangulation()
        p_nos = []
        for p in self.getPentachoronList():
            p_nos.append(p.getPentNo())
            nTri.addPentachoron(p.getPentNo())
        nTri.setNoCusps(self.getNoCusps())
        
        for p_no in p_nos:
            nTet = nTri.getPentachoronByNum(p_no)
            oldTet = self.getPentachoronByNum(p_no)
            nTet.setGluings(oldTet.getGluings())
            nTet.setCuspNos(oldTet.getCuspNos())
                
        for p_no in p_nos:
            nTet = nTri.getPentachoronByNum(p_no)
            oldTet = self.getPentachoronByNum(p_no)
            for j in range(5):
                if oldTet.getNeighbour(j) is not None:
                    neighbour_no = oldTet.getNeighbour(j).getPentNo()
                    nTet.setNeighbour(j, nTri.getPentachoronByNum(neighbour_no))
        return nTri
        
    def getPentachoronByNum(self, n):
        return self.pent_by_num[n]
        #for pent in self.getPentachoronList():
        #    if pent.getPentNo() == n:
        #        return pent
                
    def getPentachoronList(self):
        return self.pentachora[:]
        
    # return i-th Pentachoron in list of pentachora
    # this may be different to Pentachoron.getPentNo() == i
    def getPentachoron(self, i):
        return self.pentachora[i]
        
    def numPentachora(self):
        return len(self.pentachora)
        
    def setNoCusps(self, n):
        self.cusps = n
        
    def getNoCusps(self):
        return self.cusps
        
    @staticmethod
    def oppositeVertex(pent):
        for i in range(5):
            if i not in pent:
                return i
                
    def identifyCusps(self, c1, c2):
        # changes all cusps numbered either c1 or c2 to min(c1, c2)
        for pent in self.pentachora:
            cusps = pent.getCuspNos()
            for i in range(5):
                if cusps[i] in [c1, c2]:
                    pent.setCuspNo(i, min(c1,c2))
                    
    def getCuspNums(self):
        cusp_nums = []
        for pen in self.pentachora:
            for i in range(5):
                if pen.cusp_no[i] not in cusp_nums:
                    cusp_nums.append(pen.cusp_no[i])
        return cusp_nums
                
    @staticmethod
    def loadFromFile(fname):
        f = open(fname, 'r')
        lines = f.readlines()
        size = int(lines[2])
        tri = FourTriangulation(size)
        pents = tri.getPentachoronList()
        
        for i in range(size):
            nbrs = lines[4*i+4][1:].strip().split('\t')
            gluings = lines[4*i+5][1:].strip().split('\t')
            for j in range(5):
                if nbrs[j] != '-':
                    nbr = pents[int(nbrs[j])]
                    g = map(int, gluings[j])
                    pents[i].joinTo(j, nbr, g)
        return tri
            
                
    def writeToFile(self, foutname):
        tri = self.copy()
        tri.renumberTriangulation()
        
        f = open(foutname, 'w')
        f.write('% Triangulation\n')
        if tri.tri_title is not None:
            f.write(tri.tri_title + '\n')
        else:
            f.write('Triangulation\n')
        f.write(str(tri.numPentachora()) + '\n')
        
        for i in range(tri.numPentachora()):
            tet = tri.getPentachoron(i)
            f.write('\n\t')
            for j in range(5):
                if tet.getNeighbour(j) is not None:
                    f.write('{0}\t'.format(tet.getNeighbour(j).getPentNo()))
                else:
                    f.write('-\t')
            f.write('\n')

            for j in range(5):
                g = tet.getGluing(j)
                if len(g) != 5:
                    f.write('\t-----')
                else:
                    f.write('\t{0}{1}{2}{3}{4}'.format(g[0],g[1],g[2],g[3],g[4]))
            f.write('\n')
            f.write('\t{0}\t{1}\t{2}\t{3}\t{4}\n'.format(tet.getCuspNo(0),tet.getCuspNo(1),tet.getCuspNo(2),
                                                        tet.getCuspNo(3),tet.getCuspNo(4)))
        f.close()
         
    @staticmethod
    def coneThreeManifold(tri):
        """returns the FourTriangulation given by coning the three dimensional triangulation
           tri. Note that the resulting FourTriangulation has some free faces."""
        ftri = FourTriangulation(tri.numTetrahedra())
        T = tri.tetrahedra
        P = ftri.pentachora
        for i in range(tri.numTetrahedra()):
            tet = T[i]
            pent = P[i]
            for myFace in range(4):
                neighbour_num = tet.getNeighbour(myFace).n
                pent.joinTo(myFace, P[neighbour_num], tet.getGluing(myFace) + [4])
        return ftri
        
    @staticmethod
    def fourTriFrom3SphereBoundary(tri1, tri2, pairings):
        """returns a FourTriangulation given by coning the 3 dimensional triangulations tri1 and tri2
           then identifying tetrahedra in tri1 with those in tri2 given by pairings.
           pairings is a list of [tet1_num, tet2_num, [0,1,2,3]] where [a,b,c,d] = [0,1,2,3] here
           means that vertex 0 of tri1's tetrahedron tet1_num is mapped to tet2's tetrahedron tet2_num
        """
        ftri = FourTriangulation(tri1.numTetrahedra() + tri2.numTetrahedra())
        P = ftri.pentachora
        for i in range(tri1.numTetrahedra()*2):
            if i >= tri1.numTetrahedra():
                tet = tri2.tetrahedra[i - tri1.numTetrahedra()]
            else:
                tet = tri1.tetrahedra[i]
            pent = P[i]
            for myFace in range(4):
                neighbour_num = tet.getNeighbour(myFace).n
                if i >= tri1.numTetrahedra():
                    pent.joinTo(myFace, P[tri1.numTetrahedra() + neighbour_num], tet.getGluing(myFace) + [4])
                else:
                    pent.joinTo(myFace, P[neighbour_num], tet.getGluing(myFace) + [4])
                
        for tet1_num, tet2_num, gluing in pairings:
            pent1 = ftri.pentachora[tet1_num]
            pent2 = ftri.pentachora[tri1.numTetrahedra() + tet2_num]
            pent1.joinTo(4, pent2, gluing + [4])
        return ftri
        
    def glueAlongBoundary(self, tri2, bdry_iso):
        tri1 = self
        n1 = tri1.numPentachora()
        n2 = tri2.numPentachora()
        print 'n1',n1,'n2',n2
        ntri = FourTriangulation()
        
        tri1_pents = tri1.getPentachoronList()
        tri2_pents = tri2.getPentachoronList()
        
        for p in tri1_pents:
            ntri.addPentachoron(p.getPentNo())
            
        next_valid = tri1.nextPentIndex() # next valid index
        
        
        for p in tri2_pents:
            ntri.addPentachoron(next_valid + p.getPentNo())

        ntri_pents = ntri.getPentachoronList()
        
        for i in range(n1+n2):
            if i >= n1:
                p = tri2_pents[i - n1]
            else:
                p = tri1_pents[i]
            pent = ntri_pents[i]
            for myFace in range(5):
                if p.getNeighbour(myFace) is not None:
                    neighbour_num = p.getNeighbour(myFace).getPentNo()
                    if i >= n1:
                        pent.joinTo(myFace, ntri.getPentachoronByNum(next_valid + neighbour_num), p.getGluing(myFace))
                    else:
                        pent.joinTo(myFace, ntri.getPentachoronByNum(neighbour_num), p.getGluing(myFace))
                        
        for t1_no, t2_no, g in bdry_iso:
            p1_no, v1 = t1_no
            p2_no, v2 = t2_no
            ntri.getPentachoronByNum(p1_no).joinTo(v1, ntri.getPentachoronByNum(next_valid + p2_no), g)
            
        return ntri

    def glueAlongBoundaries(self, tri1, debug=False):
        tri = self
    
        all_isos = tri.boundaryIsomorphisms(tri1, debug=debug)
        if debug:
            print 'number of isos', len(all_isos)
            print 'all_isos', all_isos

        assert len(all_isos) > 0, "No combinatorial isomorphism between boundary components."
        if len(all_isos) == 0:
            return False # failed
        else:
            iso = all_isos[0]
            c_tri = tri.copy()
            #n_tri = c_tri.glueAlongBoundary(tri, [[[138,2],[138,2],[0,1,2,3,4]], 
            #                            [[139,0],[139,0],[0,1,2,3,4]]])
            n_tri = c_tri.glueAlongBoundary(tri1, iso)
            n_tri = n_tri.fastSimplify()
            return n_tri
    
    def random23(self):
        tri = self
        lst = []
        for t_no in tri.boundaryTetrahedra():
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
        
    def boundaryValidMoves(self,a,b):
        # a-b pachner move
        tri = self
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
        
    def validMoves(self, a, b):
        tri = self
        lst = []
        valid_moves = []
        c_tri = tri.copy()
        for p in tri.getPentachoronList():
            p_no = p.getPentNo()
            for face in map(list, list(itertools.combinations(range(5),b))):
                if [p_no, face] not in lst:
                    c = tri.cycle(p_no, face)
                    if len(c) == a: # found possible a-b move
                        p_nos = [n for n, f in c]
                        if len(set(p_nos)) == a: # a distinct pentachora
                            faces = [f for n, f in c]
                            ntri = c_tri.pachner(p_nos, faces)
                            if ntri != tri:
                                c_tri = tri.copy()
                                valid_moves.append([p_nos, faces])
                    lst += c
        return valid_moves
        
    def random32(self):
        tri = self
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
        
    def random41(self, debug=False):
        tri = self
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
                            #print '4 distinct', t_nos
                            ntri = tri.boundaryPachner(t_nos, faces, debug)
                            if ntri != tri:
                                return ntri
                    lst += c
        return tri
        
    def anyMove(self, a, b):
        # tries to perform some pachner a-b move
        # a-b pachner move
        tri = self
        lst = []
        pents = tri.getPentachoronList()
        r = random.randint(0, len(pents)-1)
        pents = pents[r:] + pents[:r]
        for p in pents:
            p_no = p.getPentNo()
            for face in map(list, list(itertools.combinations(range(5),b))):
                if [p_no, set(face)] not in lst:
                    c = tri.cycle(p_no, face)
                    if len(c) == a: # found possible a-b move
                        p_nos = [n for n, f in c]
                        if len(set(p_nos)) == a: # a distinct pentachora
                            faces = [f for n, f in c]
                            ntri = tri.pachner(p_nos, faces)
                            if ntri != False:
                                ntri.renumberTriangulation()
                                return ntri
                    lst += [[x[0], set(x[1])] for x in c]
        return tri
        
    def fastSimplify(self, allow51 = True):
        # modifies the triangulation
        # try 5-1 and 4-2 moves
        tri = self
        size = tri.numPentachora()
        if allow51:
            tri = tri.anyMove(5,1)
        tri = tri.anyMove(4,2)
        while tri.numPentachora() < size:
            size = tri.numPentachora()
            if allow51:
                tri = tri.anyMove(5,1)
            tri = tri.anyMove(4,2)
        return tri
    
    def boundaryFastSimplify(self):
        tri = self
        ntri = tri.random41()
        ntri = ntri.random32()
        if ntri != tri:
            return ntri.boundaryFastSimplify()
        else:
            return tri
            
    def reduceWithDepth(self, depth=2, size=None):
        tri = self
        if depth == 0:
            return tri.fastSimplify()
        if size is None:
            size = tri.numPentachora()
        moves = self.validMoves(2,4)
    
        for move in moves:
            ntri = tri.copy().pachner(move[0], move[1])
            stri = ntri.reduceWithDepth(depth-1, size)
            if stri.numPentachora() < size:
                return stri
        return tri

    def boundaryReduceWithDepth(self, depth=2, size=None):
        tri = self
        if depth == 0:
            return tri.boundaryFastSimplify()
        if size is None:
            size = tri.numBoundaryTet()
        moves = self.boundaryValidMoves(2,3)
    
        for move in moves:
            ntri = tri.boundaryPachner(move[0], move[1])
            stri = ntri.boundaryReduceWithDepth(depth-1, size)
            if stri.numBoundaryTet() < size:
                return stri
        return tri
    
    def deepReduce(self):
        tri = self
        moves1 = tri.validMoves(2,3)
        size = len(tri.boundaryTetrahedra())
        for move in moves1:
            ntri = tri.boundaryPachner(move[0], move[1])
            moves2 = validMoves(ntri,2,3)
            for move2 in moves2:
                nntri = fastSimplify(ntri.boundaryPachner(move2[0], move2[1]))
                #nntri = random41(nntri)
                if len(nntri.boundaryTetrahedra()) < size:
                    return nntri
        return tri
        
    def simplifyWithDepth(self, depth=2, size=None):
        tri = self
        if depth == 0:
            return tri.fastSimplify()
        if size is None:
            size = tri.numPentachora()
        moves = tri.validMoves(3,3)

        c_tri = tri.copy()
        for move in moves:
            ntri = c_tri.pachner(move[0], move[1])
            if ntri != False:
                stri = ntri.simplifyWithDepth(depth-1, size)
                c_tri = tri.copy()
                if stri.numPentachora() < size:
                    return stri
        return tri
            
    def simplify(self, allow51=True):
        tri = self
        moves = tri.validMoves(3,3)
        #print 'valid 3-3 moves:', len(moves)
        c_tri = tri.copy()
        for move in moves:
            ntri = c_tri.pachner(move[0], move[1])
            if ntri == False:
                continue
            c_tri = tri.copy()
            ntri = ntri.fastSimplify(allow51)
            if ntri.numPentachora() < tri.numPentachora():
                tri.pachner(move[0], move[1])
                tri.fastSimplify()
                #print 'simplified'
                return tri
        return tri
        
    def simplify2(self):
        tri = self
        moves = tri.validMoves(2,4)
        print 'valid 2-4 moves:', len(moves)
        c_tri = tri.copy()
        for move in moves:
            ntri = c_tri.pachner(move[0], move[1])
            if ntri == False:
                continue
            c_tri = tri.copy()
            ntri = ntri.fastSimplify()
            if ntri.numPentachora() < tri.numPentachora():
                tri.pachner(move[0], move[1])
                tri.fastSimplify()
                print 'simplified'
                return tri
        return tri
        
    def boundaryReduce(self):
        tri = self
        moves = tri.boundaryValidMoves(2,3)
        size = len(tri.boundaryTetrahedra())
        best_tri = tri
        best_size = size
        for move in moves:
            ntri = tri.boundaryPachner(move[0], move[1]).boundaryFastSimplify()
            if len(ntri.boundaryTetrahedra()) < best_size:
                best_size = len(ntri.boundaryTetrahedra())
                best_tri = ntri
                return ntri
        return best_tri
        
    @staticmethod
    def signDifference(p1, p2):
        # signDifference([2,1,0], [0,1,2]) returns -1
        # since an odd number of transpositions required..
        if len(p1) == 0:
            return 1
        i = p2.index(p1[0])
        return (-1)**i * FourTriangulation.signDifference(p1[1:],p2[:i]+p2[i+1:])
        
    @staticmethod
    def boundarySimplices(simp):
        # simp = [8, [1,3,2]]
        # 8 is pentachoron number
        p_no, s = simp
        lst = []
        for i in range(len(s)):
            bd = [p_no, s[:i]+s[i+1:]]
            if i % 2 == 1 and len(bd[1]) > 1:
                bd[1][0], bd[1][1] = bd[1][1], bd[1][0]
            lst.append(bd)
        return lst
        
    def H2(self):
        tets = self.cycleByDim(3)
        faces = self.cycleByDim(2)
        edges = self.cycleByDim(1)
        # compute boundary map from C_2 to C_1
        bd2 = []
        
        def search(lst, elt):
            for i in range(len(lst)):
                cycle = lst[i]
                for c in cycle:
                    if c[0] == elt[0] and set(elt[1]) == set(c[1]):
                        return (i, self.signDifference(elt[1], c[1]))
        
        for j in range(len(faces)):
            f_cycle = faces[j]
            bd2.append([0]*len(edges))
            f = f_cycle[0]
            for e in self.boundarySimplices(f):
                i, sign = search(edges, e)
                bd2[j][i] += sign
        # bd2 is the transpose of what it should be
        
        bd3 = []
        
        for j in range(len(tets)):
            tet_cycle = tets[j]
            bd3.append([0]*len(faces))
            tet = tet_cycle[0]
            for f in self.boundarySimplices(tet):
                i, sign = search(faces, f)
                bd3[j][i] += sign
        # bd3 is the transpose of what it should be
        
        #########################################
        """
        We know B_2 B_3 = 0 (chain complex property)
        We change the basis of C_2
        B_2 -> B_2*A =: D_2
        B_3 -> A^-1*B_3 =: D_3
        now the kernel of B_2*A corresponds to the 0-columns of the matrix
        and can be read off easily.
        """
        #########################################
        
        bd2_mat = Matrix(bd2)
        (R, U) = bd2_mat.RRDecomposition()
        D_2 = R.transpose()
        D_3 = U.transpose().inverse().multiply(Matrix(bd3).transpose())
        E_3, P, Q = D_3.smithNormalForm()
        # E_3 is the smith normal form of D_3
        # E_3 = P * D_3 * Q
        # Let E_2 = D_2 * P^-1
        # Read off the non-zero columns of E_2
        # let d1,...,dk be the pivots of the corresponding rows of E_3 
        # then H_2 is Z/d1 + ... + Z/dk
        E_2 = D_2.multiply(P.inverse())
        E_2_trans = E_2.transpose()
        inv_factors = []
        for i in range(len(E_2_trans.mat)):
            if E_2_trans.mat[i] == [0]*len(E_2_trans.mat[0]):
                factor = E_3.mat[i][i]
                if factor not in [-1,1]:
                    inv_factors.append(E_3.mat[i][i])
        print 'inv_factors:', inv_factors
        return inv_factors
        
    def intersectionForm(self):
        tets = self.cycleByDim(3)
        faces = self.cycleByDim(2)
        edges = self.cycleByDim(1)
        # compute boundary map from C_2 to C_1
        bd2 = []
        
        def search(lst, elt):
            for i in range(len(lst)):
                cycle = lst[i]
                for c in cycle:
                    if c[0] == elt[0] and set(elt[1]) == set(c[1]):
                        return (i, self.signDifference(elt[1], c[1]))
        
        for j in range(len(faces)):
            f_cycle = faces[j]
            bd2.append([0]*len(edges))
            f = f_cycle[0]
            for e in self.boundarySimplices(f):
                i, sign = search(edges, e)
                bd2[j][i] += sign
        # bd2 is the transpose of what it should be
        #print bd2
        
        bd3 = []
        
        for j in range(len(tets)):
            tet_cycle = tets[j]
            bd3.append([0]*len(faces))
            tet = tet_cycle[0]
            for f in self.boundarySimplices(tet):
                i, sign = search(faces, f)
                bd3[j][i] += sign
        # bd3 is the transpose of what it should be
        #print bd3
        
        #########################################
        """
        We know B_2 B_3 = 0 (chain complex property)
        We change the basis of C_2
        B_2 -> B_2*A =: D_2
        B_3 -> A^-1*B_3 =: D_3
        now the kernel of B_2*A corresponds to the 0-columns of the matrix
        and can be read off easily.
        """
        #########################################
        
        bd2_mat = Matrix(bd2)
        
        (R, U) = bd2_mat.RRDecomposition()
        D_2 = R.transpose()
        U_trans = U.transpose()
        U_trans_inv = U_trans.inverse()
        D_3 = U_trans_inv.multiply(Matrix(bd3).transpose())
        # U.transpose() here is A in the above discussion
        E_3, P, Q = D_3.smithNormalForm()
        # E_3 is the smith normal form of D_3
        # E_3 = P * D_3 * Q
        # Let E_2 = D_2 * P^-1
        # So E_2 = B_2*A*P^-1
        # Read off the non-zero columns of E_2
        # let d1,...,dk be the pivots of the corresponding rows of E_3 
        # then H_2 is Z/d1 + ... + Z/dk
        P_inv = P.inverse()
        E_2 = D_2.multiply(P_inv)

        E_2_trans = E_2.transpose()

        inv_factors = []
        for i in range(len(E_2_trans.mat)):
            if E_2_trans.mat[i] == [0]*len(E_2_trans.mat[0]):
                factor = E_3.mat[i][i]
                if factor not in [-1,1]:
                    inv_factors.append(E_3.mat[i][i])
        print 'inv_factors:', inv_factors
        
        # has no torsion free part
        if len([x for x in inv_factors if x == 0]) == 0:
            return []
        
        self.orient()
        cycles = [self.cycleByDim(i) for i in range(5)]
        
        def search(lst, elt):
            for i in range(len(lst)):
                cycle = lst[i]
                for c in cycle:
                    if c[0] == elt[0] and set(elt[1]) == set(c[1]):
                        return (i, self.signDifference(elt[1], c[1]))
        
        def mappedTo(p_no, simp):
            # the centroid of simp (of p_no) is mapped
            # to which vertex of p_no?
            cycle = cycles[len(simp)]
            i, sign = search(cycles, [p_no, simp])
            matching_elt = cycle[i]
            return matching_elt[1][0] # picks out first vertex
            
        def dualBitMappedTo(p_no, face):
            # the dual bit of face inside p_no is mapped
            # to what 2-chain in the original complex?
            
            # face = [3,1,2] then c_face = [3,1,2,0,4]
            c_face = face[:]
            for i in range(5):
                if i not in c_face:
                    c_face.append(i)
            if self.signDifference(c_face, [0,1,2,3,4])*self.getPentachoronByNum(p_no).orientation == -1:
                c_face[-1], c_face[-2] = c_face[-2], c_face[-1]
                
            tet1 = face + [c_face[-2]]
            tet2 = face + [c_face[-1]]
            
            def getVert(p_no, simp):
                cycle = cycles[len(simp)-1]
                for i in range(len(cycle)):
                    equiv_class = cycle[i]
                    for j in range(len(equiv_class)):
                        no, simp2 = equiv_class[j]
                        if no == p_no and set(simp2) == set(simp):
                            return simp2[0]
            
            a = getVert(p_no, face)
            b = getVert(p_no, tet1)
            c = getVert(p_no, range(5))
            d = getVert(p_no, tet2)
            # dual face gets homotoped to [a,b,c], [c,d,a]
            
            dual1 = [a,b,c]
            dual2 = [c,d,a]
            return [[dual1, dual2], [tet1, tet2]]
                
        mat = []
        dualBd2 = []
        for cyc in cycles[2]:
            mat.append([0]*len(cycles[2]))
            dualBd2.append([0]*len(cycles[3]))
            for p_no, simp in cyc:
                fs, bd_lines = dualBitMappedTo(p_no, simp)
                for f in fs:
                    if len(set(f)) == 3:
                        i, sign = search(cycles[2], [p_no, f])
                        mat[-1][i] += sign
                tet1, tet2 = bd_lines
                i = search(cycles[3], [p_no, tet1])[0]
                first_tet = cycles[3][i][0]

                if first_tet[0] == p_no and set(first_tet[1]) == set(tet1):
                    dualBd2[-1][i] += -1
                    
                i = search(cycles[3], [p_no, tet2])[0]
                first_tet = cycles[3][i][0]
                if first_tet[0] == p_no and set(first_tet[1]) == set(tet2):
                    dualBd2[-1][i] += 1

        # dualBd2 is the transpose of what it should be
        # change basis of domain of dualBd2 so it's easy to read off kernel
        dualBd2Reduced, basisChange = Matrix(dualBd2).RRDecomposition()
        # dualBd2Reduced = basisChange * dualBd2
        
        # we want cycles in the dual complex
        # they are spanned by (in the original basis) the rows of basisChange
        # corresponding to the 0 rows of dualBd2Reduced
        cycleBasis = []
        for i in range(len(dualBd2Reduced.mat)):
            if dualBd2Reduced.mat[i] == [0]*len(dualBd2Reduced.mat[i]):
                cycleBasis.append(basisChange.mat[i])
        
        # image of the cycles spanned by columns of imageOfCycles
        imageOfCycles = Matrix(cycleBasis).multiply(Matrix(mat)).transpose()
        
        
        # choose a nice basis for the image
        G = P.multiply(U_trans_inv).multiply(imageOfCycles)
        H = []
        
        orig = U_trans.multiply(P_inv).transpose()
        Z = []
        
        for i in range(len(E_2_trans.mat)):
            if E_2_trans.mat[i] == [0]*len(E_2_trans.mat[0]):
                factor = E_3.mat[i][i]
                if factor not in [-1,1]:
                    H.append(G.mat[i])
                    Z.append(orig.mat[i])

        J, K = Matrix(H).transpose().RRDecomposition()
        K = K.transpose() # column operations

        # first two columns of L map to the basis cycles of H_2
        L = Matrix(cycleBasis).transpose().multiply(K)
        L_trans = L.transpose()
        intersection_form = Matrix(L_trans.mat[:len(H)]).multiply(Matrix(Z).transpose())
        return intersection_form
        
        
    def orient(self, pent=None):
        if pent is None:
            for p in self.pentachora:
                p.orientation = None
            if len(self.pentachora) > 0:
                pent = self.pentachora[0]
            else:
                return
        if pent.orientation is None:
            pent.orientation = 1
            
        for i in range(5):
            nbr = pent.getNeighbour(i)
            g = pent.getGluing(i)
            if nbr.orientation is None:
                nbr.orientation = -1*self.signDifference(g, [0,1,2,3,4])*pent.orientation
                self.orient(nbr)


    def randomize(self):
        """ Randomizes the triangulation (mutates triangulation). """
        tri = self
        moves = tri.validMoves(3,3)
        if len(moves) == 0:
            return tri
        c = random.randint(0,len(moves)-1)
        moves = moves[c:] + moves[:c]
        used = []
        for m in moves:
            if m[0][0] not in used and m[0][1] not in used and m[0][2] not in used:
                if random.randint(0,1) == 0:
                    tri = tri.pachner(m[0], m[1])
                    used.extend(m[0])
        return tri
                
    def simplifyForTime(self, time_in_millis):
        """
        Attempts to simplify the triangulation for a time period of `time_in_millis`
        milliseconds. Returns a new triangulation (leaves current triangulation unmodified).
        """
        """
        import time
        init_time = int(round(time.time() * 1000))
        
        tri = self.copy()
        curr_time = int(round(time.time() * 1000))
        print "num pent: ", tri.numPentachora()
        while curr_time - init_time < time_in_millis:
            for i in range(30):
                tri = tri.randomize()
                tri = tri.fastSimplify()
                print "num pent: ", tri.numPentachora()
                
            for i in range(10):
                tri = tri.reduceWithDepth(1)
                for j in range(30):
                  tri = tri.randomize()
                print "depth 1, num pent: ", tri.numPentachora()
                
            
            for i in range(2):
                tri = tri.reduceWithDepth(2)
                curr_time = int(round(time.time() * 1000))
                if curr_time - init_time >= time_in_millis:
                    return tri
                print "depth 2, num pent: ", tri.numPentachora()
            
            for i in range(2):
                tri = tri.reduceWithDepth(3)
                curr_time = int(round(time.time() * 1000))
                if curr_time - init_time >= time_in_millis:
                    return tri
                print "num pent: ", tri.numPentachora()
            
            for i in range(20):
                curr_time = int(round(time.time() * 1000))
                if curr_time - init_time >= time_in_millis:
                    return tri
                for j in range(3):
                    tri = tri.anyMove(2, 4)
                    for k in range(5):
                        tri = tri.randomize()
                for j in range(3):
                    tri = tri.reduceWithDepth(1)
                    print "num pent: ", tri.numPentachora()
            curr_time = int(round(time.time() * 1000))
        return tri
        """
        import simplify
        return simplify.simplifyForTime(self, time_in_millis, allow51=False)
    
    def simplifyBoundaryUntil(self, goal_num_bdry_tet, max_num_pent=None, max_rounds=100):
        """
        Layers 4-simplices along boundary to try to simplify the boundary 3-manifold
        triangulation so that it has at most `goal_num_bdry_tet` tetrahedra.
        The method aims to keep the number of pentachora in the triangulation
        as at most `max_num_pent`.
        Attempts at most `max_rounds` tries starting with current triangulation.
        Does not modify triangulation. Returns new triangulation.
        """
        if max_num_pent is None:
            max_num_pent = self.numPentachora()*2 + self.numBoundaryTet()
            
        from s2xd2 import *
        tri = self.copy()
        print "num tet in bdry: ", tri.numBoundaryTet()
        for rounds in range(max_rounds):
            while tri.numPentachora() < max_num_pent:
                tri = tri.boundaryFastSimplify()
                while True:
                    if tri.numBoundaryTet() <= goal_num_bdry_tet:
                        return tri
                    n = tri.numBoundaryTet()
                    tri = reduce(tri)
                    print "reduce, num tet in bdry: ", tri.numBoundaryTet()
                    if tri.numBoundaryTet() >= n:
                        break
                    
                while True:
                    if tri.numBoundaryTet() <= goal_num_bdry_tet:
                        return tri
                    n = tri.numBoundaryTet()
                    tri = deepReduce(tri)
                    print "deepReduce, num tet in bdry: ", tri.numBoundaryTet()
                    if tri.numBoundaryTet() >= n:
                        break

                """
                while tri.numBoundaryTet() <= 6:
                    if tri.numBoundaryTet() <= goal_num_bdry_tet:
                        return tri
                    n = tri.numBoundaryTet()
                    tri = reduceWithDepth(tri, 3)
                    print "depth 3, num tet in bdry: ", tri.numBoundaryTet()
                    if tri.numBoundaryTet() >= n:
                        break
                """
                
                if tri.numBoundaryTet() <= goal_num_bdry_tet:
                    return tri
                tri = boundaryRandomize(tri)
                """
                for i in range(20):
                    for j in range(10):
                        tri = random23(tri)
                    for j in range(10):
                        tri = reduce(tri)
                """
                print "num tet in bdry: ", tri.numBoundaryTet()
                if tri.numBoundaryTet() <= goal_num_bdry_tet:
                    return tri
            tri = self.copy()
        return tri

    def isoSig(self):
        self.renumberTriangulation()

        min_encoding = None
        for p_no in range(self.numPentachora()):
            for perm_ind in range(120):
                enc = self.isoSigFrom(p_no, Perm5.fromIndex(perm_ind))
                if min_encoding is None or enc < min_encoding:
                    min_encoding = enc
        return min_encoding

    def isoSigFrom(self, init_pent_no, perm_to_relab):
        """
        Returns isomorphism signature of a particular canonical relabelling.
        A canonical relabelling is determined by a single "initial" pentachora
        (self.pentachora[init_pent_no], and a permutation `perm_to_relab` which
        maps a vertex of self.pentachora[init_pent_no] to its relabelling.

        NOTE: Requires that self.renumberTriangulation() has been called.
        Port of Regina code isosig-impl.h (isoSigFrom function).
        """
        # renumbers triangulation (modifies)
        # based on Regina code
        
        nSimp = self.numPentachora()
        nBoundaryTet = len(self.boundaryTetrahedra())

        # nFaces is number of facets we record for type-sequence.
        # For each pair of facets glued only the first encounter is recorded.
        # Boundary facets are all recorded.
        # Options: 0 - boundary
        #          1 - joined to a simplex not yet seen
        #          2 - joined to a simplex already seen (or joined to simplex 0).
        nFacets = (5*nSimp + nBoundaryTet) / 2

        # Full destination sequence
        joinDest = [None]*nFacets
        # Permutation sequence (with permutation as integer index)
        joinGluing = [None]*nFacets
        # joinPos is index into joinGluing and joinDest used in main loop
        joinPos = 0
        
        # Consider a map f : self -> a_canonical_labelling (image is f, preImage is f^-1)
        # image[i] = j means that pentachora[i] is the j-th pentachora in relabelling
        # preImage is the inverse map
        # vertexMap is the way f maps vertices
        image = [-1]*nSimp
        preImage = [-1]*nSimp
        image[init_pent_no] = 0
        preImage[0] = init_pent_no

        vertexMap = [Perm5([0,1,2,3,4])]*nSimp
        vertexMap[init_pent_no] = perm_to_relab
        nextUnusedSimp = 1

        # Type sequence
        facetAction = [0]*nFacets
        # facetPos is index into facetAction used in main loop
        facetPos = 0
        
        for simpImg in range(nSimp):
            if preImage[simpImg] < 0:
                break
            simpSrc = preImage[simpImg]
            s = self.pentachora[simpSrc]
            
            for facetImg in range(5):
                facetSrc = vertexMap[simpSrc].preImageOf(facetImg)
                # is it a boundary facet?
                if s.neighbour[facetSrc] is None:
                    facetAction[facetPos] = 0 # boundary facet
                    facetPos += 1
                    continue
                
                dest = s.neighbour[facetSrc].n
                # Have we already glued this facet from the other side?
                if image[dest] >= 0:
                    if (image[dest] < image[simpSrc] or
                        (dest == simpSrc and
                         vertexMap[simpSrc].perm[s.gluing[facetSrc][facetSrc]]
                         < vertexMap[simpSrc].perm[facetSrc])):
                        # Yes, skip gluing.
                        continue

                if image[dest] < 0: # new simplex
                    image[dest] = nextUnusedSimp
                    nextUnusedSimp += 1
                    preImage[image[dest]] = dest
                    vertexMap[dest] = vertexMap[simpSrc].compose(
                        Perm5(s.gluing[facetSrc]).inverse())
                    facetAction[facetPos] = 1 # completely new simplex
                    facetPos += 1
                    continue

                # Simplex that we've seen before, record gluing.
                joinDest[joinPos] = image[dest]                
                joinGluing[joinPos] = (vertexMap[dest].compose(
                    Perm5(s.gluing[facetSrc]).compose(vertexMap[simpSrc].inverse()))).index()
                joinPos += 1
                facetAction[facetPos] = 2 # old gluing
                facetPos += 1

        # encode information into string
        ans = ""
        if nSimp < 63:
            nChars = 1
        else:
            nChars = 0
            tmp = nSimp
            while tmp > 0:
                tmp >>= 6
                nChars += 1
            ans = IsoSigHelper.sChar(63)
            ans += IsoSigHelper.sChar(nChars)

        ans += IsoSigHelper.toStr(nSimp, nChars)
        i = 0
        while i < facetPos:
            if facetPos >= i + 3:
                nTrits = 3
            else:
                nTrits = facetPos - i
            ans += IsoSigHelper.toTritsStr(facetAction, i, nTrits)
            i += 3
        #print "After trits:", ans
        for i in range(joinPos):
            ans += IsoSigHelper.toStr(joinDest[i], nChars)
        for i in range(joinPos):
            ans += IsoSigHelper.toStr(joinGluing[i], 2) # need 2 chars for 120 permutations

        return ans
    
class IsoSigHelper:
    @staticmethod
    def sChar(c):
        # returns character representing integer value
        if c < 26:
            return chr(c + ord('a'))
        if c < 52:
            return chr(c - 26 + ord('A'))
        if c < 62:
            return chr(c - 52 + ord('0'))
        if c == 62:
            return '+'
        return '-'

    @staticmethod
    def toStr(val, nChars):
        s = ""
        while nChars > 0:
            s += IsoSigHelper.sChar(val & 0x3F)
            val = val >> 6
            nChars -= 1
        return s

    @staticmethod
    def toTritsStr(trits, i_start, nTrits):
        ans = 0
        if nTrits >= 1:
            ans |= trits[i_start]
        if nTrits >= 2:
            ans |= (trits[i_start + 1] << 2)
        if nTrits >= 3:
            ans |= (trits[i_start + 2] << 4)
        return IsoSigHelper.sChar(ans)
    
class Perm5:
    index_cache = {}
    def __init__(self, perm):
        self.perm = perm

    def preImageOf(self, i):
        for j in range(5):
            if self.perm[j] == i:
                return j

    def inverse(self):
        inv = [0]*5
        for i in range(5):
            inv[self.perm[i]] = i
        return Perm5(inv)
    
    def compose(self, g):
        return Perm5([self.perm[g.perm[i]] for i in range(5)])

    
    def computeIndex(self):
        m = 24 # (5-1)!
        ind = 0
        for i in range(len(self.perm)):
            c = 0
            for j in range(i+1, len(self.perm)):
                if self.perm[j] < self.perm[i]:
                    c = c + 1
            ind = ind + m*c
            m = m / max(1, (4-i))
        return ind

    def index(self):
        if tuple(self.perm) in Perm5.index_cache:
            return Perm5.index_cache[tuple(self.perm)]
        else:
            Perm5.index_cache[tuple(self.perm)] = self.computeIndex()
            return Perm5.index_cache[tuple(self.perm)]
        
    @staticmethod
    def fromIndex(i):
        assert i >= 0 and i < 120, "Invalid index of permutation " + str(i)
        perm = []
        labels = range(5)
        m = 24
        d = 4
        while labels:
            perm.append(labels[i / m])
            del labels[i / m]
            i %= m
            m /= max(d, 1)
            d -= 1
        return Perm5(perm)

    def __str__(self):
        return str(self.perm)
        
class Pentachoron:
    # Pentachoron neighbour[0-4] 
    #   where neighbour[i] reference to Pentachoron glued to tetrahedron opposite vertex i.
    # Gluings gluing[0-4][0-4]
    #   where gluing[i][j] pairs face i and vertex j to pentachoron neighbour[i]
    #   vertex gluing[i][j].
    
    def __init__(self, n=0):
        self.neighbour = [None]*5
        self.gluing = [[]]*5
        self.cusp_no = [-1]*5
        self.n = n
        self.orientation = None
        
    def joinTo(self, myFace, you, gluing, tri=None):
        # tri is the FourTriangulation self is part of, this is used to update cusp numbers
        self.setNeighbour(myFace, you)
        self.setGluing(myFace, gluing)
        
        you.setNeighbour(gluing[myFace], self)
        
        you_gluing = [0]*5
        for i in range(5):
            you_gluing[gluing[i]] = i
            
        you.setGluing(gluing[myFace], you_gluing)
        
        #if tri is not None:
        #    for i in range(5):
        #        if i != myFace:
        #            tri.identifyCusps(self.cusp_no[i], you.cusp_no[gluing[i]])
        
    def setNeighbour(self, myFace, you):
        self.neighbour[myFace] = you
        
    def setGluings(self, gluings):
        self.gluing = [x[:] for x in gluings]
        
    def setGluing(self, myFace, gluing):
        self.gluing[myFace] = gluing
        
    def setCuspNo(self, vert, n):
        self.cusp_no[vert] = n
        
    def setCuspNos(self, cusps):
        self.cusp_no = cusps[:]
        
    def getNeighbour(self, myFace):
        return self.neighbour[myFace]
        
    def getGluing(self, myFace):
        return self.gluing[myFace][:]
        
    def getGluings(self):
        return [x[:] for x in self.gluing]
        
    def getCuspNos(self):
        return self.cusp_no[:]
        
    def getCuspNo(self, vert):
        return self.cusp_no[vert]
        
    def getPentNo(self):
        return self.n
        
def oppositeTetrahedron(m):
    return filter(lambda x: x != m, range(5))
    
def oppositeTriangle(m):
    return filter(lambda x: x != m, range(4))
                
def enum_bij(lst):
    bijs = []
    
    if len(lst) == 0:
        return [[]]
        
    for i in range(len(lst)):
        for bij in enum_bij(lst[:i] + lst[(i+1):]):
            bijs.append([lst[i]] + bij)
    return bijs
    
def perm_sign(perm):
    sgn = 1
    for i in range(len(perm)):
        for j in range(i):
            if perm[i]-perm[j] < 0:
                sgn *= -1
    return sgn

v_perms = {}
def valid_perms(m, n):
    # face m joins to face n
    if (m,n) in v_perms:
        return v_perms[(m,n)]
    bijs = enum_bij(oppositeTetrahedron(n))
    valid_bijs = []
    for bij in bijs:
        perm = bij[:m] + [n] + bij[m:]
        if perm_sign(perm) == -1:
            valid_bijs.append(perm)
    v_perms[(m,n)] = valid_bijs
    return valid_bijs
    
def validFaceCycle(lst):
        # a face cycle is valid if the same face doesn't get identified to itself
        # in more than one way
    if len(set([tuple(x[:-1]) for x in lst])) == len(lst):
        return True
    else:
        return False

count = 0
def enumerate(tri, free_faces):
    cusps = tri.getCuspNums()
    if len(cusps) < 3:
        return
        
    if not free_faces:
        global count
        if len(cusps) == 3:
            count += 1
            if count % 1000 == 0:
                print count
                #print cusps
                sys.stdout.flush()
        # tri is done! what do we want to do with it?
        #print tri.triFaceCycle([[0, 0, 1, [2,3,4]]])

        if len(cusps) == 3:
            valid = True
            for p_no in [0,1]:
                for tet_no in range(5):
                    for tri_no in oppositeTetrahedron(tet_no):
                        triangle = filter(lambda x: x != tri_no, oppositeTetrahedron(tet_no))
                        if not validFaceCycle(tri.triFaceCycle([[p_no, tet_no, tri_no, triangle]])):
                            valid = False
                        
            if valid:
                #tri.link().writeSnapPea('C:/pr2011/4dim_triangulation/lnk' + str(count) + '.tri')
                dir = '/mount/autofs/home_stude/aiissa/4dim_triangulation/'
                tri.link().writeSnapPea(dir + 'lnk' + str(count) + '.tri')
                reg_tri = regina.readSnapPea(dir + 'lnk' + str(count) + '.tri')
                comps = reg_tri.splitIntoComponents()
                correct_link = True
                if comps != 3:
                    correct_link = False
                for i in range(comps):
                    comp_tri = reg_tri.nextTreePacket()
                    if not comp_tri.isThreeSphere():
                        correct_link = False
                if not correct_link:
                    os.remove(dir + 'lnk' + str(count) + '.tri')
                else:
                    tri.writeToFile(dir + 'mfld' + str(count) + '.tri')
                #sys.exit(0)
                pass
            else:
                count -= 1
        return
        
    c_free_faces = copy.deepcopy(free_faces)
    # pair the first free face with all others
    first_face = c_free_faces[0]
    for i in range(1, len(c_free_faces)):
        second_face = c_free_faces[i]
        for perm in valid_perms(first_face[1], second_face[1]):
            c_tri = tri.copy()
            pent1 = c_tri.getPentachoron(first_face[0])
            pent2 = c_tri.getPentachoron(second_face[0])
            pent1.joinTo(first_face[1], pent2, perm, c_tri)
            enumerate(c_tri, free_faces[1:i] + free_faces[i+1:])

def find_pillow():
    tri = FourTriangulation(2, tri_title="Potential pillow with tunnel")
    free_faces = [[0,i] for i in range(1,5)] + [[1,i] for i in range(1,5)]
    tri.getPentachoron(0).joinTo(0, tri.getPentachoron(1), [0,1,2,4,3])
    enumerate(tri, free_faces)
    print count
    
