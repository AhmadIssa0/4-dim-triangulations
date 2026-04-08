
## deals with formatting a triangulation for snappea.
import sys
import copy

"""
In SnapPy a positively oriented tetrahedra is one where the edges (30, 31, 32)
form a right handed coordinate system. The tetrahedra shape parameter in SnapPy
corresponds to the edge 23 (and 01).
23 <-> z
02 <-> 1/(1-z)
12 <-> (z-1)/z
"""

class Triangulation:
    # tetrahedra = list of Tetrahedron
    # cusps = number of cusps
    
    def __init__(self, n=0, tri_title=None):
        # tri_title is the title of the triangulation as written in the triangulation file
        self.tetrahedra = []
        self.cusps = 0
        for i in range(n):
            self.addTetrahedron()
        self.tri_title = tri_title
            
    def addTetrahedron(self, num=None):
        if num == None:
            num = self.nextTetIndex()
        tet = Tetrahedron(num)
        self.tetrahedra.append(tet)
        return tet
        
    def removeTetrahedron(self, num):
        # Removes Tetrahedron with getTetNo() == num
        for i in range(self.numTetrahedra()):
            if self.tetrahedra[i].getTetNo() == num:
                del self.tetrahedra[i]
                break
                
    def renumberTriangulation(self):
        for i in range(self.numTetrahedra()):
            self.tetrahedra[i].n = i
                
    def renumberTetrahedron(self, old, new):
        # Since tet.getNeighbour(i) returns a reference to a tetrahedron,
        # rather than a tetrahedron number, this is simple.
        self.getTetrahedronByNum(old).n = new
        
        
    def threeToTwo(self, tet_nos, t1_faces):
        # Assumes that the three tetrahedra do not externally glue to each other.
        # Sets all the new cusp nos to 0. This may not work well if there are multiple cusps.
        # t1_faces 
        nTri = self.copy()
        
        tet_nos[1], face = nTri.gluedFace(tet_nos[0], [t1_faces[0][0], t1_faces[1][0], t1_faces[1][1]])
        opVert = nTri.oppositeVertex(face)
        t2_faces = [[face[0], face[2], opVert], [face[1], opVert, face[2]]]
        print 't2_faces', t2_faces
        
        tet_nos[2], face = nTri.gluedFace(tet_nos[1], [t2_faces[0][0], t2_faces[1][0], t2_faces[1][1]])
        opVert = nTri.oppositeVertex(face)
        t3_faces = [[face[0], face[2], opVert], [face[1], opVert, face[2]]]
        print 't3_faces', t3_faces
        
        tet1 = nTri.getTetrahedronByNum(tet_nos[0])
        tet2 = nTri.getTetrahedronByNum(tet_nos[1])
        tet3 = nTri.getTetrahedronByNum(tet_nos[2])
        
        n1 = nTri.nextTetIndex()
        ntet1 = nTri.addTetrahedron(n1)
        
        n2 = nTri.nextTetIndex()
        ntet2 = nTri.addTetrahedron(n2)  
        
        nTri.switchFace(face)
        ntet1.joinTo(0, ntet2, [0,1,3,2])
        ntet2.joinTo(0, ntet1, [0,1,3,2])
        # def switchFace(self, tri):
        # def gluedFace(self, tet_no, tri):
        
        ntet1.setCuspNos([0]*4)
        ntet2.setCuspNos([0]*4)
        
        print 'z_new_tet1(01) = ', 'z_' + str(tet_nos[1]), str(t2_faces[0][0]) + str(t2_faces[0][2]), '* z_' + str(tet_nos[2]), str(t3_faces[0][0]) + str(t3_faces[0][1])                
        print 'z_new_tet2(01) = ', 'z_' + str(tet_nos[1]), str(t2_faces[1][0]) + str(t2_faces[1][1]), '* z_' + str(tet_nos[2]), str(t3_faces[1][0]) + str(t3_faces[1][2])                  
        opVert = nTri.oppositeVertex(t1_faces[0])
        g = tet1.getGluing(opVert)
        print 'neighbour:', tet1.getNeighbour(opVert).n
        ntet1.joinTo(1, tet1.getNeighbour(opVert),
                           [g[t1_faces[0][0]],
                           g[opVert],
                           g[t1_faces[0][1]],
                           g[t1_faces[0][2]]])
                           
        opVert = nTri.oppositeVertex(t2_faces[0])
        g = tet2.getGluing(opVert)
        print 'neighbour:', tet2.getNeighbour(opVert).n
        ntet1.joinTo(2, tet2.getNeighbour(opVert),
                           [g[t2_faces[0][0]],
                           g[t2_faces[0][2]],
                           g[opVert],
                           g[t2_faces[0][1]]])        

        opVert = nTri.oppositeVertex(t3_faces[0])
        g = tet3.getGluing(opVert)
        print 'neighbour:', tet3.getNeighbour(opVert).n
        ntet1.joinTo(3, tet3.getNeighbour(opVert),
                           [g[t3_faces[0][0]],
                           g[t3_faces[0][1]],
                           g[t3_faces[0][2]],
                           g[opVert]])
                           
        # tet 2                   
        opVert = nTri.oppositeVertex(t1_faces[1])
        g = tet1.getGluing(opVert)
        print 'neighbour:', tet1.getNeighbour(opVert).n
        ntet2.joinTo(1, tet1.getNeighbour(opVert),
                           [g[t1_faces[1][0]],
                           g[opVert],
                           g[t1_faces[1][1]],
                           g[t1_faces[1][2]]])
                                   
        opVert = nTri.oppositeVertex(t2_faces[1])
        g = tet2.getGluing(opVert)
        print 'neighbour:', tet2.getNeighbour(opVert).n
        ntet2.joinTo(3, tet2.getNeighbour(opVert),
                           [g[t2_faces[1][0]],
                           g[t2_faces[1][1]],
                           g[t2_faces[1][2]],
                           g[opVert]])        

        opVert = nTri.oppositeVertex(t3_faces[1])
        g = tet3.getGluing(opVert)
        print 'neighbour:', tet3.getNeighbour(opVert).n
        ntet2.joinTo(2, tet3.getNeighbour(opVert),
                           [g[t3_faces[1][0]],
                           g[t3_faces[1][2]],
                           g[opVert],
                           g[t3_faces[1][1]]])
               
        print tet_nos
        for i in range(3):
            nTri.removeTetrahedron(tet_nos[i])

        nTri.renumberTriangulation()
        return nTri
        
    def two_to_three(self, tet1_no, tet2_no, vert1, vert2, face1, face2):
        nTri = self.copy()
        
        tet = []
        tet.append(nTri.getTetrahedronByNum(tet1_no))
        tet.append(nTri.getTetrahedronByNum(tet2_no))
        
        ntet = []
        
        n1 = nTri.nextTetIndex()
        ntet.append(nTri.addTetrahedron(n1))
        
        n2 = nTri.nextTetIndex()
        ntet.append(nTri.addTetrahedron(n2))

        n3 = nTri.nextTetIndex()
        ntet.append(nTri.addTetrahedron(n3))
        
        # def joinTo(self, myFace, you, gluing):
        # set cusp numbers
        for i in range(3):
            ntet[i].setCuspNos([tet[1].getCuspNo(vert2), tet[0].getCuspNo(vert1), tet[0].getCuspNo(face1[(i+1)%3]), tet[0].getCuspNo(face1[i%3])])
        
        # 1,2,3,4,5 are the vertices of a 4-dimensional simplex
        two_tets = [[1,2,3,4], [1,2,3,5]]
        two_tet_free = [[0,[1,2,4]], [0,[1,3,4]], [0,[2,3,4]], [1,[1,2,5]], [1,[1,3,5]], [1,[2,3,5]]]
        three_tets = [[2,3,4,5], [1,3,4,5], [1,2,4,5]]
        three_tet_free = [[2,[1,2,4]], [1,[1,3,4]], [0,[2,3,4]], [2,[1,2,5]], [1,[1,3,5]], [0,[2,3,5]]]
        three_tet_internal = [[0,[3,4,5]], [1,[3,4,5]], [0,[2,4,5]], [2,[2,4,5]], [1,[1,4,5]], [2,[1,4,5]]]
        free_faces = [[1,2,4], [1,3,4], [2,3,4], [1,2,5], [1,3,5], [2,3,5]]
        
        # need a map from free_faces to faces of tet1 and tet2
        
        # common verts in two_tets is [1,2,3]
        # common face of two_tets must map to common face of tet1 and tet2
        # then we can extend to the rest in whatever way we want.
        # ASSUME: face1[i] pairs to face2[i] for i=0,1,2.
        # two_tets[0] = [1,2,3,4] <-> face1 + [vert1] of tet1
        # two_tets[1] = [1,2,3,5] <-> face2 + [vert2] of tet2
        two_tet_pairing = [face1 + [vert1], face2 + [vert2]]
        # two_tet_pairing[0][i] is the vert of tet1 that two_tets[i] corresponds to
        free_abs_to_actual_two = [[t_no, [two_tet_pairing[t_no][two_tets[t_no].index(v)] for v in face]] for t_no, face in two_tet_free]
        
        # [ntet1,ntet2,ntet3] correspond to three_tets
        # so ntet1's vertices [0,1,2,3] <-> [2,3,4,5]
        # induces a map from free_faces to faces of ntets
        three_tet_pairing = [[0,1,2,3], [0,1,2,3], [0,1,2,3]]
        free_abs_to_actual_three = [[t_no, [three_tet_pairing[t_no][three_tets[t_no].index(v)] for v in face]] for t_no, face in three_tet_free]
        internal_abs_to_actual_three = [[t_no, [three_tet_pairing[t_no][three_tets[t_no].index(v)] for v in face]] for t_no, face in three_tet_internal]
        
        def glue_faces(t1, t2, f1, f2):
            # t1 tetrahedron face f1 gets glued to f2 of t2
            v1 = Triangulation.oppositeVertex(f1)
            v2 = Triangulation.oppositeVertex(f2)
            g = {}
            g[v1] = v2
            for i in range(len(f1)):
                g[f1[i]] = f2[i]
            t1.joinTo(v1, t2, [g[0],g[1],g[2],g[3]])
            
        def two_tet_to_three_tet(n_two_tet, face):
            for i in range(len(free_abs_to_actual_two)):
                t_no, f = free_abs_to_actual_two[i]
                if t_no == n_two_tet and set(face).issubset(f):
                    new_t_no, new_f = free_abs_to_actual_three[i]
                    dic = {f[0]:new_f[0], f[1]:new_f[1], f[2]:new_f[2]}
                    return new_t_no, [dic[face[0]], dic[face[1]], dic[face[2]]]
        
        for i in range(len(internal_abs_to_actual_three)/2):
            nt1, f1 = internal_abs_to_actual_three[2*i]
            nt2, f2 = internal_abs_to_actual_three[2*i+1]
            glue_faces(ntet[nt1], ntet[nt2], f1, f2)
            
        for i in range(len(free_abs_to_actual_two)):
            nt1, f1 = free_abs_to_actual_two[i]
            v = Triangulation.oppositeVertex(f1)
            nbr = tet[nt1].getNeighbour(v)
            g = tet[nt1].getGluing(v)
            outer_face = [g[j] for j in f1]
            # deal with the case where external face is of an internal tetrahedron
            if nbr == tet[0]:
                nt, nf = two_tet_to_three_tet(0, outer_face)
                nbr = ntet[nt]
                g = nbr.getGluing(Triangulation.oppositeVertex(nf))
                outer_face = nf
            elif nbr == tet[1]:
                nt, nf = two_tet_to_three_tet(1, outer_face)
                nbr = ntet[nt]
                g = nbr.getGluing(Triangulation.oppositeVertex(nf))
                outer_face = nf
            nt2, f2 = free_abs_to_actual_three[i]
            glue_faces(nbr, ntet[nt2], outer_face, f2)
            
        nTri.removeTetrahedron(tet1_no)
        nTri.removeTetrahedron(tet2_no)
        nTri.renumberTriangulation()
        return nTri
            
    def twoToThree(self, tet1_no, tet2_no, vert1, vert2, face1, face2):
        # Assumes that the pair of tetrahedra only glue along one face, i.e. no external gluings.
        # Returns a new triangulation object
        nTri = self.copy()
        
        tet1 = nTri.getTetrahedronByNum(tet1_no)
        tet2 = nTri.getTetrahedronByNum(tet2_no)
        
        n1 = nTri.nextTetIndex()
        ntet1 = nTri.addTetrahedron(n1)
        
        n2 = nTri.nextTetIndex()
        ntet2 = nTri.addTetrahedron(n2)     

        n3 = nTri.nextTetIndex()
        ntet3 = nTri.addTetrahedron(n3)
        
        # def joinTo(self, myFace, you, gluing):
        
        # set cusp numbers
        ntet1.setCuspNos([tet2.getCuspNo(vert2), tet1.getCuspNo(vert1), tet1.getCuspNo(face1[1]), tet1.getCuspNo(face1[0])])
        ntet2.setCuspNos([tet2.getCuspNo(vert2), tet1.getCuspNo(vert1), tet1.getCuspNo(face1[2]), tet1.getCuspNo(face1[1])])
        ntet3.setCuspNos([tet2.getCuspNo(vert2), tet1.getCuspNo(vert1), tet1.getCuspNo(face1[0]), tet1.getCuspNo(face1[2])])
        
        
        print 'z_new_tet1(23) = ', 'z_' + str(tet1_no), str(face1[0]) + str(face1[1]), '* z_' + str(tet2_no), str(face2[0]) + str(face2[1])                
        print 'z_new_tet2(23) = ', 'z_' + str(tet1_no), str(face1[1]) + str(face1[2]), '* z_' + str(tet2_no), str(face2[1]) + str(face2[2])     
        print 'z_new_tet3(23) = ', 'z_' + str(tet1_no), str(face1[2]) + str(face1[0]), '* z_' + str(tet2_no), str(face2[2]) + str(face2[0])
        
        # face pairings of first new tetrahedron
        g = tet1.getGluing(face1[2])
        ntet1.joinTo(0, tet1.getNeighbour(face1[2]), 
                    [g[face1[2]], 
                    g[vert1], 
                    g[face1[1]], 
                    g[face1[0]]]
                    )
                    
        g = tet2.getGluing(face2[2]) 
        ntet1.joinTo(1, tet2.getNeighbour(face2[2]), 
                    [g[vert2], 
                    g[face2[2]], 
                    g[face2[1]], 
                    g[face2[0]]]
                    )
        ntet1.joinTo(2, ntet3, [0,1,3,2])
        ntet1.joinTo(3, ntet2, [0,1,3,2])
        
        
        # face pairings of second new tetrahedron
        g = tet1.getGluing(face1[0])

        ntet2.joinTo(0, tet1.getNeighbour(face1[0]), 
                    [g[face1[0]], 
                    g[vert1], 
                    g[face1[2]], 
                    g[face1[1]]]
                    )
        g = tet2.getGluing(face2[0])     

        ntet2.joinTo(1, tet2.getNeighbour(face2[0]), 
                    [g[vert2], 
                    g[face2[0]], 
                    g[face2[2]], 
                    g[face2[1]]]
                    )
        ntet2.joinTo(2, ntet1, [0,1,3,2])
        ntet2.joinTo(3, ntet3, [0,1,3,2])
        
        # face pairings of third new tetrahedron
        g = tet1.getGluing(face1[1])

        ntet3.joinTo(0, tet1.getNeighbour(face1[1]), 
                    [g[face1[1]], 
                    g[vert1], 
                    g[face1[0]], 
                    g[face1[2]]]
                    )
        g = tet2.getGluing(face2[1])     

        ntet3.joinTo(1, tet2.getNeighbour(face2[1]), 
                    [g[vert2], 
                    g[face2[1]], 
                    g[face2[0]], 
                    g[face2[2]]]
                    )
        ntet3.joinTo(2, ntet2, [0,1,3,2])
        ntet3.joinTo(3, ntet1, [0,1,3,2])

        nTri.removeTetrahedron(tet1_no)
        nTri.removeTetrahedron(tet2_no)
        nTri.renumberTriangulation()
        return nTri

        
    def idealToFinite(self):
        # Idea for this was taken from Regina source code
        # WARNING: Calls self.renumberTriangulation()
        
        def perm4(i,j):
            p = range(4)
            p[i], p[j] = j, i
            return p
            
        self.renumberTriangulation()
        numOldTet = self.numTetrahedra()
        oldTets = []
        for i in range(numOldTet):
            oldTets.append(self.getTetrahedronByNum(i))
            
        tri = Triangulation()
        
        interior = []
        edge = []
        vertex = []

        for k in range(numOldTet):
            interior.append([])
            edge.append([])
            vertex.append([])
            for i in range(4):
                edge[k].append([None]*4)
                vertex[k].append([None]*4)
        
            for i in range(4):
                interior[k].append(tri.addTetrahedron())
                for j in range(4):
                    if i != j:
                        edge[k][i][j] = tri.addTetrahedron()
                        vertex[k][i][j] = tri.addTetrahedron()
                    
        # glue tetrahedra inside the same old tetrahedron together
        for i in range(numOldTet):
            for j in range(4):
                for k in range(4):
                    if j != k:
                        interior[i][j].joinTo(k, vertex[i][k][j], range(4))
                        edge[i][j][k].joinTo(j, edge[i][k][j], perm4(j,k))
                        for m in range(4):
                            if m != j and m != k:
                                edge[i][j][k].joinTo(m, vertex[i][j][m], perm4(k,m))
                                
        # gluings between adjacent big tetrahedra
        for i in range(numOldTet):
            oldTet = oldTets[i]
            for j in range(4):
                if oldTet.getNeighbour(j) is not None:
                    oppTet = oldTet.getNeighbour(j).getTetNo()
                    p = oldTet.getGluing(j)
                    
                    for k in range(4):
                        if j != k:
                            edge[i][j][k].joinTo(k, edge[oppTet][p[j]][p[k]], p)
                            vertex[i][j][k].joinTo(k, vertex[oppTet][p[j]][p[k]], p)
                       
   
    def nextTetIndex(self):
        if len(self.tetrahedra) == 0:
            return 0
        return max(map(lambda tet: tet.getTetNo(), self.tetrahedra))+1
     
    def copy(self):
        # Requires the tetrahedra in self to be numbered 0..(n-1)
        nTri = Triangulation(self.numTetrahedra())
        nTri.setNoCusps(self.getNoCusps())
        
        for i in range(self.numTetrahedra()):
            nTet = nTri.getTetrahedron(i)
            oldTet = self.getTetrahedron(i)
            nTet.setGluings(oldTet.getGluings())
            nTet.setCuspNos(oldTet.getCuspNos())
                
        for i in range(self.numTetrahedra()):
            nTet = nTri.getTetrahedron(i)
            oldTet = self.getTetrahedron(i)
            for j in range(4):
                neighbour_no = oldTet.getNeighbour(j).getTetNo()
                nTet.setNeighbour(j, nTri.getTetrahedronByNum(neighbour_no))
        return nTri
            
    def getTetrahedronByNum(self, n):
        for tet in self.getTetrahedraList():
            if tet.getTetNo() == n:
                return tet
                
    def getTetrahedraList(self):
        return self.tetrahedra
        
    # return i-th Tetrahedron in list of tetrahedra
    # this may be different to Tetrahedron.getTetNo() == i
    def getTetrahedron(self, i):
        return self.tetrahedra[i]
        
    def numTetrahedra(self):
        return len(self.tetrahedra)
        
    def setNoCusps(self, n):
        self.cusps = n
        
    def getNoCusps(self):
        return self.cusps
        
    @staticmethod
    def oppositeVertex(tri):
        for i in range(4):
            if i not in tri:
                return i
        
    def switchFace(self, tri):
        """ if tri = [a,b,c], then want [a,b,d] """
        return tri[:2] + [self.oppositeVertex(tri)]
                
    def gluedFace(self, tet_no, tri):
        tet = self.tetrahedra[tet_no]
        face = self.oppositeVertex(tri)
        g = tet.getGluing(face)
        return (tet.getNeighbour(face).n, [g[tri[0]], g[tri[1]], g[tri[2]]])
    
    def rotateFace(self, tet_no, tri):
        """ tri = [a,b,c], then a-b is the special edge """
        return self.gluedFace(tet_no, self.switchFace(tri))
        
    def writeSnapPea(self, foutname):
        f = open(foutname, 'w')
        f.write('% Triangulation\n')
        if self.tri_title is not None:
            f.write(self.tri_title + '\n')
        else:
            f.write('tt_triangulation\n')
        f.write('not_attempted 0.0\n')
        f.write('unknown_orientability\n')
        f.write('CS_unknown\n')
        f.write('{0} 0\n'.format(self.cusps))
        for i in range(self.cusps):
            f.write('\ttorus\t0.0\t0.0\n')
        f.write(str(self.numTetrahedra()) + '\n')
        
        for i in range(self.numTetrahedra()):
            tet = self.getTetrahedron(i)
            for j in range(4):
                f.write(' ' * max(0, 4 - len(str(tet.getNeighbour(j).getTetNo()))))
                f.write('{0}'.format(tet.getNeighbour(j).getTetNo()))
                if j < 4:
                    f.write(' ')
            f.write('\n')
            #f.write('  {0}   {1}   {2}   {3}\n'.format(tet.getNeighbour(0).getTetNo(),tet.getNeighbour(1).getTetNo(),
                                         #           tet.getNeighbour(2).getTetNo(),tet.getNeighbour(3).getTetNo()))

            for j in range(4):
                g = tet.getGluing(j)
                f.write(' {0}{1}{2}{3}'.format(g[0],g[1],g[2],g[3]))
            f.write('\n')
            f.write('  {0}   {1}   {2}   {3} \n'.format(tet.getCuspNo(0),tet.getCuspNo(1),tet.getCuspNo(2),tet.getCuspNo(3)))
            for k in range(4):
                f.write('  0'*16 + '\n')
            f.write('0.0 0.0\n\n')
            


class Tetrahedron:
    # Tetrahedron neighbour[0-3]
    # Gluings gluing[0-3][4] where gluing[i][j] pairs face i and vertex j to tet neighbour[i] vertex gluing[i][j]
    # Cusp numbers corresponding to vertices cusp_no[0-3]
    
    def __init__(self, n=0):
        self.neighbour = [0]*4
        self.gluing = [[]]*4
        self.cusp_no = [-1]*4
        self.n = n
        
    def joinTo(self, myFace, you, gluing):
        self.setNeighbour(myFace, you)
        self.setGluing(myFace, gluing)
        
        you.setNeighbour(gluing[myFace], self)
        
        you_gluing = [0]*4
        for i in range(4):
            you_gluing[gluing[i]] = i
            
        you.setGluing(gluing[myFace], you_gluing)
        
    def setNeighbour(self, myFace, you):
        self.neighbour[myFace] = you
        
    def setGluings(self, gluings):
        self.gluing = copy.deepcopy(gluings)
        
    def setGluing(self, myFace, gluing):
        self.gluing[myFace] = gluing
        
    def setCuspNo(self, vert, n):
        self.cusp_no[vert] = n
        
    def setCuspNos(self, cusps):
        self.cusp_no = copy.deepcopy(cusps)
        
    def getNeighbour(self, myFace):
        if type(self.neighbour[myFace]) == int:
            print 'returning int, tet:', self.n
        return self.neighbour[myFace]
        
    def getGluing(self, myFace):
        return self.gluing[myFace]
        
    def getGluings(self):
        return copy.deepcopy(self.gluing)
        
    def getCuspNos(self):
        return self.cusp_no[:]
        
    def getCuspNo(self, vert):
        return self.cusp_no[vert]
        
    def getTetNo(self):
        return self.n

def SnappyToEquation(vecs):
    """ takes a snappy gluing equation (form="rect") and writes it out as an actual equation."""
    eqns = []
    for vec in vecs:
        A, B, c = vec
        eqn = ""
        for i in range(len(A)):
            if A[i] != 0:
                eqn += "* z{0}^{1}".format(i,A[i],B[i])
            if B[i] != 0:
                eqn += "* (1 - z{0})^{2}".format(i,A[i],B[i])
        eqn = eqn[2:] # remove leading '* '
        eqn += " = " + str(c)
        eqns.append(eqn)
    return eqns
    
#print SnappyToEquation([([2, -1], [-1, 2], 1), ([-2, 1], [1, -2], 1), ([1, 0], [0, 1], 1), ([0, -2], [0, 4], 1)])

"""
tri = Triangulation(2, 'fig8')
T0 = tri.getTetrahedron(0)
T1 = tri.getTetrahedron(1)
T0.joinTo(0, T1, [0,1,3,2])
T0.joinTo(1, T1, [1,2,3,0])
T0.joinTo(2, T1, [2,3,1,0])
T0.joinTo(3, T1, [2,1,0,3])
ntri = tri.two_to_three(0, 1, 0, 0, [1,2,3], [1,3,2])
ntri.writeSnapPea('C:/pr2011/4dim_triangulation/fig8_pachner.tri')
"""

"""
rect = [([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, -1], 1), ([0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1], [-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0], 1), ([-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1], [1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 1], -1), ([-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, -2, -1], [1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 1, 0, 0, 0, 0, 2, 1], 1), ([0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-1, 1, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 1), ([1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0], 1), ([0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0], 1), ([0, 0, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, -1, 0, 1, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 1), ([0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0], [0, 0, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0], 1), ([0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1, 1, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 1), ([0, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 1), ([0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0], [0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0], 1), ([0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, -1, 1, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 1), ([0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 1), ([0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 1, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 1), ([0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0], 1), ([0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 1), ([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, 0], 1), ([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 1), ([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0], 1), ([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0], 1), ([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, -1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, -1], 1), ([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 1, -1, 0, 0, 0, 0, 0, 0, 0], 1), ([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 1, 1, -1, 0, 0, 0, 0], -1), ([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 1, -1, 0, 0], 1), ([1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1], -1), ([-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, -1, 0, 0], [1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], -1), ([-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, -1, -1, -1, -1, -1, -2, -2], [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 1, 0, 0, 0, 0, 2, 1], -1), ([-1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 2, 0, 0, 1, 0, 1, 0, 1], [1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0, 0, -1], -1)]
"""

"""rect = [([2, 2], [-1, -1], 1),
 ([-2, -2], [1, 1], 1),
 ([0, 4], [0, -2], 1),
 ([-1, -1], [0, 1], -1)]"""

"""
for eq in SnappyToEquation(rect):
    print eq
"""
    
def friendly_print(vec):
    eq = ''
    for i in range(len(vec)/3):
        var = 'z' + str(i)
        if vec[3*i] != 0:
            eq += var + '^' + str(vec[3*i]) + '*'
        if vec[3*i+1] != 0:
            eq += '(1/(1-' + var + '))^' + str(vec[3*i+1]) + '*'
        if vec[3*i+2] != 0:
            eq += '((' + var + '-1)/' + var + ')^' + str(vec[3*i+2]) + '*'
    return eq
"""  
print friendly_print([1, 0, 2, 3, 0, 0, 0, 2, -1, 2, 0, 0, 1, 0, 0, 2, 0, 0, 1, 0, 0, 2, 0
, 0, 2, 0, 0, 3, 0, 0, 1, 0, 0, 2, 0, 0, 1, 0, 0, 2, 0, 0, 1, 0, 0, 2, 0, 0, 2,
0, 0, 3, 0, 0, -1, -1, 0, 1, 0, 0, 2, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 1, 1, 0, 0])
"""
