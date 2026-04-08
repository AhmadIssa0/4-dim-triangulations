
import copy 
from triangulation import FourTriangulation

# Triangulates punctured Cappell-Shaneson bundles

class Point:
    def __init__(self, x=0, y=0, z=0):
        self.x = x
        self.y = y
        self.z = z
        
    def __copy__(self):
        return Point(self.x, self.y, self.z)
        
    def __add__(self, pt):
        return Point(self.x + pt.x, self.y + pt.y, self.z + pt.z)
        
    def __eq__(self, pt):
        return self.x == pt.x and self.y == pt.y and self.z == pt.z
        
    def __neg__(self):
        return Point(-self.x, -self.y, -self.z)
        
    def __repr__(self):
        return str((self.x, self.y, self.z))
        
    def toTuple(pt):
        return (pt.x, pt.y, pt.z)
        
class Layer:
    def __init__(self, e1=Point(1,0,0), e2=Point(0,1,0), e3=Point(0,0,1), tetrahedra=None):
        self.e1 = e1
        self.e2 = e2
        self.e3 = e3
        
        if tetrahedra is None:
            """ # this is the Budney/Hillman/Burton starting triangulation
            self.tetrahedra = [[Point(0,0,0), e3, e1, e1+e2+e3],
                               [Point(0,0,0), e3, e2+e3, e1+e2+e3],
                               [e1, e3, e1+e3, e1+e2+e3],
                               [Point(0,0,0), e1+e2, e1+e2+e3, e1],
                               [Point(0,0,0), e2+e3, e1+e2+e3, e1+e2],
                               [Point(0,0,0), e1+e2, e2, e2+e3]]
            """
            self.tetrahedra = [[Point(0,0,0), e1+e2+e3, e3, e1+e3],
                               [Point(0,0,0), e1+e2+e3, e1, e1+e3],
                               [Point(0,0,0), e1+e2+e3, e1, e1+e2],
                               [Point(0,0,0), e1+e2+e3, e3, e2+e3],
                               [Point(0,0,0), e1+e2+e3, e2, e2+e3],
                               [Point(0,0,0), e1+e2+e3, e2, e1+e2]]
                               
        else:
            self.tetrahedra = tetrahedra[:]
            
        self.above = [] # associations above this layer
        self.below = [] # associations below this layer
        
    def printTetrahedra(self):
        for tet in self.tetrahedra:
            print map(Point.toTuple, tet)
    
    @staticmethod
    def printTetList(tets):
        for tet in tets:
            print map(Point.toTuple, tet)
            
    def copy(self):
        return Layer(e1=self.e1, e2=self.e2, e3=self.e3, 
                        tetrahedra=copy.deepcopy(self.tetrahedra))

    def eqInTorus(self, pt1, pt2):
        return (pt1 == pt2 or pt1 + self.e1 == pt2 or pt1 + self.e2 == pt2
                or pt1 + self.e3 == pt2 or pt2 + self.e1 == pt1 or pt2 + self.e2 == pt1
                or pt2 + self.e3 == pt1)
        
    def tetClosure(self, subsimplex, delete=False):
        """ subsimplex: list of points
        delete: bool, if true we delete matching tetrahedra.
        returns a list of tetrahedra which contain the vertices of subsimplex
        here we identify across a FD face if necessary.
        """
        
        res = []
        indices = [] # indices of matched tetrahedra
       
        for i in range(len(self.tetrahedra)):
            tet = self.tetrahedra[i]
            matched_pts = 0
            ordered_tet = []
            for pt in subsimplex:
                for pt2 in tet:
                    if pt == pt2:
                        ordered_tet.append(pt2)
                        matched_pts += 1

            if matched_pts == len(subsimplex):
                indices.append(i)
                for pt in tet:
                    if pt not in ordered_tet:
                        ordered_tet.append(pt)
                res.append(ordered_tet)

        if delete:
            for i in reversed(indices):
                del self.tetrahedra[i]
        return res
        
    @staticmethod
    def eqTets(tet1, tet2):
        return set(map(Point.toTuple, tet1)) == set(map(Point.toTuple, tet2))
        
    def assocAbove(self, tet, pent_verts, pent=False):
        self.above.append((tet, pent_verts, pent))
        
    def getBelow(self, tet):
        for (t, pts, is_pent) in self.below:
            if Layer.eqTets(t, tet):
                res = []
                for pt in tet:
                    res.append(pts[t.index(pt)])
                return (res, is_pent)
        return None
        
    def getAbove(self, tet):
        for (t, pts, is_pent) in self.above:
            if Layer.eqTets(t, tet):
                res = []
                for pt in tet:
                    res.append(pts[t.index(pt)])
                return (res, is_pent)
        return None

    def printAbove(self):
        for (x,y,p) in self.above:
            print (map(Point.toTuple, x), y)
        
    def assocBelow(self, tet, pent_verts, pent=False):
        self.below.append((tet, pent_verts, pent))
        
    def printBelow(self):
        for (x,y,p) in self.below:
            print (map(Point.toTuple, x), y)
            
    def clearAssoc(self):
        self.above = []
        self.below = []
            
    @staticmethod
    def tetInTets(tet, tets):
        return reduce(lambda x, y: x or y, map(lambda x: Layer.eqTets(tet, x), tets))
        
    @staticmethod
    def translate(pts, e):
        return map(lambda pt: pt+e, pts)
        
    def pachner(self, subsimplex, onFace=False, e=None):
        # if onFace is True then subsimplex is an edge,
        # and the edge lies on a face of the FD.
        # In this case subsimplex should be on the face where the flat
        # tetrahedron lies.
        
        # e is a Point which tells that the opposite face is given by going
        # in the e direction.
        # WARNING: subsimplex must be specified on face, such that
        # the opposite face is given by ADDING e. The flat tetrahedron
        # must be on that opposite face.
        
        # returns a new Layer with new triangulation given by
        # performing a pachner move on common face subsimplex
        layer = self.copy()
        layer.clearAssoc()
        tets = layer.tetClosure(subsimplex, delete=True)
        
        if onFace: # extra tets in closure given by translating face by (-e)
            # this will be the unflattened tetrahedron
            flat_tet = layer.tetClosure(map(lambda t: t + e, subsimplex), delete=True)[0]
        
        # update trivial associations
        if onFace:
            other_tets = filter(lambda t: not Layer.tetInTets(t, tets+[flat_tet]), self.tetrahedra)
        else:
            other_tets = filter(lambda t: not Layer.tetInTets(t, tets), self.tetrahedra)
        for t in other_tets:
            self.assocAbove(t, t, False)
            layer.assocBelow(t, t, False)
        
        if len(subsimplex) == 3: # we're doing a 2-3 pachner move
            pt1, pt2 = tets[0][-1], tets[1][-1]
            
            # 5 points of 4-simplex: pt1, pt2, subsimplex[0-2]
            dic = {}
            dic[pt1.toTuple()] = 0
            dic[pt2.toTuple()] = 1
            for i in range(len(subsimplex)):
                dic[subsimplex[i].toTuple()] = i+2
                
            for t in tets:
                self.assocAbove(t, map(lambda x: dic[x.toTuple()], t), pent=True)
            
            for i in range(len(subsimplex)):
                n_tet = [pt1, pt2, subsimplex[i%3], subsimplex[(i+1)%3]]
                layer.tetrahedra.append([pt1, pt2, subsimplex[i%3], subsimplex[(i+1)%3]])
                print 'new tet:', [pt1, pt2, subsimplex[i%3], subsimplex[(i+1)%3]]
                layer.assocBelow(n_tet, map(lambda x: dic[x.toTuple()], n_tet), pent=True)
        elif len(subsimplex) == 2: # we're doing a 3-2 pachner move
            tri = []
            for x in tets:
                for pt in x[2:]:
                    if pt not in tri:
                        tri.append(pt)
                        
            n_tet1 = [subsimplex[0]] + tri
            n_tet2 = [subsimplex[1]] + tri
            layer.tetrahedra.append(n_tet1)
            layer.tetrahedra.append(n_tet2)
                        
            # 5 points of 4-simplex: subsimplex + tri
            dic = {}
            dic[subsimplex[0].toTuple()] = 0
            dic[subsimplex[1].toTuple()] = 1
            for i in range(len(tri)):
                dic[tri[i].toTuple()] = i+2
                
            for t in tets:
                self.assocAbove(t, map(lambda x: dic[x.toTuple()], t), pent=True)
                
            if onFace:
                print 'flat_tet', flat_tet
                print 'assoc to', map(lambda pt: dic[pt.toTuple()], Layer.translate(flat_tet,-e))
                self.assocAbove(flat_tet, map(lambda pt: dic[pt.toTuple()], Layer.translate(flat_tet,-e)), pent=True)

            layer.assocBelow(n_tet1, map(lambda x: dic[x.toTuple()], n_tet1), pent=True)
            layer.assocBelow(n_tet2, map(lambda x: dic[x.toTuple()], n_tet2), pent=True)
        return layer
        
    def shearFD(self, subsimp, Z):
        # still need to manually update new e1, e2, e3
        layer = self.copy()
        layer.clearAssoc()
        
        for tet in layer.tetrahedra:
            if subsimp[0] in tet or subsimp[1] in tet:
                ctet = tet[:]
                for i in range(len(tet)):
                    tet[i] = tet[i] + Z
                self.assocAbove(ctet, tet)
                layer.assocBelow(tet, ctet)
            else:
                self.assocAbove(tet, tet, pent=False)
                layer.assocBelow(tet, tet, pent=False)
        return layer
        
class TriangulationBuilder:

    def __init__(self):
        self.layers = [Layer()]
        self.layer_to_pno = {} # layer_to_pno[layer] = p_no is above layer
        self.n_pents = 0
        
    def pachner(self, subsimplex, onFace=False, e=None):
        self.layer_to_pno[len(self.layers)-1] = self.n_pents
        self.n_pents += 1
        self.layers.append(self.layers[-1].pachner(subsimplex, onFace, e))
        
    def setIso(self, fn):
        # fn is a function from R^3 to R^3 (Point -> Point)
        # which is the isomorphism from layers[0] to layers[-1]
        for tet in self.layers[0].tetrahedra:
            assoc_tet = map(fn, tet)
            self.layers[0].assocBelow(tet, assoc_tet, pent=False)
            self.layers[-1].assocAbove(assoc_tet, tet, pent=False)
            
    def goDown(self, n_layer, tet):
        res = self.layers[n_layer].getBelow(tet)
        #print 'bad tet', tet
        #print self.layers[n_layer].tetrahedra
        #print 'searching in layer', n_layer, 'below', tet
        #print 'tets below this layer', self.layers[n_layer].below
        if res[1]: # found a pentachoron face
            n = (n_layer-1)%len(self.layers)
            return (self.layer_to_pno[n], res[0])
        else:
            return self.goDown((n_layer-1)%len(self.layers), res[0])
        
    def goUp(self, n_layer, tet):
        res = self.layers[n_layer].getAbove(tet)
        #print 'searching in layer', n_layer, 'above', tet, 'got', res
        
        if res[1]: # found a pentachoron face
            return (self.layer_to_pno[n_layer], res[0])
        else:
            return self.goUp((n_layer+1)%len(self.layers), res[0])
        
    def shearFD(self, subsimp, Z, e1, e2, e3):
        # simplices containing subsimp get shifted by Z
        # this updates e1, e2, e3 for the new fundamental domain
        self.layers.append(self.layers[-1].shearFD(subsimp, Z))
        layer = self.layers[-1]
        layer.e1, layer.e2, layer.e3 = e1, e2, e3
        
    def shearXZ(self):
        # shears X axis in direction of Z axis
        # the FD x coord vector e1 now goes to e1+e3
        # Assumes FD triangulation is in standard form.
        # WARNING: assumes shearFD doesn't change e1,e2,e3
        layer = self.layers[-1]
        e1, e2, e3 = layer.e1, layer.e2, layer.e3
        
        self.shearFD([e1, e1+e2], e3, e1+e3, e2, e3)
        
        self.pachner([e1+e3, e3, e1+e2+e3])
        self.pachner([e3, e1+e2+e3])
        
        self.pachner([Point(0,0,0), e1+e2+e3, e2+e3])
        self.pachner([e3, e1+e3], onFace=True, e=e2)

    def shearYZ(self):
        layer = self.layers[-1]
        e1, e2, e3 = layer.e1, layer.e2, layer.e3
        
        self.shearFD([e2, e1+e2], e3, e1, e2+e3, e3)
        
        self.pachner([e2+e3, e3, e1+e2+e3])
        self.pachner([e3, e1+e2+e3])
        
        self.pachner([Point(0,0,0), e1+e2+e3, e1+e3])
        self.pachner([e3, e2+e3], onFace=True, e=e1)
        
    def shearYX(self):
        layer = self.layers[-1]
        e1, e2, e3 = layer.e1, layer.e2, layer.e3
        
        self.shearFD([e2, e3+e2], e1, e1, e2+e1, e3)
        
        self.pachner([e2+e1, e1, e1+e2+e3])
        self.pachner([e1, e1+e2+e3])
        
        self.pachner([Point(0,0,0), e1+e2+e3, e1+e3])
        self.pachner([e1, e2+e1], onFace=True, e=e3)
        
    def shearZX(self):
        layer = self.layers[-1]
        e1, e2, e3 = layer.e1, layer.e2, layer.e3
        
        self.shearFD([e3, e3+e2], e1, e1, e2, e3+e1)
        
        self.pachner([e3+e1, e1, e1+e2+e3])
        self.pachner([e1, e1+e2+e3])
        
        self.pachner([Point(0,0,0), e1+e2+e3, e1+e2])
        self.pachner([e1, e3+e1], onFace=True, e=e2)
        
    def shearXY(self):
        layer = self.layers[-1]
        e1, e2, e3 = layer.e1, layer.e2, layer.e3
        
        self.shearFD([e1, e1+e3], e2, e1+e2, e2, e3)
        
        self.pachner([e1+e2, e2, e1+e2+e3])
        self.pachner([e2, e1+e2+e3])
        
        self.pachner([Point(0,0,0), e1+e2+e3, e2+e3])
        self.pachner([e2, e1+e2], onFace=True, e=e3)
        
    def shearZY(self):
        layer = self.layers[-1]
        e1, e2, e3 = layer.e1, layer.e2, layer.e3
        
        self.shearFD([e3, e1+e3], e2, e1, e2, e3+e2)
        
        self.pachner([e3+e2, e2, e1+e2+e3])
        self.pachner([e2, e1+e2+e3])
        
        self.pachner([Point(0,0,0), e1+e2+e3, e2+e1])
        self.pachner([e2, e3+e2], onFace=True, e=e1)
        
    def shearXZinv(self):
        layer = self.layers[-1]
        e1, e2, e3 = layer.e1, layer.e2, layer.e3
        
        self.pachner([Point(0,0,0), e1+e2+e3, e2])
        self.pachner([Point(0,0,0), e1+e3], onFace=True, e=e2)
        
        self.pachner([Point(0,0,0), e1+e2+e3, e1])
        self.pachner([Point(0,0,0), e1+e2+e3])
        
        self.shearFD([e1+e3, e1+e2+e3], -e3, e1+(-e3), e2, e3)
        
    def shearXYinv(self):
        layer = self.layers[-1]
        e1, e2, e3 = layer.e1, layer.e2, layer.e3
        
        self.pachner([Point(0,0,0), e1+e2+e3, e3])
        self.pachner([Point(0,0,0), e1+e2], onFace=True, e=e3)
        
        self.pachner([Point(0,0,0), e1+e2+e3, e2+e3])
        self.pachner([Point(0,0,0), e1+e2+e3])
        
        self.shearFD([e1+e2, e1+e2+e3], -e2, e1+(-e2), e2, e3)
        
    def shearYZinv(self):
        layer = self.layers[-1]
        e1, e2, e3 = layer.e1, layer.e2, layer.e3
        
        self.pachner([Point(0,0,0), e1+e2+e3, e1])
        self.pachner([Point(0,0,0), e2+e3], onFace=True, e=e1)
        
        self.pachner([Point(0,0,0), e1+e2+e3, e2])
        self.pachner([Point(0,0,0), e1+e2+e3])
        
        self.shearFD([e2+e3, e1+e2+e3], -e3, e1, e2+(-e3), e3)
        
    def shearYXinv(self):
        layer = self.layers[-1]
        e1, e2, e3 = layer.e1, layer.e2, layer.e3
        
        self.pachner([Point(0,0,0), e1+e2+e3, e3])
        self.pachner([Point(0,0,0), e2+e1], onFace=True, e=e3)
        
        self.pachner([Point(0,0,0), e1+e2+e3, e2])
        self.pachner([Point(0,0,0), e1+e2+e3])
        
        self.shearFD([e2+e1, e1+e2+e3], -e1, e1, e2+(-e1), e3)
        
    def shearZXinv(self):
        layer = self.layers[-1]
        e1, e2, e3 = layer.e1, layer.e2, layer.e3
        
        self.pachner([Point(0,0,0), e1+e2+e3, e2])
        self.pachner([Point(0,0,0), e3+e1], onFace=True, e=e2)
        
        self.pachner([Point(0,0,0), e1+e2+e3, e1+e2])
        self.pachner([Point(0,0,0), e1+e2+e3])
        
        self.shearFD([e3+e1, e1+e2+e3], -e1, e1, e2, e3+(-e1))
        
    def shearZYinv(self):
        layer = self.layers[-1]
        e1, e2, e3 = layer.e1, layer.e2, layer.e3
        
        self.pachner([Point(0,0,0), e1+e2+e3, e1])
        self.pachner([Point(0,0,0), e3+e2], onFace=True, e=e1)
        
        self.pachner([Point(0,0,0), e1+e2+e3, e1+e2])
        self.pachner([Point(0,0,0), e1+e2+e3])
        
        self.shearFD([e3+e2, e1+e2+e3], -e2, e1, e2, e3+(-e2))
        
    def build(self):
        triang = FourTriangulation(self.n_pents)
        for i in range(len(self.layers)):
            layer = self.layers[i]
            
            for tet in layer.tetrahedra:
                p1_no, vs1 = self.goUp(i, tet)
                p2_no, vs2 = self.goDown(i, tet)
                op1_vert = [x for x in range(5) if x not in vs1][0]
                op2_vert = [x for x in range(5) if x not in vs2][0]
                gluing = [0]*5
                gluing[op1_vert] = op2_vert
                for j in range(len(vs1)):
                    gluing[vs1[j]] = vs2[j]
                p1 = triang.getPentachoronByNum(p1_no)
                p2 = triang.getPentachoronByNum(p2_no)
                p1.joinTo(op1_vert, p2, gluing)
        return triang

    @staticmethod
    def isInGL3Z(a):
        return (a[0][0] * (a[1][1] * a[2][2] - a[2][1] * a[1][2])
                -a[1][0] * (a[0][1] * a[2][2] - a[2][1] * a[0][2])
                +a[2][0] * (a[0][1] * a[1][2] - a[1][1] * a[0][2])) in [-1, 1]
    @staticmethod
    def buildBundle(monod):
        assert TriangulationBuilder.isInGL3Z(monod), "Matrix given is not in GL_3(Z)."
        # monod given as a matrix
        tri_b = TriangulationBuilder()
        # triangulate T^3 x [0,1] first by layering on tetrahedra.
        """
        tri_b.shearXY()
        tri_b.shearXYinv()
        tri_b.shearYX()
        tri_b.shearYXinv()
        tri_b.shearXZ()
        tri_b.shearXZinv()
        tri_b.shearYZ()
        tri_b.shearYZinv()
        tri_b.shearZX()
        tri_b.shearZXinv()
        tri_b.shearZY()
        tri_b.shearZYinv()
        """
        
        tri_b.shearXZ()
        tri_b.shearXZinv()
        tri_b.shearXY()
        tri_b.shearXYinv()
        tri_b.shearYZ()
        tri_b.shearYZinv()
        
        prod = TriangulationBuilder.normalForm(monod)
        # the permutation matrix changes what we consider our axes
        # so our shearing matrices don't represent shears 
        # in terms of the xyz coordinates, but a permuted version
        perm = [0,1,2]
        for P in prod:
            if P[0] == 'SwitchRows':
                a, b = P[1:]
                perm[a], perm[b] = perm[b], perm[a]
        
        for move in prod:
            if move[0] == 'Shear':
                n, i, j = move[1:]
                i, j = perm[i], perm[j] # we act on the permuted axes
                if n > 0:
                    for k in range(n):
                        if i == 1 and j == 0:
                            tri_b.shearXY()
                        elif i == 2 and j == 0:
                            tri_b.shearXZ()
                        elif i == 0 and j == 1:
                            tri_b.shearYX()
                        elif i == 2 and j == 1:
                            tri_b.shearYZ()
                        elif i == 0 and j == 2:
                            tri_b.shearZX()
                        elif i == 1 and j == 2:
                            tri_b.shearZY()
                elif n < 0:
                    for k in range(-n):
                        if i == 1 and j == 0:
                            tri_b.shearXYinv()
                        elif i == 2 and j == 0:
                            tri_b.shearXZinv()
                        elif i == 0 and j == 1:
                            tri_b.shearYXinv()
                        elif i == 2 and j == 1:
                            tri_b.shearYZinv()
                        elif i == 0 and j == 2:
                            tri_b.shearZXinv()
                        elif i == 1 and j == 2:
                            tri_b.shearZYinv()
                print tri_b.layers[-1].e1, tri_b.layers[-1].e2, tri_b.layers[-1].e3
                
        def iso(p):
            A = monod
            return Point(p.x*A[0][0] + p.y*A[0][1] + p.z*A[0][2], 
                         p.x*A[1][0] + p.y*A[1][1] + p.z*A[1][2], 
                         p.x*A[2][0] + p.y*A[2][1] + p.z*A[2][2])
        tri_b.setIso(iso)
        # do the right sequence of shears...
        return tri_b.build()
        
    @staticmethod
    def normalAsMatrices(prod):
        # takes normalForm()'s given as prod output and writes as matrices
        mats = []
        for x in prod:
            if x[0] == 'Shear':
                id = [[1,0,0],[0,1,0],[0,0,1]]
                n, i, j = x[1:]
                id[i][j] = n
                mats.append(id)
            elif x[0] == 'SwitchRows':
                mats.append([[1,0,0],[0,1,0],[0,0,1]])
                i,j = x[1],x[2]
                M = mats[-1]
                M[i],M[j] = M[j],M[i]
        return mats
        
    @staticmethod
    def normalForm(A):
        # writes the matrix A in GL_3(Z) as a product:
        # A = P A_1 ... A_n where
        # 1) P is a permutation matrix
        # 2) A_i is obtained from the identity matrix by changing a non-diagonal entry to 1 or -1.
        
        # ('Shear', n, i, j) is obtained from identity by changing (i,j) entry to n
        
        n_cols = 3
        n_rows = 3
        
        prod = []
        A = [x[:] for x in A] # make a copy so we don't edit the input matrix
        
        def switch(i,j):
            A[i], A[j] = A[j], A[i]
            prod.append(['SwitchRows', i, j])
            
        def negateRow(i):
            A[i] = [-x for x in A[i]]
            # negation of a row can be written as a product of shears and a permutation
            if i == 0:
                prod.extend([['Shear',1,1,0], ['Shear',-1,0,1], ['Shear',1,1,0], ['SwitchRows',0,1]])
            elif i == 1:
                prod.extend([['Shear',1,2,1], ['Shear',-1,1,2], ['Shear',1,2,1], ['SwitchRows',1,2]])
            elif i == 2:
                prod.extend([['Shear',1,1,2], ['Shear',-1,2,1], ['Shear',1,1,2], ['SwitchRows',1,2]])
            
        def addMultiple(c, i, j):
            for k in range(3):
                A[j][k] += c*A[i][k]
            prod.append(['Shear', -1*c, j, i])
        
        p_i = 0 # pivot column index
        for p_j in range(n_cols):
            # switch rows to make pivot non-zero
            for i in range(p_i,n_rows):
                if A[i][p_j] != 0:
                    if i > p_i:
                        switch(i,p_i)
                    break
            if A[p_i][p_j] == 0: # no pivot found..
                raise Exception('Matrix input: ' + str(A) + ' not in GL_3(Z)', 'Pivot attempt:' + str(p_i))
            if A[p_i][p_j] < 0:
                negateRow(p_i)
            for i in range(p_i+1,n_rows):
                while A[i][p_j] != 0:
                    if A[i][p_j] < 0:
                        negateRow(i)
                    if A[i][p_j] < A[p_i][p_j]:
                        switch(i,p_i)
                    c = A[i][p_j] / A[p_i][p_j]
                    addMultiple(-1*c, p_i, i)
            ## eliminate entries above ##
            for i in range(0, p_i):
                if A[i][p_j] % A[p_i][p_j] == 0:
                    c = A[i][p_j] / A[p_i][p_j]
                    if c != 0:
                        addMultiple(-1*c, p_i, i)
            ####
            p_i += 1
            
        
        # if P is a permutation matrix and A_ij is a shear matrix
        # then P A_ij P^-1 = A_{P(i),P(j)}
        # If P switches 2 rows then P = P^-1.
        
        # bring all the permutation matrices to the front
        for i in range(len(prod)):
            if prod[i][0] == 'SwitchRows':
                j = i
                a,b = prod[i][1], prod[i][2]
                while j-1 >= 0 and prod[j-1][0] == 'Shear':
                    prod[j], prod[j-1] = prod[j-1], prod[j]
                    
                    for k in range(2,4):
                        if prod[j][k] == a:
                            prod[j][k] = b
                        elif prod[j][k] == b:
                            prod[j][k] = a
                    j -= 1
                
        return prod
        

"""
# BBH triangulation
tri_b = TriangulationBuilder()
tri_b.pachner([Point(0,0,0), Point(1,1,1), Point(1,1,0)])
tri_b.pachner([Point(0,0,0), Point(1,1,1)])
tri_b.shearFD([Point(1,0,1), Point(1,1,1)], Point(0,0,-1))

def iso(p):
    return Point(p.z, p.x, p.y-p.z)
tri_b.setIso(iso)
"""
"""
tri_b = TriangulationBuilder()
tri_b.shearXY()

tri_b.shearYX()

tri_b.shearXZ()
tri_b.shearXZ()
tri_b.shearYZ()
tri_b.shearZX()
tri_b.shearZY()

tri_b.shearZYinv()
tri_b.shearZXinv()
tri_b.shearYZinv()
tri_b.shearXZinv()
tri_b.shearXZinv()

tri_b.shearYXinv()

tri_b.shearXYinv()



#print tri_b.layers[-2].getAbove([Point(0, 1, 1), Point(1, 1, 1), Point(0, 1, 0), Point(1, 1, 2)])


def iso(p):
    return Point(p.y, p.z, p.x)
    #return Point(p.y, p.y+p.z, p.x+p.z)
tri_b.setIso(iso)
"""
#tri_b.build().writeToFile("C:/pr2011/cs_bundle.tri")
#tri_b.build().link().writeSnapPea("C:/pr2011/4dim_triangulation/cs0_bundle_link.tri")
#print tri_b.goUp(0, [Point(0,0,0), Point(1,0,0), Point(0,0,1), Point(1,1,1)])


#print TriangulationBuilder.normalForm([[0,1,0],[0,1,1],[1,0,-5]])

"""
monod = [[0,1,0],[0,1,1],[1,0,-5]]
monad = [[1,0,0],[0,1,0],[0,0,1]]
filepath = "C:/Users/Ahmad/Dropbox/4dim_triangulation/cs_test_link.tri"
TriangulationBuilder.buildBundle(monad).link().writeSnapPea(filepath)
"""

"""
filepath = "C:/Users/Ahmad/Dropbox/4dim_triangulation/cs_test_link.tri"
tri = FourTriangulation.loadFromFile("C:/Users/Ahmad/Dropbox/4dim_triangulation/cappell.tri")
tri.link().writeSnapPea(filepath)
"""

"""
layer = Layer()
e1, e2, e3 = layer.e1, layer.e2, layer.e3
layer.printTetrahedra()
new_layer = layer.pachner([Point(0,0,0), Point(1,1,1), Point(1,1,0)])

print 'new:'
new_layer.printTetrahedra()

old_layer = new_layer.pachner([Point(0,0,0), Point(1,1,1)])
print 'old:'
old_layer.printTetrahedra()

print 'sheared'
old_layer.shearFD([Point(1,0,1), Point(1,1,1)], Point(0,0,-1)).printTetrahedra()



old_layer.printBelow()
new_layer.printAbove()

new_layer.printAbove()
print 'below:'
old_layer.printBelow()
"""

