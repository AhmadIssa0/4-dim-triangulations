
# two dimensional triangulations

class TwoTriangulation:
    # triangles = list of Triangles
    # cusps = number of cusps
    
    def __init__(self, n=0, tri_title=None):
        # tri_title is the title of the triangulation as written in the triangulation file
        self.triangles = []
        self.cusps = 0
        for i in range(n):
            self.addTriangle()
        self.tri_title = tri_title
            
    def addTriangle(self, num=None):
        if num == None:
            num = self.nextTriIndex()
        tri = Triangle(num)
        self.triangles.append(tri)
        return tri
        
    def removeTriangle(self, num):
        # Removes Triangle with getTriNo() == num
        for i in range(self.numtriangles()):
            if self.triangles[i].getTriNo() == num:
                del self.triangles[i]
                break
                
    def renumberTriangulation(self):
        for i in range(self.numtriangles()):
            self.triangles[i].n = i
                
    def renumberTriangle(self, old, new):
        # Since tet.getNeighbour(i) returns a reference to a Triangle,
        # rather than a Triangle number, this is simple.
        self.getTriangleByNum(old).n = new


    def nextTriIndex(self):
        if len(self.triangles) == 0:
            return 0
        return max(map(lambda tet: tet.getTriNo(), self.triangles))+1
     
    def copy(self):
        # Requires the triangles in self to be numbered 0..(n-1)
        nTri = Triangulation(self.numtriangles())
        nTri.setNoCusps(self.getNoCusps())
        
        for i in range(self.numtriangles()):
            nTet = nTri.getTriangle(i)
            oldTet = self.getTriangle(i)
            nTet.setGluings(oldTet.getGluings())
            nTet.setCuspNos(oldTet.getCuspNos())
                
        for i in range(self.numtriangles()):
            nTet = nTri.getTriangle(i)
            oldTet = self.getTriangle(i)
            for j in range(3):
                neighbour_no = oldTet.getNeighbour(j).getTriNo()
                nTet.setNeighbour(j, nTri.getTriangleByNum(neighbour_no))
        return nTri
            
    def getTriangleByNum(self, n):
        for tet in self.gettrianglesList():
            if tet.getTriNo() == n:
                return tet
                
    def gettrianglesList(self):
        return self.triangles
        
    # return i-th Triangle in list of triangles
    # this may be different to Triangle.getTriNo() == i
    def getTriangle(self, i):
        return self.triangles[i]
        
    def numtriangles(self):
        return len(self.triangles)
        
    def setNoCusps(self, n):
        self.cusps = n
        
    def getNoCusps(self):
        return self.cusps
        
    @staticmethod
    def oppositeVertex(edge):
        for i in range(2):
            if i not in edge:
                return i
        
    def switchFace(self, tri):
        """ if tri = [a,c], then want [a,d] """
        return tri[:1] + [self.oppositeVertex(tri)]
                
    def gluedFace(self, tri_no, edge):
        tri = self.triangles[tri_no]
        face = self.oppositeVertex(edge)
        g = tri.getGluing(face)
        return (tri.getNeighbour(face).n, [g[edge[0]], g[edge[1]]])
    
    def writeToFile(self, foutname):
        tri = self.copy()
        tri.renumberTriangulation()
        
        f = open(foutname, 'w')
        f.write('% Triangulation\n')
        if tri.tri_title is not None:
            f.write(tri.tri_title + '\n')
        else:
            f.write('Triangulation\n')

        f.write(str(tri.numTriangles()) + '\n')
        
        for i in range(tri.numTriangles()):
            tet = tri.getTriangles(i)
            f.write('\n\t{0}\t{1}\t{2}\n'.format(tet.getNeighbour(0).getTriNo(),tet.getNeighbour(1).getTriNo(),
                                                    tet.getNeighbour(2).getTriNo()))
            
            for j in range(5):
                g = tet.getGluing(j)
                f.write('\t{0}{1}{2}'.format(g[0],g[1],g[2]))
            f.write('\n')
            f.write('\t{0}\t{1}\t{2}\n'.format(tet.getCuspNo(0),tet.getCuspNo(1),tet.getCuspNo(2),))
         
			
class Triangle:
	# Triangle neighbour[0-2]
	# Gluings gluing[0-2][3] where gluing[i][j] pairs face i and
	#	vertex j to tri neighbour[i] vertex gluing[i][j]
	# Cusp numbers corresponding to vertices cusp_no[0-3]
	
	def __init__(self, n=0):
		self.neighbour = [0]*3
		self.gluing = [[]]*3
		self.cusp_no = [0]*3
		self.n = n
		
	def joinTo(self, myFace, you, gluing):
		self.setNeighbour(myFace, you)
		self.setGluing(myFace, gluing)
		you.setNeighbour(gluing[myFace], self)
		
		you_gluing = [0]*3
		for i in range(3):
			you_gluing[gluing[i]] = i
			
		you.setGluing(gluing[myFace], you_gluing)
		
	def setNeighbour(self, myFace, you):
		self.neighbour[myFace] = you
		
	def setGluings(self, gluings):
		self.gluing = [x[:] for x in gluings]
		
	def setGluing(self, myFace, gluing):
		self.gluing[myFace] = gluing[:]
		
	def setCuspNo(self, vert, n):
		self.cusp_no[vert] = n
		
	def setCuspNos(self, cusps):
		self.cusp_no = cusps[:]
		
    def getNeighbour(self, myFace):
        return self.neighbour[myFace]
        
    def getGluing(self, myFace):
        return self.gluing[myFace][:]
        
    def getGluings(self):
        return [x[:] for x in gluings]
        
    def getCuspNos(self):
        return self.cusp_no[:]
        
    def getCuspNo(self, vert):
        return self.cusp_no[vert]
        
    def getTriNo(self):
        return self.n