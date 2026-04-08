
# Triangulation of Budney-Burton-Hillman Cappell-Shaneson knot complement

from triangulation import *
import itertools

tri = FourTriangulation(2)
P = tri.pentachora
P[0].joinTo(3, P[0], [1,2,3,0,4])
P[1].joinTo(4, P[1], [0,2,3,4,1])
P[1].joinTo(2, P[0], [0,1,4,2,3])
P[1].joinTo(3, P[0], [0,2,3,1,4])
P[1].joinTo(0, P[0], [2,0,1,3,4])

#NTri = ftri2.pachner([0,1], [[0,1,2,3], [0,1,3,4]])
#NTri = ftri2.pachner([0], [[]])
#ftri2.writeToFile('C:/pr2011/4dim_triangulation/cappell.tri')
#NTri.writeToFile('C:/pr2011/4dim_triangulation/cappell_pachner.tri')
#NTri.link().writeSnapPea('C:/pr2011/4dim_triangulation/cappell_link_pachner.tri')
    
tri = tri.idealToFinite()
print "Size of finite vertex triangulation:", tri.numPentachora()
for i in range(3):
    tri = tri.boundaryReduce()

tri.boundaryTriangulation().writeSnapPea('C:/pr2011/4dim_triangulation/cappell_boundary.tri')
bTet = tri.boundaryTetrahedra()
print bTet
for t_no in bTet:
    p_no, v = t_no
    for j in [x for x in range(5) if x != v]:
        print t_no, 'vert', j, 'glues to', tri.boundaryGluing(t_no, j)
