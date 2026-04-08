
from triangulation import *
from s2xd2 import *

"""
tri1 = FourTriangulation.loadFromFile('C:/pr2011/4dim_triangulation/s2xd2_simp_86.tri')
tri = FourTriangulation.loadFromFile('C:/pr2011/4dim_triangulation/s2xd2_gluck_simp.tri')

#tri = FourTriangulation.loadFromFile("C:/pr2011/4dim_triangulation/cappell.tri")
glueAlongBoundaries(tri1, tri, 'C:/pr2011/4dim_triangulation/s2xs2.tri', debug=False)
tri = FourTriangulation.loadFromFile("C:/pr2011/4dim_triangulation/s2xs2.tri")
"""

"""
tri1 = FourTriangulation.loadFromFile('C:/pr2011/4dim_triangulation/triv_sph_complement_simp.tri')
tri1 = tri1.idealToFinite()
print tri1.numBoundaryTet(), tri1.numPentachora()
tri1.boundaryTriangulation().writeSnapPea('C:/pr2011/4dim_triangulation/s2xd2_boundary.tri')
"""

"""
for j in range(5):
    tri1 = randomBoundaryMove(tri1,2,3)
    tri1 = reduceVerts(tri1)
    tri1 = randomBoundaryMove(tri1,3,2)
    tri1 = reduceVerts(tri1)
"""

"""
for i in range(9):
    print 'deepReduce', i, tri1.numBoundaryTet()
    tri1 = deepReduce(tri1)
"""
    
"""
for i in range(2):
    print 'depth 3:', i, tri1.numBoundaryTet()
    tri1 = reduceWithDepth(tri1, 3)
    if tri1.numBoundaryTet() == 2:
        break
for i in range(20):
    for j in range(5):
        tri1 = randomBoundaryMove(tri1,2,3)
        tri1 = reduceVerts(tri1)
    for j in range(10):
        tri1 = reduce(tri1)
    if tri1.numBoundaryTet() == 2:
        break
"""

"""
print tri1.numBoundaryTet(), tri1.numPentachora()
if tri1.numBoundaryTet() == 4:
    tri1.writeToFile('C:/pr2011/4dim_triangulation/unknot_comp.tri')
tri = FourTriangulation.loadFromFile('C:/pr2011/4dim_triangulation/s2xd2_gluck_simp.tri')
"""

"""
tri1 = FourTriangulation.loadFromFile("C:/pr2011/4dim_triangulation/unknot_comp_capped.tri")
tri = FourTriangulation.loadFromFile('C:/pr2011/4dim_triangulation/s2xd2_gluck_simp.tri')
glueAlongBoundaries(tri1, tri, 'C:/pr2011/4dim_triangulation/hsphere.tri', debug=False)



tri = FourTriangulation.loadFromFile("C:/pr2011/4dim_triangulation/triv_sph_complement_simp.tri")
print len(tri.cycleByDim(0))
print tri.intersectionForm()
"""

"""
tri1 = FourTriangulation.loadFromFile('C:/pr2011/4dim_triangulation/unknot_comp.tri')
print tri1.boundaryTetrahedra()
print tri1.boundaryCycle([297,1], [0,2,3])
print tri1.boundaryCycle([297,1], [0,2,4])
print tri1.boundaryCycle([297,1], [2,3,4])
print tri1.boundaryCycle([297,1], [0,3,4])
p1 = tri1.addPentachoron()
p2 = tri1.addPentachoron()
for i in range(4):
    p1.joinTo(i, p2, range(5))
p1.joinTo(4, tri1.getPentachoronByNum(297), [0,2,3,4,1])
p2.joinTo(4, tri1.getPentachoronByNum(300), [3,4,2,1,0])
tri1.writeToFile('C:/pr2011/4dim_triangulation/unknot_comp_capped.tri')
"""

tri = FourTriangulation.loadFromFile("C:/Users/Ahmad/Dropbox/4dim_triangulation/cappell.tri")
tri = tri.idealToFinite()
print "boundary tet: ", tri.numBoundaryTet()
import cProfile as cp
cp.run("tri = tri.simplifyBoundaryUntil(2)")
print "boundary tet: ", tri.numBoundaryTet()
