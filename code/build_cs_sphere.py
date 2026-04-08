
from punctured_cs_bundle import TriangulationBuilder
from triangulation import FourTriangulation, Perm5

def csSpheres(monod):
    bundle = TriangulationBuilder.buildBundle(monod)
    bundle = bundle.simplifyForTime(1200000) # simplify for 3 seconds
    print "Bundle simplified to", bundle.numPentachora(), "pentachora"    
    finite_tri = bundle.idealToFinite()
    finite_tri = finite_tri.simplifyBoundaryUntil(2)
    assert finite_tri.numBoundaryTet() == 2, "Could not simplify boundary to 2 tetrahedra."
    filepath = "C:/Users/Ahmad/Dropbox/4dim_triangulation/"
    spheres = []
    for filename in ["s2xd2_gluck_piece1.tri", "s2xd2_gluck_piece2.tri"]:
        s2xd2 = FourTriangulation.loadFromFile(filepath + filename)
        spheres.append(finite_tri.glueAlongBoundaries(s2xd2))
    return spheres

def simplifyTriangulation(filepath, time_millis=3600000):
    tri = FourTriangulation.loadFromFile(filepath)
    tri = tri.simplifyForTime(time_millis)
    tri.writeToFile(filepath)

def createTriangulations():
    import pickle
    conj_reps = pickle.load(open("C:/Users/Ahmad/Dropbox/mphil/sage/conj_class_reps.p", "rb"))
    for tr in sorted(conj_reps.keys()):
        for monod in conj_reps[tr]:
            bundle = TriangulationBuilder.buildBundle(monod)
            bundle = bundle.simplifyForTime(240000)
            bundle.writeToFile("C:/Users/Ahmad/Dropbox/4dim_triangulation/data/bundles_simp2/tr" + str(tr) + "_" + str(monod) + ".tri")

def glueS2xD2(idealTri,
              s2xd2Filepath="C:/Users/Ahmad/Dropbox/4dim_triangulation/",
              s2xd2Filename1="s2xd2_gluck_piece1.tri",
              s2xd2Filename2="s2xd2_gluck_piece2.tri"):
    """
    Truncate ideal vertex then glue in S^2 x D^2 in two ways.
    Returns pair of resulting triangulations.
    """
    finiteTri = idealTri.idealToFinite().simplifyBoundaryUntil(2)
    assert finiteTri.numBoundaryTet() == 2, "Could not simplify boundary to 2 tetrahedra."
    spheres = []
    for filename in [s2xd2Filename1, s2xd2Filename2]:
        s2xd2 = FourTriangulation.loadFromFile(s2xd2Filepath + filename)
        spheres.append(finiteTri.glueAlongBoundaries(s2xd2))
    return spheres

"""
tri = FourTriangulation.loadFromFile("C:/Users/Ahmad/Dropbox/4dim_triangulation/cappell.tri")
spheres = glueS2xD2(tri)
for csSphereIndex in range(2):
    spheres[csSphereIndex].writeToFile("C:/Users/Ahmad/Dropbox/4dim_triangulation/bbh_cs_" + str(csSphereIndex) + ".tri")
"""

def batchPuncturedBundleToSurgered(
        bundleFilepath="C:/Users/Ahmad/Dropbox/4dim_triangulation/data/bundles_simp/",
        outputFilepath="C:/Users/Ahmad/Dropbox/4dim_triangulation/data/cs_spheres_tris/"
        ):
    """
    For every file in bundleFilepath, construct two Cappell-Shaneson spheres by gluing.
    Save these to files in outputFilepath.
    """
    import os
    filepath = bundleFilepath
    #filepath = "/home/aissa/Documents/MPhil_6_21_17/4dim_triangulation/data/bundles_simp/"
    csSpheresFilepath = "C:/Users/Ahmad/Dropbox/4dim_triangulation/data/cs_spheres_tris/"
    tris = [(FourTriangulation.loadFromFile(filepath + filename), filename) for filename in os.listdir(filepath)]
    tris.sort(key=lambda p: p[0].numPentachora())
    for tri, filename  in tris:
        print filename, tri.numPentachora()
        spheres = glueS2xD2(tri)
        filenameWithNoExt = filename[:-4] # remove '.tri'
        for csSphereIndex in range(2):
            spheres[csSphereIndex].writeToFile(csSpheresFilepath + filenameWithNoExt + "_cs" + str(csSphereIndex + 1) + '.tri')

            
"""
monod =[[0, 1, -12], [0, 0, 1], [1, 0, 13]]
tr = 13
#bundle = TriangulationBuilder.buildBundle(monod)
#bundle = bundle.simplifyForTime(240000)
filepath = "/home/aissa/Documents/MPhil_6_21_17/4dim_triangulation/data/bundles_simp/"
#filename = "C:/Users/Ahmad/Dropbox/4dim_triangulation/data/bundles_simp/tr" + str(tr) + "_" + str(monod) + ".tri"
#bundle.writeToFile(filepath + "tr" + str(tr) + "_" + str(monod) + ".tri")

##simplifyTriangulation(filename)
simplifyTriangulation(filepath + "tr" + str(tr) + "_" + str(monod) + ".tri")
"""

def filenameToMonodromy(filename):
    matrixStr = filename.split('_')[1][:-4]
    matrixEntries = [int(entry) for entry in matrixStr.replace('[', '').replace(']', '').split(',')]
    return [matrixEntries[0:3], matrixEntries[3:6], matrixEntries[6:]]

def trace(matrix):
    return sum(matrix[i][i] for i in range(3))

def matrixToLatex(matrix):
    body = ' \\\\ '.join([' & '.join(str(n) for n in matrix[row]) for row in range(3)])
    return "\\begin{pmatrix} " + body + " \\end{pmatrix}"
    
def tableOfCS(
        bundleFilepath="C:/Users/Ahmad/Dropbox/4dim_triangulation/data/bundles_simp/",
        csFilepath="C:/Users/Ahmad/Dropbox/4dim_triangulation/data/cs_spheres_tris/"
        ):
    import os
    rowData = []
    for filename in os.listdir(bundleFilepath):
        bundle = FourTriangulation.loadFromFile(bundleFilepath + filename)
        monodromy = filenameToMonodromy(filename)
        tr = trace(monodromy)
        cs1 = FourTriangulation.loadFromFile(csFilepath + filename[:-4] + '_cs1.tri')
        isoSig = bundle.isoSig()
        rowData.append((tr, monodromy, bundle.numPentachora(), cs1.numPentachora(), isoSig))
        
    rowData.sort(key=lambda row: row[0]) # sort by trace
    
    for row in rowData:
        tr, monodromy, bundlePentCount, cs1PentCount, isoSig = row
        latexRow = "{0} & ${1}$ & {2} & {3} & {{\\footnotesize \\texttt{{\\seqsplit{{{4}}}}}}} \\\\ \\hline".format(
            tr, matrixToLatex(monodromy), bundlePentCount, cs1PentCount, isoSig)
        print latexRow

def simplifyCS(inputFilename, outputFilename, time_in_millis=9*3600*1000):
    print "Simplifying:", inputFilename
    import time
    start_time = int(round(time.time())) # in seconds
    tri = FourTriangulation.loadFromFile(inputFilename)
    import simplify
    tri, isS4 = simplify.simplifyUntilS4(tri, time_in_millis=time_in_millis)
    #tri = simplify.simplifyForTime(tri, time_in_millis=time_in_millis)
    end_time = int(round(time.time())) # in seconds
    print "Simplified to S4?", isS4, "Time taken:", (end_time-start_time), inputFilename
    tri.writeToFile(outputFilename)

#filename = "C:/Users/Ahmad/Dropbox/4dim_triangulation/data/cs_spheres_tris_simp_backup/tr0_[[0, 1, 1], [0, 0, 1], [1, 0, 0]]_cs2.tri"

"""
inputFilename="C:/Users/Ahmad/Dropbox/4dim_triangulation/data/cs_spheres_tris_simp/tr1_[[0, 1, 0], [0, 0, 1], [1, 0, 1]]_cs2.tri"
outputFilename="C:/Users/Ahmad/Dropbox/4dim_triangulation/data/cs_spheres_tris_simp/tr1_[[0, 1, 0], [0, 0, 1], [1, 0, 1]]_cs2.tri"
simplifyCS(inputFilename=inputFilename, outputFilename=outputFilename)
FourTriangulation.debugFile.flush()
FourTriangulation.debugFile.close()
"""

#tri = FourTriangulation.loadFromFile(filename)
#print tri.isoSig()
filename = "C:/Users/Ahmad/Dropbox/4dim_triangulation/tr1_cs1.txt"
max_size = 0
min_size = 1000
min_sizes = []
max_sizes = []
with open(filename, "r") as fp:
    for line in fp:
        pachner, size = [int(x) for x in line.split(" ")]
        max_size = max(size, max_size)
        if size < min_size:
            min_size = size
            min_sizes.append((pachner, size))
            max_sizes.append(max_size)
            max_size = 0
print 'max size:', max_size, pachner
print 'max_sizes:', max_sizes
print 'min_sizes:', min_sizes



#simplifyTriangulation("C:/Users/Ahmad/Dropbox/4dim_triangulation/data/bundles_simp/tr13_[[0, 1, -12], [0, 0, 1], [1, 0, 13]].tri")
#spheres = csSpheres([[0,0,1], [-5, 2, 0], [-8, 3, -7]])
#spheres = csSpheres([[0,0,1], [5, 2, 0], [7, 3, 8]])
#print "CS Sphere sizes:", spheres[0].numPentachora(), spheres[1].numPentachora()
#trans = [[0, 5, 17], [0, 2, 7], [1, 0, 18]]
#orig = [[0,0,1], [5, 2, 0], [17, 7, 18]]
#trans2 = [[0, -11, -48], [0, 3, 13], [1, 0, -23]]
#print TriangulationBuilder.normalForm(trans2)
#spheres = csSpheres(trans2)
