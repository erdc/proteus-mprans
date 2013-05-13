"""
convenience script to generate xdmf and vtk output from tetgen files
"""
from proteus import Comm
Comm.init()
comm = Comm.get()
def refine_uniform_from_tetgen(filebase,refinementLevels,index_base=0,EB=False):
    #import pdb
    #pdb.set_trace()
    from proteus import MeshTools,Archiver
    mesh = MeshTools.TetrahedralMesh()
    mesh.generateFromTetgenFiles(filebase,index_base,skipGeometricInit=False)

    MLMesh = MeshTools.MultilevelTetrahedralMesh(0,0,0,skipInit=True)

    MLMesh.generateFromExistingCoarseMesh(mesh,refinementLevels)

    MLMesh.meshList[-1].writeTetgenFiles(filebase+'_out',index_base)

    ar = Archiver.XdmfArchive('.',filebase+'_out')
    import xml.etree.ElementTree as ElementTree
    ar.domain = ElementTree.SubElement(ar.tree.getroot(),"Domain")
    mesh.writeMeshXdmf(ar,'mesh_coarse'+'_out',init=True,EB=EB,tCount=0)
    MLMesh.meshList[-1].writeMeshXdmf(ar,'mesh_fine'+'_out',init=True,EB=EB,tCount=1)
    ar.close()
    
def refine_uniform_from_triangle(filebase,refinementLevels,index_base=0,EB=False):
    #import pdb
    #pdb.set_trace()
    from proteus import MeshTools,Archiver
    mesh = MeshTools.TriangularMesh()
    mesh.generateFromTriangleFiles(filebase,index_base)

    MLMesh = MeshTools.MultilevelTriangularMesh(0,0,0,skipInit=True)

    MLMesh.generateFromExistingCoarseMesh(mesh,refinementLevels)

    MLMesh.meshList[-1].writeTriangleFiles(filebase+'_out',index_base)

    ar = Archiver.XdmfArchive('.',filebase+'_out')
    import xml.etree.ElementTree as ElementTree
    ar.domain = ElementTree.SubElement(ar.tree.getroot(),"Domain")
    mesh.writeMeshXdmf(ar,'mesh_coarse'+'_out',init=True,EB=EB,tCount=0)
    MLMesh.meshList[-1].writeMeshXdmf(ar,'mesh_fine'+'_out',init=True,EB=EB,tCount=1)
    ar.close()
    
if __name__ == '__main__':
    import sys
    if len(sys.argv) < 4 or len(sys.argv) > 5:
        print "error! usage %s filebase spaceDim nrefine [index_base] " % sys.argv[0]
    filebase = sys.argv[1]
    spaceDim = int(sys.argv[2])
    assert spaceDim in [2,3], "only 2d or 3d for now, got dim= %s" % spaceDim
    nrefine = int(sys.argv[3])
    index_base = 0
    if len(sys.argv) > 4:
        index_base = int(sys.argv[4])
    if spaceDim == 3:
        refine_uniform_from_tetgen(filebase,nrefine,index_base,EB=False)
    else:
        refine_uniform_from_triangle(filebase,nrefine,index_base,EB=False)
