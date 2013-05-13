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
    
if __name__ == '__main__':
    import sys
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        print "error! usage %s filebase nrefine [index_base] " % sys.argv[0]
    filebase = sys.argv[1]
    nrefine = int(sys.argv[2])
    index_base = 0
    if len(sys.argv) > 3:
        index_base = int(sys.argv[3])
    refine_uniform_from_tetgen(filebase,nrefine,index_base,EB=False)
