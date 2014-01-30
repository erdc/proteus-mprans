#! /usr/bin/env python
"""
define the basic geometry for the jet erosion test

see jet_domain.tex for description
"""

#define the lengths [cm]
#radius of nozel 
ld=0.32
#length of velocity inlet
lin=4
#length of separation between inlet region and wall
lp=1
#distance from nozzel to bottom
hb=10.
#vertical length of nozzel at restriction
hd=.127
#distance from nozzel to inflow
hin=lin
#length of tank
lT=15
#height of tank
hT=hb+hd+hin
#length of pressure outlet (in vertical)
hout=2

openfoam_header = """
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.2.2                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters  1.0;

"""
openfoam_footer = """
mergePatchPairs
(
);

// ************************************************************************* //
"""

def build_vertices_and_segments():
    """
    create vertices and the segments connecting them
    """
    vertices = [ (0,0),
                 (hb,0),
                 (hb+hd,0),
                 (hT,0),
                 (hT,ld),
                 (hT,lin),
                 (hb+hd,lin),
                 (hb+hd,ld),
                 (hb,ld),
                 (hb,lin+lp),
                 (hT-hout,lin+lp),
                 (hT,lin+lp),
                 (hT,lT),
                 (hT-hout,lT),
                 (hb,lT),
                 (0,lT),
                 (0,lin+lp),
                 (0,ld)
                 ]
    nv = len(vertices)

    segments = []
    for iv in range(nv-1):
        segments.append((iv,iv+1))
    #
    segments.append((nv-1,0))
    ## include flags as well
    wall_id = 1
    inflow_id  = 2
    outlet_id  = 3
    symmetry_id= 4
    boundaryTags = {'wall':wall_id,
                    'inflow':inflow_id,
                    'outlet':outlet_id,
                    'symmetry':symmetry_id}

    vertexFlags = [wall_id]*nv
    #inflow
    for iv in [3,4,5]:
        vertexFlags[iv] = inflow_id
    #outflow
    for iv in [12,13]:
        vertexFlags[iv] = outlet_id
    #symmetry boundary
    for iv in [0,1,2]:
        vertexFlags[iv] = symmetry_id

    ns = len(segments)
    segmentFlags = [0]*ns
    #inflow
    for ii in [3,4]:
        segmentFlags[ii] = inflow_id
    #outflow
    for ii in [12]:
        segmentFlags[ii] = outlet_id
    #symmetry boundary
    for ii in [0,1,2]:
        segmentFlags[ii] = symmetry_id

    
    return vertices,segments,vertexFlags,segmentFlags,boundaryTags


###convenience functions
def build_bottom_faces_2d():
    """
    cell connectivity for 2d geometry (bottom of domain)
    """
    faces_bottom = [(0,1,8,17),
                    (1,2,7,8),
                    (2,3,4,7),
                    (7,4,5,6),
                    (17,8,9,16),
                    (16,9,14,15),
                    (9,10,13,14),
                    (10,11,12,13)]
    return faces_bottom
def build_inlet_faces_2d(offset):
    inlet = [(3, 4, 4+offset, 3+offset),
             (4, 5, 5+offset, 4+offset)]
    return inlet
def build_outlet_faces_2d(offset):
    p_outer = [(12, 13, 13+offset, 12+offset)]
    return p_outer
def build_symmetry_faces_2d(offset):
    symmetry = [(0, 1, 1+offset, 0+offset),
                (1, 2, 2+offset, 1+offset),
                (2, 3, 3+offset, 2+offset)]
    return symmetry
def build_walls_2d(offset): 
    walls_bottom = [(0, 0+offset, 17+offset, 17),
                    (17, 17+offset, 16+offset, 16),
                    (16, 16+offset, 15+offset, 15)]

    walls_outer  = [(15, 15+offset, 14+offset, 14),
                    (14, 14+offset, 13+offset, 13)]

    walls_top = [(11, 12, 12+offset, 11+offset),
                 (10, 11, 11+offset, 10+offset),
                 (9, 10, 10+offset, 9+offset),
                 (9, 9+offset, 8+offset, 8),
                 (7, 7+offset, 6+offset, 6),
                 (6, 5, 5+offset, 6+offset)]

    walls = walls_bottom+walls_outer+walls_top
    return walls

def write_vertices(f,vertices):
    f.write('vertices \n (\n ')
    for v in vertices:
        f.write('\t(%12.5e %12.5e %12.5e)\n' % (v[0],v[1],v[2]))
    f.write(');\n\n')
    
def write_faces(f,faces):
    for face in faces:
        f.write('\t\t(')
        for iv in face:
            f.write(' %d ' % iv)
        f.write(') \n')
    f.write('\t\t);\n')
    f.write('\t} \n')

def write_blocks(f,blocks,nc,grading):
    f.write('blocks\n (\n')
    for b in blocks:
        f.write('  hex (')
        for iv in b:
            f.write(' %d ' % iv)
        f.write(') ')
        #discretization details
        f.write(' (%d %d %d) simpleGrading (%12.5e %12.5e %12.5e)\n' % (nc[0],nc[1],nc[2],grading[0],grading[1],grading[2]))
    f.write(');\n\n')
    
def generate_blockMeshDict_2d(width=1.0,filename='blockMeshDict',header=openfoam_header,footer=openfoam_footer):
    """
    Try to generate an OpenFOAM blockMeshDict for 2d simulation using vertices and segments
    """

    vertices2d,segments2d,vflags2d,sflags2d,btags2d = build_vertices_and_segments()

    nv2d = len(vertices2d)
    #vertices in 3d
    z_bottom=0.0; z_top = z_bottom+width
    vertices = []
    
    for v in vertices2d:
        vertices.append((v[0],v[1],z_bottom))
    for v in vertices2d:
        vertices.append((v[0],v[1],z_top))

    nv = len(vertices)

    faces_bottom = build_bottom_faces_2d()

    faces_top = []
    for face in faces_bottom:
        top = tuple([nv2d+iv for iv in face])
        faces_top.append(top)

    blocks = []
    for i,face in enumerate(faces_bottom):
        hex = face+faces_top[i]
        blocks.append(hex)

    #boundary faces
    inlet = build_inlet_faces_2d(nv2d)

    p_outer = build_outlet_faces_2d(nv2d)

    symmetry = build_symmetry_faces_2d(nv2d)


    walls = build_walls_2d(nv2d)

    #reverse faces on the bottom for boundary conditions
    faces_bottom_bc = []
    for face in faces_bottom:
        face_rev = [face[0],face[3],face[2],face[1]]
        faces_bottom_bc.append(face_rev)
    faces_empty = faces_bottom_bc+faces_top

    #####
    f = open(filename,'w')
    f.write(header)
    
    ##vertices
    write_vertices(f,vertices)

    ##blocks
    write_blocks(f,blocks,[20,20,1],[1,1,1])
    ##edges
    f.write('edges \n( \n); \n\n')
    ##boundary
    f.write('boundary \n( \n')
    #
    f.write('\tinlet \n\t{ \n')
    f.write('\t\ttype patch;\n')
    f.write('\t\tfaces \n\t\t( \n')
    write_faces(f,inlet)
    #
    f.write('\tp_outer \n\t{ \n')
    f.write('\t\ttype patch;\n')
    f.write('\t\tfaces \n\t\t( \n')
    write_faces(f,p_outer)
    #
    f.write('\twalls \n\t{ \n')
    f.write('\t\ttype wall;\n')
    f.write('\t\tfaces \n\t\t( \n')
    write_faces(f,walls)
    #
    f.write('\tjet-symm \n\t{ \n')
    f.write('\t\ttype patch;\n')
    f.write('\t\tfaces \n\t\t( \n')
    write_faces(f,symmetry)
    #
    f.write('\tfrontAndBack\n\t{ \n')
    f.write('\t\ttype empty;\n')
    f.write('\t\tfaces \n\t\t( \n')
    write_faces(f,faces_empty)
    #
    #end boundary
    f.write(');\n\n')

    f.write(footer)
    f.close()
#end generate_blockMeshDict

def build_PSLG():
    """
    build PlanarStraightLineGraph and return
    """
    v,s,vf,sf,btags = build_vertices_and_segments()

    domain = Domain.PlanarStraightLineGraphDomain(vertices=v,segments=s,vertexFlags=vf,segmentFlags=sf)
    domain.boundaryTags = btags
    
    return domain

if __name__ == "__main__":
    from proteus import Domain
    v,s,vf,sf,btags = build_vertices_and_segments()

    domain = Domain.PlanarStraightLineGraphDomain(vertices=v,segments=s,vertexFlags=vf,segmentFlags=sf)
    fileprefix='jet_domain'
    domain.writePoly(fileprefix)
    domain.writeGeo(fileprefix)

    generate_blockMeshDict_2d(width=1.0,filename='test_blockMeshDict_2d')
