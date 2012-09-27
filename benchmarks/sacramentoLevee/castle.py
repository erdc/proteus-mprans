from proteus import Domain

outline_left = [#(45.11,15.14),
    (15.14,15.14),
    (15.14,22.76),
    (18.13,24.82),
    (18.13,46.97),
    (15.14,49.03),
    (15.14,56.03),
    (19.05,56.03),
    (19.05,52.01),
    (22.14,52.01),
    (22.14,56.03),
    (26.06,56.03),
    (26.06,52.01),
    (29.15,52.01),
    (29.15,56.03),
    (33.06,56.03),
    (33.06,49.03),
    (29.97,46.97),
    (29.97,37.80),
    (38.01,37.80),
    (38.01,40.89),
    (36.15,42.02),
    (36.15,49.95),
    (40.07,49.95),
    (40.07,46.04),
    (43.15,46.04),
    (43.15,49.95),
    (45.11,49.95)]

outline_left_new = [(x-45.11,y-15.14) for x,y in outline_left]
outline = outline_left_new[:]
for x,y in outline_left_new[-2::-1]:
    outline.append((-x,y))

    #from matplotlib import pyplot

ll_window = [(22.04,21.94),(22.04,31.93),(26.06,31.93),(26.06,21.94)]
ul_window = [(x,y+(36.98-21.94)) for x,y in ll_window]
ml_window  = [(33.27,24.20),(33.27,30.69),(36.77,30.69),(36.69,24.20)]
door_left = [#(45.11,19.98),
             (40.17,19.98),(40.17,29.66),(41.09,30.8),(42.12,31.62),(43.04,31.93),(44.08,32.34),(45.11,33.06)]

ll_window_new = [(x-45.11,y-15.14) for x,y in ll_window]
ul_window_new = [(x-45.11,y-15.14) for x,y in ul_window]
ml_window_new = [(x-45.11,y-15.14) for x,y in ml_window]
door_left_new = [(x-45.11,y-15.14) for x,y in door_left]
door = door_left_new[:]
for x,y in door_left_new[-2::-1]:
    door.append((-x,y))
lr_window_new = [(-x,y) for x,y in ll_window_new]
ur_window_new = [(-x,y) for x,y in ul_window_new]
mr_window_new = [(-x,y) for x,y in ml_window_new]

#pyplot.plot([p[0] for p in outline],[p[1] for p in outline])
#pyplot.plot([p[0] for p in ll_window_new],[p[1] for p in ll_window_new])
#pyplot.plot([p[0] for p in ul_window_new],[p[1] for p in ul_window_new])
#pyplot.plot([p[0] for p in ml_window_new],[p[1] for p in ml_window_new])
#pyplot.plot([p[0] for p in mr_window_new],[p[1] for p in mr_window_new])
#pyplot.plot([p[0] for p in ur_window_new],[p[1] for p in ur_window_new])
#pyplot.plot([p[0] for p in lr_window_new],[p[1] for p in lr_window_new])
#pyplot.plot([p[0] for p in door],[p[1] for p in door])
#pyplot.show()


xmin = 0
xmax = 0
zmin = 0
zmax = 0
for x,z in outline:
    xmin = min(xmin,x)
    xmax = max(xmax,x)
    zmin = min(zmin,z)
    zmax = max(zmax,z)

Lx = 2*(xmax-xmin)
Ly = Lx
Lz = 1.25*(zmax-zmin)

boundaries=['left','right','bottom','top','front','back','obstacle']
boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
bt = boundaryTags
#tank
holes = []
vertices=[[-0.5*Lx-Lx,-0.5*Ly,zmin],#0
          [-0.5*Lx-Lx,0.5*Ly,zmin],#1
          [0.5*Lx,0.5*Ly,zmin],#2
          [0.5*Lx,-0.5*Ly,zmin],#3
          [-0.5*Lx-Lx,-0.5*Ly,zmin+Lz],#4
          [-0.5*Lx-Lx,0.5*Ly,zmin+Lz],#5
          [0.5*Lx,0.5*Ly,zmin+Lz],#6
          [0.5*Lx,-0.5*Ly,zmin+Lz]]#7
vertexFlags=[boundaryTags['left'],
             boundaryTags['left'],
             boundaryTags['right'],
             boundaryTags['right'],
             boundaryTags['left'],
             boundaryTags['left'],
             boundaryTags['right'],
             boundaryTags['right']]
facets=[[[0,1,2,3]],
        [[4,5,6,7]],
        [[0,1,5,4]],
        [[1,2,6,5]],
        [[2,3,7,6]],
        [[3,0,4,7]]]
facetHoles=[[],
            [],
            [],
            [],
            [],
            []]
facetFlags=[boundaryTags['bottom'],
            boundaryTags['top'],
            boundaryTags['left'],
            boundaryTags['back'],
            boundaryTags['right'],
            boundaryTags['front']]
#castle bottom
bl = outline[0]
br = outline[-1]
#castle front and back
front_facet = []
vn = 7
f_outline_start = vn+1
for p in outline:
    vertices.append([p[0],-0.3*Ly,p[1]])
    vertexFlags.append(boundaryTags['obstacle'])
    vn = vn+1
    front_facet.append(vn)
f_ll_window_start = vn+1
front_facet_ll_window = []
for p in ll_window_new:
    vertices.append([p[0],-0.3*Ly,p[1]])
    vertexFlags.append(boundaryTags['obstacle'])
    vn = vn+1
    front_facet_ll_window.append(vn)
f_ul_window_start = vn+1
front_facet_ul_window = []
for p in ul_window_new:
    vertices.append([p[0],-0.3*Ly,p[1]])
    vertexFlags.append(boundaryTags['obstacle'])
    vn = vn+1
    front_facet_ul_window.append(vn)
f_lr_window_start = vn+1
front_facet_lr_window = []
for p in lr_window_new:
    vertices.append([p[0],-0.3*Ly,p[1]])
    vertexFlags.append(boundaryTags['obstacle'])
    vn = vn+1
    front_facet_lr_window.append(vn)
f_ur_window_start = vn+1
front_facet_ur_window = []
for p in ur_window_new:
    vertices.append([p[0],-0.3*Ly,p[1]])
    vertexFlags.append(boundaryTags['obstacle'])
    vn = vn+1
    front_facet_ur_window.append(vn)
f_ml_window_start = vn+1
front_facet_ml_window = []
for p in ml_window_new:
    vertices.append([p[0],-0.3*Ly,p[1]])
    vertexFlags.append(boundaryTags['obstacle'])
    vn = vn+1
    front_facet_ml_window.append(vn)
f_mr_window_start = vn+1
front_facet_mr_window = []
for p in mr_window_new:
    vertices.append([p[0],-0.3*Ly,p[1]])
    vertexFlags.append(boundaryTags['obstacle'])
    vn = vn+1
    front_facet_mr_window.append(vn)
f_door_start = vn+1
front_facet_door = []
for p in door:
    vertices.append([p[0],-0.3*Ly,p[1]])
    vertexFlags.append(boundaryTags['obstacle'])
    vn = vn+1
    front_facet_door.append(vn)
facets.append([front_facet,
               front_facet_ll_window,
               front_facet_ul_window,
               front_facet_ml_window,
               front_facet_door,
               front_facet_mr_window,
               front_facet_lr_window,
               front_facet_ur_window])
facetFlags.append(boundaryTags['obstacle'])
facetHoles.append([])
#back
b_outline_start = vn+1
back_facet = []
for p in outline:
    vertices.append([p[0],0.3*Ly,p[1]])
    vertexFlags.append(boundaryTags['obstacle'])
    vn = vn+1
    back_facet.append(vn)
b_ll_window_start = vn+1
back_facet_ll_window = []
for p in ll_window_new:
    vertices.append([p[0],0.3*Ly,p[1]])
    vertexFlags.append(boundaryTags['obstacle'])
    vn = vn+1
    back_facet_ll_window.append(vn)
b_ul_window_start = vn+1
back_facet_ul_window = []
for p in ul_window_new:
    vertices.append([p[0],0.3*Ly,p[1]])
    vertexFlags.append(boundaryTags['obstacle'])
    vn = vn+1
    back_facet_ul_window.append(vn)
b_lr_window_start = vn+1
back_facet_lr_window = []
for p in lr_window_new:
    vertices.append([p[0],0.3*Ly,p[1]])
    vertexFlags.append(boundaryTags['obstacle'])
    vn = vn+1
    back_facet_lr_window.append(vn)
b_ur_window_start = vn+1
back_facet_ur_window = []
for p in ur_window_new:
    vertices.append([p[0],0.3*Ly,p[1]])
    vertexFlags.append(boundaryTags['obstacle'])
    vn = vn+1
    back_facet_ur_window.append(vn)
b_ml_window_start = vn+1
back_facet_ml_window = []
for p in ml_window_new:
    vertices.append([p[0],0.3*Ly,p[1]])
    vertexFlags.append(boundaryTags['obstacle'])
    vn = vn+1
    back_facet_ml_window.append(vn)
b_mr_window_start = vn+1
back_facet_mr_window = []
for p in mr_window_new:
    vertices.append([p[0],0.3*Ly,p[1]])
    vertexFlags.append(boundaryTags['obstacle'])
    vn = vn+1
    back_facet_mr_window.append(vn)
b_door_start = vn+1
back_facet_door = []
for p in door:
    vertices.append([p[0],0.3*Ly,p[1]])
    vertexFlags.append(boundaryTags['obstacle'])
    vn = vn+1
    back_facet_door.append(vn)
facets.append([back_facet,
               back_facet_ll_window,
               back_facet_ul_window,
               back_facet_ml_window,
               back_facet_door,
               back_facet_mr_window,
               back_facet_lr_window,
               back_facet_ur_window])
facetFlags.append(boundaryTags['obstacle'])
facetHoles.append([])
#cross
for fN,bN in zip(range(f_outline_start,f_outline_start+len(outline)-1),range(b_outline_start,b_outline_start+len(outline)-1)):
    facets.append([[fN,bN,bN+1,fN+1]])
    facetFlags.append(boundaryTags['obstacle'])
    facetHoles.append([])
for fN,bN in zip(range(f_ll_window_start,f_ll_window_start+len(ll_window)-1),range(b_ll_window_start,b_ll_window_start+len(ll_window)-1)):
    facets.append([[fN,bN,bN+1,fN+1]])
    facetFlags.append(boundaryTags['obstacle'])
    facetHoles.append([])
facets.append([[bN+1,fN+1,f_ll_window_start,b_ll_window_start]])
facetFlags.append(boundaryTags['obstacle'])
facetHoles.append([])
for fN,bN in zip(range(f_ul_window_start,f_ul_window_start+len(ul_window)-1),range(b_ul_window_start,b_ul_window_start+len(ul_window)-1)):
    facets.append([[fN,bN,bN+1,fN+1]])
    facetFlags.append(boundaryTags['obstacle'])
    facetHoles.append([])
facets.append([[bN+1,fN+1,f_ul_window_start,b_ul_window_start]])
facetFlags.append(boundaryTags['obstacle'])
facetHoles.append([])
for fN,bN in zip(range(f_ml_window_start,f_ml_window_start+len(ml_window)-1),range(b_ml_window_start,b_ml_window_start+len(ml_window)-1)):
    facets.append([[fN,bN,bN+1,fN+1]])
    facetFlags.append(boundaryTags['obstacle'])
    facetHoles.append([])
facets.append([[bN+1,fN+1,f_ml_window_start,b_ml_window_start]])
facetFlags.append(boundaryTags['obstacle'])
facetHoles.append([])
for fN,bN in zip(range(f_mr_window_start,f_mr_window_start+len(mr_window_new)-1),range(b_mr_window_start,b_mr_window_start+len(mr_window_new)-1)):
    facets.append([[fN,bN,bN+1,fN+1]])
    facetFlags.append(boundaryTags['obstacle'])
    facetHoles.append([])
facets.append([[bN+1,fN+1,f_mr_window_start,b_mr_window_start]])
facetFlags.append(boundaryTags['obstacle'])
facetHoles.append([])
for fN,bN in zip(range(f_lr_window_start,f_lr_window_start+len(lr_window_new)-1),range(b_lr_window_start,b_lr_window_start+len(lr_window_new)-1)):
    facets.append([[fN,bN,bN+1,fN+1]])
    facetFlags.append(boundaryTags['obstacle'])
    facetHoles.append([])
facets.append([[bN+1,fN+1,f_lr_window_start,b_lr_window_start]])
facetFlags.append(boundaryTags['obstacle'])
facetHoles.append([])
for fN,bN in zip(range(f_ur_window_start,f_ur_window_start+len(ur_window_new)-1),range(b_ur_window_start,b_ur_window_start+len(ur_window_new)-1)):
    facets.append([[fN,bN,bN+1,fN+1]])
    facetFlags.append(boundaryTags['obstacle'])
    facetHoles.append([])
facets.append([[bN+1,fN+1,f_ur_window_start,b_ur_window_start]])
facetFlags.append(boundaryTags['obstacle'])
facetHoles.append([])
for fN,bN in zip(range(f_door_start,f_door_start+len(door)-1),range(b_door_start,b_door_start+len(door)-1)):
    facets.append([[fN,bN,bN+1,fN+1]])
    facetFlags.append(boundaryTags['obstacle'])
    facetHoles.append([])
facets.append([[bN+1,fN+1,f_door_start,b_door_start]])
facetFlags.append(boundaryTags['obstacle'])
facetHoles.append([])
holes.append([vertices[8][0]+1.0e-3,vertices[8][1]+1.0e-3,vertices[8][2]+1.0e-3])
print len(facetHoles),len(facetFlags),len(facets)
print facetHoles
facets[0].append([f_outline_start,b_outline_start,b_ll_window_start-1,f_ll_window_start-1])
facetHoles[0].append((vertices[f_outline_start][0]+1.0e-3,
                      vertices[f_outline_start][1]+1.0e-3,
                      zmin))
domain = Domain.PiecewiseLinearComplexDomain(vertices=vertices,
                                             vertexFlags=vertexFlags,
                                             facets=facets,
                                             facetFlags=facetFlags,
                                             facetHoles=facetHoles,
                                             holes=holes)
#domain.writePoly("castle")
#domain.writePLY("castle")
