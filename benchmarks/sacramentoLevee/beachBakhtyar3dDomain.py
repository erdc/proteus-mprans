import math
from proteus import Domain

def beachBakhtyar3d(L=[8.5,1.0,0.8],
                    Lb=3.5):
    boundaries=['left','right','bottom','top','front','back','obstacle']
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
    vertices=[[0.0,0.0,0.0],#0
              [Lb,0.0,0.0],#1
              [L[0],0.0,0.5],#2
              [L[0],0.0,L[2]],#3
              [0.0,0.0,L[2]],#4
              [0.0,L[1],0.0],#5
              [Lb,L[1],0.0],#6
              [L[0],L[1],0.5],#7
              [L[0],L[1],L[2]],#8
              [0.0,L[1],L[2]]]#9
    vertexFlags=[boundaryTags['left'],#0
                 boundaryTags['bottom'],#1
                 boundaryTags['right'],#2
                 boundaryTags['right'],#3
                 boundaryTags['left'],#4
                 boundaryTags['left'],#5
                 boundaryTags['bottom'],#6
                 boundaryTags['right'],#7
                 boundaryTags['right'],#8
                 boundaryTags['left']]#9
    facets=[[[0,1,2,3,4]],
            [[0,4,9,5]],
            [[2,3,8,7]],
            [[3,4,9,8]],
            [[0,1,6,5]],
            [[1,2,7,6]],
            [[5,6,7,8,9]]]
    facetFlags=[boundaryTags['front'],
                boundaryTags['left'],
                boundaryTags['right'],
                boundaryTags['top'],
                boundaryTags['bottom'],
                boundaryTags['bottom'],
                boundaryTags['back']]
    domain = Domain.PiecewiseLinearComplexDomain(vertices=vertices,
                                                 vertexFlags=vertexFlags,
                                                 facets=facets,
                                                 facetFlags=facetFlags)
    #go ahead and add a boundary tags member 
    domain.boundaryTags = boundaryTags
    return domain

if __name__=='__main__':
    import os
    domain = beachBakhtyar3d()
    #domain.writeAsymptote("beachBakhtyar3d")
    domain.writePoly("beachBakhtyar3d")
    domain.writePLY("beachBakhtyar3d")
    #os.system("asy -V beachBakhtyar3d")
