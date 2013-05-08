try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

#wigley_all4_xmf = XDMFReader( FileName='/Users/cekees/mprans/benchmarks/wigley/debug_4pe/wigley_all4.xmf' )

#AnimationScene1 = GetAnimationScene()
#wigley_all4_xmf.PointArrays = ['nodeMaterialTypes','phi']

#AnimationScene1.EndTime = 6.2603025177899996
#AnimationScene1.PlayMode = 'Snap To TimeSteps'

RenderView1 = GetRenderView()
# a1_NodeMapL2G_PVLookupTable = GetLookupTableForArray( "NodeMapL2G", 1, RGBPoints=[0.0, 0.0, 0.0, 1.0, 10849.0, 1.0, 0.0, 0.0], VectorMode='Magnitude', NanColor=[0.49803921568627502, 0.49803921568627502, 0.49803921568627502], ScalarOpacityFunction=[], ColorSpace='HSV', ScalarRangeInitialized=1.0, AllowDuplicateScalars=1, Annotations=['2.1390766798348069 0.019832821328072146 -0.045358380840157572', '2.1390766798348069 0.019832821328072146 -0.045358380840157572', '2.1296637604868858 0.01145792042302154 -0.027697776661718167', '2.1296637604868858 0.01145792042302154 -0.027697776661718167', '2.126570461187808 -0.015070373011410895 -0.019931675178417988', '2.126570461187808 -0.015070373011410895 -0.019931675178417988', '2.1206453220756747 0.037365161731833234 -0.02989327769398073', '2.1206453220756747 0.037365161731833234 -0.02989327769398073', '2.1189493128396233 -0.081537783383271362 0.073401323404416033', '2.1189493128396233 -0.081537783383271362 0.073401323404416033', '2.1134575016305939 0.0136999031006529 -0.0069069716848219921', '2.1134575016305939 0.0136999031006529 -0.0069069716848219921', '2.1095600431708448 -0.080725114211490229 -0.037390823198444077', '2.1095600431708448 -0.080725114211490229 -0.037390823198444077', '2.1073859384208999 -0.01123890133428845 0.0032381956681317306', '2.1073859384208999 -0.01123890133428845 0.0032381956681317306', '2.1027959246473777 -0.0094526635782381952 0.0023222551003503128', '2.1027959246473777 -0.0094526635782381952 0.0023222551003503128', '2.0980602170603895 0.00026811714991760265 -5.8104048155014279e-06', '2.0980602170603895 0.00026811714991760265 -5.8104048155014279e-06', '2.0980117818974491 -0.00043556199942712179 0.00021900470031305423', '2.0980117818974491 -0.00043556199942712179 0.00021900470031305423', '2.0979089790705414 -0.00021044555330085882 0.00011997670691345165', '2.0979089790705414 -0.00021044555330085882 0.00011997670691345165', '2.0978012786625624 -9.8152368748262778e-05 8.7332086881109359e-05', '2.0978012786625624 -9.8152368748262778e-05 8.7332086881109359e-05', '2.0975370915829084 -9.6023415144978627e-05 1.4989249898214919e-05', '2.0975370915829084 -9.6023415144978627e-05 1.4989249898214919e-05', '2.0975152722046126 0.00046440848773055736 0.00013800935336849598', '2.0975152722046126 0.00046440848773055736 0.00013800935336849598', '2.0974798056875557 0.00032923845172051888 -7.6371479679501568e-05', '2.0974798056875557 0.00032923845172051888 -7.6371479679501568e-05', '2.0974612304546989 -0.00043882876398771228 0.00011915709376935022', '2.0974612304546989 -0.00043882876398771228 0.00011915709376935022', '2.0974448455203776 0 0', '2.0974448455203776 0 0', '2.0973766268958438 -6.6880795038571757e-05 -1.538720721463146e-07', '2.0973766268958438 -6.6880795038571757e-05 -1.538720721463146e-07', '2.0973260866109062 -0.00050745132467635201 0.00022404414153184165', '2.0973260866109062 -0.00050745132467635201 0.00022404414153184165', '2.097312382183282 0.0001781389750294839 -7.4141951191850362e-05', '2.097312382183282 0.0001781389750294839 -7.4141951191850362e-05', '2.0973122090868874 0.00017246793888824768 -6.7524157154809957e-05', '2.0973122090868874 0.00017246793888824768 -6.7524157154809957e-05', '2.0972884121709341 -0.00054101388336901695 0.00022449936685750346', '2.0972884121709341 -0.00054101388336901695 0.00022449936685750346', '2.0971529745473223 -0.00035407028128771009 -0.0002629667660479576', '2.0971529745473223 -0.00035407028128771009 -0.0002629667660479576', '2.0955606534603222 0.0050896381321391126 0.00089255132896870753', '2.0955606534603222 0.0050896381321391126 0.00089255132896870753', '2.0945511822832708 -0.00046601306145391578 -0.00062413067202789858', '2.0945511822832708 -0.00046601306145391578 -0.00062413067202789858', '2.0933754222401375 -0.0044940766587617366 0.00071013057158392885', '2.0933754222401375 -0.0044940766587617366 0.00071013057158392885', '2.0924729676601022 -0.00031526966938574916 -0.00087922026771666644', '2.0924729676601022 -0.00031526966938574916 -0.00087922026771666644', '2.0923896043975598 -0.0045399114511279594 0.00078401982616313698', '2.0923896043975598 -0.0045399114511279594 0.00078401982616313698', '2.0901127575644294 0.026888515488696967 -0.00063866913260711411', '2.0901127575644294 0.026888515488696967 -0.00063866913260711411', '2.070639316980881 -0.063865977140751345 -0.025763060062894548', '2.070639316980881 -0.063865977140751345 -0.025763060062894548', '2.0619751938824713 -0.043504821196674634 0.035398496507998242', '2.0619751938824713 -0.043504821196674634 0.035398496507998242'] )

# a1_NodeMapL2G_PiecewiseFunction = CreatePiecewiseFunction( Points=[2.0112321896875098, 0.0, 0.5, 0.0, 2.1540010523897899, 1.0, 0.5, 0.0] )

# DataRepresentation13 = Show()
# DataRepresentation13.EdgeColor = [0.0, 0.0, 0.50000762951094835]
# DataRepresentation13.SelectionPointFieldDataArrayName = 'NodeMapL2G'
# DataRepresentation13.SelectionCellFieldDataArrayName = 'CellMapL2G'
# DataRepresentation13.ScalarOpacityFunction = []
# DataRepresentation13.ColorArrayName = 'NodeMapL2G'
# DataRepresentation13.ScalarOpacityUnitDistance = 0.041662962906450239
# DataRepresentation13.LookupTable = a1_NodeMapL2G_PVLookupTable
# DataRepresentation13.ExtractedBlockIndex = 1
# DataRepresentation13.ScaleFactor = 0.15000000000000002

#wigley_all4_xmf.Grids = ['Mesh Spatial_Domain']

RenderView1.CenterOfRotation = [0.0, 0.0, 0.09375]

Transform3 = Transform( Transform="Transform" )

RenderView1.CameraPosition = [0.0, -2.9772812767548129, 0.09375]
RenderView1.CameraFocalPoint = [0.0, 0.0, 0.09375]
RenderView1.CameraClippingRange = [2.6490084639872649, 3.3976904959061351]
RenderView1.CameraParallelScale = 0.7705770970512944

Transform3.Transform = "Transform"

active_objects.source.SMProxy.InvokeEvent('UserEvent', 'ShowWidget')


# DataRepresentation14 = Show()
# DataRepresentation14.EdgeColor = [0.0, 0.0, 0.50000762951094835]
# DataRepresentation14.SelectionPointFieldDataArrayName = 'NodeMapL2G'
# DataRepresentation14.SelectionCellFieldDataArrayName = 'CellMapL2G'
# DataRepresentation14.ScalarOpacityFunction = []
# DataRepresentation14.ColorArrayName = 'NodeMapL2G'
# DataRepresentation14.ScalarOpacityUnitDistance = 0.041662962969182267
# DataRepresentation14.LookupTable = a1_NodeMapL2G_PVLookupTable
# DataRepresentation14.ExtractedBlockIndex = 1
# DataRepresentation14.ScaleFactor = 0.15000000000000002

# DataRepresentation13.Visibility = 0

Transform3.Transform.Translate = [0.0, 0.0, -0.9375]

RenderView1.CameraClippingRange = [2.6490084521259405, 3.3976905108370987]

Transform3.Transform.Translate = [0.0, 0.0, -0.09375]

MergeBlocks2 = MergeBlocks()

# DataRepresentation15 = Show()
# DataRepresentation15.EdgeColor = [0.0, 0.0, 0.50000762951094835]
# DataRepresentation15.SelectionPointFieldDataArrayName = 'NodeMapL2G'
# DataRepresentation15.SelectionCellFieldDataArrayName = 'CellMapL2G'
# DataRepresentation15.ScalarOpacityFunction = []
# DataRepresentation15.ColorArrayName = 'NodeMapL2G'
# DataRepresentation15.ScalarOpacityUnitDistance = 0.041662962969182267
# DataRepresentation15.LookupTable = a1_NodeMapL2G_PVLookupTable
# DataRepresentation15.ExtractedBlockIndex = 0
# DataRepresentation15.ScaleFactor = 0.15000000000000002

# DataRepresentation14.Visibility = 0

Clip3 = Clip( ClipType="Scalar" )

Clip3.Scalars = ['POINTS', 'nodeMaterialTypes']
Clip3.ClipType = "Scalar"
Clip3.Value = 7.0

# DataRepresentation16 = Show()
# DataRepresentation16.EdgeColor = [0.0, 0.0, 0.50000762951094835]
# DataRepresentation16.SelectionPointFieldDataArrayName = 'NodeMapL2G'
# DataRepresentation16.SelectionCellFieldDataArrayName = 'CellMapL2G'
# DataRepresentation16.ColorArrayName = 'NodeMapL2G'
# DataRepresentation16.ScalarOpacityUnitDistance = 0.040349781414118609
# DataRepresentation16.ScaleFactor = 0.10000000000000001

# DataRepresentation15.Visibility = 0

# DataRepresentation16.ScalarOpacityFunction = []
# DataRepresentation16.LookupTable = a1_NodeMapL2G_PVLookupTable

ExtractSurface3 = ExtractSurface()

RenderView1.CameraClippingRange = [2.8480084625045992, 3.147190497772506]

# DataRepresentation17 = Show()
# DataRepresentation17.EdgeColor = [0.0, 0.0, 0.50000762951094835]
# DataRepresentation17.SelectionPointFieldDataArrayName = 'NodeMapL2G'
# DataRepresentation17.SelectionCellFieldDataArrayName = 'CellMapL2G'
# DataRepresentation17.ColorArrayName = 'NodeMapL2G'
# DataRepresentation17.LookupTable = a1_NodeMapL2G_PVLookupTable
# DataRepresentation17.ScaleFactor = 0.10000000000000001

# DataRepresentation16.Visibility = 0

Contour4 = Contour( PointMergeMethod="Uniform Binning" )

# Contour4.PointMergeMethod = "Uniform Binning"
# Contour4.ContourBy = ['POINTS', 'NodeMapL2G']
# Contour4.Isosurfaces = [4683.5]

Contour4.ContourBy = ['POINTS', 'phi']
Contour4.ComputeNormals = 0
Contour4.Isosurfaces = [0.0]

# DataRepresentation18 = Show()
# DataRepresentation18.EdgeColor = [0.0, 0.0, 0.50000762951094835]
# DataRepresentation18.SelectionPointFieldDataArrayName = 'nodeMaterialTypes'
# DataRepresentation18.SelectionCellFieldDataArrayName = 'CellMapL2G'
# DataRepresentation18.ColorArrayName = 'NodeMapL2G'
# DataRepresentation18.ScaleFactor = 0.10000000000000001

# DataRepresentation17.Visibility = 0

# DataRepresentation18.ColorArrayName = ''

#AnimationScene1.AnimationTime = 6.2603025177899996

#RenderView1.ViewTime = 6.2603025177899996
#RenderView1.CacheKey = 6.2603025177899996
#RenderView1.UseCache = 1

Transform4 = Transform( Transform="Transform" )

RenderView1.UseCache = 0

Transform4.Transform = "Transform"

RenderView1.CameraClippingRange = [2.8504913999672272, 3.1441450536378825]

#active_objects.source.SMProxy.InvokeEvent('UserEvent', 'ShowWidget')


DataRepresentation19 = Show()
DataRepresentation19.ScaleFactor = 0.10000000000000001
DataRepresentation19.EdgeColor = [0.0, 0.0, 0.50000762951094835]
DataRepresentation19.SelectionPointFieldDataArrayName = 'nodeMaterialTypes'
#DataRepresentation19.SelectionCellFieldDataArrayName = 'CellMapL2G'

# DataRepresentation18.Visibility = 0

Transform4.Transform.Scale = [1.0, 1.0, 10.0]

RenderView1.CameraClippingRange = [2.8504913999672272, 3.1441450536378825]

DataRepresentation19.CustomBounds = [-0.5, 0.5, -0.0487345, 0.048805300000000003, -0.101045, 0.350885]
DataRepresentation19.CubeAxesYAxisTickVisibility = 0
DataRepresentation19.CubeAxesYAxisMinorTickVisibility = 0
DataRepresentation19.CubeAxesUseDefaultYTitle = 0
DataRepresentation19.CubeAxesColor = [0.0, 0.0, 0.0]
DataRepresentation19.CustomRange = [-0.5, 0.5, -0.0487345, 0.048805300000000003, -0.101045, 0.350885]
DataRepresentation19.CubeAxesYAxisVisibility = 0
DataRepresentation19.CubeAxesVisibility = 1

#DataRepresentation18.DiffuseColor = [0.0, 0.0, 0.0]

RenderView1.CameraPosition = [0.0, -1.3889236885663774, 0.09375]
RenderView1.Background2 = [0.0, 0.0, 0.16470588235294117]
RenderView1.Background = [1.0, 1.0, 1.0]
RenderView1.CameraClippingRange = [1.278017387660676, 1.5319621016266205]

DataRepresentation19.LineWidth = 3.0
DataRepresentation19.DiffuseColor = [0.0, 0.0, 0.0]

RenderView1.CameraParallelScale = 0.35947990279502062
RenderView1.InteractionMode = '2D'

Render()
