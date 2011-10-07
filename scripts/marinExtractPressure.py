#! /usr/bin/env pvpython

from paraview import servermanager
from paraview.simple import *
from optparse import OptionParser

usage = ""
parser = OptionParser(usage=usage)
parser.add_option("-f","--filebase",
                    help="filename",
                    action="store",
                    type="string",
                    dest="filename",
                    default="simulation")

parser.add_option("-e","--eps",
                    help="tolerance",
                    action="store",
                    type="float",
                    dest="eps",
                    default=0.005)
		    
(opts,args) = parser.parse_args()	

if not paraview.servermanager.ActiveConnection:
	connection = paraview.servermanager.Connect()

reader = servermanager.sources.XdmfReader(FileName=opts.filename)
reader.UpdatePipeline()
timesteps = reader.TimestepValues

points=[]

points.append(PointSource(Center=[2.3995-opts.eps ,0.5255, 0.021         ],NumberOfPoints=1))
points.append(PointSource(Center=[2.3995-opts.eps ,0.5255, 0.061         ],NumberOfPoints=1))
points.append(PointSource(Center=[2.4955          ,0.4745, 0.161+opts.eps],NumberOfPoints=1))
points.append(PointSource(Center=[2.5355          ,0.4745, 0.161+opts.eps],NumberOfPoints=1))

probes=[]
for point in points:
  probes.append(ProbePoint(Source=point,Input=reader))

outfile = open("pressure.txt",'w')
for time in timesteps:
  outfile.write(str(time))
  for probe in probes:
     probe.UpdatePipeline (time)

     fp = servermanager.Fetch(probe)
     pdata= fp.GetPointData()

     pressure = pdata.GetArray("p").GetTuple1(0)

     outfile.write("  " + str(pressure))
  outfile.write("\n")
outfile.close()
