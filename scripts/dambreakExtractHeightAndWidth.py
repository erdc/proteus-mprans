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

parser.add_option("-r","--resolution",
                    help="filename",
                    action="store",
                    type="int",
                    dest="resolution",
                    default=100)

parser.add_option("-a","--accuracy",
                    help="filename",
                    action="store",
                    type="float",
                    dest="accuracy",
                    default=0.01)

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

lines=[]
lines.append(LineSource(Point1=[opts.eps, 0.073, 0.0     ], Point2=[opts.eps, 0.073, 0.35    ], Resolution=opts.resolution))
lines.append(LineSource(Point1=[0.0,      0.073, opts.eps], Point2=[0.584,    0.073, opts.eps], Resolution=opts.resolution))

probes=[]
for line in lines:
  probes.append(ProbePoint(Source=line,Input=reader))
  
point=PointSource(Center=[0.0,0.5,0.0],NumberOfPoints=1)  
pprobe=ProbePoint(Source=point,Input=reader)

outfile = open("height.txt",'w')
for time in timesteps:
  print "Time =" + str(time)
  outfile.write(str(time))
  
  phi_old = 99.99
  height  = 0.0

  # Column height
  probes[0].UpdatePipeline (time)

  fp	= servermanager.Fetch(probes[0])
  pdata = fp.GetPointData()
  for i in  range(opts.resolution+1):
    phi = pdata.GetArray("phid").GetTuple1(i)
    print phi	     
    if (phi > 0.0) and (phi_old < 0.0): 	 
       height = 0.35* (float(i-1) + (phi_old/(phi_old-phi)))/float(opts.resolution)	       
    phi_old=phi
    
  if (height > 0.0): 
     phi = 111.0
     while abs(phi) > opts.accuracy: 
     	 point.Center=[opts.eps, 0.073, height]
  	 fp2 = servermanager.Fetch(pprobe)
  	 pdata2= fp2.GetPointData()
     	 phi  = pdata2.GetArray("phid").GetTuple1(0)	       
     	 height = height - phi
     	 print time,height, phi
  else:      
     height = 0.0    
  outfile.write("  " + str(height))
     
  
  phi_old = 99.99
  height  = 0.0

  # Column width
  probes[1].UpdatePipeline (time)

  fp	= servermanager.Fetch(probes[1])
  pdata = fp.GetPointData()
  for i in  range(opts.resolution+1):
    phi = pdata.GetArray("phid").GetTuple1(i)
    print phi	     
    if (phi > 0.0) and (phi_old < 0.0): 	 
       height = 0.584*(float(i-1) + (phi_old/(phi_old-phi)))/float(opts.resolution)	       
    phi_old=phi
    
  if (height > 0.0): 
     phi = 111.0
     while abs(phi) > opts.accuracy: 
     	 point.Center=[height,0.073, opts.eps]
  	 fp2 = servermanager.Fetch(pprobe)
  	 pdata2= fp2.GetPointData()
     	 phi  = pdata2.GetArray("phid").GetTuple1(0)	       
     	 height = height - phi
     	 print time,height, phi
  else:      
     height = 0.584  
  outfile.write("  " + str(height))
          
     
     
  outfile.write("\n")
outfile.close()
