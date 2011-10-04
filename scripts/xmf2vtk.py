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
		    
(opts,args) = parser.parse_args()	

if not paraview.servermanager.ActiveConnection:
	connection = paraview.servermanager.Connect()

reader = servermanager.sources.XdmfReader(FileName=opts.filename)
#reader.UpdatePipeline()
#timesteps = reader.TimestepValues

i=0

#for time in timesteps:
#print "Time = " + str(time)

# create VTU writer and connect it to the reader
writer = servermanager.writers.XMLUnstructuredGridWriter(Input=reader,
							 FileName="sol"+str(i)+".vtu")

# Trigger execution of pipeline
writer.UpdatePipeline()

# 
i = i + 1 
  
  
