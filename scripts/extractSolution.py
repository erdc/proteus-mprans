#! /usr/bin/env python

import os
import tables
import sys

def splitH5(NavStoBaseName,PhiBaseName,CorrBaseName,size,start,finaltime,stride):

    for proc in range(0,size):
    	     print "Processor", proc
	     # Loop over entries and put in appropriate file
     	     filename = NavStoBaseName+str(proc)+".h5"
	     print " Open:",filename
     	     nsFile = tables.openFile(filename)  

       	     filename = PhiBaseName+str(proc)+".h5"
	     print " Open:",filename
     	     phiFile = tables.openFile(filename) 
	     
       	     filename = CorrBaseName+str(proc)+".h5"
	     print " Open:",filename
     	     corrFile = tables.openFile(filename) 	     
	     
	     print "   Step:",
	     
     	     for step in range(start,finaltime+1,stride):
                print  step,
		sys.stdout.flush()
		
      	    	filename="sol.p"+str(proc)+"."+str(step)+".h5"
     	     	hdfFile = tables.openFile(filename,
				    mode = "w",
				    title = filename+" Data")
		
		if nsFile.__contains__("/elements_c0q2_Lagrange"+str(step)):
		        name =  "elements_c0q2_Lagrange"+str(step)
		elif nsFile.__contains__("/elements_c0p2_Lagrange"+str(step)):
		        name =  "elements_c0p2_Lagrange"+str(step)
		elif nsFile.__contains__("/elementsSpatial_Domain"+str(step)):
		        name =  "elementsSpatial_Domain"+str(step)
		else:
			print nsFile.__contains__("/elements_c0q2_Lagrange"+str(step))
			print nsFile.__contains__("/elements_c0p2_Lagrange"+str(step))
			print nsFile.__contains__("/elementsSpatial_Domain"+str(step))
			sys.exit("No element connectivity")
					    
	        hdfFile.createArray("/","elements",nsFile.getNode("/",name)[:])			

		if nsFile.__contains__("/nodes_c0q2_Lagrange"+str(step)):
		        name =  "nodes_c0q2_Lagrange"+str(step)
		elif nsFile.__contains__("/nodes_c0p2_Lagrange"+str(step)):
		        name =  "nodes_c0p2_Lagrange"+str(step)
		elif nsFile.__contains__("/nodesSpatial_Domain"+str(step)):
		        name =  "nodesSpatial_Domain"+str(step)
		else:
			print nsFile.__contains__("/nodes_c0q2_Lagrange"+str(step))
			print nsFile.__contains__("/nodes_c0p2_Lagrange"+str(step))
			print nsFile.__contains__("/nodesSpatial_Domain"+str(step))
			sys.exit("No node coordinates")
     		
		hdfFile.createArray("/","nodes",nsFile.getNode("/",name)[:]) 

     		name =  "u"+str(step)
     		hdfFile.createArray("/","u",nsFile.getNode("/",name)[:])

     		name =  "v"+str(step)
     		hdfFile.createArray("/","v",nsFile.getNode("/",name)[:])
		     
     		name =  "w"+str(step)
     		hdfFile.createArray("/","w",nsFile.getNode("/",name)[:])
		     
     		name =  "p"+str(step)
     		hdfFile.createArray("/","p",nsFile.getNode("/",name)[:])
		
     		name =  "phid"+str(step)
     		hdfFile.createArray("/","phi",phiFile.getNode("/",name)[:])	
						
     		name =  "phiCorr"+str(step)
     		hdfFile.createArray("/","phic",corrFile.getNode("/",name)[:])	
		
                hdfFile.close()

             nsFile.close()	            	 
             phiFile.close()			    			     
             corrFile.close()
	     
	     print "finished"
		     	    
def H5toXMF(basename,size,start,finaltime,stride):

# Open XMF files

     t1="  "
     t2=t1+t1
     t3=t2+t1
     t4=t3+t1
     t5=t4+t1

     XMFfile1 = open(basename+".xmf","w")
     XMFfile1.write('<?xml version="1.0" ?>'+"\n")
     XMFfile1.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'+"\n")
     XMFfile1.write('<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">'+"\n")
     XMFfile1.write(t1 + '<Domain>'+"\n")
     XMFfile1.write(t1 + '<Grid GridType="Collection"   CollectionType="Temporal">'+"\n")
	
     string="" 
     print "   Step:",
     for step in range(start,finaltime+1,stride):
                print step,		
		sys.stdout.flush()
		
                filename = basename+"."+str(step)+".h5"
     	     	hdfFile=  tables.openFile(filename,
				    mode = "w",
				    title = filename+" Data")     
     
                XMFfile2 = open(basename+"."+str(step)+".xmf","w")
                XMFfile2.write('<?xml version="1.0" ?>'+"\n")
                XMFfile2.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'+"\n")
                XMFfile2.write('<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">'+"\n")
                XMFfile2.write(t1 + '<Domain>'+"\n")
	
		string = t2 + '<Grid GridType="Collection" CollectionType="Spatial">'+"\n"	
	        string = string + t3 + '<Time Value="'+str(step)+'" />'+"\n"

                for proc in range(0,size):
                        group = hdfFile.createGroup(hdfFile.root, 'p'+str(proc))
			
                        solname="sol.p"+str(proc)+"."+str(step)+".h5"
     	        	f1 = tables.openFile(solname)  
			
		        string = string + t3+'<Grid GridType="Uniform">'+"\n"	

     			data=f1.getNode("/","elements")[:]
     			hdfFile.createArray(group,"elements",data)
		
                  	string = string + t4 + '<Topology NumberOfElements="' +str(len(data))+ '" Type="Tetrahedron">'+"\n"
         		string = string + t5 + '<DataItem DataType="Int" Dimensions="' +str(len(data))+ ' 4" Format="HDF">'+"\n"
			string = string + t5 + filename + ':/p'+str(proc)+'/elements'+"\n"
			string = string + t5 +'</DataItem>'+"\n"
        		string = string + t4 + '</Topology>'+"\n"

     			data=f1.getNode("/","nodes")[:]
     			hdfFile.createArray(group,"nodes",data)
					     
                   	string = string + t4 + '<Geometry Type="XYZ">'+"\n"
                  	string = string + t5 + '<DataItem DataType="Float" Dimensions="' +str(len(data))+ ' 3" Format="HDF" Precision="8">' + "\n"
			string = string + t5 + filename + ':/p'+str(proc)+'/nodes'+"\n"
			string = string + t5 + '</DataItem>'+"\n"
                        string = string + t4 + '</Geometry>'+"\n"

     			data=f1.getNode("/","u")[:]
     			hdfFile.createArray(group,"u",data)
			  			
                  	string = string + t4 + '<Attribute AttributeType="Scalar" Center="Node" Name="u">'+"\n"
                  	string = string + t5 + '<DataItem DataType="Float" Dimensions="' +str(len(data))+ '" Format="HDF" Precision="8">' + "\n"
			string = string + t5 + filename + ':/p'+str(proc)+'/u'+"\n"
			string = string + t5 + '</DataItem>'+"\n"
                  	string = string + t4 + '</Attribute>'+"\n"

     			data=f1.getNode("/","v")[:]
     			hdfFile.createArray(group,"v",data)
			
     			string = string + t4 + '<Attribute AttributeType="Scalar" Center="Node" Name="v">'+"\n"
                  	string = string + t5 +'<DataItem DataType="Float" Dimensions="' +str(len(data))+ '" Format="HDF" Precision="8">' + "\n"
			string = string + t5 + filename + ':/p'+str(proc)+'/v'+"\n"
			string = string + t5 + '</DataItem>'+"\n"
                  	string = string + t4 + '</Attribute>'+"\n"

     			data=f1.getNode("/","w")[:]
     			hdfFile.createArray(group,"w",data)
						
                  	string = string + t4 + '<Attribute AttributeType="Scalar" Center="Node" Name="w">'+"\n"
                  	string = string + t5 + '<DataItem DataType="Float" Dimensions="' +str(len(data))+ '" Format="HDF" Precision="8">' + "\n"
			string = string + t5 + filename + ':/p'+str(proc)+'/w'+"\n"
			string = string + t5 + '</DataItem>'+"\n"
                  	string = string + t4 + '</Attribute>'+"\n"

     			data=f1.getNode("/","p")[:]
     			hdfFile.createArray(group,"p",data)
						
                  	string = string + t4 + '<Attribute AttributeType="Scalar" Center="Node" Name="p">'+"\n"
                  	string = string + t5 + '<DataItem DataType="Float" Dimensions="' +str(len(data))+ '" Format="HDF" Precision="8">' + "\n"
			string = string + t5 + filename + ':/p'+str(proc)+'/p'+"\n"
			string = string + t5 + '</DataItem>'+"\n"
                  	string = string + t4 + '</Attribute>'+"\n"

     			data=f1.getNode("/","phi")[:]
     			hdfFile.createArray(group,"phi",data)
			
                  	string = string + t4 + '<Attribute AttributeType="Scalar" Center="Node" Name="phid">'+"\n"
                        string = string + t5 + '<DataItem DataType="Float" Dimensions="' +str(len(data))+ '" Format="HDF" Precision="8">' + "\n"
			string = string + t5 + filename + ':/p'+str(proc)+'/phid'+"\n"
			#string = string + t5 + '</DataItem>'+"\n"
                  	string = string + t4 + '</Attribute>'+"\n"
			
			string = string + t3+'</Grid>'+"\n"
			     
	     		f1.close()
			os.remove(solname)

        	string = string + t2 + '</Grid>'+"\n"  
		     
                XMFfile1.write(string)	
		
                XMFfile2.write(string)			
		XMFfile2.write(t1 + '</Domain>'+"\n")		
                XMFfile2.write('</Xdmf>'+"\n")		
		XMFfile2.close()
		
		hdfFile.close()
		
		
		
     XMFfile1.write(t1 + '</Grid>'+"\n")			
     XMFfile1.write(t1 + '</Domain>'+"\n")			
     XMFfile1.write('</Xdmf>'+"\n")		
     XMFfile1.close()

def H5toVTK(basename,size,start,finaltime,stride):

     if not paraview.servermanager.ActiveConnection:
		connection = paraview.servermanager.Connect()
	
     for step in range(start,finaltime+1,stride):
                print step,		
		sys.stdout.flush()
		
                reader = servermanager.sources.XdmfReader(basename+"."+str(step)+".xmf")
     		writer = servermanager.writers.XMLUnstructuredGridWriter(Input=reader,
							 FileName=basename +"."+str(step) + ".vtu")

		writer.UpdatePipeline()


def H5toVTK_old(basename,size,start,finaltime,stride):

# Open XMF files

     t1="  "
     t2=t1+t1
     t3=t2+t1
     t4=t3+t1
     t5=t4+t1


#	
     nodes = 0
     elems = 0	
     
     offsets=range(size+1)
     
     for proc in range(0,size):
     	     solname="sol.p"+str(proc)+"."+str(start)+".h5"
     	     f1 = tables.openFile(solname)   
	     data=f1.getNode("/","elements")[:]
	     elems = elems + len(data)
	     nshl  = len(data[0])
	     
	     data=f1.getNode("/","nodes")[:]
	     nodes = nodes + len(data)	
	     offsets[proc+1] = nodes
	     f1.close()
	          
     print "   Nodes  = " + str(nodes)
     print "   Elems  = " + str(elems)     	
     print "   Nshl   = " + str( nshl)    	

     if (nshl==4):
     	etype = "10\n"
     elif (nshl == 8):
     	etype = "12\n"
     elif (nshl == 10):
     	etype = "24\n"
     elif (nshl == 27):
     	etype = "27\n"
     else:	
	sys.exit("No element type for NHSL: "+ str(nshl))
	
     string="" 
     print "   Step:",
     
     try: 
     	os.mkdir(basename)
     except:
     	print "Directory already exists"
    
     for step in range(start,finaltime+1,stride):
                print step,		
		sys.stdout.flush()
		
                filename = "solution."+str(step)+".vtk"
     	     	vtkFile =  open(os.path.join(basename,filename),"w")     

     	     	vtkFile.write('# vtk DataFile Version 3.0\n')
     	     	vtkFile.write('vtk output\n')
     	     	vtkFile.write('ASCII\n')
     	     	vtkFile.write('DATASET UNSTRUCTURED_GRID\n')

     	     	vtkFile.write('POINTS ' + str(nodes) +  ' double\n')
                for proc in range(0,size):
      
                        solname="sol.p"+str(proc)+"."+str(step)+".h5"
     	        	f1 = tables.openFile(solname)  

     			data=f1.getNode("/","nodes")			
			for d in data:
				vtkFile.write(str(d[0]) + " " + str(d[1])  +  " "   + str(d[2]) +  "\n" )
			f1.close()
			
     	     	vtkFile.write('CELLS '+ str(elems) + " " + str(elems*(nshl+1)) + "\n")			     
                for proc in range(0,size):
      
                        solname="sol.p"+str(proc)+"."+str(step)+".h5"
     	        	f1 = tables.openFile(solname)  
			
     			data=f1.getNode("/","elements")
			for d in data:
			        vtkFile.write(str(nshl) + " ")
				for i in range(nshl):
					vtkFile.write(str(d[i]+offsets[proc]) + " ")
				vtkFile.write( "\n" )
					
			f1.close()
						
     	     	vtkFile.write('CELL_TYPES '+ str(elems) + "\n")
		for proc in range(0,elems):
		        vtkFile.write(etype)

     	     	vtkFile.write('POINT_DATA '+ str(nodes) + "\n")	
		
     	     	vtkFile.write('VECTORS u double \n')					     
                for proc in range(0,size):
      
                        solname="sol.p"+str(proc)+"."+str(step)+".h5"
     	        	f1 = tables.openFile(solname)  	
								
     			udata=f1.getNode("/","u")			  			
     			vdata=f1.getNode("/","v")
			wdata=f1.getNode("/","w")
			
			for i in range(0,len(udata)):
				vtkFile.write(str(udata[i]) + " " + str(vdata[i])  +  " "   + str(wdata[i]) +  "\n" )
			       			     
	     		f1.close()

     	     	vtkFile.write('SCALARS p double \n')
     	     	vtkFile.write('LOOKUP_TABLE default \n')
                for proc in range(0,size):
      
                        solname="sol.p"+str(proc)+"."+str(step)+".h5"
     	        	f1 = tables.openFile(solname)  	
								
     			data=f1.getNode("/","p")
			for d in data:				
				vtkFile.write(str(d) + "\n" )		
			       			     
	     		f1.close()

     	     	vtkFile.write('SCALARS phi double \n')
     	     	vtkFile.write('LOOKUP_TABLE default \n')
                for proc in range(0,size):
      
                        solname="sol.p"+str(proc)+"."+str(step)+".h5"
     	        	f1 = tables.openFile(solname)  	
								
     			data=f1.getNode("/","phi")
     			data2=f1.getNode("/","phic")
			for d in data:				
				vtkFile.write(str(d+data2[i]) + "\n" )		
			       			     
	     		f1.close()
								
		vtkFile.close()
		for proc in range(0,size):
		        solname="sol.p"+str(proc)+"."+str(step)+".h5"
			os.remove(solname)
				

if __name__ == '__main__':
    from optparse import OptionParser
    
    usage = ""
    parser = OptionParser(usage=usage)

    parser.add_option("-n","--size",
                      help="number of processors for run",
                      action="store",
                      type="int",
                      dest="size",
                      default=1)
		      
    parser.add_option("-s","--start",
                      help="start time step",
                      action="store",
                      type="int",
                      dest="start",
                      default=0)

    parser.add_option("-e","--end",
                      help="end time step",
                      action="store",
                      type="int",
                      dest="end",
                      default=10000)

    parser.add_option("-i","--increment",
                      help="increment",
                      action="store",
                      type="int",
                      dest="increment",
                      default=0)
	    	      
    parser.add_option("-f","--filebase_flow",
                      help="base filename for the flow/ns data",
                      action="store",
                      type="string",
                      dest="ns_filebase",
                      default="twp_navier_stokes_p")

    (opts,args) = parser.parse_args()
    
    if opts.ns_filebase == "twp_navier_stokes_p" :
      rd_default="redist_p"
      corr_default="ls_consrv_p"
      sol_default="solution"
    else:
      rd_default=opts.ns_filebase
      corr_default=opts.ns_filebase
      sol_default=opts.ns_filebase.rstrip("_p")

          		      
    parser.add_option("-p","--filebase_phi",
                      help="base filename for the phi data",
                      action="store",
                      type="string",
                      dest="phi_filebase",
                      default=rd_default)

    parser.add_option("-c","--filebase_corr",
                      help="base filename for the phi-correction data",
                      action="store",
                      type="string",
                      dest="corr_filebase",
                      default=corr_default)

    parser.add_option("-o","--filebase_out",
                      help="base filename for solution files",
                      action="store",
                      type="string",
                      dest="sol_filebase",
                      default=sol_default)
		      
#    parser.add_option("-x","--xmf",
#                      help="flag for indicating xmf output",
#                      action="store",
#                      type="boolean",
#                      dest="usexmf",
#                      default=False)		      

#    parser.add_option("-k","--vtk",
#                      help="flag for inicating vtk output",
#                      action="store",
#                      type="boolean",
#                      dest="usevtk",
#                      default=True)		      
    
    
    (opts,args) = parser.parse_args()
    
    if opts.increment == 0:
       opts.increment = 1
       opts.end = opts.start
       
    print "=================="
    print "   Extracting" 
    print "=================="     
    splitH5(opts.ns_filebase,opts.phi_filebase,opts.corr_filebase,opts.size,opts.start,opts.end,opts.increment)

    print "=================="
    print "   Composing" 
    print "=================="       
    H5toXMF(opts.sol_filebase,opts.size,opts.start,opts.end,opts.increment)
    H5toVTK(opts.sol_filebase,opts.size,opts.start,opts.finaltime,opts.stride)
