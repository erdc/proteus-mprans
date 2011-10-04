default: marinHexMesh gmsh2tetgen gmsh2tetgen_withtags

marinHexMesh: marinHexMesh.f 
	${FC}  marinHexMesh.f -o marinHexMesh

gmsh2tetgen:gmsh2tetgen.f 	
	${FC}  gmsh2tetgen.f -o gmsh2tetgen

gmsh2tetgen_withtags:gmsh2tetgen_withtags.f 	
	${FC}  gmsh2tetgen_withtags.f -o gmsh2tetgen_withtags
			
clean:
	rm -f  marinHexMesh gmsh2tetgen gmsh2tetgen_withtags
