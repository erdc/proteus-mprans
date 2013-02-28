#!/usr/bin/env python
def tetgen23dm(filebase,append_default_mat=True):
    #won't work with comments
    felem = open(filebase+'.ele')
    fnode = open(filebase+'.node')

    line = felem.readline()
    while line and line.split()[0] == "#":
        line = felem.readline()
    nelem = int(line.split()[0])
    
    fout = open(filebase+'.3dm','w')
    fout.write("MESH3D\n")
    while line:
        line = felem.readline()
        if line and line.split()[0] != "#":
            if append_default_mat:
                line = line.replace('\n','  1 \n')
                fout.write("E4T "+line)
            else:
                fout.write("E4T "+line)
    felem.close()
    line = fnode.readline()
    while line and line.split()[0] == "#":
        line = fnode.readline()
    nnode = int(line.split()[0])
    while line:
        line = fnode.readline()
        if line and line.split()[0] != "#":
            fout.write("ND  "+line)
    fnode.close()
    fout.close()
if __name__ == "__main__":

    import sys

    name = "open-book"
    if len(sys.argv) > 1:
        name = sys.argv[1]

    
    tetgen23dm(name)

