Include "addbox.geo";
Include "hull.geo";

//Mesh.RandomFactor = 1e-10;
Mesh.SaveElementTagType= 2;
Mesh.QualityType = 1;
Mesh.Optimize = 1;

//Mesh.Algorithm=6;
//Mesh.Algorithm3D=1;

xdim = 6;
x[0] =  -1.25;      hx[0] = 1.00;
x[1] =  -0.90;      hx[1] = 0.60;
x[2] =  -0.65;      hx[2] = 0.50;
x[3] =   0.75;      hx[3] = 0.40;
x[4] =   1.25;      hx[4] = 0.50;
x[5] =   1.50;      hx[5] = 1.00;

ydim = 6;
y[0+ydim*0] =  -1.0;  hy[0+ydim*0] = 0.80;
y[1+ydim*0] =  -0.25; hy[1+ydim*0] = 0.50;
y[2+ydim*0] =  -0.15; hy[2+ydim*0] = 0.30;
y[3+ydim*0] =   0.15; hy[3+ydim*0] = 0.30;
y[4+ydim*0] =   0.25; hy[4+ydim*0] = 0.50;
y[5+ydim*0] =   1.0;  hy[5+ydim*0] = 0.80;

y[0+ydim*1] =  -1.0;  hy[0+ydim*1] = 0.80;
y[1+ydim*1] =  -0.2;  hy[1+ydim*1] = 0.60;
y[2+ydim*1] =  -0.1;  hy[2+ydim*1] = 0.30;
y[3+ydim*1] =   0.1;  hy[3+ydim*1] = 0.30;
y[4+ydim*1] =   0.2;  hy[4+ydim*1] = 0.60;
y[5+ydim*1] =   1.0;  hy[5+ydim*1] = 0.80;

y[0+ydim*2] =  -1.0;  hy[0+ydim*2] = 0.80;
y[1+ydim*2] =  -0.2;  hy[1+ydim*2] = 0.40;
y[2+ydim*2] =  -0.1;  hy[2+ydim*2] = 0.30;
y[3+ydim*2] =   0.1;  hy[3+ydim*2] = 0.30;
y[4+ydim*2] =   0.2;  hy[4+ydim*2] = 0.40;
y[5+ydim*2] =   1.0;  hy[5+ydim*2] = 0.80;

y[0+ydim*3] =  -1.0;  hy[0+ydim*3] = 0.80;
y[1+ydim*3] =  -0.5;  hy[1+ydim*3] = 0.35;
y[2+ydim*3] =  -0.2;  hy[2+ydim*3] = 0.30;
y[3+ydim*3] =   0.2;  hy[3+ydim*3] = 0.30;
y[4+ydim*3] =   0.5;  hy[4+ydim*3] = 0.35;
y[5+ydim*3] =   1.0;  hy[5+ydim*3] = 0.80;

y[0+ydim*4] =  -1.0;  hy[0+ydim*4] = 0.80;
y[1+ydim*4] =  -0.6;  hy[1+ydim*4] = 0.30;
y[2+ydim*4] =  -0.15; hy[2+ydim*4] = 0.30;
y[3+ydim*4] =   0.15; hy[3+ydim*4] = 0.30;
y[4+ydim*4] =   0.6;  hy[4+ydim*4] = 0.30;
y[5+ydim*4] =   1.0;  hy[5+ydim*4] = 0.80;

y[0+ydim*5] =  -1.0;  hy[0+ydim*5] = 0.80;
y[1+ydim*5] =  -0.75; hy[1+ydim*5] = 0.30;
y[2+ydim*5] =  -0.15; hy[2+ydim*5] = 0.30;
y[3+ydim*5] =   0.15; hy[3+ydim*5] = 0.30;
y[4+ydim*5] =   0.75; hy[4+ydim*5] = 0.30;
y[5+ydim*5] =   1.0;  hy[5+ydim*5] = 0.80;


zdim = 5;
z[0] = 0.00;      hz[0] = 1.00;
z[1] = 0.35;      hz[1] = 0.40;
z[2] = 0.50;      hz[2] = 0.25;
z[3] = 0.70;      hz[3] = 0.60;
z[4] = 0.80;      hz[4] = 1.00;


hcl = 0.8*(hz[2]*(hx[2]+hx[3])*(hy[2]+hy[3])/4.0);

Mesh.CharacteristicLengthMax=0.2; 
 
Call addHull;
Call addBoxes;
Call addRef;
