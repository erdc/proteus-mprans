Include "addbox.geo";
Include "hull.geo";

Mesh.RandomFactor = 1e-4;
Mesh.SaveElementTagType= 2;
Mesh.QualityType = 1;
Mesh.Optimize = 1;

//Mesh.Algorithm=6;
//Mesh.Algorithm3D=1;

xdim = 6;
x[0] =  -5.0;      hx[0] = 1.00;
x[1] =  -1.25;     hx[1] = 0.50;
x[2] =  -0.75;     hx[2] = 0.40;
x[3] =   6.25;     hx[3] = 0.30;
x[4] =   8.0;      hx[4] = 0.40;
x[5] =  15.0;      hx[5] = 1.00;

ydim = 6;
y[0+ydim*0] =  -7.5; hy[0+ydim*0] = 1.50;
y[1+ydim*0] =  -1.2; hy[1+ydim*0] = 0.40;
y[2+ydim*0] =  -0.8; hy[2+ydim*0] = 0.20;
y[3+ydim*0] =   0.8; hy[3+ydim*0] = 0.20;
y[4+ydim*0] =   1.2; hy[4+ydim*0] = 0.40;
y[5+ydim*0] =   7.5; hy[5+ydim*0] = 1.50;

y[0+ydim*1] =  -7.5; hy[0+ydim*1] = 1.50;
y[1+ydim*1] =  -1.0; hy[1+ydim*1] = 0.40;
y[2+ydim*1] =  -0.6; hy[2+ydim*1] = 0.20;
y[3+ydim*1] =   0.6; hy[3+ydim*1] = 0.20;
y[4+ydim*1] =   1.0; hy[4+ydim*1] = 0.40;
y[5+ydim*1] =   7.5; hy[5+ydim*1] = 1.50;

y[0+ydim*2] =  -7.5; hy[0+ydim*2] = 1.50;
y[1+ydim*2] =  -1.0; hy[1+ydim*2] = 0.30;
y[2+ydim*2] =  -0.6; hy[2+ydim*2] = 0.20;
y[3+ydim*2] =   0.6; hy[3+ydim*2] = 0.20;
y[4+ydim*2] =   1.0; hy[4+ydim*2] = 0.30;
y[5+ydim*2] =   7.5; hy[5+ydim*2] = 1.50;

y[0+ydim*3] =  -7.5; hy[0+ydim*3] = 1.50;
y[1+ydim*3] =  -2.0; hy[1+ydim*3] = 0.25;
y[2+ydim*3] =  -0.8; hy[2+ydim*3] = 0.20;
y[3+ydim*3] =   0.8; hy[3+ydim*3] = 0.20;
y[4+ydim*3] =   2.0; hy[4+ydim*3] = 0.25;
y[5+ydim*3] =   7.5; hy[5+ydim*3] = 1.50;

y[0+ydim*4] =  -7.5; hy[0+ydim*4] = 1.50;
y[1+ydim*4] =  -2.2; hy[1+ydim*4] = 0.20;
y[2+ydim*4] =  -0.8; hy[2+ydim*4] = 0.20;
y[3+ydim*4] =   0.8; hy[3+ydim*4] = 0.20;
y[4+ydim*4] =   2.2; hy[4+ydim*4] = 0.20;
y[5+ydim*4] =   7.5; hy[5+ydim*4] = 1.50;

y[0+ydim*5] =  -7.5; hy[0+ydim*5] = 1.50;
y[1+ydim*5] =  -2.8; hy[1+ydim*5] = 0.20;
y[2+ydim*5] =  -0.8; hy[2+ydim*5] = 0.20;
y[3+ydim*5] =   0.8; hy[3+ydim*5] = 0.20;
y[4+ydim*5] =   2.8; hy[4+ydim*5] = 0.20;
y[5+ydim*5] =   7.5; hy[5+ydim*5] = 1.50;


zdim = 5;
z[0] = -3.25;      hz[0] = 1.00;
z[1] = -0.25;      hz[1] = 0.30;
z[2] =  0.241984;  hz[2] = 0.25;
z[3] =  0.80;      hz[3] = 0.40;
z[4] =  2.50;      hz[4] = 1.00;


hcl = 0.75*(hz[2]*(hx[2]+hx[3])*(hy[2]+hy[3])/4.0);

Mesh.CharacteristicLengthMax=0.25; 
 
Call addHull;
Call addOneBox;
Call addRef;
