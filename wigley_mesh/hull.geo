/*******************************************

Specifies the hull

********************************************/

// Bounding box hull  
xhull[0] = -0.5;
xhull[1] =  0.5;

yhull[0] = -0.09;
yhull[1] =  0.09;

zhull[0] =  0.4375;
zhull[1] =  0.6550;

hsl = -1;

/*******************************************
Add hull
********************************************/
Function addHull
  Geometry.Tolerance = 1e-5;
  Geometry.OCCSewFaces = 1;

  Merge 'wigley.igs';

  // Surfaces & Volume
  hsl = newsl; Surface Loop(hsl) = {1,2,3,4,5,6}; 
  Physical Surface(11) = {1,2,3,4,5,6}; 

  // Set Characteristic lengths on hull
  Characteristic Length {5,6} = hcl;    // Keel
  Characteristic Length {4} = hcl;      // Waterline stern
  Characteristic Length {2} = hcl;      // Waterline bow
  Characteristic Length {1,3} = hcl;    // Deck lower
  Characteristic Length {7,8} = hcl;    // Deck upper
   
Return

/*******************************************
Specifies the hull
********************************************/
Function addRef

  // Bow boxes
  Field[11] = Box;

  Field[11].XMin =  -0.55;
  Field[11].XMax =  -0.45;

  Field[11].YMin = -0.05;
  Field[11].YMax =  0.05;

  Field[11].ZMin =  0.40;
  Field[11].ZMax =  0.60;

  Field[11].VIn  =  0.7*hcl;
  Field[11].VOut = 99.9;

  Field[12] = Box;

  Field[12].XMin = -0.70;
  Field[12].XMax = -0.30;

  Field[12].YMin = -0.10;
  Field[12].YMax =  0.10;

  Field[12].ZMin =  0.30;
  Field[12].ZMax =  0.70;

  Field[12].VIn  =  0.9*hcl;
  Field[12].VOut = 99.9;

  
  // Stern boxes 
  Field[31] = Box;

  Field[31].XMin =   0.40;
  Field[31].XMax =   0.60;

  Field[31].YMin = -0.15;
  Field[31].YMax =  0.15;

  Field[31].ZMin =  0.40;
  Field[31].ZMax =  0.60;

  Field[31].VIn  =  0.6*hcl;
  Field[31].VOut = 99.9;

  Field[32] = Box;

  Field[32].XMin =   0.30;
  Field[32].XMax =   0.70;

  Field[32].YMin = -0.25;
  Field[32].YMax =  0.25;

  Field[32].ZMin =   0.30;
  Field[32].ZMax =   0.70;

  Field[32].VIn  =  0.7*hcl;
  Field[32].VOut = 99.9;
    
  // Combine
  Field[100] = Min;
  Field[100].FieldsList = {11,12,31,32};
  Background Field = 100;

Return
