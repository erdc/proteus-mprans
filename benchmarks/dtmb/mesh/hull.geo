/*******************************************

Specifies the hull

********************************************/

// Bounding box hull  
xhull[0] = -0.5;
xhull[1] =  5.8;

yhull[0] = -0.5;
yhull[1] =  0.5;

zhull[0] = -0.10;
zhull[1] =  0.70;

hsl = -1;

/*******************************************
Add hull
********************************************/
Function addHull
  Geometry.Tolerance = 1e-2;
  Geometry.OCCSewFaces = 1;

   Merge 'hull.igs';

  // Surfaces & Volume
  hsl = newsl; Surface Loop(hsl) = {1,2,3,4,5,6,7,8,9,10,11,12,13}; 
  Physical Surface(7) = {1,2,3,4,5,6,7,8,9,10,11,12,13}; 

  // Set Characteristic lengths on hull
  Characteristic Length {3}  = hcl;    // Sonar Dome
  Characteristic Length {1}  = hcl;      // Keel
  Characteristic Length {13} = hcl;      // Keel

  Characteristic Length {4}    = hcl;    // Bow Deck
  Characteristic Length {5,7}  = hcl;    // Bow Deck
  Characteristic Length {2,6}  = hcl;    // Bow Deck
  Characteristic Length {8,12} = hcl;    // Bow Deck

  Characteristic Length {10}    = hcl;   // Stern upper   keel
  Characteristic Length {16}    = hcl;   // Stern medium  keel
  Characteristic Length {14}    = hcl;   // Stern lower   keel
  Characteristic Length {9,11}  = hcl;   // Stern upper   sides
  Characteristic Length {15,17} = hcl;   // Stern upper   sides
   
Return

/*******************************************
Specifies the hull
********************************************/
Function addRef

  // Dome boxes
  Field[1] = Box;

  Field[1].XMin =  0.55;
  Field[1].XMax =  0.70;

  Field[1].YMin = -0.1;
  Field[1].YMax =  0.1;

  Field[1].ZMin = -0.1;
  Field[1].ZMax =  0.05;

  Field[1].VIn  =  0.66*hcl;
  Field[1].VOut = 99.9;

  Field[2] = Box;

  Field[2].XMin =  0.40;
  Field[2].XMax =  0.90;

  Field[2].YMin = -0.2;
  Field[2].YMax =  0.2;

  Field[2].ZMin = -0.3;
  Field[2].ZMax =  0.15;

  Field[2].VIn  =  hcl;
  Field[2].VOut =  99.9;


  // Stern box 1
  Field[11] = Box;

  Field[11].XMin =  5.00;
  Field[11].XMax =  10.0;

  Field[11].YMin = -1.0;
  Field[11].YMax =  1.0;

  Field[11].ZMin =  0.00;
  Field[11].ZMax =  0.60;

  Field[11].VIn  =  1.5*hcl;
  Field[11].VOut = 99.9;

  // Stern box 2
  Field[12] = Box;

  Field[12].XMin =  5.25;
  Field[12].XMax =  7.50;

  Field[12].YMin = -0.85;
  Field[12].YMax =  0.85;

  Field[12].ZMin =  0.10;
  Field[12].ZMax =  0.40;

  Field[12].VIn  =  hcl;
  Field[12].VOut = 99.9;

  // Stern box 3
  Field[13] = Box;

  Field[13].XMin =  5.60;
  Field[13].XMax =  6.50;

  Field[13].YMin = -0.75;
  Field[13].YMax =  0.75;

  Field[13].ZMin =  0.20;
  Field[13].ZMax =  0.30;

  Field[13].VIn  = 0.5*hcl;
  Field[13].VOut = 99.9;
  
  // Stern box 2
  Field[31] = Box;

  Field[31].XMin =  -0.25;
  Field[31].XMax =   0.50;

  Field[31].YMin = -0.45;
  Field[31].YMax =  0.45;

  Field[31].ZMin =  0.20;
  Field[31].ZMax =  0.40;

  Field[31].VIn  =  0.7*hcl;
  Field[31].VOut = 99.9;

  Field[32] = Box;

  Field[32].XMin =  -0.25;
  Field[32].XMax =   0.25;

  Field[32].YMin = -0.25;
  Field[32].YMax =  0.25;

  Field[32].ZMin =  -0.20;
  Field[32].ZMax =   0.30;

  Field[32].VIn  =  0.8*hcl;
  Field[32].VOut = 99.9;
    
  // Combine
  Field[100] = Min;
  Field[100].FieldsList = {1,2,11,12,13,31,32};
  Background Field = 100;

Return
