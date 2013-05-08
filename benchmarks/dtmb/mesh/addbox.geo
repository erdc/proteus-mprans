


Function addBox
  // Points 
  p1  = newp+99; Point(p1) = {x1, y1, z1, hx1*hy1*hz1};
  p2  = newp+99; Point(p2) = {x2, y3, z1, hx2*hy3*hz1};
  p3  = newp+99; Point(p3) = {x2, y4, z1, hx2*hy4*hz1};
  p4  = newp+99; Point(p4) = {x1, y2, z1, hx1*hy2*hz1};

  p11 = newp+99; Point(p11) = {x1, y1, z2, hx1*hy1*hz2};
  p12 = newp+99; Point(p12) = {x2, y3, z2, hx2*hy3*hz2};
  p13 = newp+99; Point(p13) = {x2, y4, z2, hx2*hy4*hz2};
  p14 = newp+99; Point(p14) = {x1, y2, z2, hx1*hy2*hz2};

  // Lines 
  l1  = newl+99; Line(l1 ) = {p1,p2};
  l2  = newl+99; Line(l2 ) = {p2,p3};
  l3  = newl+99; Line(l3 ) = {p3,p4};
  l4  = newl+99; Line(l4 ) = {p4,p1};

  l11 = newl+99; Line(l11) = {p11,p12};
  l12 = newl+99; Line(l12) = {p12,p13};
  l13 = newl+99; Line(l13) = {p13,p14};
  l14 = newl+99; Line(l14) = {p14,p11};

  l21 = newl+99; Line(l21) = {p1,p11};
  l22 = newl+99; Line(l22) = {p2,p12};
  l23 = newl+99; Line(l23) = {p3,p13};
  l24 = newl+99; Line(l24) = {p4,p14};

  // Surfaces
  ll1 = newll+99; Line Loop(ll1) = {l1,l2,l3,l4};     

  ll2 = newll+99; Line Loop(ll2) = {l1,l22,-l11,-l21};
  ll3 = newll+99; Line Loop(ll3) = {l2,l23,-l12,-l22};
  ll4 = newll+99; Line Loop(ll4) = {l3,l24,-l13,-l23};
  ll5 = newll+99; Line Loop(ll5) = {l4,l21,-l14,-l24};

  ll6 = newll+99; Line Loop(ll6) = {l11,l12,l13,l14}; 

  s1 = news+99; Plane Surface(s1) = {ll1};

  s2 = news+99; Plane Surface(s2) = {ll2};
  s3 = news+99; Plane Surface(s3) = {ll3};
  s4 = news+99; Plane Surface(s4) = {ll4};
  s5 = news+99; Plane Surface(s5) = {ll5};

  s6 = news+99; Plane Surface(s6) = {ll6};

  // Volume
  sl = newsl+99; Surface Loop (sl) = {s1,s2,s3,s4,s5,s6}; 
  v  = newv+99; Volume(v) = {sl};
    
Return

Function addPoints
  For i In {0:xdim-1}
    For j In {0:ydim-1}
      For k In {0:zdim-1}
        p  = newp+99; Point(p) = {x[i], y[j+ydim*i], z[k], hx[i]*hy[j]*hz[k]};
      EndFor
    EndFor
  EndFor
Return


Function addOneBox 
  count = 0; 


  //Call addPoints;

  boc = 0;
  frc = 0;
  ric = 0;
  bac = 0;
  lec = 0;
  toc = 0;

        x1=x[0];       hx1=hx[0];
        x2=x[5];       hx2=hx[5];

        y1=y[0];    hy1=hy[0];
        y2=y[5];    hy2=hy[5];
        
        y3=y[0];  hy3=hy[0];
        y4=y[5];  hy4=hy[5];
              
        z1=z[0];       hz1=hz[0];
        z2=z[4];       hz2=hz[4];
          
        Call addBox;

        // Add volume and surfaces to global lists       
        If (v > 0)
	  vol[count] = v;
          count = count+1;
	EndIf
		

  // Build volume around hull
  bsl = newsl; Surface Loop(bsl) = {lowerSurf[],upperSurf[]}; 
  If (hsl > 0)
    v   = newv;  Volume(v) = {bsl,hsl};
  EndIf
  If (hsl == -1) // No hull !!!
    v   = newv;  Volume(v) = {bsl};
  EndIf
  
  vol[count] = v;

	  left[lec] = s5;
          lec = lec  + 1;
       
	  right[ric] = s3;
          ric = ric + 1;

	  front[frc] = s2;
          frc = frc + 1;
       
	  back[bac] = s4;
          bac = bac  + 1;

	  bottom[boc] = s1;
          boc = boc  + 1;

	  top[toc] = s6;
          toc = toc + 1;
  //
  Printf("bottom = %g",boc);
  Printf("top    = %g",toc);
  Printf("right  = %g",ric);
  Printf("left   = %g",lec);
  Printf("front  = %g",frc);
  Printf("back   = %g",bac);
  
  // Define physical entities
  Physical Surface(1) = {bottom[]};

  Physical Surface(2) = {front[]};
  Physical Surface(3) = {right[]};
  Physical Surface(4) = {back[]};
  Physical Surface(5) = {left[]};

  Physical Surface(6) = {top[]};
 
  Physical Volume(2) = {vol[]};

  // Delete redundant entities
  Coherence; 
Return

Function addBoxes 
  count = 0; 


  //Call addPoints;

  boc = 0;
  frc = 0;
  ric = 0;
  bac = 0;
  lec = 0;
  toc = 0;

  For i In {0:xdim-2}
    For j In {0:ydim-2}
      For k In {0:zdim-2}

      
        x1=x[i];         hx1=hx[i];
        x2=x[i+1];       hx2=hx[i+1];

        y1=y[j  +ydim*(i)];    hy1=hy[j   + ydim*(i)];
        y2=y[j+1+ydim*(i)];    hy2=hy[j+1 + ydim*(i)];
        
        y3=y[j  +ydim*(i+1)];  hy3=hy[j   +ydim*(i+1)];
        y4=y[j+1+ydim*(i+1)];  hy4=hy[j+1 +ydim*(i+1)];
              
        z1=z[k];         hz1=hz[k];
        z2=z[k+1];       hz2=hz[k+1];
          
        Call addBox;

        // Delete volume and approproriate surface when hull is present	
	If ((x1 <= xhull[0]) && (x2 >= xhull[0]) || (x1 <= xhull[1]) && (x2 >= xhull[1]))
  	  If ((y1 <= yhull[0]) && (y2 >= yhull[0]) || (y1 <= yhull[1]) && (y2 >= yhull[1]))
  	    If ((z1 <= zhull[0]) && (z2 >= zhull[0]) || (z1 <= zhull[1]) && (z2 >= zhull[1]))
              Printf("Hull volume at [%g,%g,%g]",i,j,k);              
	      Printf("x = [%g,%g]",x1,x2);
	      Printf("y = [%g,%g]",y1,y2);
	      Printf("z = [%g,%g]",z1,z2);

              Printf("Deleting volume %g",v);
	      Delete {Volume{v};}

              If (z2 > zhull[1])      // Upper volume
                Printf("Deleting surface %g",s1);
                Delete {Surface{s1};}
		lowerSurf[] = {s2,s3,s4,s5,s6};
		
	      EndIf	
	      
              If (z1 < zhull[0])      // Lower volume      
	        Printf("Deleting surface %g",s6);
                Delete {Surface{s6};}
		upperSurf[] = {s1,s2,s3,s4,s5};
              EndIf
	
	      v = -1;
	    EndIf	
	  EndIf 	
	EndIf 

        // Add volume and surfaces to global lists       
        If (v > 0)
	  vol[count] = v;
          count = count+1;
	EndIf
		
	If (i == 0)    // Left
	  left[lec] = s5;
          lec = lec  + 1;
	EndIf
       
	If (i == xdim-2)   // Right
	  right[ric] = s3;
          ric = ric + 1;
	EndIf

	If (j == 0)  // Front
	  front[frc] = s2;
          frc = frc + 1;
	EndIf
       
	If (j == ydim-2)  // Back
	  back[bac] = s4;
          bac = bac  + 1;
	EndIf       

	If (k == 0)  // Bottom
	  bottom[boc] = s1;
          boc = boc  + 1;
	EndIf
       
	If (k == zdim-2)   // Top
	  top[toc] = s6;
          toc = toc + 1;
	EndIf      
	 
        
      EndFor
    EndFor
  EndFor

  // Build volume around hull
  bsl = newsl; Surface Loop(bsl) = {lowerSurf[],upperSurf[]}; 
  If (hsl > 0)
    v   = newv;  Volume(v) = {bsl,hsl};
  EndIf
  If (hsl == -1) // No hull !!!
    v   = newv;  Volume(v) = {bsl};
  EndIf
  
  vol[count] = v;

  //
  Printf("bottom = %g",boc);
  Printf("top    = %g",toc);
  Printf("right  = %g",ric);
  Printf("left   = %g",lec);
  Printf("front  = %g",frc);
  Printf("back   = %g",bac);
  
  // Define physical entities
  Physical Surface(1) = {bottom[]};

  Physical Surface(2) = {front[]};
  Physical Surface(3) = {right[]};
  Physical Surface(4) = {back[]};
  Physical Surface(5) = {left[]};

  Physical Surface(6) = {top[]};
 
  Physical Volume(2) = {vol[]};

  // Delete redundant entities
  Coherence; 
Return
