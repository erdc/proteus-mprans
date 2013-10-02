
unitsize(4.0 inches / 270.000000);
size(11 inches);
real Lx=270.000000;
real Ly=18.190403;
real offset=0.0125Lx;
real x=0.000000;
real y=-11.938669;
string strx="$270.00\mbox{m}$";
string stry="$18.19\mbox{m}$";
draw(strx,(x,y-offset)--(x+Lx,y-offset),S,black,Bars,Arrows,PenMargins);
draw(stry,(x-offset,y)--(x-offset,y+Ly),W,black,Bars,Arrows,PenMargins);
import graph;
import palette;
pen[] allPens=Wheel();
pen[] myPens = new pen[3+1];
for(int i=0;i<3+1;++i)
  {
   int iPen = round(i*allPens.length/(3+1));
   myPens[i] = allPens[iPen];
  }
draw((270.000000,6.251734)--(0.000000,6.251734),myPens[3]+linewidth(0.01));
draw((0.000000,6.251734)--(0.000000,-11.938669),myPens[0]+linewidth(0.01));
draw((0.000000,-11.938669)--(101.710000,-8.530000),myPens[2]+linewidth(0.01));
draw((101.710000,-8.530000)--(130.000000,-7.581900),myPens[2]+linewidth(0.01));
draw((130.000000,-7.581900)--(135.000000,-7.428100),myPens[2]+linewidth(0.01));
draw((135.000000,-7.428100)--(140.000000,-7.273100),myPens[2]+linewidth(0.01));
draw((140.000000,-7.273100)--(145.000000,-7.122700),myPens[2]+linewidth(0.01));
draw((145.000000,-7.122700)--(150.000000,-6.964500),myPens[2]+linewidth(0.01));
draw((150.000000,-6.964500)--(155.000000,-6.841200),myPens[2]+linewidth(0.01));
draw((155.000000,-6.841200)--(160.000000,-6.677300),myPens[2]+linewidth(0.01));
draw((160.000000,-6.677300)--(165.000000,-6.502900),myPens[2]+linewidth(0.01));
draw((165.000000,-6.502900)--(170.000000,-6.399100),myPens[2]+linewidth(0.01));
draw((170.000000,-6.399100)--(175.000000,-6.252100),myPens[2]+linewidth(0.01));
draw((175.000000,-6.252100)--(180.000000,-6.137300),myPens[2]+linewidth(0.01));
draw((180.000000,-6.137300)--(185.000000,-6.068900),myPens[2]+linewidth(0.01));
draw((185.000000,-6.068900)--(190.000000,-2.633900),myPens[2]+linewidth(0.01));
draw((190.000000,-2.633900)--(195.000000,-1.360000),myPens[2]+linewidth(0.01));
draw((195.000000,-1.360000)--(215.000000,-1.360000),myPens[2]+linewidth(0.01));
draw((215.000000,-1.360000)--(215.000000,3.220000),myPens[2]+linewidth(0.01));
draw((215.000000,3.220000)--(215.300000,3.220000),myPens[2]+linewidth(0.01));
draw((215.300000,3.220000)--(215.300000,-1.360000),myPens[2]+linewidth(0.01));
draw((215.300000,-1.360000)--(235.000000,-1.360000),myPens[2]+linewidth(0.01));
draw((235.000000,-1.360000)--(235.000000,-7.920000),myPens[2]+linewidth(0.01));
draw((235.000000,-7.920000)--(270.000000,-7.920000),myPens[2]+linewidth(0.01));
draw((270.000000,-7.920000)--(270.000000,6.251734),myPens[2]+linewidth(0.01));
