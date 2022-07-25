ls = 5;

Xi = 100; 
L = Xi*0.2;
R = 50.0;
A = R * 0.05;
 
Point(1) = {0, 0, 0, ls};
Point(2) = {Xi, 0, 0, ls};
Point(3) = {Xi, R, 0, ls};
Point(4) = {0, R, 0, ls};
Point(5) = {Xi + L, 0, 0, ls};
Point(8) = {Xi + L, R, 0, ls};
Point(9) = {0,2*R,0,ls};
Point(10) = {Xi,2*R,0,ls};
Point(11) = {Xi+L,2*R,0,ls};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(8) = {8, 5};
Line(9) = {2, 5};
Line(10) = {3,10};
Line(11) = {10,9};
Line(12) = {9,4};
Line(13) = {8,11};
Line(14) = {11,10};


pList[0] = 3; // First point label
nPoints = 21; // Number of discretization points (top-right point of the inlet region)
r = Xi+L; 
l = Xi; 
n = 3;
For i In {1 : nPoints}
  x = Xi + L*i/(nPoints + 1);
  pList[i] = newp;
  Point(pList[i]) = {x,
                ( R + A*Sin((2*Pi/(r-l))*n*(x-l))),
                0,
                ls};
EndFor

pList[nPoints+1] = 8; // Last point label (top-left point of the outlet region)
Spline(newl) = pList[];

Line Loop(11) = {4, 1, 2, 3};
Plane Surface(12) = {11};
Line Loop(12) = {-3, 10,11,12};
Plane Surface(13) = {12};

Line Loop(13) = {-9,2,15,8};
Plane Surface(14) = {13};
Line Loop(14) = {-15,10,-14,-13};
Plane Surface(15) = {14};

Physical Surface("top", 16) = {13, 15};
Physical Surface("bottom", 17) = {12, 14};

Transfinite Curve {15} = 50 Using Bump 1;

Physical Surface("top", 16) += {13, 15};
Physical Surface("bottom", 17) += {14, 12};
