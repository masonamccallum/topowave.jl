ls = 5;

Xi = 100; // um
Xo = 100; // um
L = 100.0; // um

x0 = Xi + L/2.0;
R = 50.0;   // um
f0 = 0.5;   // 0--1

Z = 5;

Point(1) = {0, 0, 0, ls};
#Point(2) = {Xi, 0, 0, ls};
Point(3) = {Xi, R, 0, ls};
Point(4) = {0, R, 0, ls};

Point(5) = {Xi + L, 0, 0, ls};
Point(8) = {Xi + L, R, 0, ls};

Line(1) = {1, 5};
Line(2) = {5, 8};
Line(3) = {3, 4};
Line(4) = {4, 1};

pList[0] = 3; // First point label
nPoints = 21; // Number of discretization points (top-right point of the inlet region)
r = Xi+L; 
l = Xi; 
n = 3;
A = R * 0.05;
For i In {1 : nPoints}
  x = Xi + L*i/(nPoints + 1);
  pList[i] = newp;
  Point(pList[i]) = {x,
                ( R + A*Sin((4*Pi/r)*n*(x-l))),
                0,
                ls};
EndFor
pList[nPoints+1] = 8; // Last point label (top-left point of the outlet region)

Spline(newl) = pList[];


Transfinite Line {9, 10} = Ceil(L/ls) Using Progression 1;
Transfinite Line {4, -2, 8, -6} = Ceil(R/ls) Using Progression 1.1;
Transfinite Line {1, 3} = Ceil(Xi/ls) Using Progression 1;
Transfinite Line {5, 7} = Ceil(Xo/ls) Using Progression 1;

Line Loop(11) = {4, 1, 8, 10, 3};
Plane Surface(12) = {11};
