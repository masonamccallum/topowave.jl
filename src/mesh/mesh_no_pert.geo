ls = 5;

Xi = 120; 
R = 50.0;
A = R * 0.05;
 
Point(1) = {0, 0, 0, ls};
Point(2) = {Xi, 0, 0, ls};
Point(3) = {Xi, R, 0, ls};
Point(4) = {0, R, 0, ls};
Point(9) = {0,2*R,0,ls};
Point(10) = {Xi,2*R,0,ls};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(10) = {3,10};
Line(11) = {10,9};
Line(12) = {9,4};

Line Loop(11) = {4, 1, 2, 3};
Plane Surface(12) = {11};
Line Loop(12) = {-3, 10,11,12};
Plane Surface(13) = {12};

Physical Surface("top", 16) = {13};
Physical Surface("bottom", 17) = {12};

