//+
SetFactory("OpenCASCADE");
Rectangle(1) = {-4, 0, 3, 17, 5, 0};
//+
Circle(9) = {2, 2.5, 3, 1.25, 0, 2*Pi};
//+
Rectangle(3) = {3.75, 1, 3, 1.5, 3, 0};
//+
Point(15) = {6, 2.5, 3, 1.0};
//+
Point(16) = {8.5, 3.5, 3, 1.0};
//+
Point(17) = {6.5, 4.5, 3, 1.0};
//+
Line(15) = {15, 16};
//+
Line(16) = {16, 17};
//+
Line(17) = {17, 15};
//+
Ellipse(18) = {7.5, 1.6, 3, 1, 0.5, 0, 2*Pi};
//+
Ellipse(19) = {7.5, 1.5, 3, 1.6, 0.9, 0, 2*Pi};
//+
Curve Loop(4) = {9};
//+
Surface(4) = {4};
//+
Curve Loop(8) = {17, 15, 16};
//+
Surface(6) = {8};
//+
Curve Loop(12) = {18};
//+
Surface(7) = {12};
//+
Curve Loop(14) = {19};
//+
Surface(8) = {14};
//+
BooleanDifference{ Surface{8}; Delete; }{ Surface{7}; }
//+
BooleanDifference{ Surface{1}; Delete; }{ Surface{4}; Surface{3}; Surface{8}; Surface{7}; Surface{6}; }
//+
Point(24) = {-4, 0, 0, 1.0};
//+
Point(25) = {13, 0, 0, 1.0};
//+
Point(26) = {13, 5, 0, 1.0};
//+
Point(27) = {-4, 5, 0, 1.0};
//+
Line(24) = {24, 25};
//+
Line(25) = {25, 26};
//+
Line(26) = {26, 27};
//+
Line(27) = {27, 24};
//+
Line(28) = {24, 20};
//+
Line(29) = {27, 22};
//+
Line(30) = {25, 21};
//+
Line(31) = {26, 23};
//+
Transfinite Curve {15, 17, 16} = 20 Using Progression 1;
//+
Transfinite Curve {10, 12} = 12 Using Progression 1;
//+
Transfinite Curve {9, 19} = 60 Using Progression 1;
//+
Transfinite Curve {18} = 40 Using Progression 1;
//+
Curve Loop(15) = {24, 25, 26, 27};
//+
Surface(9) = {15};
//+
Curve Loop(17) = {28, 20, -30, -24};
//+
Surface(10) = {17};
//+
Curve Loop(19) = {30, 22, -31, -25};
//+
Surface(11) = {19};
//+
Curve Loop(21) = {31, 23, -29, -26};
//+
Surface(12) = {21};
//+
Curve Loop(23) = {21, -28, -27, 29};
//+
Surface(13) = {23};
//+
Surface Loop(1) = {10, 9, 11, 12, 1, 13, 3, 8, 7, 4, 6};
//+
Volume(1) = {1};
//+
Transfinite Curve {30, 31, 29, 28} = 8 Using Progression 0.7;
//+
Transfinite Curve {21, 22} = 40 Using Progression 1;
//+
Transfinite Curve {27, 25} = 8 Using Progression 1;
//+
Transfinite Curve {20, 23} = 80 Using Progression 1;
//+
Transfinite Curve {24, 26} = 20 Using Progression 1;
//+
Transfinite Curve {13, 11} = 25 Using Progression 1;
//+
Physical Surface("plate") = {4, 3, 6, 8};
//+
Physical Surface("freeSurface") = {1, 7};
//+
Physical Surface("inlet") = {13};
//+
Physical Surface("outlet") = {11};
//+
Physical Surface("sides") = {10, 12};
//+
Physical Surface("bottom") = {9};
//+
Physical Volume("water") = {1};
