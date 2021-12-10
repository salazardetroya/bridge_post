size = 0.03;
//+
Point(1) = {0, 0, 0, size};
//+
Point(2) = {0.5, 0, 0, size};
//+
Point(3) = {0.5, 1.0, 0, size};
//+
Point(4) = {0, 1.0, 0, size};
//+
Point(5) = {0, 0, 0.9, size};
//+
Point(6) = {0.5, 0, 0.9, size};
//+
Point(7) = {0.5, 1, 0.9, size};
//+
Point(8) = {0, 1, 0.9, size};
//+
Point(9) = {0, 0, 1, size};
//+
Point(10) = {0.5, 0, 1, size};
//+
Point(11) = {0.5, 1, 1, size};
//+
Point(12) = {0, 1, 1, size};
//+
//+
Point(13) = {0.2, 0.3, 0, size};
//+
Point(14) = {0, 0.3, 0, size};
//+
Point(15) = {0.2, 0, 0, size};
//+
Line(1) = {1, 15};
//+
Line(2) = {15, 2};
//+
Line(3) = {2, 3};
//+
Line(4) = {3, 4};
//+
Line(5) = {4, 14};
//+
Line(6) = {14, 1};
//+
Line(7) = {15, 13};
//+
Line(8) = {13, 14};
//+
Line(9) = {1, 5};
//+
Line(10) = {5, 9};
//+
Line(11) = {9, 12};
//+
Line(12) = {12, 11};
//+
Line(13) = {11, 10};
//+
Line(14) = {10, 6};
//+
Line(15) = {6, 5};
//+
Line(16) = {9, 10};
//+
Line(17) = {7, 6};
//+
Line(18) = {5, 8};
//+
Line(19) = {8, 7};
//+
Line(20) = {12, 8};
//+
Line(21) = {11, 7};
//+
Line(22) = {7, 3};
//+
Line(23) = {4, 8};
//+
Line(24) = {6, 6};
//+
Line(25) = {2, 6};
//+
Line Loop(1) = {7, 8, 6, 1};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {2, 25, 15, -9, 1};
//+
Plane Surface(2) = {2};
//+
Line Loop(3) = {3, 4, 5, -8, -7, 2};
//+
Plane Surface(3) = {3};
//+
Line Loop(4) = {23, 19, 22, 4};
//+
Plane Surface(4) = {4};
//+
Line Loop(5) = {23, -18, -9, -6, -5};
//+
Plane Surface(5) = {5};
//+
Line Loop(6) = {18, -20, -11, -10};
//+
Plane Surface(6) = {6};
//+
Line Loop(7) = {19, -21, -12, 20};
//+
Plane Surface(7) = {7};
//+
Line Loop(8) = {13, 14, -17, -21};
//+
Plane Surface(8) = {8};
//+
Line Loop(9) = {15, 10, 16, 14};
//+
Plane Surface(9) = {9};
//+
Line Loop(10) = {16, -13, -12, -11};
//+
Plane Surface(10) = {10};
//+
Line Loop(11) = {17, 15, 18, 19};
//+
Plane Surface(11) = {11};
//+
Surface Loop(1) = {10, 9, 6, 7, 8, 11};
//+
Volume(1) = {1};
//+
Line Loop(12) = {22, -3, 25, -17};
//+
Plane Surface(12) = {12};
//+
Surface Loop(2) = {12, 4, 5, 2, 3, 1, 11};
//+
Volume(2) = {2};
//+
Physical Surface(2) = {1};
//+
Physical Surface(3) = {10};
//+
Physical Volume(0) = {1};
//+
Physical Volume(1) = {2};
//+
Physical Surface(4) = {5};
//+
Physical Surface(5) = {4};
