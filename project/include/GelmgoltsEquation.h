#pragma once

#include <iostream>

//#define ITERAT 11000
//#define ITERJAC ITERAT
//#define ITERZEID ITERAT
//#define STEP 0.01

using namespace std;

// Main
void Zeidel(double** _mesh, int _rows, int _cols, double _k, double _step, int ITERAT);
void Jacobi(double** _mesh, int _rows, int _cols, double _k, double _step, int ITERAT);
void ZeidelParal(double** _mesh, int _rows, int _cols, double _k, double _step, int ITERAT);
void JacobiParal(double** _mesh, int _rows, int _cols, double _k, double _step, int ITERAT);
double f(double _x, double _y, double _k);
double exactSolution(double _x, double _y);
double** rightPart(double _step, double _rows, double _cols, double _k);

// Additional
double** createMesh(double _xBorder, double _yBorder, double _step, double* _rows, double* _cols);  // Step must be valid
void zeroLayer(double** _mesh, int _rows, int _cols);
double** copyMesh(double** _mesh, int _rows, int _cols);
bool checkResult(double** _result, int _rows, int _cols, double _step);

// Boundaries
void leftBoundary(double** _mesh, int _rows, int _cols, double _boundValue);
void rightBoundary(double** _mesh, int _rows, int _cols, double _boundValue);
void topBoundary(double** _mesh, int _rows, int _cols, double _boundValue);
void bottomBoundary(double** _mesh, int _rows, int _cols, double _boundValue);
