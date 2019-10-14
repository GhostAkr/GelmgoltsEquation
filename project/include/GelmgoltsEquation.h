#pragma once

#include <iostream>

using namespace std;

// Main
double** Zeidel(double** _mesh, int _rows, int _cols);  // TODO: make Zeidel
double** Jacobi(double** _mesh, int _rows, int _cols);  // TODO: make Jacobi

// Additional
double** createMesh(double _xBorder, double _yBorder, double _step, double* _rows, double* _cols);  // Step must be valid
void zeroLayer(double** _mesh, int _rows, int _cols);

// Boundaries
void leftBoundary(double** _mesh, int _rows, int _cols, double _boundValue = 0.0);
void rightBoundary(double** _mesh, int _rows, int _cols, double _boundValue = 0.0);
void topBoundary(double** _mesh, int _rows, int _cols, double _boundValue = 0.0);
void bottomBoundary(double** _mesh, int _rows, int _cols, double _boundValue = 0.0);
