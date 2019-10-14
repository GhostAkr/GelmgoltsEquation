#pragma once

#include <iostream>

using namespace std;

// Main
void Zeidel(double** _mesh, int _rows, int _cols, double _k);  // TODO: make Zeidel
void Jacobi(double** _mesh, int _rows, int _cols, double _k);  // TODO: make Jacobi
double rightPart(double _x, double _y, double _k);
double exactSolution(double _x, double _y);

// Additional
double** createMesh(double _xBorder, double _yBorder, double _step, double* _rows, double* _cols);  // Step must be valid
void zeroLayer(double** _mesh, int _rows, int _cols);

// Boundaries
void leftBoundary(double** _mesh, int _rows, int _cols, double _boundValue = 0.0);
void rightBoundary(double** _mesh, int _rows, int _cols, double _boundValue = 0.0);
void topBoundary(double** _mesh, int _rows, int _cols, double _boundValue = 0.0);
void bottomBoundary(double** _mesh, int _rows, int _cols, double _boundValue = 0.0);
