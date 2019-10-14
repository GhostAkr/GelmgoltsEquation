#include "../include/GelmgoltsEquation.h"
#include "../include/Matr.h"

using namespace std;

// Additional

double** createMesh(double _xBorder, double _yBorder, double _step, double* _rows, double* _cols) {
	double zeroPoint = 0.0;
	// Step must be valid. Here is no exception for it!
	int rows = (_yBorder - zeroPoint) / _step + 1;
	int cols = (_xBorder - zeroPoint) / _step + 1;
	*_rows = rows;
	*_cols = cols;
	double** resMesh = createMatr(rows, cols);
	return resMesh;
}

void zeroLayer(double** _mesh, int _rows, int _cols) {
	double zeroValue = 0.0;
	for (int i = 0; i < _rows; ++i) {
		for (int j = 0; j < _cols; ++j) {
			_mesh[i][j] = zeroValue;
		}
	}
	leftBoundary(_mesh, _rows, _cols);
	rightBoundary(_mesh, _rows, _cols);
	topBoundary(_mesh, _rows, _cols);
	bottomBoundary(_mesh, _rows, _cols);
}

// Boundaries

void leftBoundary(double** _mesh, int _rows, int _cols, double _boundValue = 0.0) {
	for (int i = 0; i < _rows; ++i) {
		_mesh[i][0] = _boundValue;
	}
}

void rightBoundary(double** _mesh, int _rows, int _cols, double _boundValue = 0.0) {
	for (int i = 0; i < _rows; ++i) {
		_mesh[i][_cols - 1] = _boundValue;
	}
}

void topBoundary(double** _mesh, int _rows, int _cols, double _boundValue = 0.0) {
	for (int i = 0; i < _cols; ++i) {
		_mesh[_rows - 1][i] = _boundValue;
	}
}

void bottomBoundary(double** _mesh, int _rows, double _cols, double _boundValue = 0.0) {
	for (int i = 0; i < _cols; ++i) {
		_mesh[0][i] = _boundValue;
	}
}
