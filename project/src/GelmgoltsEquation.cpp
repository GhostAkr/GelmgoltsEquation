#include "../include/GelmgoltsEquation.h"
#include "../include/Matr.h"
#include <cmath>

#define Pi 3.14

using namespace std;

// Main

double f(double _x, double _y, double _k) {
	return 2.0 * sin(Pi * _y) + _k * _k * (1 - _x) * _x * sin(Pi * _y) + \
		Pi * Pi * (1 - _x) * _x * sin(Pi * _y);
}

double exactSolution(double _x, double _y) {
	return (1 - _x) * _x * sin(Pi * _y);
}

double** rightPart(double _step, double _rows, double _cols, double _k) {
	double** result = createMatr(_rows, _cols);
	for (int i = 1; i < _rows - 1; ++i) {
		for (int j = 1; j < _cols - 1; ++j) {
			result[i][j] = f(i * _step, j * _step, _k);
		}
	}
	double boundaryValue = 0.0;
	/*leftBoundary(result, _rows, _cols, boundaryValue);
	rightBoundary(result, _rows, _cols, boundaryValue);
	topBoundary(result, _rows, _cols, boundaryValue);
	bottomBoundary(result, _rows, _cols, boundaryValue);*/
	return result;
}

void Jacobi(double** _mesh, int _rows, int _cols, double _k, double _step) {
	double** rPart = rightPart(_step, _rows, _cols, _k);
	//double** previousLayer = nullptr;
	for (int s = 0; s < 42; ++s) {
		cout << "_mesh[1][1] = " << _mesh[1][1] << endl;
		cout << "_mesh[2][1] = " << _mesh[2][1] << endl;
		double** previousLayer = copyMesh(_mesh, _rows, _cols);
		for (int i = 1; i < _rows - 1; ++i) {
			for (int j = 1; j < _cols - 1; ++j) {
				_mesh[i][j] = 0.25 * (previousLayer[i - 1][j] + previousLayer[i + 1][j] + previousLayer[i][j - 1] + \
					previousLayer[i][j + 1] - _k * _k * _step * _step * previousLayer[i][j]) + rPart[i][j];
			}
		}
		deleteMatr(previousLayer, _rows);
	}
	if (checkResult(_mesh, _rows, _cols, _step)) {
		cout << "Answer is correct" << endl;
	}
	else {
		cout << "Answer is INcorrect" << endl;
	}
	deleteMatr(rPart, _rows);
	//deleteMatr(previousLayer, _rows);
}



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
	double boundaryValue = 0.0;
	/*leftBoundary(_mesh, _rows, _cols, boundaryValue);
	rightBoundary(_mesh, _rows, _cols, boundaryValue);
	topBoundary(_mesh, _rows, _cols, boundaryValue);
	bottomBoundary(_mesh, _rows, _cols, boundaryValue);*/
}

double** copyMesh(double** _mesh, int _rows, int _cols) {
	double** result = createMatr(_rows, _cols);
	for (int i = 0; i < _rows; ++i) {
		for (int j = 0; j < _cols; ++j) {
			_mesh[i][j] = result[i][j];
		}
	}
	return result;
}

bool checkResult(double** _result, int _rows, int _cols, double _step) {
	double epsNull = 1e-5;
	for (int i = 0; i < _rows; ++i) {
		for (int j = 0; j < _cols; ++j) {
			if (fabs(_result[i][j] - exactSolution(i * _step, j * _step)) > epsNull) {
				cout << "_result[i][j] = " << _result[i][j] << endl;
				cout << "Exact = " << exactSolution(i * _step, j * _step) << endl;
				cout << "i = " << i << endl;
				cout << "j = " << j << endl;
				return false;
			}
		}
	}
	return true;
}

// Boundaries

void leftBoundary(double** _mesh, int _rows, int _cols, double _boundValue) {
	for (int i = 0; i < _rows; ++i) {
		_mesh[i][0] = _boundValue;
	}
}

void rightBoundary(double** _mesh, int _rows, int _cols, double _boundValue) {
	for (int i = 0; i < _rows; ++i) {
		_mesh[i][_cols - 1] = _boundValue;
	}
}

void topBoundary(double** _mesh, int _rows, int _cols, double _boundValue) {
	for (int i = 0; i < _cols; ++i) {
		_mesh[_rows - 1][i] = _boundValue;
	}
}

void bottomBoundary(double** _mesh, int _rows, double _cols, double _boundValue) {
	for (int i = 0; i < _cols; ++i) {
		_mesh[0][i] = _boundValue;
	}
}
