#include "../include/GelmgoltsEquation.h"
#include "../include/Matr.h"
#include <cmath>
#include <omp.h>

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
	double c = 1 / (4 + _k * _k *_step*_step);
	//double** previousLayer = nullptr;
	double tt = 0.0;
	double tt1 = 0.0;
	double tt2 = 0.0;
	for (int s = 0; s < 9000; ++s) {
		//cout << "_mesh[1][1] = " << _mesh[1][1] << endl;
		//cout << "_mesh[2][1] = " << _mesh[2][1] << endl;
		double t1 = omp_get_wtime();
		double** previousLayer = copyMesh(_mesh, _rows, _cols);
		double t2 = omp_get_wtime();
		tt += t2 - t1;
		//cout << "previousLayer = " << previousLayer[1][1] << endl;
		t1 = omp_get_wtime();
		for (int i = 1; i < _rows - 1; ++i) {
			for (int j = 1; j < _cols - 1; ++j) {
				_mesh[i][j] = c * (previousLayer[i - 1][j] + previousLayer[i + 1][j] + previousLayer[i][j - 1] + \
					previousLayer[i][j + 1] + _step * _step*rPart[i][j]);
			}
		}
		t2 = omp_get_wtime();
		tt1 += t2 - t1;
		t1 = omp_get_wtime();
		deleteMatr(previousLayer, _rows);
		t2 = omp_get_wtime();
		tt2 += t2 - t1;
	}
	cout << "Time of deleting pointers: " << tt2 << endl;
	cout << "Time of cycle: " << tt1 << endl;
	cout << "Time of copying matrices: " << tt << endl;
	if (checkResult(_mesh, _rows, _cols, _step)) {
		cout << "Answer is correct" << endl;
	}
	else {
		cout << "Answer is INcorrect" << endl;
	}
	deleteMatr(rPart, _rows);
	//deleteMatr(previousLayer, _rows);
}

//void Jacobi(double** _mesh, int _rows, int _cols, double _k, double _step) {
//	double** rPart = rightPart(_step, _rows, _cols, _k);
//	double c = 1 / (4 + _k * _k *_step*_step);
//	double temp = 0;
//	//double** previousLayer = nullptr;
//	for (int s = 0; s < 100000; ++s) {
//		//cout << "_mesh[1][1] = " << _mesh[1][1] << endl;
//		//cout << "_mesh[99][99] = " << _mesh[99][99] << endl;
//		double** previousLayer = copyMesh(_mesh, _rows, _cols);
//		//cout << "previousLayer = " << previousLayer[1][1] << endl;
//		for (int i = 1; i < _rows / 2; ++i) {
//			for (int j = 1; j < _cols - 1; ++j) {
//				temp = c * (previousLayer[i - 1][j] + previousLayer[i + 1][j] + previousLayer[i][j - 1] + \
//					previousLayer[i][j + 1] + _step * _step*rPart[i][j]);
//				_mesh[i][j] = temp;
//				_mesh[_rows - 1 - i][_cols - 1 - j] = temp;
//				//cout << i << " " << _rows - 1 - i << endl;
//			}
//			for (int i = _rows / 2; i < _rows / 2 + 1; ++i) {
//				for (int j = 1; j < _cols - 1; ++j) {
//					_mesh[i][j] = c * (previousLayer[i - 1][j] + previousLayer[i + 1][j] + previousLayer[i][j - 1] + \
//						previousLayer[i][j + 1] + _step * _step*rPart[i][j]);
//				}
//			}
//
//		}
//		deleteMatr(previousLayer, _rows);
//	}
//	if (checkResult(_mesh, _rows, _cols, _step)) {
//		cout << "Answer is correct" << endl;
//	}
//	else {
//		cout << "Answer is INcorrect" << endl;
//	}
//	deleteMatr(rPart, _rows);
//	//deleteMatr(previousLayer, _rows);
//}

void Zeidel(double** _mesh, int _rows, int _cols, double _k, double _step) {
	double** rPart = rightPart(_step, _rows, _cols, _k);
	double c = 1.0 / (4.0 + _k * _k * _step * _step);
	for (int s = 0; s < 4800; ++s) {
		for (int i = 1; i < _rows - 1; ++i) {
			for (int j = 1; j < _cols - 1; ++j) {
				_mesh[i][j] = c * (_mesh[i - 1][j] + _mesh[i + 1][j] + _mesh[i][j - 1] + _mesh[i][j + 1] + _step * _step * rPart[i][j]);
			}
		}
	}
	if (checkResult(_mesh, _rows, _cols, _step)) {
		cout << "Answer is correct" << endl;
	}
	else {
		cout << "Answer is INcorrect" << endl;
	}
	deleteMatr(rPart, _rows);
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
			result[i][j] = _mesh[i][j];
		}
	}
	return result;
}

bool checkResult(double** _result, int _rows, int _cols, double _step) {
	double epsNull = 1e-2;
	for (int i = 0; i < _rows; ++i) {
		for (int j = 0; j < _cols; ++j) {
			if (fabs(_result[i][j] - exactSolution(i * _step, j * _step)) > epsNull) {
				cout << endl;
				cout << "_result[i][j] = " << _result[i][j] << endl;
				cout << "Exact = " << exactSolution(i * _step, j * _step) << endl;
				cout << "i = " << i << endl;
				cout << "j = " << j << endl;
				cout << endl;
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
