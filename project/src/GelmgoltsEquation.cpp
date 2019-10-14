#include "../include/GelmgoltsEquation.h"
#include "../include/Matr.h"

using namespace std;

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
