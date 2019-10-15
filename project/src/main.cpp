#include "../include/Matr.h"
#include "../include/GelmgoltsEquation.h"

using namespace std;

int main() {
	double rows = 0.0;
	double cols = 0.0;
	double step = 0.01;
	double** mesh = createMesh(1, 1, step, &rows, &cols);
	Jacobi(mesh, rows, cols, 1, step);
	//printMatr(mesh, rows, cols);
	deleteMatr(mesh, rows);
	system("pause");
}