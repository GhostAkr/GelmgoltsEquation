#include "../include/Matr.h"
#include "../include/GelmgoltsEquation.h"
#include <omp.h>

using namespace std;

int main() {
	double rows = 0.0;
	double cols = 0.0;
	double step = 0.01;
	double** mesh = createMesh(1, 1, step, &rows, &cols);
	cout << "Jacobi" << endl;
	double t1 = omp_get_wtime();
	Jacobi(mesh, rows, cols, 1, step);
	double t2 = omp_get_wtime();
	cout << "Timesa spent (Jacobi): " << t2 - t1 << endl;
	cout << "Zeidel" << endl;
	t1 = omp_get_wtime();
	Zeidel(mesh, rows, cols, 1, step);
	t2 = omp_get_wtime();
	cout << "Time spent (Zeidel):" << t2 - t1 << endl;
	//printMatr(mesh, rows, cols);
	deleteMatr(mesh, rows);
	system("pause");
}