#include "../include/Matr.h"
#include "../include/GelmgoltsEquation.h"
#include <omp.h>

using namespace std;

int main() {
	// Init block
	double rows = 0.0;
	double cols = 0.0;
	double step = 0.01;
	double** mesh = createMesh(1, 1, step, &rows, &cols);
	double** mesh1 = copyMesh(mesh, rows, cols);
	double t1 = 0.0, t2 = 0.0;

	// Jacobi
	cout << "Jacobi" << endl;
	t1 = omp_get_wtime();
	Jacobi(mesh, rows, cols, 0, step);
	t2 = omp_get_wtime();
	cout << "Time spent (Jacobi): " << t2 - t1 << endl;

	// Zeidel
	cout << "Zeidel" << endl;
	t1 = omp_get_wtime();
	Zeidel(mesh1, rows, cols, 1, step);
	t2 = omp_get_wtime();
	cout << "Time spent (Zeidel):" << t2 - t1 << endl;
	
	// Cleanup
	deleteMatr(mesh, rows);
	deleteMatr(mesh1, rows);
	system("pause");
}