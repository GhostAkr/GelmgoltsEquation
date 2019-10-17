#include "../include/Matr.h"
#include "../include/GelmgoltsEquation.h"
#include <omp.h>

using namespace std;

int main() {
	// Init block
	double rows = 0.0;
	double cols = 0.0;
	double step = 0.01;
	double** mesh1 = createMesh(1, 1, step, &rows, &cols);
	double** mesh2 = copyMesh(mesh1, rows, cols);
	double** mesh3 = copyMesh(mesh1, rows, cols);
	double** mesh4 = copyMesh(mesh1, rows, cols);
	double t1 = 0.0, t2 = 0.0;

	// Jacobi
	cout << "Jacobi" << endl;
	t1 = omp_get_wtime();
	Jacobi(mesh1, rows, cols, 0, step);
	t2 = omp_get_wtime();
	cout << "Time spent (Jacobi): " << t2 - t1 << endl << endl;

	// Parallel Jacobi
	/*cout << "Parallel Jacobi" << endl;
	t1 = omp_get_wtime();
	JacobiParal(mesh2, rows, cols, 0, step);
	t2 = omp_get_wtime();
	cout << "Time spent (Parallel Jacobi): " << t2 - t1 << endl << endl;*/

	// Zeidel
	cout << "Zeidel" << endl;
	t1 = omp_get_wtime();
	Zeidel(mesh3, rows, cols, 1, step);
	t2 = omp_get_wtime();
	cout << "Time spent (Zeidel):" << t2 - t1 << endl << endl;

	// Parallel Zeidel
	/*cout << "Parallel Zeidel" << endl;
	t1 = omp_get_wtime();
	ZeidelParal(mesh4, rows, cols, 1, step);
	t2 = omp_get_wtime();
	cout << "Time spent (Parallel Zeidel):" << t2 - t1 << endl << endl;*/
	
	// Cleanup
	deleteMatr(mesh1, rows);
	deleteMatr(mesh2, rows);
	deleteMatr(mesh3, rows);
	deleteMatr(mesh4, rows);
	system("pause");
}