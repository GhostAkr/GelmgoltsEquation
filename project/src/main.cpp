#include "../include/Matr.h"
#include "../include/GelmgoltsEquation.h"

using namespace std;

int main() {
	double rows = 0.0;
	double cols = 0.0;
	double** testMesh = createMesh(1, 1, 0.5, &rows, &cols);
	printMatr(testMesh, rows, cols);
	deleteMatr(testMesh, rows);
	system("pause");
}