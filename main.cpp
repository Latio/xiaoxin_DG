#include"NdgPhysMat.h"
#include<time.h>

using namespace std;

MeshUnion mesh;
const MeshUnion *meshunion = &mesh;

int main()
{
	clock_t begintime, endtime;
	begintime = clock();

	NdgPhysMat Solver;
	Solver.matSolver();

	endtime = clock();

	cout << "\n\nRunning Time : " << endtime - begintime << " ms\n" << endl;
	cout << "hello";
	system("pause");
	return 0;

}

