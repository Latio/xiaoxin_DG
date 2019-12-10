#pragma once
#include"NdgQuadFreeStrongFormAdvSolver2d.h"
#include"AbstractOutputFile.h"

extern "C" {

}
class NdgPhysMat
{
public:
	NdgPhysMat();
	~NdgPhysMat();

	void matSolver();
	void matEvaluateSSPRK22();
	void UpdateExternalField(double tloc);
	void EvaluateRHS(double *fphys, double *frhs);
	void UpdateOutputResult(double& time, double *fphys,int Nvar);


	//double* EvaluatePostFunc(double *fphys);
	//double UpdateTimeInterval(double *fphys);
	//void EvaluateSourceTerm(double *fphys);
	
protected:

	SWEAbstract2d sweabstract2d;
	SWEConventional2d sweconventional2d;
	NdgQuadFreeStrongFormAdvSolver2d ndgquadfreestrongformadvsolver2d;
	AbstractOutputFile abstractoutputfile;
	
	//double *fext;
	//int *outputfile;
	//double *limiter;
	//double *davectionSolver;
	//double *viscositySolver;
	//double *NonhydrostaticSolver;
	//double ftime;
	std::string casename;

	double *frhs;
	double *fext;
	double *fphys0;
	double *fphys;
	double *zGrad;

	double ftime;
	int outputIntervalNum;
	int *Np;
	int *K;
	int Nfield;
	int *Nv;
	int Nvar;
	int *boundarydge_Nfp;
	int *boundarydge_Ne;
	//double gra;
	//double hmin;
	//int Np;// dimension 
	//int K;// dimension 
	//int Nfield;// dimension 
	//int Nvar;// dimension 

};

