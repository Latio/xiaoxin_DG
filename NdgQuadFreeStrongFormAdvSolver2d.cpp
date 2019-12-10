#include "NdgQuadFreeStrongFormAdvSolver2d.h"

NdgQuadFreeStrongFormAdvSolver2d::NdgQuadFreeStrongFormAdvSolver2d()
{
}

NdgQuadFreeStrongFormAdvSolver2d::~NdgQuadFreeStrongFormAdvSolver2d()
{
}

void dotmul(const int num, double *const matrix1, double *matrix2)
{

	for (int i = 0; i < num; i++)
	{
		*(matrix2 + i) = (*(matrix1 + i))*(*(matrix2 + i));
	}
	return;
}

void multiply(double *const matrix1, double *const matrix2, double *const result)
{
	const enum CBLAS_ORDER Order = CblasColMajor;
	const enum CBLAS_TRANSPOSE TransA = CblasNoTrans;
	const enum CBLAS_TRANSPOSE TransB = CblasNoTrans;
	const int M = *meshunion->cell_p->Np;//A��������C������
	const int N = *meshunion->K;//B��������C������
	const int L = *meshunion->cell_p->Np;//A��������B������
	const double alpha = 1.0;
	const double beta = 0.0;
	const int lda = M;//A����        
	const int ldb = L;//B����
	const int ldc = M;//C����   //��������ȣ��ֱ�дABC����

	cblas_dgemm(Order, TransA, TransB, M, N, L, alpha, matrix1, lda, matrix2, ldb, beta, result, ldc);
	//std::cout << "NdgQuadFreeStrongFormAdvSolver2d.cpp" << std::endl;
}


void NdgQuadFreeStrongFormAdvSolver2d::evaluateAdvectionRHS(double *fphys, double *frhs, double *fext)
{
	int Nfield = meshunion->Nfield;
	int *Np = meshunion->cell_p->Np;
	int *K = meshunion->K;
	double *invM = meshunion->cell_p->invM;
	double *J = meshunion->J;
	int *Nfp = meshunion->inneredge_p->Nfp;
	int *Ne = meshunion->inneredge_p->Ne;
	int *Nfp_b = meshunion->boundarydge_p->Nfp;
	int *Ne_b = meshunion->boundarydge_p->Ne;
	const int NVAR = 3;

	double *fm, *fp, *fluxM, *fluxP, *fluxS;
	requestmemory(&fm, Nfp, Ne, Nfield);
	requestmemory(&fp, Nfp, Ne, Nfield);
	requestmemory(&fluxM, Nfp, Ne, NVAR);
	requestmemory(&fluxP, Nfp, Ne, NVAR);
	requestmemory(&fluxS, Nfp, Ne, NVAR);

	// evaluate inner edge
	double *nx = meshunion->inneredge_p->nx;
	double *ny = meshunion->inneredge_p->ny;
	mesh.inneredge.EvaluateSurfValue(fphys, fm, fp, Np, K, Nfield);//return fm,fp
	sweabstract2d.EvaluateSurfFlux(nx, ny, fm, fluxM, Nfp, Ne);//return fluxM
	sweabstract2d.EvaluateSurfFlux(nx, ny, fp, fluxP, Nfp, Ne);//return fluxP
	sweabstract2d.EvaluateSurfNumFlux(nx, ny, fm, fp, fluxS, Nfp, Ne);//retuen fluxS
	mesh.inneredge.EvaluateStrongFromEdgeRHS(fluxM, fluxP, fluxS, frhs, invM, J, Np, K, NVAR);

	freememory(&fm);
	freememory(&fp);
	freememory(&fluxM);
	freememory(&fluxP);
	freememory(&fluxS);

	requestmemory(&fm, Nfp_b, Ne_b, Nfield);
	requestmemory(&fp, Nfp_b, Ne_b, Nfield);
	requestmemory(&fluxM, Nfp_b, Ne_b, NVAR);
	//requestmemory(&fluxP, Nfp_b, Ne_b, NVAR);
	requestmemory(&fluxS, Nfp_b, Ne_b, NVAR);

	// evaluate boundary edge
	nx = meshunion->boundarydge_p->nx;
	ny = meshunion->boundarydge_p->ny;
	mesh.boundarydge.EvaluateSurfValue(fphys, fm, fp, Np, K, Nfield);
	sweabstract2d.ImposeBoundaryCondition(nx, ny, fm, fp, fext);
	sweabstract2d.EvaluateSurfFlux(nx, ny, fm, fluxM, Nfp_b, Ne_b);//return fluxM
	sweabstract2d.EvaluateSurfNumFlux(nx, ny, fm, fp, fluxS, Nfp_b, Ne_b);

	double *frhs_temp;
	requestmemory(&frhs_temp, Np, K, NVAR);

	mesh.boundarydge.EvaluateStrongFromEdgeRHS(invM, J, fluxM, fluxS, Np, K, NVAR, frhs_temp);

	const int num = (*Np)*(*K)*NVAR;
	const double alpha = 1;
	const int dis1 = 1;
	const int dis2 = 1;
	cblas_daxpy(num, alpha, frhs_temp, dis1, frhs, dis2);


	freememory(&frhs_temp);
	freememory(&fm);
	freememory(&fp);
	freememory(&fluxM);
	//freememory(&fluxP);
	freememory(&fluxS);

	double *E;
	double *G;
	requestmemory(&E, Np, K, NVAR);
	requestmemory(&G, Np, K, NVAR);

	swepreblanaced2d.EvaluateFlux(fphys, E, G);

	double *const Dr = meshunion->cell_p->Dr;
	double *const Ds = meshunion->cell_p->Ds;
	double *const rx = meshunion->rx;
	double *const ry = meshunion->ry;
	double *const sx = meshunion->sx;
	double *const sy = meshunion->sy;

	double *rx_dr_e, *sx_ds_e, *ry_dr_g, *sy_ds_g;


	//cblas_dgemm(chn, chn, &np, &oneI, &np, &one, invM, &np, rhs_, &np, &zero, temp, &np);

	//for i = 1:phys.Nvar
	//	phys.frhs{ m }(:, : , i) = ...
	//	phys.frhs{ m }(:, : , i) + ...
	//	- obj.rx{ m }.*(obj.Dr{ m } *E(:, : , i)) ...
	//	- obj.sx{ m }.*(obj.Ds{ m } *E(:, : , i)) ...
	//	- obj.ry{ m }.*(obj.Dr{ m } *G(:, : , i)) ...
	//	- obj.sy{ m }.*(obj.Ds{ m } *G(:, : , i)); ...

	int dis = (*Np)*(*K);
	double alpha_ = -1.0;
	for (int i = 0; i < NVAR; i++)
	{
		requestmemory(&rx_dr_e, Np, K);
		requestmemory(&sx_ds_e, Np, K);
		requestmemory(&ry_dr_g, Np, K);
		requestmemory(&sy_ds_g, Np, K);

		double *frhs_ = frhs + i * dis;
		double *E_ = E + i * dis;
		double *G_ = G + i * dis;

		multiply(Dr, E_, rx_dr_e);
		multiply(Ds, E_, sx_ds_e);
		multiply(Dr, G_, ry_dr_g);
		multiply(Ds, G_, sy_ds_g);

		dotmul(dis, rx, rx_dr_e);
		dotmul(dis, sx, sx_ds_e);
		dotmul(dis, ry, ry_dr_g);
		dotmul(dis, sy, sy_ds_g);

		//	cblas_daxpy(num,alpha,a,1,a,1);
		//cblas_dscal(dis, alpha, rx_dr_e, 1);
		//cblas_dscal(dis, alpha, sx_ds_e, 1);
		//cblas_dscal(dis, alpha, ry_dr_g, 1);
		//cblas_dscal(dis, alpha_, sy_ds_g, 1);

		cblas_daxpy(dis, alpha, rx_dr_e, 1, sx_ds_e, 1);
		cblas_daxpy(dis, alpha, sx_ds_e, 1, ry_dr_g, 1);
		cblas_daxpy(dis, alpha, ry_dr_g, 1, sy_ds_g, 1);
		cblas_daxpy(dis, alpha_, sy_ds_g, 1, frhs_, 1);

		freememory(&rx_dr_e);
		freememory(&sx_ds_e);
		freememory(&ry_dr_g);
		freememory(&sy_ds_g);
	}
	freememory(&E);
	freememory(&G);


	//for (int i = 0; i < 1080 * 3; i++)
	//{
	//	std::cout << " frhs[" << i << "];  " << frhs[i] << std::endl;
	//}//test

};


