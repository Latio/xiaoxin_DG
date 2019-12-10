#ifdef _OPENMP
#include <omp.h>
#endif

void myfree(double **arr);
//#include "mex.h"
#include "cblas.h"
#include"stdlib.h"

// #if !defined(_WIN32)
// #define dgemm dgemm_
// #endif
//
//#define DEBUG 0

//#define NRHS 8
//#define NLHS 1

void c_boundary_EvaluateStrongFormEdgeRHS(double *invM_, double *M_, double *FToE_, double *FToN1_, double *Js_, double *J_, double *fluxM__, double *fluxS__, int *Np_, int *K_, int *Nfp_, int *Ne_, int Nfield_,double *const frhs_temp_)
{

	//void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	//	/* check input & output */
	//	if (nrhs != NRHS) {
	//		mexPrintf("Matlab:%s:InvalidNumberInput,\n", __FILE__);
	//		mexPrintf("%d inputs required.\n", NRHS);
	//	}

	//	if (nlhs != NLHS) {
	//		mexPrintf("Matlab:%s:InvalidNumberOutput,\n", __FILE__);
	//		mexPrintf("%d inputs required.\n", NLHS);
	//	}

	double *invM = invM_;
	double *Mb = M_;
	double *FToE = FToE_;
	double *FToN1 = FToN1_;
	double *Js = Js_;
	double *J = J_;
	double *fluxM = fluxM__;
	double *fluxS = fluxS__;

	// dims = mxGetDimensions(prhs[6]);
	const int Np = *Np_;  // num of interp nodes
	const int K = *K_;   // num of elements
	//const mwSize *dims = mxGetDimensions(prhs[6]);
	const int Nfp = *Nfp_;
	const int Ne = *Ne_;  // num of edges
	int Nfield = Nfield_;

	//if (mxGetNumberOfDimensions(prhs[6]) > 2) {
	//	Nfield = dims[2];
	//}
	//else {
	//	Nfield = 1;  // fluxM is a 2D matrix
	//}

	/*const size_t ndimOut = 3;
	const mwSize dimOut[3] = { Np, K, Nfield };

	plhs[0] = mxCreateNumericArray(ndimOut, dimOut, mxDOUBLE_CLASS, mxREAL);*/
	double *frhs = frhs_temp_;

	//char *chn = "N";
	//double one = 1.0, zero = 0.0;
	int oneI = 1;
	//int np = Np;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int fld = 0; fld < Nfield; fld++) {
		double *rhs = frhs + Np * K * fld;
		double *fluxM_ = fluxM + Nfp * Ne * fld;
		// double *fluxP_ = fluxP + Nfp * Ne * fld;
		double *fluxS_ = fluxS + Nfp * Ne * fld;

		double *rhsM = (double*)malloc(sizeof(double)*Nfp);
		
		for (int k = 0; k < Ne; k++) {  // evaluate rhs on each edge
			const int e1 = (int)FToE[2 * k] - 1;
			// const int e2 = (int)FToE[2 * k + 1] - 1;
			const int ind1 = e1 * Np - 1;
			// const int ind2 = e2 * Np - 1;
			const int ind = k * Nfp;

			//double rhsM[Nfp];  // rhsP[Nfp];//

			for (int n = 0; n < Nfp; n++) {
				rhsM[n] = 0;
				// rhsP[n] = 0;
			}

			for (int n = 0; n < Nfp; n++) {
				const int sk = n + ind;
				double dfM = fluxM_[sk] - fluxS_[sk];
				// double dfP = fluxP_[sk] - fluxS_[sk];
				double j = Js[sk];

				double *mb = Mb + n * Nfp;

				for (int m = 0; m < Nfp; m++) {
					rhsM[m] += mb[m] * j * dfM;
					// rhsP[m] -= mb[m] * j * dfP;
				}
			}

			for (int n = 0; n < Nfp; n++) {
				const int sk = n + ind;
				const int m1 = (int)FToN1[sk] + ind1;
				// const int m2 = (int)FToN2[sk] + ind2;
				rhs[m1] += rhsM[n];
				// rhs[m2] += rhsP[n];
			}

		}


		//double temp[Np];
		//double *temp = (double*)malloc(sizeof(double)*Np);
		double *temp = (double*)malloc(sizeof(double)*Np);
		for (int k = 0; k < K; k++) {
			double *rhs_ = rhs + k * Np;
			double *j = J + k * Np;


			const enum CBLAS_ORDER Order = CblasColMajor;
			const enum CBLAS_TRANSPOSE TransA = CblasNoTrans;
			const enum CBLAS_TRANSPOSE TransB = CblasNoTrans;
			const int M = Np;//A��������C������
			const int N = oneI;//B��������C������
			const int K = Np;//A��������B������
			const double alpha = 1.0;
			const float beta = 0.0;
			const int lda = M;//A����        
			const int ldb = K;//B����
			const int ldc = M;//C����   //��������ȣ��ֱ�дABC����

			//cblas_dgemm(chn, chn, &np, &oneI, &np, &one, invM, &np, rhs_, &np, &zero, temp, &np);
			cblas_dgemm(Order, TransA, TransB, M, N, K, alpha, invM, lda, rhs_, ldb, beta, temp, ldc);

			//printf( "c_boundary_EvaluateStrongFormEdgeRHS.c\n");

			//dgemm(chn, chn, &np, &oneI, &np, &one, invM, &np, rhs_, &np, &zero, temp,
			//	&np);

			// copy rhs
			for (int n = 0; n < Np; n++) {
				rhs_[n] = temp[n] / j[n];
			}

			
			
		}
		myfree(&temp);
		myfree(&rhsM);
	}
	return;
}


//void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
//  /* check input & output */
//  if (nrhs != NRHS) {
//    mexPrintf("Matlab:%s:InvalidNumberInput,\n", __FILE__);
//    mexPrintf("%d inputs required.\n", NRHS);
//  }
//
//  if (nlhs != NLHS) {
//    mexPrintf("Matlab:%s:InvalidNumberOutput,\n", __FILE__);
//    mexPrintf("%d inputs required.\n", NLHS);
//  }
//
//  double *invM = mxGetPr(prhs[0]);
//  double *Mb = mxGetPr(prhs[1]);
//  double *FToE = mxGetPr(prhs[2]);
//  double *FToN1 = mxGetPr(prhs[3]);
//  double *Js = mxGetPr(prhs[4]);
//  double *J = mxGetPr(prhs[5]);
//  double *fluxM = mxGetPr(prhs[6]);
//  double *fluxS = mxGetPr(prhs[7]);
//
//  // dims = mxGetDimensions(prhs[6]);
//  const int Np = mxGetM(prhs[5]);  // num of interp nodes
//  const int K = mxGetN(prhs[5]);   // num of elements
//  const mwSize *dims = mxGetDimensions(prhs[6]);
//  const int Nfp = dims[0];
//  const int Ne = dims[1];  // num of edges
//  int Nfield;
//
//  if (mxGetNumberOfDimensions(prhs[6]) > 2) {
//    Nfield = dims[2];
//  } else {
//    Nfield = 1;  // fluxM is a 2D matrix
//  }
//
//  const size_t ndimOut = 3;
//  const mwSize dimOut[3] = {Np, K, Nfield};
//
//  plhs[0] = mxCreateNumericArray(ndimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
//  double *frhs = mxGetPr(plhs[0]);
//
//  char *chn = "N";
//  double one = 1.0, zero = 0.0;
//  ptrdiff_t oneI = 1;
//  ptrdiff_t np = Np;
//
//#ifdef _OPENMP
//#pragma omp parallel for num_threads(DG_THREADS)
//#endif
//  for (int fld = 0; fld < Nfield; fld++) {
//    double *rhs = frhs + Np * K * fld;
//    double *fluxM_ = fluxM + Nfp * Ne * fld;
//    // double *fluxP_ = fluxP + Nfp * Ne * fld;
//    double *fluxS_ = fluxS + Nfp * Ne * fld;
//
//    for (int k = 0; k < Ne; k++) {  // evaluate rhs on each edge
//      const int e1 = (int)FToE[2 * k] - 1;
//      // const int e2 = (int)FToE[2 * k + 1] - 1;
//      const int ind1 = e1 * Np - 1;
//      // const int ind2 = e2 * Np - 1;
//
//      const int ind = k * Nfp;
//      double rhsM[Nfp];  // rhsP[Nfp];
//      for (int n = 0; n < Nfp; n++) {
//        rhsM[n] = 0;
//        // rhsP[n] = 0;
//      }
//
//      for (int n = 0; n < Nfp; n++) {
//        const int sk = n + ind;
//        double dfM = fluxM_[sk] - fluxS_[sk];
//        // double dfP = fluxP_[sk] - fluxS_[sk];
//        double j = Js[sk];
//
//        double *mb = Mb + n * Nfp;
//
//        for (int m = 0; m < Nfp; m++) {
//          rhsM[m] += mb[m] * j * dfM;
//          // rhsP[m] -= mb[m] * j * dfP;
//        }
//      }
//
//      for (int n = 0; n < Nfp; n++) {
//        const int sk = n + ind;
//        const int m1 = (int)FToN1[sk] + ind1;
//        // const int m2 = (int)FToN2[sk] + ind2;
//        rhs[m1] += rhsM[n];
//        // rhs[m2] += rhsP[n];
//      }
//    }
//    double temp[Np];
//    for (int k = 0; k < K; k++) {
//      double *rhs_ = rhs + k * Np;
//      double *j = J + k * Np;
//
//      dgemm(chn, chn, &np, &oneI, &np, &one, invM, &np, rhs_, &np, &zero, temp,
//            &np);
//      // copy rhs
//      for (int n = 0; n < Np; n++) {
//        rhs_[n] = temp[n] / j[n];
//      }
//    }
//  }
//  return;
//}