//#include<iostream>
//#include"cblas.h"
//
//using namespace std;
//
//int main()
//{
//	double dt[5] = { 1,3,5,7,9 };
//	double *p = dt;
//	const double *pg = p;
//	p++;
//
//	cout << pg << endl << p << endl;
//	system("pause");
//	return 0;
//}
//
//class test
//{
//public:
//	static const double pi ;
//	test();
//	~test();
//
//private:
//
//};
//
//
//test::test()
//{
//}
//
//test::~test()
//{
//}

//#include<cblas.h>   // <strong>����cblas.h�ļ��Ѿ�����������Ŀ¼�У�ֻ����˫���� </strong>
//#include<stdio.h>
//int main(void) {
//	const enum CBLAS_ORDER Order = CblasRowMajor;
//	const enum CBLAS_TRANSPOSE TransA = CblasNoTrans;
//	const enum CBLAS_TRANSPOSE TransB = CblasNoTrans;
//	const int M = 4;//A��������C������
//	const int N = 2;//B��������C������
//	const int K = 3;//A��������B������
//	const float alpha = 1;
//	const float beta = 0;
//	const int lda = K;//A����
//	const int ldb = N;//B����
//	const int ldc = N;//C����
//	const float A[12] = { 1,2,3,4,5,6,7,8,9,8,7,6 };
//	const float B[6] = { 5,4,3,2,1,0 };
//	float C[8];
//
//	cblas_sgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
//
//	for (int i = 0; i < M; i++)
//	{
//		for (int j = 0; j < N; j++)
//		{
//			printf("c:%f",C[i*N + j]);
//			
//		}
//	}
//
//	return 0;
//}


//	
//	//{
//	//	
//	//	double	*ptr;
//	//	ptr = x.EToV;
//	//	for (int i = 0; i < 160; i++)
//	//	{
//	//		cout << *ptr++ << endl;
//	//	}
//
//	//}
//
////	cout << "This is the mesh file test" << endl;
////
////	cout << *(x.BoundaryEdge_FToE + 6) << endl;
////
////	cout << "This is the object x" << endl;
////
////
////	MeshUnion * testMesh;
////
////	cout << *(testMesh->BoundaryEdge_FToE + 6) << endl;
////	cout << "Test finished" << endl;
////
////	cout << x.K <<",,"<<x.Np<< endl;
////	cout << x.cell_Dr << endl;
////
////	//x.ncvar_read(x.BoundaryEdge_FToE,"BoundaryEdge_FToE",x.Ne_boundary,x.two);
////	/*x.ncvar_read(x.BoundaryEdge_FToE,"BoundaryEdge_FToE");*/
////	for (int i = 0; i < x.Np*x.Np; i++)
////	{
////		cout << *((x.cell_Dr)++)<< endl;
////	}
////	cout << x.cell_Dr << endl;
//
//
//#include <iostream>
//#include <netcdfcpp.h>
//#include <vector>
//using namespace std;
//using namespace netCDF;
//using namespace netCDF::exceptions;
//
//static const int NX = 6;
//static const int NY = 12;
//
//static const int NC_ERR = 2;
//
//int main()
//{
//	int dataOut[NX][NY];
//
//	for (int i = 0; i < NX; i++)
//		for (int j = 0; j < NY; j++)
//			dataOut[i][j] = i * NY + j;
//
//	try
//	{
//		NcFile dataFile("simple_xy.nc", NcFile::replace);
//
//		NcDim xDim = dataFile.addDim("x", NX);
//		NcDim yDim = dataFile.addDim("y", NY);
//
//		vector<NcDim> dims;
//		dims.push_back(xDim);
//		dims.push_back(yDim);
//		NcVar data = dataFile.addVar("data", ncInt, dims);
//
//		data.putVar(dataOut);
//		system("pause");
//		return 0;
//	}
//	catch (NcException& e)
//	{
//		e.what();
//		return NC_ERR;
//	}
//}
////
////#include<iostream>
////#include <netcdfcpp.h>
////using namespace std;
////using namespace netCDF;
////using namespace netCDF::exceptions;
////
////static const int NX = 6;
////static const int NY = 12;
////
////static const int NC_ERR = 2;
////
////int main()
////{
////	try
////	{
////		int dataIn[NX][NY];
////
////		NcFile dataFile("simple_xy.nc", NcFile::read);
////
////		NcVar data = dataFile.getVar("data");
////		if (data.isNull()) return NC_ERR;
////		data.getVar(dataIn);
////
////		for (int i = 0; i < NX; i++)
////			for (int j = 0; j < NY; j++)
////			{
////				cout << dataIn[i][j] << endl;
////
////			}
////
////		system("pause");
////		return 0;
////	}
////	catch (NcException& e)
////	{
////		e.what();
////		cout << "FAILURE*************************************" << endl;
////		return NC_ERR;
////	}
////}
////
////
////
////
////#include<cblas.h>   // <strong>����cblas.h�ļ��Ѿ�����������Ŀ¼�У�ֻ����˫���� </strong>
////#include<iostream>
////using namespace std;
////int main(void) {
////	const enum CBLAS_ORDER Order = CblasRowMajor;
////	const enum CBLAS_TRANSPOSE TransA = CblasNoTrans;
////	const enum CBLAS_TRANSPOSE TransB = CblasNoTrans;
////	const int M = 4;//A��������C������
////	const int N = 2;//B��������C������
////	const int K = 3;//A��������B������
////	const float alpha = 1;
////	const float beta = 0;
////	const int lda = K;//A����
////	const int ldb = N;//B����
////	const int ldc = N;//C����
////	const float A[K*M] = { 1,2,3,4,5,6,7,8,9,8,7,6 };
////	const float B[K*N] = { 5,4,3,2,1,0 };
////	float C[M*N];
////
////	cblas_sgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
////
////	for (int i = 0; i < M; i++)
////	{
////		for (int j = 0; j < N; j++)
////		{
////			cout << C[i*N + j] << "  " ;
////			
////		}
////		cout << endl;
////	}
////	system("pause");
////
////	return EXIT_SUCCESS;
////}
////
//#include<cblas.h>   // <strong>����cblas.h�ļ��Ѿ�����������Ŀ¼�У�ֻ����˫���� </strong>
//#include<iostream>
//using namespace std;
//int main(void) {
//	const enum CBLAS_ORDER Order = CblasRowMajor;
//	const enum CBLAS_TRANSPOSE TransA = CblasNoTrans;
//	const enum CBLAS_TRANSPOSE TransB = CblasNoTrans;
//	const int M = 4;//A��������C������
//	const int N = 2;//B��������C������
//	const int K = 3;//A��������B������
//	const float alpha = 1;
//	const float beta = 0;
//	const int lda = M;//A����        //��������ȣ��ֱ�дABC����
//	const int ldb = K;//B����
//	const int ldc = N;//C����
//	const double A[K*M] = { 1,4,7,8,2,5,8,7,3,6,9,6 };
//	const double B[K*N] = { 5,3,1,4,2,0 };
//	double C[M*N];
//
//	cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);




//	int a[5] = { 1,2,3,4,5 };
//	int b[8] = { 8,7,6,5,4,3,2,1 };
//	int *pa = a;
//	int *pb = b;
//	double *c = C;
//
//	double *testcopy1 = new double[10]();
//	int t1 = 2;
//	int t2 = 4;
//	cblas_dcopy(t1*t2, c, 1, testcopy1, 1);
//
//
//	for (int i = 0; i < 10; i++)
//	{
//		cout << *(testcopy1 + i) << endl;
//	}
//	delete[]testcopy1;
//	//for (int i = 0; i < 5; i++)
//	//{
//	//	*(a + i) = *(a + i) + (*(b + i));
//	//	cout << *(a + i) << endl;
//	//}
//
//
//	//for (int i = 0; i < M; i++)
//	//{
//	//	for (int j = 0; j < N; j++)
//	//	{
//	//		cout << C[i*N + j] << "  ";
//
//	//	}
//	//	cout << endl;
//	//}
//	system("pause");
//
//	return EXIT_SUCCESS;
//}
//
//
////#include"cblas.h"
////#include <stdio.h>
////#include<iostream>
////extern "C"
////using namespace std;
////
////void main()
////{
////	int i = 0;
////	double A[6] = { 1.0, 2.0, 1.0, -3.0, 4.0, -1.0 };
////	double B[6] = { 1.0, 2.0, 1.0, -3.0, 4.0, -1.0 };
////	double C[9] = { .5, .5, .5, .5, .5, .5, .5, .5, .5 };
////	double D[9] = { .5, .5, .5, .5, .5, .5, .5, .5, .5 };
////	double E[9] = { .5, .5, .5, .5, .5, .5, .5, .5, .5 };
////	double F[9] = { .5, .5, .5, .5, .5, .5, .5, .5, .5 };
////
////	/*��������չ��*/
////	//1������ת��
////	cout << "��������չ��,����ת�ã�" << endl;
////	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 3, 3, 2, 1, A, 3, B, 2, 0, C, 3);
////	for (i = 0; i < 9; i++)
////		printf("%lf ", C[i]);
////	printf("\n");
////
////	//2\����Bת��
////	cout << "��������չ��������Bת�ã�" << endl;
////	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, 3, 3, 2, 1, A, 3, B, 3, 0, D, 3);
////	for (i = 0; i < 9; i++)
////		printf("%lf ", D[i]);
////	printf("\n");
////
////	/*��������չ��*/
////	//1������ת��
////	cout << "��������չ��,����ת��:" << endl;
////	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 2, 1, A, 2, B, 3, 0, E, 3);
////	for (i = 0; i < 9; i++)
////		printf("%lf ", E[i]);
////	printf("\n");
////
////	//2������Bת��
////	cout << "��������չ��,����Bת��:" << endl;
////	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 3, 3, 2, 1, A, 2, B, 2, 0, F, 3);
////	for (i = 0; i < 9; i++)
////		printf("%lf ", F[i]);
////	printf("\n");
////
////	system("pause");
////}
//
////
////#include <iostream> 
////using namespace std;
////
////class Shape {
////protected:
////	int width, height;
////public:
////	Shape(int a = 0, int b = 0)
////	{
////		width = a;
////		height = b;
////	}
////	virtual int area()
////	{
////		cout << "Parent class area :" << endl;
////		return 0;
////	}
////};
////class Rectangle : public Shape {
////public:
////	Rectangle(int a = 0, int b = 0) :Shape(a, b) { }
////	int area()
////	{
////		cout << "Rectangle class area :" << endl;
////		return (width * height);
////	}
////};
////class Triangle : public Shape {
////public:
////	Triangle(int a = 0, int b = 0) :Shape(a, b) { }
////	int area()
////	{
////		cout << "Triangle class area :" << endl;
////		return (width * height / 2);
////	}
////};
////
////class Cuboid:public Rectangle
////{
////public:
////	Cuboid(int a = 0, int b = 0, int c = 1) :Rectangle(a,b), length(c) {};
////	~Cuboid() {};
////	int area()
////	{
////		cout << "Cuboid class area:" << endl;
////		return (width*height*length);
////	}
////protected:
////	int length;
////};
////
////
////
////// �����������
////int main()
////{
////	Shape *shape;
////	Rectangle rec(10, 7);
////	Triangle  tri(10, 5);
////	Cuboid cub(2, 3, 4);
////	// �洢���εĵ�ַ
////	shape = &rec;
////	// ���þ��ε���������� area
////
////	cout << shape->area()<<endl;
////
////	// �洢�����εĵ�ַ
////	shape = &tri;
////	// ���������ε���������� area
////	
////	cout << shape->area();
////
////	shape = &cub;
////	cout << shape->area();
////
////	system("pause");
////	return 0;
////
////
////#include <iostream>
////
////using namespace std;
////const int MAX = 3;
////
////int main()
////{
////	int  var[MAX] = { 10, 100, 200 };
////	int  *ptr;
////
////	// ָ���е������ַ
////	ptr = var;
////	for (int i = 0; i < 5; i++)
////	{
////		cout << "Address of var[" << i << "] = ";
////		cout << ptr << endl;
////
////		cout << "Value of var[" << i << "] = ";
////		cout << *ptr << endl;
////
////		// �ƶ�����һ��λ��
////		ptr++;
////	}
////	cout << var << endl;
////	system("pause");
////	return 0;
////}
//
//
//
//

//#include<cblas.h>   // <strong>����cblas.h�ļ��Ѿ�����������Ŀ¼�У�ֻ����˫���� </strong>
//#include<iostream>
//using namespace std;
//int main(void) {
//	const enum CBLAS_ORDER Order = CblasRowMajor;
//	const enum CBLAS_TRANSPOSE TransA = CblasNoTrans;
//	const enum CBLAS_TRANSPOSE TransB = CblasNoTrans;
//	const int M = 4;//A��������C������
//	const int N = 2;//B��������C������
//	const int K = 3;//A��������B������
//	const float alpha = 1;
//	const float beta = 0;
//	const int lda = K;//A����
//	const int ldb = N;//B����
//	const int ldc = N;//C����
//	const float A[K*M] = { 1,2,3,4,5,6,7,8,9,8,7,6 };
//	const float B[K*N] = { 5,4,3,2,1,0 };
//	float C[M*N];
//
//	cblas_sgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
//
//	for (int i = 0; i < M; i++)
//	{
//		for (int j = 0; j < N; j++)
//		{
//			cout << C[i*N + j] << "  ";
//
//		}
//		cout << endl;
//	}
//	system("pause");
//
//	return EXIT_SUCCESS;
//}
//
//#include<cblas.h>   
//#include<iostream>
//using namespace std;
//int main(void) {
//	const enum CBLAS_ORDER Order = CblasRowMajor;
//	const enum CBLAS_TRANSPOSE TransA = CblasNoTrans;
//	const enum CBLAS_TRANSPOSE TransB = CblasNoTrans;
//	const int M = 4;//A��������C������
//	const int N = 2;//B��������C������
//	const int K = 3;//A��������B������
//	const float alpha = 1;
//	const float beta = 0;
//	const int lda = M;//A����        //��������ȣ��ֱ�дABC����
//	const int ldb = K;//B����
//	const int ldc = M;//C����
//	const float A[K*M] = { 1,4,7,8,2,5,8,7,3,6,9,6 };
//	const float B[K*N] = { 5,3,1,4,2,0 };
//	float C[M*N];
//
//	cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
//		M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
//
//	for (int i = 0; i < M; i++)
//	{
//		for (int j = 0; j < N; j++)
//		{
//			cout << C[i*N + j] << "  ";
//
//		}
//		cout << endl;
//	}
//	system("pause");
//
//	return EXIT_SUCCESS;
//}

//#include"cblas.h"
//#include<stdio.h>
//#include<stdlib.h>
//
//
//int main() {
//	double a[5] = { 1.1,3.3,5.5,7.8,9.1 };
//	double b[5] = { 0 };
//	//for (int i = 0; i < 5; i++)
//	//{
//	//	*(b + i) = a[i] + 10;
//	//}
//	int dis = 3;
//	int disss = 2;
//	double *temp_a = a + disss;
//	double *temp_b = b + disss;
//	 
//	int num = 5;
//	double alpha = 2.0;
//	
//	cblas_daxpy(num,alpha,a,1,a,1);
//	//cblas_dcopy(dis, temp_a, 1, temp_b, 1);
//	for (int i = 0; i < 5; i++)
//	{
//		printf("a:%f,b:%f", a[i], b[i]);
//	}
//	system("pause");
//}
//


//
//void printdata()
//{
//	for (int i = 0; i < 1080; i++)
//	{
//		std::cout << *(meshunion->EToV+i)<<std::endl;
//	}
//	


//#include<stdlib.h>
//#include<stdio.h>


//int main()
//{
//	const int a = 5;
//	typedef int(myarry)[a];
//	myarry xxx;
//
//	for (int i = 0; i < 5; i++)
//	{
//		xxx[i] = i * i;
//	};
//	for (int j = 0; j < 5; j++)
//	{
//		printf("xxx[%d]:%d", j, xxx);
//	};
//
//	system("pause");
//	return 0;
//}
