#ifndef _MATRIX
#define _MATRIX

#include <iostream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#define THRESHOLD 0.000001
#define MAXITER 2000

using namespace std;

class Matrix
{
	public:
		Matrix();
		Matrix(int ,int);
		~Matrix();
		int init(int,int);
		int initAsIdentity();
		Matrix* addMatrix(Matrix*);
		Matrix* subtractMatrix(Matrix*);
		Matrix* multiplyMatrix(Matrix*);

		int addWithMatrix(Matrix*);
		int subtractWithMatrix(Matrix*);
		int multiplyWithMatrix(Matrix*);
		
		int addScalar(double);
		int subtractScalar(double);
		int multiplyScalar(double);
		int divideScalar(double);

		int setValue(double,int,int);
		int setAllValues(double);
		double getValue(int,int);

		Matrix* invMatrix(); 		
		Matrix* transMatrix();

		bool dimequal(Matrix*);
		int getRowCnt();
		int getColCnt();

		double detMatrix();
		int showMatrix(ostream& o=cout);
		int normalize();
		int normalizeVector();
		double getMax();
		int makeUncorrelated();
		bool rowZero();
		bool colZero();
		int makePositive();
		Matrix* copyMe();
		Matrix* findClosest();
		Matrix* getNextClosest(Matrix*,double);
		double getDistance(Matrix*);
	private:
		int convertToFloat(gsl_matrix_float*,const gsl_matrix*,int,int);
		int convertFromFloat(const gsl_matrix_float*, gsl_matrix*,int,int);
		gsl_matrix* matrix;
		int row;
		int col;
};
#endif
