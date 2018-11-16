#ifndef LSCM
#define LSCM

#include <vector>
#include <set>
using namespace std;


#include "Primitives.h"

/*
extern "C"
{
	#include "taucs.h"

	extern void taucs_ccs_free    (taucs_ccs_matrix*);
	extern int taucs_linsolve(taucs_ccs_matrix* A, void** F, int nrhs, void* X, void* B, char* options[], void* opt_arg[]);
	extern taucs_ccs_matrix* taucs_ccs_create               (int m, int n, int nnz, int flags);
}
*/


#include <mkl_pardiso.h>


/* PARDISO prototype. */
//extern int omp_get_max_threads();
//extern int PARDISO(void *, int *, int *, int *, int *, int *, double *, int *, int *, int*, int *, int *, int *, double *, double *, int*);



#define MAX2(x,y)		( (x)>(y) ? (x) : (y) )
#define MAX3(x,y,z)		MAX2( MAX2((x),(y)) , (z) )
#define SHIFT3(type, a, b, c) { type tmp; tmp = a; a = c; c = b; b = tmp; }


class LSCM_Variable
{
public:
	double m_value;
	int m_index;
	bool m_pinned;
};


class MatrixElement
{
public:
	double m_value;
	int m_index;
};


class MatrixRow
{
public:
	MatrixRow()
	{
		size = 0;
		capacity = 0;
		al = 0.0f;
		m_elements = NULL;
	}

	~MatrixRow()
	{
		size = 0;
		capacity = 0;
		al = 0.0f;
		if(m_elements != NULL) delete[] m_elements;
		m_elements = NULL;
	}

	void reset()
	{
		size = 0;
		capacity = 0;
		al = 0.0f;
		if(m_elements != NULL) delete[] m_elements;
		m_elements = NULL;
	}

	float al;
	int size;					// current number of elements
	int capacity;				// current capacity for elements
	MatrixElement* m_elements;
	void addElement(int Index, double Value);
	void AddRowColumn(int ColIndex, double Value);
	int sizeToDiag(int Index);
	void Sort();

private:
	void expand();
};






class LSCMap
{
public:
	LSCMap()
	{
		m_variables = NULL;
		MatrixBElements = NULL;
		m_matrixRows = NULL;
		m_MatrixXEntries = NULL;
	};

	~LSCMap()
	{
		if(m_variables != NULL) delete[] m_variables;
		if(MatrixBElements != NULL) delete[] MatrixBElements;
		if(m_MatrixXEntries != NULL) delete[] m_MatrixXEntries;
		if(m_matrixRows != NULL) delete[] m_matrixRows;

		m_variables = NULL;
		MatrixBElements = NULL;
		m_MatrixXEntries = NULL;
		m_matrixRows = NULL;
	};

	bool Unwrap(vector<RKVertex*> &rVertices, vector<RKFace*> &rFaces, int Pin1, int Pin2, bool LSCM);

	void initVariables(int NumberOfVariables);
	void lockSetVariable(int Index, double Value);
	void setUpMatrix();

	void addRowElement(int Index, double Value);
	void endRow();

//	bool SolveMatrixTaucs();
	bool SolveMatrixMKL();
	double GetVariable(int Index);


	static void InitSolver();
	static bool m_solverSetUp;
	static void *pt[64];					// MKL solver memory management..  set up once!!
	static int iparm[64];

private:

	int m_numberOfVariables;
	LSCM_Variable* m_variables;					// where the results go

	int m_matrixSize;							// rows == columns == square matrix
	MatrixRow* m_matrixRows;

	double* MatrixBElements;
	MatrixRow m_newLine;

	double* m_MatrixXEntries;

	void MatrixAdd(int RowIndex, int ColIndex, double Value);
};


#endif
