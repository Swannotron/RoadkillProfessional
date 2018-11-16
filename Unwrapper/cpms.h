#ifndef CPMSCALING
#define CPMSCALING

#include <vector>
#include <set>
using namespace std;

//#include <SuperLu/ssp_defs.h>
//#include <SuperLu/util.h>
/*
extern "C"
{
	#include "taucs.h"

	extern void taucs_ccs_free    (taucs_ccs_matrix*);
	extern int taucs_linsolve(taucs_ccs_matrix* A, void** F, int nrhs, void* X, void* B, char* options[], void* opt_arg[]);
	extern taucs_ccs_matrix* taucs_ccs_create               (int m, int n, int nnz, int flags);
}
*/

#include "Primitives.h"

#include "lscm.h"				// move the matrix stuff to a separate place





// Curvature Prescription and Metric Scaling

class CPMS
{
public:
	CPMS()
	{
//		m_MatrixAEntries = NULL;
		m_MatrixXEntries = NULL;
		m_pEdgeLengths = NULL;
		m_pAngles = NULL;
		m_pkOrig = NULL;
		m_pBoundaryVerts = NULL;
		m_pBoundaryReOrder = NULL;
		m_matrixRows = NULL;
	};


	~CPMS()
	{
		FreeMemory();
	};

	bool SolveCPMS(vector<RKVertex*> &rVertices, vector<RKFace*> &rFaces);


private:

	void EdgeLengthsAndAngles(vector<RKFace*> &rFaces);
	void FindKorig(vector<RKVertex*> &rVertices, vector<RKFace*> &rFaces);
	void FindLaplacian(vector<RKFace*> &rFaces);
	bool SolveMatrix();
	bool SolveMatrixTaucs();
	bool SolveMatrixMKL();
	void GetNewEdgeLengths(vector<RKFace*> &rFaces, vector<RKVertex*> &rVertices);


	double ScalePHI(int Index1, int Index2);
	void MatrixAdd(int RowIndex, int ColIndex, double Value);
	double EdgeLength(RKVertex* pVert1, RKVertex* pVert2);
	double CoTan(double x);
	void FreeMemory();

//	void TestTaucs();


	double* m_pEdgeLengths;
	double* m_pAngles;
	double* m_pkOrig;

	int* m_pBoundaryVerts;
	int* m_pBoundaryReOrder;

	int m_numberInterior;


//	int m_MatANumColumns;
//	int m_MatANumRows;
//	int m_MatANumEntries;
//	double* m_MatrixAEntries;

//	int* m_MatAColumns;
//	int* m_MatARows;

//	double* m_MatrixBEntries;
	double* m_MatrixXEntries;

//    SuperMatrix m_A, m_B;						// SuperLU's matrices
//	SuperMatrix L, U;
//	int* perm_r;
//	int* perm_c;


	int m_matrixSize;							// rows == columns == square matrix
	MatrixRow* m_matrixRows;

};



#endif		// CPMSCALING
