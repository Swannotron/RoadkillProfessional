#include "cpms.h"


#include <math.h>




bool CPMS::SolveCPMS(vector<RKVertex*> &rVertices, vector<RKFace*> &rFaces)
{
	m_matrixSize = rVertices.size();
	m_pBoundaryVerts = new int[m_matrixSize];
	m_pBoundaryReOrder = new int[m_matrixSize];
	m_matrixRows = new MatrixRow[m_matrixSize];

	m_numberInterior = 0;
	int VertIndex = 0;
	vector<RKVertex*>::iterator iter = rVertices.begin();

	while(iter != rVertices.end())
	{
		RKVertex* pVert = *iter;
		if(pVert->m_OnEdge == true)
		{
			m_pBoundaryVerts[VertIndex] = 0;
			m_pBoundaryReOrder[VertIndex] = -1;

		}
		else
		{
			m_pBoundaryVerts[VertIndex] = 1;
			m_pBoundaryReOrder[VertIndex] = m_numberInterior++;
		}

		VertIndex++;
		iter++;
	}

	if(m_numberInterior == 0)
	{
		FreeMemory();
		return(false);
	}


	EdgeLengthsAndAngles(rFaces);
	FindKorig(rVertices, rFaces);
	FindLaplacian(rFaces);
	SolveMatrixMKL();			//Taucs();

	GetNewEdgeLengths(rFaces, rVertices);
	FreeMemory();

	return(true);
}




void CPMS::EdgeLengthsAndAngles(vector<RKFace*> &rFaces)
{
//	m_pEdgeLengths = new float[rFaces.size()*3];
	m_pAngles = new double[rFaces.size()*3];

	for(int TriIndex = 0; TriIndex < rFaces.size(); TriIndex++)
	{
		RKFace* pFace = rFaces[TriIndex];
		RKVertex* pVert1 = pFace->pVert1;
		RKVertex* pVert2 = pFace->pVert2;
		RKVertex* pVert3 = pFace->pVert3;

		double Length1 = EdgeLength(pVert2, pVert3);
		double Length2 = EdgeLength(pVert1, pVert3);
		double Length3 = EdgeLength(pVert1, pVert2);
		double Angle1 = (Length2 * Length2 + Length3 * Length3 - Length1 * Length1) / (2.0 * Length2 * Length3);
		double Angle2 = (Length1 * Length1 + Length3 * Length3 - Length2 * Length2) / (2.0 * Length1 * Length3);
		double Angle3 = (Length1 * Length1 + Length2 * Length2 - Length3 * Length3) / (2.0 * Length1 * Length2);

		if(pVert1->UnwrapIndex < pVert2->UnwrapIndex) pFace->pVert1->getEdge(pVert2->UnwrapIndex, Length3);
		if(pVert2->UnwrapIndex < pVert1->UnwrapIndex) pFace->pVert2->getEdge(pVert1->UnwrapIndex, Length3);
		if(pVert2->UnwrapIndex < pVert3->UnwrapIndex) pFace->pVert2->getEdge(pVert3->UnwrapIndex, Length1);
		if(pVert3->UnwrapIndex < pVert2->UnwrapIndex) pFace->pVert3->getEdge(pVert2->UnwrapIndex, Length1);
		if(pVert3->UnwrapIndex < pVert1->UnwrapIndex) pFace->pVert3->getEdge(pVert1->UnwrapIndex, Length2);
		if(pVert1->UnwrapIndex < pVert3->UnwrapIndex) pFace->pVert1->getEdge(pVert3->UnwrapIndex, Length2);


//		m_pEdgeLengths[TriIndex*3] = Length1;
//		m_pEdgeLengths[TriIndex*3+1] = Length2;
//		m_pEdgeLengths[TriIndex*3+2] = Length3;

		m_pAngles[TriIndex*3] = acos(Angle1);
		m_pAngles[TriIndex*3+1] = acos(Angle2);
		m_pAngles[TriIndex*3+2] = acos(Angle3);
	}
}


void CPMS::FindKorig(vector<RKVertex*> &rVertices, vector<RKFace*> &rFaces)
{
	int numVertices = rVertices.size();
	m_pkOrig = new double[numVertices];

	for(int Index = 0; Index < numVertices; Index++) m_pkOrig[Index] = 0.0;

	for(int TriIndex = 0; TriIndex < rFaces.size(); TriIndex++)
	{
		int VertIndex1 = rFaces[TriIndex]->pVert1->UnwrapIndex;
		int VertIndex2 = rFaces[TriIndex]->pVert2->UnwrapIndex;
		int VertIndex3 = rFaces[TriIndex]->pVert3->UnwrapIndex;
		double Angle1 = m_pAngles[TriIndex*3];
		double Angle2 = m_pAngles[TriIndex*3+1];
		double Angle3 = m_pAngles[TriIndex*3+2];

		m_pkOrig[VertIndex1] += Angle1;
		m_pkOrig[VertIndex2] += Angle2;
		m_pkOrig[VertIndex3] += Angle3;
	}

	vector<RKVertex*>::iterator iter = rVertices.begin();
	int VertIndex = 0;

	while(iter != rVertices.end())
	{
		RKVertex* pVert = *iter;
		double Value = m_pkOrig[VertIndex];
		double BoundaryValue = 0.0;

		if(pVert->m_OnEdge) BoundaryValue = 1.0;
		m_pkOrig[VertIndex] = (2.0 * M_PI) - (M_PI * BoundaryValue) - Value;
		VertIndex++;
		iter++;
	}
}



void CPMS::FindLaplacian(vector<RKFace*> &rFaces)
{
	for(int TriIndex = 0; TriIndex < rFaces.size(); TriIndex++)
	{
		int VertIndex1 = rFaces[TriIndex]->pVert1->UnwrapIndex;
		int VertIndex2 = rFaces[TriIndex]->pVert2->UnwrapIndex;
		int VertIndex3 = rFaces[TriIndex]->pVert3->UnwrapIndex;
		double Angle1 = 0.5 * CoTan(m_pAngles[TriIndex*3+2]);
		double Angle2 = 0.5 * CoTan(m_pAngles[TriIndex*3]);
		double Angle3 = 0.5 * CoTan(m_pAngles[TriIndex*3+1]);

		MatrixAdd(VertIndex1, VertIndex2, -Angle1);
		MatrixAdd(VertIndex2, VertIndex3, -Angle2);
		MatrixAdd(VertIndex3, VertIndex1, -Angle3);

		MatrixAdd(VertIndex2, VertIndex1, -Angle1);
		MatrixAdd(VertIndex3, VertIndex2, -Angle2);
		MatrixAdd(VertIndex1, VertIndex3, -Angle3);

		MatrixAdd(VertIndex1, VertIndex1, Angle1);
		MatrixAdd(VertIndex2, VertIndex2, Angle2);
		MatrixAdd(VertIndex3, VertIndex3, Angle3);

		MatrixAdd(VertIndex2, VertIndex2, Angle1);
		MatrixAdd(VertIndex3, VertIndex3, Angle2);
		MatrixAdd(VertIndex1, VertIndex1, Angle3);
	}
}



bool CPMS::SolveMatrixMKL()
{
	int MatANumEntries = 0;

	for(int RowIndex = 0; RowIndex < m_matrixSize; RowIndex++)
	{
		if(m_pBoundaryReOrder[RowIndex] != -1)							// boundary vertex
		{
			m_matrixRows[RowIndex].Sort();
			MatANumEntries += m_matrixRows[RowIndex].sizeToDiag(RowIndex);				// factor in removed rows
		}
	}

	double* MatrixBElements = new double[m_numberInterior];		//(double*)malloc(sizeof(double) * m_numberInterior);
	m_MatrixXEntries = new double[m_numberInterior];
	int* ColumnPtr = new int[m_numberInterior+1];		//m_matrixSize+1];
	int* RowIndices = new int[MatANumEntries];
	double* MatrixAEntries = new double[MatANumEntries];


//	taucs_ccs_matrix* matrixA = taucs_ccs_create(m_numberInterior, m_numberInterior, m_MatANumEntries, TAUCS_DOUBLE | TAUCS_SYMMETRIC | TAUCS_LOWER);

	int count = 0;
	int RowCount = 0;
	for(int RowIndex = 0; RowIndex < m_matrixSize; RowIndex++)
	{
		if(m_pBoundaryReOrder[RowIndex] != -1)							// boundary vertex
		{
			ColumnPtr[RowCount++] = count+1;
//			matrixA->colptr[RowCount++] = count;

			int newRowIndex = m_pBoundaryReOrder[RowIndex];
			MatrixBElements[newRowIndex] = m_pkOrig[RowIndex];

			m_matrixRows[RowIndex].Sort();
			for(int ColIndex = 0; ColIndex < m_matrixRows[RowIndex].size; ColIndex++)
			{
				int oldIndex = m_matrixRows[RowIndex].m_elements[ColIndex].m_index;
				int newIndex = m_pBoundaryReOrder[oldIndex];

				if(newIndex != -1)										// interior vertex
				{
					if(newIndex >= newRowIndex)							// lower matrix only
					{
						MatrixAEntries[count] = m_matrixRows[RowIndex].m_elements[ColIndex].m_value;
						RowIndices[count] = newIndex+1;
						//matrixA->values.d[count] = m_matrixRows[RowIndex].m_elements[ColIndex].m_value;
						//matrixA->rowind[count] = newIndex;				//m_matrixRows[RowIndex].m_elements[ColIndex].m_index;
						count++;
					}
				}
			}
		}
	}
	ColumnPtr[m_numberInterior] = count+1;
//	matrixA->colptr[m_numberInterior] = count;



	if(LSCMap::m_solverSetUp == false) LSCMap::InitSolver();

	int mtype = -2;				/* Real symmetric matrix */
	int nrhs = 1;				/* Number of right hand sides. */
	double ddum;				/* Double dummy*/
	int idum;					/* Integer dummy.*/

	int maxfct = 1;				/* Maximum number of numerical factorizations. */
	int mnum = 1;				/* Which factorization to use. */
	int msglvl = 0;				/* Don't print statistical information in file */
	int error = 0;				/* Initialize error flag */


	LSCMap::iparm[0] = 1;								/* No solver default*/
	LSCMap::iparm[1] = 2;								/* Fill-in reordering from METIS */
	LSCMap::iparm[2] = 0;		// mkl_get_max_threads();			/* Numbers of processors, value of MKL_NUM_THREADS */
	LSCMap::iparm[3] = 0;								/* No iterative-direct algorithm */
	LSCMap::iparm[4] = 0;								/* No user fill-in reducing permutation */
	LSCMap::iparm[5] = 0;								/* Write solution into x */
	LSCMap::iparm[6] = 16;								/* Default logical fortran unit number for output */
	LSCMap::iparm[7] = 2;								/* Max numbers of iterative refinement steps */
	LSCMap::iparm[8] = 0;								/* Not in use*/
	LSCMap::iparm[9] = 13;								/* Perturb the pivot elements with 1E-13 */
	LSCMap::iparm[10] = 1;								/* Use nonsymmetric permutation and scaling MPS */
	LSCMap::iparm[11] = 0;								/* Not in use*/
	LSCMap::iparm[12] = 0;								/* Not in use*/
	LSCMap::iparm[13] = 0;								/* Output: Number of perturbed pivots */
	LSCMap::iparm[14] = 0;								/* Not in use*/
	LSCMap::iparm[15] = 0;								/* Not in use*/
	LSCMap::iparm[16] = 0;								/* Not in use*/
	LSCMap::iparm[17] = -1;								/* Output: Number of nonzeros in the factor LU */
	LSCMap::iparm[18] = -1;								/* Output: Mflops for LU factorization */
	LSCMap::iparm[19] = 0;								/* Output: Numbers of CG Iterations */


/* --------------------------------------------------------------------*/
/* .. Reordering and Symbolic Factorization. This step also allocates */
/* all memory that is necessary for the factorization. */
/* --------------------------------------------------------------------*/
	int phase = 11;
	PARDISO (LSCMap::pt, &maxfct, &mnum, &mtype, &phase,
	&m_numberInterior, MatrixAEntries, ColumnPtr, RowIndices, &idum, &nrhs,
	LSCMap::iparm, &msglvl, &ddum, &ddum, &error);

	if (error != 0)
	{
		//printf("\nERROR during symbolic factorization: %d", error);
		return(false);			//exit(1);
	}

//	printf("\nReordering completed ... ");
//	printf("\nNumber of nonzeros in factors = %d", iparm[17]);
//	printf("\nNumber of factorization MFLOPS = %d", iparm[18]);


/* --------------------------------------------------------------------*/
/* .. Numerical factorization.*/
/* --------------------------------------------------------------------*/
	phase = 22;
	PARDISO (LSCMap::pt, &maxfct, &mnum, &mtype, &phase,
	&m_numberInterior, MatrixAEntries, ColumnPtr, RowIndices, &idum, &nrhs,
	LSCMap::iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0)
	{
		//printf("\nERROR during numerical factorization: %d", error);
		return(false);				//exit(2);
	}
//	printf("\nFactorization completed ... ");


/* --------------------------------------------------------------------*/
/* .. Back substitution and iterative refinement. */
/* --------------------------------------------------------------------*/
	phase = 33;
	LSCMap::iparm[7] = 2; /* Max numbers of iterative refinement steps. */
	/* Set right hand side to one.*/


	PARDISO (LSCMap::pt, &maxfct, &mnum, &mtype, &phase,
	&m_numberInterior, MatrixAEntries, ColumnPtr, RowIndices, &idum, &nrhs,
	LSCMap::iparm, &msglvl, MatrixBElements, m_MatrixXEntries, &error);
	if (error != 0)
	{
//		printf("\nERROR during solution: %d", error);
		return(false);			//exit(3);
	}
//	printf("\nSolve completed ... ");


	/* --------------------------------------------------------------------*/
	/* .. Termination and release of memory. */
	/* --------------------------------------------------------------------*/
	phase = -1; /* Release internal memory. */
	PARDISO (LSCMap::pt, &maxfct, &mnum, &mtype, &phase,
	&m_numberInterior, &ddum, ColumnPtr, RowIndices, &idum, &nrhs,
	LSCMap::iparm, &msglvl, &ddum, &ddum, &error);

	delete[] MatrixBElements;
	delete[] ColumnPtr;
	delete[] RowIndices;
	delete[] MatrixAEntries;

//	MKL_FreeBuffers();
	return(true);

}



bool CPMS::SolveMatrixTaucs()
{
/*	m_MatANumEntries = 0;

	for(int RowIndex = 0; RowIndex < m_matrixSize; RowIndex++)
	{
		if(m_pBoundaryReOrder[RowIndex] != -1)							// boundary vertex
		{
			m_matrixRows[RowIndex].Sort();
			m_MatANumEntries += m_matrixRows[RowIndex].sizeToDiag(RowIndex);				// factor in removed rows
		}
	}

	m_MatrixBEntries = new double[m_numberInterior];		//(double*)malloc(sizeof(double) * m_numberInterior);
	m_MatrixXEntries = new double[m_numberInterior];

	//m_MatANumEntries = 336;			// temp!!
	taucs_ccs_matrix* matrixA = taucs_ccs_create(m_numberInterior, m_numberInterior, m_MatANumEntries, TAUCS_DOUBLE | TAUCS_SYMMETRIC | TAUCS_LOWER);

	int count = 0;
	int RowCount = 0;
	for(int RowIndex = 0; RowIndex < m_matrixSize; RowIndex++)
	{
		if(m_pBoundaryReOrder[RowIndex] != -1)							// boundary vertex
		{
			matrixA->colptr[RowCount++] = count;

			int newRowIndex = m_pBoundaryReOrder[RowIndex];
			m_MatrixBEntries[newRowIndex] = m_pkOrig[RowIndex];

			m_matrixRows[RowIndex].Sort();
			for(int ColIndex = 0; ColIndex < m_matrixRows[RowIndex].size; ColIndex++)
			{
				int oldIndex = m_matrixRows[RowIndex].m_elements[ColIndex].m_index;
				int newIndex = m_pBoundaryReOrder[oldIndex];

				if(newIndex != -1)							// interior vertex
				{
					if(newIndex <= newRowIndex)				// lower matrix only
					{
						matrixA->values.d[count] = m_matrixRows[RowIndex].m_elements[ColIndex].m_value;
						matrixA->rowind[count] = newIndex;				//m_matrixRows[RowIndex].m_elements[ColIndex].m_index;
						count++;
					}
				}
			}
		}
	}
	matrixA->colptr[m_numberInterior] = count;


//	char* options[2] = { "taucs.factor.LLT=true", NULL };
//	int result = taucs_linsolve(matrixA, NULL, 1, m_MatrixXEntries, m_MatrixBEntries, options, NULL);
//	taucs_ccs_free(matrixA);


//	char* orderOptions[2] = { "genmmd", NULL };


	int* perm;
	int* invperm;

	taucs_ccs_order(matrixA, &perm, &invperm, "metis");	//genmmd");		// "metis"
	taucs_ccs_matrix* matrixAA = taucs_ccs_permute_symmetrically(matrixA, perm, invperm);

	double* m_MatrixBEntries2 = new double[m_numberInterior];
	taucs_vec_permute(m_numberInterior, matrixA->flags, m_MatrixBEntries, m_MatrixBEntries2, perm);

	void* L = taucs_ccs_factor_llt_ll(matrixAA);
//	void* L = taucs_ccs_factor_llt_mf(matrixAA);

	if(L == NULL)
	{
		taucs_ccs_free(matrixA);
		taucs_ccs_free(matrixAA);

		delete[] m_MatrixBEntries2;
		delete[] m_MatrixBEntries;
		return(false);
	}

	taucs_supernodal_solve_llt(L, m_MatrixBEntries, m_MatrixBEntries2); // direct solver
	taucs_vec_ipermute(m_numberInterior, matrixA->flags, m_MatrixBEntries, m_MatrixXEntries, perm);

	taucs_supernodal_factor_free(L);
	taucs_ccs_free(matrixA);
	taucs_ccs_free(matrixAA);

	delete[] m_MatrixBEntries2;
	delete[] m_MatrixBEntries;
	m_MatrixBEntries = NULL;
*/
	return(true);
}




void CPMS::GetNewEdgeLengths(vector<RKFace*> &rFaces, vector<RKVertex*> &rVertices)
{
	int outIndex = 0;
	for(int vertIndex = 0; vertIndex < m_matrixSize; vertIndex++)
	{
		if(m_pBoundaryReOrder[vertIndex] == -1)
		{
			m_pkOrig[vertIndex] = 0.0f;
		}
		else
		{
			m_pkOrig[vertIndex] = m_MatrixXEntries[outIndex++];
		}
	}

	for(int VertIndex = 0; VertIndex < rVertices.size(); VertIndex++)
	{
		RKVertex* pThisVert = rVertices[VertIndex];
		int VertIndex1 = pThisVert->UnwrapIndex;

		for(int EdgeIndex = 0; EdgeIndex < pThisVert->m_LinkedEdges.size(); EdgeIndex++)
		{
			RKEdge* pThisEdge = pThisVert->m_LinkedEdges[EdgeIndex];

			int VertIndex2 = pThisEdge->m_Vert2;
			double length = pThisEdge->m_length;
			double Scale = exp(-m_pkOrig[VertIndex2]) + exp(-m_pkOrig[VertIndex1]);			//ScalePHI(VertIndex1, VertIndex2);
			pThisEdge->m_length = length * Scale;
		}
	}


	bool balanced = false;

	for(int TryCount = 0; TryCount < 10 ; TryCount++)
	{
		if(balanced == true) break;

		balanced = true;

		for(int TriIndex = 0; TriIndex < rFaces.size(); TriIndex++)
		{
			RKFace* pThisface = rFaces[TriIndex];
			RKEdge *pEdge1, *pEdge2, *pEdge3;

			int VertIndex1 = pThisface->pVert1->UnwrapIndex;
			int VertIndex2 = pThisface->pVert2->UnwrapIndex;
			int VertIndex3 = pThisface->pVert3->UnwrapIndex;

			if(VertIndex1 < VertIndex2) pEdge1 = pThisface->pVert1->getEdge(VertIndex2, 0);
			else pEdge1 = pThisface->pVert2->getEdge(VertIndex1, 0);
			if(VertIndex2 < VertIndex3) pEdge2 = pThisface->pVert2->getEdge(VertIndex3, 0);
			else pEdge2 = pThisface->pVert3->getEdge(VertIndex2, 0);
			if(VertIndex3 < VertIndex1) pEdge3 = pThisface->pVert3->getEdge(VertIndex1, 0);
			else pEdge3 = pThisface->pVert1->getEdge(VertIndex3, 0);

			double EL1 = pEdge1->m_length;
			double EL2 = pEdge2->m_length;
			double EL3 = pEdge3->m_length;

			if(EL1 > EL2 + EL3)
			{
				double frac = EL1 / (EL2+EL3);
				frac *= 1.0025;
				EL2 *= frac;
				EL3 *= frac;
				pEdge2->m_length = EL2;
				pEdge3->m_length = EL3;
				balanced = false;
			}

			if(EL2 > EL1 + EL3)
			{
				double frac = EL2 / (EL1+EL3);
				frac *= 1.0025;
				EL1 *= frac;
				EL3 *= frac;
				pEdge1->m_length = EL1;
				pEdge3->m_length = EL3;
				balanced = false;
			}

			if(EL3 > EL1 + EL2)
			{
				double frac = EL3 / (EL1+EL2);
				frac *= 1.0025;
				EL1 *= frac;
				EL2 *= frac;
				pEdge1->m_length = EL1;
				pEdge2->m_length = EL2;
				balanced = false;
			}
		}
	}


	for(int TriIndex = 0; TriIndex < rFaces.size(); TriIndex++)
	{
		RKFace* pThisface = rFaces[TriIndex];

		int VertIndex1 = pThisface->pVert1->UnwrapIndex;
		int VertIndex2 = pThisface->pVert2->UnwrapIndex;
		int VertIndex3 = pThisface->pVert3->UnwrapIndex;
/*
		float Scale1 = ScalePHI(VertIndex2, VertIndex3);
		float Scale2 = ScalePHI(VertIndex3, VertIndex1);
		float Scale3 = ScalePHI(VertIndex1, VertIndex2);

		float EL1 = m_pEdgeLengths[TriIndex*3] * Scale1;
		float EL2 = m_pEdgeLengths[TriIndex*3+1] * Scale2;
		float EL3 = m_pEdgeLengths[TriIndex*3+2] * Scale3;
*/
		RKEdge *pEdge1, *pEdge2, *pEdge3;
		if(VertIndex1 < VertIndex2) pEdge1 = pThisface->pVert1->getEdge(VertIndex2, 0);
		else pEdge1 = pThisface->pVert2->getEdge(VertIndex1, 0);
		if(VertIndex2 < VertIndex3) pEdge2 = pThisface->pVert2->getEdge(VertIndex3, 0);
		else pEdge2 = pThisface->pVert3->getEdge(VertIndex2, 0);
		if(VertIndex3 < VertIndex1) pEdge3 = pThisface->pVert3->getEdge(VertIndex1, 0);
		else pEdge3 = pThisface->pVert1->getEdge(VertIndex3, 0);


		double EL1 = pEdge2->m_length;
		double EL2 = pEdge3->m_length;
		double EL3 = pEdge1->m_length;

/*
		if(EL1 < 0.0f || EL2 < 0.0f || EL3 < 0.0f)
		{
			int a = 0;
			a++;
		}
*/
		// Use law of cosines to get interior triangle angles
//		float Frac = ((-(EL3 * EL3) + (EL1 * EL1) + (EL2 * EL2) ) / (EL1 * EL2));
		double Frac1 = ( (EL2 * EL2) + (EL3 * EL3) - (EL1 * EL1) ) / (2*(EL2*EL3));
		double Frac2 = ( (EL1 * EL1) + (EL2 * EL2) - (EL3 * EL3) ) / (2*(EL1*EL2));

		double Angle1, Angle2, Angle3;

		if(Frac1 < 1.0 && Frac2 < 1.0)
		{
			Angle1 = acos(Frac1);
			Angle2 = acos(Frac2);
			//float Angle2 = asin(Frac);
			Angle3 = M_PI - (Angle1 + Angle2);
		}
		else
		{
			Angle1 = pThisface->m_InteriorAngle1;
			Angle2 = pThisface->m_InteriorAngle3;
			Angle3 = pThisface->m_InteriorAngle2;
		}

//		RKFace* test = rFaces[TriIndex];

		pThisface->m_ABFAngle1 = Angle1;
		pThisface->m_ABFAngle2 = Angle3;
		pThisface->m_ABFAngle3 = Angle2;
	}
}



double CPMS::ScalePHI(int Index1, int Index2)
{
//	float test1 = m_pkOrig[Index1];
//	float test2 = exp(test1);

	double S = exp(-m_pkOrig[Index2]) + exp(-m_pkOrig[Index1]);
	return(S);
}



void CPMS::MatrixAdd(int RowIndex, int ColIndex, double Value)
{
	m_matrixRows[RowIndex].AddRowColumn(ColIndex, Value);
}


double CPMS::CoTan(double x)
{
	double y = tan(x);
    if (y == 0)
    {
        return (double)HUGE_VAL;
    }
    return 1.0 / y;
}


double CPMS::EdgeLength(RKVertex* pVert1, RKVertex* pVert2)
{
	double dX = pVert2->X - pVert1->X;
	double dY = pVert2->Y - pVert1->Y;
	double dZ = pVert2->Z - pVert1->Z;
	double Length = sqrt((dX * dX) + (dY * dY) + (dZ * dZ));
	return(Length);
}



void CPMS::FreeMemory()
{
	if(m_pBoundaryVerts != NULL) delete[] m_pBoundaryVerts;
	if(m_pBoundaryReOrder != NULL) delete[] m_pBoundaryReOrder;
	if(m_matrixRows != NULL) delete[] m_matrixRows;
//	if(m_pEdgeLengths != NULL) delete[] m_pEdgeLengths;
	if(m_pAngles != NULL) delete[] m_pAngles;
	if(m_pkOrig != NULL) delete[] m_pkOrig;
	if(m_MatrixXEntries != NULL) delete[] m_MatrixXEntries;

	m_pBoundaryVerts = NULL;
	m_pBoundaryReOrder = NULL;
	m_matrixRows = NULL;
//	m_pEdgeLengths = NULL;
	m_MatrixXEntries = NULL;
	m_pAngles = NULL;
	m_pkOrig = NULL;
}
