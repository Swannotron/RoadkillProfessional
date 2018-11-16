#include "lscm.h"

#include <math.h>



void* LSCMap::pt[64];
bool LSCMap::m_solverSetUp = false;
int LSCMap::iparm[64];




void MatrixRow::Sort()
{
	MatrixElement* sorted = new MatrixElement[this->size];

	int min = -1;

	for(int LoopIndex = 0; LoopIndex < this->size; LoopIndex++)
	{
		sorted[LoopIndex].m_index = 0x12345678;

		for(int InnerIndex = 0; InnerIndex < this->size; InnerIndex++)
		{
			int thisIndex = m_elements[InnerIndex].m_index;
			double thisVal = m_elements[InnerIndex].m_value;

			if(m_elements[InnerIndex].m_index > min && sorted[LoopIndex].m_index > thisIndex)
			{
				sorted[LoopIndex].m_index = thisIndex;
				sorted[LoopIndex].m_value = thisVal;
			}
		}

		min = sorted[LoopIndex].m_index;
	}

	for(int LoopIndex = 0; LoopIndex < this->size; LoopIndex++)
	{
		m_elements[LoopIndex].m_index = sorted[LoopIndex].m_index;
		m_elements[LoopIndex].m_value = sorted[LoopIndex].m_value;
	}

	delete[] sorted;
}


void MatrixRow::addElement(int Index, double Value)
{
	if(this->size == this->capacity)
	{
		this->expand();
	}

	this->m_elements[size].m_index = Index;
	this->m_elements[size].m_value = Value;
	size++;
}


void MatrixRow::expand()
{
	if(capacity != 0)
	{
		capacity *= 2;							// double the capacity
		MatrixElement* pNewArray = new MatrixElement[capacity];
		for(int Index = 0; Index < this->size; Index++)
		{
			pNewArray[Index].m_index = m_elements[Index].m_index;
			pNewArray[Index].m_value = m_elements[Index].m_value;
		}
		delete[] m_elements;
		m_elements = pNewArray;
	}
	else
	{
		capacity = 4;
		m_elements = new MatrixElement[4];
	}
}


void MatrixRow::AddRowColumn(int ColIndex, double Value)
{
	for(int Index = 0; Index < this->size; Index++)
	{
		if(m_elements[Index].m_index == ColIndex)
		{
			m_elements[Index].m_value += Value;
			return;
		}
	}


	addElement(ColIndex, Value);
}




int MatrixRow::sizeToDiag(int DiagIndex)
{
	int count = 0;

	for(int Index = 0; Index < this->size; Index++)
	{
		if(m_elements[Index].m_index > DiagIndex)
		{
			return(count);
		}

		count++;
	}

	return(count);
}








void LSCMap::initVariables(int NumberOfVariables)
{
	m_numberOfVariables = NumberOfVariables;
	m_variables = new LSCM_Variable[NumberOfVariables];

	for(int Index = 0; Index < NumberOfVariables; Index++)
	{
		m_variables[Index].m_pinned = false;
		m_variables[Index].m_value = 0.0f;
	}
}


void LSCMap::lockSetVariable(int Index, double Value)
{
	m_variables[Index].m_pinned = true;
	m_variables[Index].m_value = Value;
}





void LSCMap::setUpMatrix()
{
	int varIndex = 0;

	for(int Index = 0; Index < m_numberOfVariables; Index++)
	{
		if(m_variables[Index].m_pinned == false)
		{
			m_variables[Index].m_index = varIndex++;
		}
		else
		{
			m_variables[Index].m_index = -1;
		}
	}

	m_matrixSize = varIndex;
	m_matrixRows = new MatrixRow[varIndex];
	MatrixBElements = new double[varIndex];

	for(int Index = 0; Index < varIndex; Index++) MatrixBElements[Index] = 0.0f;

	m_newLine.reset();
}




void LSCMap::addRowElement(int Index, double Value)
{
	if(m_variables[Index].m_pinned == true)
	{
		m_newLine.al += Value * m_variables[Index].m_value;
	}
	else
	{
		int finalIndex = m_variables[Index].m_index;
		m_newLine.addElement(finalIndex, Value);
	}
}



void LSCMap::endRow()
{
	int newLineSize = m_newLine.size;

	for(int Outer =0; Outer < newLineSize; Outer++)
	{
		for(int Inner = 0; Inner < newLineSize; Inner++)
		{
			MatrixAdd(m_newLine.m_elements[Outer].m_index, m_newLine.m_elements[Inner].m_index, m_newLine.m_elements[Outer].m_value * m_newLine.m_elements[Inner].m_value);
		}
	}

	double S = m_newLine.al;

	for(int Index = 0; Index < newLineSize; Index++)
	{
		MatrixBElements[m_newLine.m_elements[Index].m_index] -= m_newLine.m_elements[Index].m_value * S;
	}


	m_newLine.reset();
}



void LSCMap::MatrixAdd(int RowIndex, int ColIndex, double Value)
{
	m_matrixRows[RowIndex].AddRowColumn(ColIndex, Value);
}



double LSCMap::GetVariable(int Index)
{
	if(m_variables[Index].m_pinned == true)
	{
		return(m_variables[Index].m_value);
	}

	int resultIndex = m_variables[Index].m_index;
	return(m_MatrixXEntries[resultIndex]);
}



void LSCMap::InitSolver()
{
	/* --------------------------------------------------------------------*/
	/* .. Setup Pardiso control parameters.*/
	/* --------------------------------------------------------------------*/
	for(int i = 0; i < 64; i++)
	{
		LSCMap::iparm[i] = 0;
	}

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
	/* .. Initialize the internal solver memory pointer. This is only */
	/* necessary for the FIRST call of the PARDISO solver. */
	/* --------------------------------------------------------------------*/
	for(int i = 0; i < 64; i++)
	{
		LSCMap::pt[i] = 0;
	}

	LSCMap::m_solverSetUp = true;
}




bool LSCMap::SolveMatrixMKL()
{
	int MatANumEntries = 0;
//	bool retVal = true;


	for(int RowIndex = 0; RowIndex < m_matrixSize; RowIndex++)
	{
		m_matrixRows[RowIndex].Sort();
		MatANumEntries += m_matrixRows[RowIndex].sizeToDiag(RowIndex);				// factor in removed rows
	}

	m_MatrixXEntries = new double[m_matrixSize];
	int* ColumnPtr = new int[m_matrixSize+1];
	int* RowIndices = new int[MatANumEntries];
	double* MatrixAEntries = new double[MatANumEntries];


	int count = 0;
	int RowCount = 0;
	for(int RowIndex = 0; RowIndex < m_matrixSize; RowIndex++)
	{
		ColumnPtr[RowCount++] = count+1;

		m_matrixRows[RowIndex].Sort();
		for(int ColIndex = 0; ColIndex < m_matrixRows[RowIndex].size; ColIndex++)
		{
			int oldIndex = m_matrixRows[RowIndex].m_elements[ColIndex].m_index;

			if(oldIndex >= RowIndex)								// upper matrix only
			{
				MatrixAEntries[count] = m_matrixRows[RowIndex].m_elements[ColIndex].m_value;
				RowIndices[count] = oldIndex+1;
				count++;
			}
		}
	}

	ColumnPtr[m_matrixSize] = count+1;


	if(LSCMap::m_solverSetUp == false) InitSolver();

	int mtype = -2;				/* Real symmetric matrix */
	int nrhs = 1;				/* Number of right hand sides. */
	double ddum;				/* Double dummy*/
	int idum;					/* Integer dummy.*/

	int maxfct = 1;				/* Maximum number of numerical factorizations. */
	int mnum = 1;				/* Which factorization to use. */
	int msglvl = 0;				/* Don't print statistical information in file */
	int error = 0;				/* Initialize error flag */


	for(int i = 0; i < 64; i++)
	{
		LSCMap::iparm[i] = 0;
	}

	LSCMap::iparm[0] = 1;								/* No solver default*/
	LSCMap::iparm[1] = 2;								/* Fill-in reordering from METIS */
	LSCMap::iparm[2] = 1;			//mkl_get_max_threads();			/* Numbers of processors, value of MKL_NUM_THREADS */
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
	PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	&m_matrixSize, MatrixAEntries, ColumnPtr, RowIndices, &idum, &nrhs,
	iparm, &msglvl, &ddum, &ddum, &error);

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
	PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	&m_matrixSize, MatrixAEntries, ColumnPtr, RowIndices, &idum, &nrhs,
	iparm, &msglvl, &ddum, &ddum, &error);
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
	iparm[7] = 2; /* Max numbers of iterative refinement steps. */
	/* Set right hand side to one.*/


	PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	&m_matrixSize, MatrixAEntries, ColumnPtr, RowIndices, &idum, &nrhs,
	iparm, &msglvl, MatrixBElements, m_MatrixXEntries, &error);
	if (error != 0)
	{
//		printf("\nERROR during solution: %d", error);
		return(false);			//exit(3);
	}
//	printf("\nSolve completed ... ");
//	printf("\nThe solution of the system is: ");

//	for (i = 0; i < n; i++)
//	{
//		printf("\n x [%d] = % f", i, x[i] );
//	}
//	printf ("\n");

	/* --------------------------------------------------------------------*/
	/* .. Termination and release of memory. */
	/* --------------------------------------------------------------------*/
	phase = -1; /* Release internal memory. */
	PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	&m_matrixSize, &ddum, ColumnPtr, RowIndices, &idum, &nrhs,
	iparm, &msglvl, &ddum, &ddum, &error);

	if(ColumnPtr != NULL) delete[] ColumnPtr;
	if(RowIndices != NULL) delete[] RowIndices;
	if(MatrixAEntries != NULL) delete[] MatrixAEntries;

//	MKL_FreeBuffers();
	return(true);
}


/*
bool LSCMap::SolveMatrixTaucs()
{
	int MatANumEntries = 0;
	bool retVal = true;
//	float minValue = 100000.0f;
//	float maxValue = -100000.0f;

	for(int RowIndex = 0; RowIndex < m_matrixSize; RowIndex++)
	{
		m_matrixRows[RowIndex].Sort();
		MatANumEntries += m_matrixRows[RowIndex].sizeToDiag(RowIndex);				// factor in removed rows
	}

	m_MatrixXEntries = new double[m_matrixSize];

	taucs_ccs_matrix* matrixA = taucs_ccs_create(m_matrixSize, m_matrixSize, MatANumEntries, TAUCS_DOUBLE | TAUCS_SYMMETRIC | TAUCS_LOWER);

	int count = 0;
	int RowCount = 0;
	for(int RowIndex = 0; RowIndex < m_matrixSize; RowIndex++)
	{
		matrixA->colptr[RowCount++] = count;

		m_matrixRows[RowIndex].Sort();
		for(int ColIndex = 0; ColIndex < m_matrixRows[RowIndex].size; ColIndex++)
		{
			int oldIndex = m_matrixRows[RowIndex].m_elements[ColIndex].m_index;

			if(oldIndex <= RowIndex)								// lower matrix only
			{
//				float val = m_matrixRows[RowIndex].m_elements[ColIndex].m_value;
//				if(val < minValue) minValue = val;
//				if(val > maxValue) maxValue = val;

				matrixA->values.d[count] = m_matrixRows[RowIndex].m_elements[ColIndex].m_value;
				matrixA->rowind[count] = oldIndex;
				count++;
			}
		}
	}
	matrixA->colptr[m_matrixSize] = count;


//	char* options[2] = { "taucs.factor.LLT=true", NULL };
//	int result = taucs_linsolve(matrixA, NULL, 1, m_MatrixXEntries, MatrixBElements, options, NULL);
//	taucs_ccs_free(matrixA);
//	return(true);


//	char* options[2] = { "taucs.factor.LLT=true", NULL };
	int* perm;
	int* invperm;

	taucs_ccs_order(matrixA, &perm, &invperm, "metis");	//genmmd");		//metis");	//genmmd");		//metis");
	taucs_ccs_matrix* matrixAA = taucs_ccs_permute_symmetrically(matrixA, perm, invperm);

	double* m_MatrixBEntries2 = new double[matrixA->n];
	taucs_vec_permute(matrixA->n, matrixA->flags, MatrixBElements, m_MatrixBEntries2, perm);

//	void* L = taucs_ccs_factor_llt_mf(matrixAA);
	void* L = taucs_ccs_factor_llt_ll(matrixAA);

	if(L)
	{
		taucs_supernodal_solve_llt(L, MatrixBElements, m_MatrixBEntries2); // direct solver
		taucs_vec_ipermute(matrixA->n, matrixA->flags, MatrixBElements, m_MatrixXEntries, perm);
	}
	else
	{
		retVal = false;
	}

	taucs_ccs_free(matrixA);
	taucs_ccs_free(matrixAA);

	delete[] m_MatrixBEntries2;
	delete[] MatrixBElements;
	MatrixBElements = NULL;

	return(retVal);
}
*/



bool LSCMap::Unwrap(vector<RKVertex*> &rVertices, vector<RKFace*> &rFaces, int Pin1, int Pin2, bool UseLSCM)
{
	initVariables((int)rVertices.size() * 2);				// the new stuff!

	if(Pin1 != -1)
	{
		lockSetVariable(2 * Pin1, -0.9235f);		//0.1f);
		lockSetVariable(2 * Pin1 + 1, -0.2203f);	//0.1f);
		lockSetVariable(2 * Pin2, 0.896f);		//0.9f);
		lockSetVariable(2 * Pin2 + 1, -0.15135f);	//0.9f);
	}
	else
	{
		vector<RKVertex*>::iterator iter = rVertices.begin();

		while(iter != rVertices.end())
		{
			RKVertex* pVert = *iter;
			if(pVert->m_Pinned == true)
			{
				int UnwrapIndex = pVert->UnwrapIndex;

				lockSetVariable(2*UnwrapIndex, pVert->U);
				lockSetVariable(2*UnwrapIndex+1, pVert->V);
			}

			iter++;
		}
	}



	/* construct matrix */

	setUpMatrix();

//	int temp = rFaces.size();

	for(int FaceIndex = 0; FaceIndex < rFaces.size(); FaceIndex++)
	{
		RKVertex *pVert1, *pVert2, *pVert3;
		pVert1 = rFaces[FaceIndex]->pVert1;
		pVert2 = rFaces[FaceIndex]->pVert2;
		pVert3 = rFaces[FaceIndex]->pVert3;
		double Angle1, Angle2, Angle3, Ratio, Cosine, Sine;
		double SinAngle1, SinAngle2, SinAngle3;				// SinMax;


		if(UseLSCM)
		{
			Angle1 = rFaces[FaceIndex]->m_InteriorAngle1;
			Angle2 = rFaces[FaceIndex]->m_InteriorAngle2;
			Angle3 = rFaces[FaceIndex]->m_InteriorAngle3;
		}
		else
		{
			Angle1 = rFaces[FaceIndex]->m_ABFAngle1;
			Angle2 = rFaces[FaceIndex]->m_ABFAngle2;
			Angle3 = rFaces[FaceIndex]->m_ABFAngle3;
		}

		double sinmax;
		SinAngle1 = sin(Angle1);
		SinAngle2 = sin(Angle2);
		SinAngle3 = sin(Angle3);
/*
		if(Angle1 < 0.34f || Angle2 < 0.34f || Angle3 < 0.34f)
		{
			int a = 0;
			a++;
			// very acute angle
		}
*/

		sinmax = MAX3(SinAngle1, SinAngle2, SinAngle3);

		// shift vertices to find most stable order
		if (SinAngle3 != sinmax)
		{
			SHIFT3(RKVertex*, pVert1, pVert2, pVert3);
			SHIFT3(double, Angle1, Angle2, Angle3);
			SHIFT3(double, SinAngle1, SinAngle2, SinAngle3);

			if (SinAngle2 == sinmax)
			{
				SHIFT3(RKVertex*, pVert1, pVert2, pVert3);
				SHIFT3(double, Angle1, Angle2, Angle3);
				SHIFT3(double, SinAngle1, SinAngle2, SinAngle3);
			}
		}


		/* angle based lscm formulation */
		Ratio = (SinAngle3 == 0.0f)? 1.0f: SinAngle2/SinAngle3;
		Cosine = cos(Angle1) * Ratio;
		Sine = SinAngle1 * Ratio;

		addRowElement(2 * pVert1->UnwrapIndex,   Cosine - 1.0f);
		addRowElement(2 * pVert1->UnwrapIndex+1, -Sine);
		addRowElement(2 * pVert2->UnwrapIndex,   -Cosine);
		addRowElement(2 * pVert2->UnwrapIndex+1, Sine);
		addRowElement(2 * pVert3->UnwrapIndex,   1.0f);
		endRow();

		addRowElement(2 * pVert1->UnwrapIndex,   Sine);
		addRowElement(2 * pVert1->UnwrapIndex+1, Cosine - 1.0f);
		addRowElement(2 * pVert2->UnwrapIndex,   -Sine);
		addRowElement(2 * pVert2->UnwrapIndex+1, -Cosine);
		addRowElement(2 * pVert3->UnwrapIndex+1, 1.0f);
		endRow();
	}


	//if(SolveMatrixTaucs())
	if(SolveMatrixMKL())
	{
		vector<RKVertex*>::iterator iter = rVertices.begin();

		while(iter != rVertices.end())
		{
			RKVertex* pVert = *iter;
			pVert->U = GetVariable(2 * pVert->UnwrapIndex);				// ->U
			pVert->V = GetVariable(2 * pVert->UnwrapIndex + 1);					// ->V
			pVert->m_NewU = GetVariable(2 * pVert->UnwrapIndex);				// ->U
			pVert->m_NewV = GetVariable(2 * pVert->UnwrapIndex + 1);					// ->V
			iter++;	
		}

		if(m_MatrixXEntries != NULL) delete[] m_MatrixXEntries;
		m_MatrixXEntries = NULL;
		return(true);
	}
	else
	{
		vector<RKVertex*>::iterator iter = rVertices.begin();

		while(iter != rVertices.end())
		{
			RKVertex* pVert = *iter;
			pVert->U = 0.0f;
			pVert->V = 0.0f;
			pVert->m_NewU = 0.0f;
			pVert->m_NewV = 0.0f;
			iter++;
		}
		if(m_MatrixXEntries != NULL) delete[] m_MatrixXEntries;
		m_MatrixXEntries = NULL;
		return(false);
	}
}
