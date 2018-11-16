#include "Primitives.h"
#include "Vector4.h"

#include <math.h>


extern int gEdgeLen;
extern double gMatchAngleTolerence;
extern double gMatchEdgeTolerence;



// Link up faces connected to vertex on edges.  No connection is a split edge and vertex on the edge

void RKVertex::LinkUp()
{
	RKVertex* pThisVert1 = NULL;
	RKVertex* pThisVert2 = NULL;

	list<RKFace*>::iterator FaceIterator;
	list<RKFace*>::iterator FaceIterator2;

	this->m_OnEdge = false;

	for(FaceIterator = m_LinkedFaces.begin(); FaceIterator != m_LinkedFaces.end(); FaceIterator++)
	{
		RKFace* pFace = *FaceIterator;
		bool* pEdge1 = NULL;
		bool* pEdge2 = NULL;

		if(pFace->pVert1 == this)
		{
			pThisVert1 = pFace->pVert2;
			pThisVert2 = pFace->pVert3;
			pEdge1 = &pFace->pEdge1Cut;
			pEdge2 = &pFace->pEdge3Cut;
		}

		if(pFace->pVert2 == this)
		{
			pThisVert1 = pFace->pVert3;
			pThisVert2 = pFace->pVert1;
			pEdge1 = &pFace->pEdge2Cut;
			pEdge2 = &pFace->pEdge1Cut;
		}

		if(pFace->pVert3 == this)
		{
			pThisVert1 = pFace->pVert1;
			pThisVert2 = pFace->pVert2;
			pEdge1 = &pFace->pEdge3Cut;
			pEdge2 = &pFace->pEdge2Cut;
		}


		FaceIterator2 = FaceIterator;
		FaceIterator2++;


		for(; FaceIterator2 != m_LinkedFaces.end(); FaceIterator2++)
		{
			RKFace* pFace2 = *FaceIterator2;

			if(pFace2->pVert1 == pThisVert1 && pFace2->pVert2 == this)
			{
				pFace2->pEdge1Cut = false;
				*pEdge1 = false;
			}

			if(pFace2->pVert2 == pThisVert1 && pFace2->pVert3 == this)
			{
				pFace2->pEdge2Cut = false;
				*pEdge1 = false;
			}

			if(pFace2->pVert3 == pThisVert1 && pFace2->pVert1 == this)
			{
				pFace2->pEdge3Cut = false;
				*pEdge1 = false;
			}



			if(pFace2->pVert1 == this && pFace2->pVert2 == pThisVert2)
			{
				pFace2->pEdge1Cut = false;
				*pEdge2 = false;
			}

			if(pFace2->pVert2 == this && pFace2->pVert3 == pThisVert2)
			{
				pFace2->pEdge2Cut = false;
				*pEdge2 = false;
			}

			if(pFace2->pVert3 == this && pFace2->pVert1 == pThisVert2)
			{
				pFace2->pEdge3Cut = false;
				*pEdge2 = false;
			}
		}

		if(*pEdge1 == true || *pEdge2 == true)
		{
			this->m_OnEdge = true;
		}
	}
}


RKEdge* RKVertex::getEdge(int Vert2, float Length)
{
	for(int Index = 0; Index < m_LinkedEdges.size(); Index++)
	{
		if(m_LinkedEdges[Index]->m_Vert2 == Vert2)
		{
			return(m_LinkedEdges[Index]);
		}
	}

	RKEdge* newEdge = new RKEdge();
	newEdge->m_Vert2 = Vert2;
	newEdge->m_length = Length;
	m_LinkedEdges.push_back(newEdge);

	return(newEdge);
}



void RKVertex::getVertexNormal(double& VNX, double& VNY, double& VNZ)
{
	VNX = VNY = VNZ = 0.0f;

	int NumberOfFaces = (int)m_LinkedFaces.size();

	list<RKFace*>::iterator FaceIterator;

	for(FaceIterator = m_LinkedFaces.begin(); FaceIterator != m_LinkedFaces.end(); FaceIterator++)
	{
		VNX += (*FaceIterator)->NX;
		VNY += (*FaceIterator)->NY;
		VNZ += (*FaceIterator)->NZ;
	}

	VNX /= NumberOfFaces;
	VNY /= NumberOfFaces;
	VNZ /= NumberOfFaces;
}


void RKVertex::getFaceNormal(double& VNX, double& VNY, double& VNZ, RKVertex* pNextVertex)
{
	VNX = VNY = VNZ = 0.0;

	list<RKFace*>::iterator FaceIterator;

	for(FaceIterator = m_LinkedFaces.begin(); FaceIterator != m_LinkedFaces.end(); FaceIterator++)
	{
		if((*FaceIterator)->pVert1 == pNextVertex || (*FaceIterator)->pVert2 == pNextVertex || (*FaceIterator)->pVert3 == pNextVertex)
		{
			VNX = (*FaceIterator)->NX;
			VNY = (*FaceIterator)->NY;
			VNZ = (*FaceIterator)->NZ;
			return;
		}
	}
}




double RKFace::Normalise(double *n)
{
	double d;

	d= (n[0]*n[0]) + (n[1]*n[1]) + (n[2]*n[2]);
	/* A larger value causes normalise errors in a scaled down models with camera xtreme close */
	if(d>1.0e-35F) {
		d = sqrt(d);

		n[0] /= d;
		n[1] /= d;
		n[2] /= d;
	}
	else
	{
		n[0] = n[1] = n[2] = 0.0;
		d = 0.0;
	}
	return d;
}


double RKFace::VecAngleCos(float *pVert1, float *pVert2, float *pVert3)
{
	double d1[3], d2[3];

	d1[0] = pVert1[0] - pVert2[0];
	d1[1] = pVert1[1] - pVert2[1];
	d1[2] = pVert1[2] - pVert2[2];

	d2[0] = pVert3[0] - pVert2[0];
	d2[1] = pVert3[1] - pVert2[1];
	d2[2] = pVert3[2] - pVert2[2];

	Normalise(d1);
	Normalise(d2);

	return d1[0]*d2[0] + d1[1]*d2[1] + d1[2]*d2[2];
}



double RKFace::VecAngle(float *pVert1, float *pVert2, float *pVert3)
{
	double Dot = VecAngleCos(pVert1, pVert2, pVert3);

	if (Dot <= -1.0)
	{
		return M_PI;
	}
	else if (Dot >= 1.0)
		return 0.0;
	else
		return acos(Dot);
}


void RKFace::TriangleAngles(float *pVert1, float *pVert2, float *pVert3, double *Angle1, double *Angle2, double *Angle3)
{
	*Angle1 = VecAngle(pVert3, pVert1, pVert2);
	*Angle2 = VecAngle(pVert1, pVert2, pVert3);
	*Angle3 = (float)M_PI - *Angle2 - *Angle1;
}


void RKFace::TriangleAngles(RKVertex *pVert1, RKVertex *pVert2, RKVertex *pVert3, double *Angle1, double *Angle2, double *Angle3)
{
	float V1[3], V2[3], V3[3];
	V1[0] = pVert1->X;  V1[1] = pVert1->Y;  V1[2] = pVert1->Z;
	V2[0] = pVert2->X;  V2[1] = pVert2->Y;  V2[2] = pVert2->Z;
	V3[0] = pVert3->X;  V3[1] = pVert3->Y;  V3[2] = pVert3->Z;

	TriangleAngles(V1, V2, V3, Angle1, Angle2, Angle3);
}


void RKFace::GetUVAngles(float& Angle1, float& Angle2, float& Angle3)
{
	float pVert1[3];
	float pVert2[3];
	float pVert3[3];

	pVert1[0] = this->pVert1->U;
	pVert1[1] = this->pVert1->V;
	pVert1[2] = 0.0f;

	pVert2[0] = this->pVert2->U;
	pVert2[1] = this->pVert2->V;
	pVert2[2] = 0.0f;

	pVert3[0] = this->pVert3->U;
	pVert3[1] = this->pVert3->V;
	pVert3[2] = 0.0f;

	Angle1 = VecAngle(pVert3, pVert1, pVert2);
	Angle2 = VecAngle(pVert1, pVert2, pVert3);
	Angle3 = (float)M_PI - Angle2 - Angle1;
}



void RKFace::GetFaceAngles()
{
	float Coefficent1[3];
	float Coefficent2[3];
	float Coefficent3[3];

	Coefficent1[0] = pVert1->X;
	Coefficent1[1] = pVert1->Y;
	Coefficent1[2] = pVert1->Z;

	Coefficent2[0] = pVert2->X;
	Coefficent2[1] = pVert2->Y;
	Coefficent2[2] = pVert2->Z;

	Coefficent3[0] = pVert3->X;
	Coefficent3[1] = pVert3->Y;
	Coefficent3[2] = pVert3->Z;

	TriangleAngles(Coefficent1, Coefficent2, Coefficent3, &m_InteriorAngle1, &m_InteriorAngle2, &m_InteriorAngle3);
}




float RKFace::GetArea()
{
	Vector4 CrossResult;
//	Vector4 Vert1, Vert2, Vert3;
//	Vector4 DeltaX, DeltaY;
	float PolygonArea;

	Vector4 DeltaX(pVert2->X - pVert1->X, pVert2->Y - pVert1->Y, pVert2->Z - pVert1->Z, 1.0f);
	Vector4 DeltaY(pVert3->X - pVert1->X, pVert3->Y - pVert1->Y, pVert3->Z - pVert1->Z, 1.0f);

	float Mag = Vector4::CrossProduct(DeltaY, DeltaX);
//	float Mag = CrossResult.Magnitude();
	PolygonArea = Mag * 0.5f;

	GetFaceAngles();
	if(this->m_InteriorAngle1 == 0.0f) PolygonArea = 0.0f;
	if(this->m_InteriorAngle2 == 0.0f) PolygonArea = 0.0f;
	if(this->m_InteriorAngle3 == 0.0f) PolygonArea = 0.0f;

//	if(PolygonArea > 0.0000000000000001f)

	return(PolygonArea);
}



float RKFace::GetUVArea()
{
	float VecU1 = pVert2->U - pVert1->U;
	float VecV1 = pVert2->V - pVert1->V;
	float VecU2 = pVert3->U - pVert1->U;
	float VecV2 = pVert3->V - pVert1->V;

	float UVArea = ((VecU1 * VecV2) - (VecV1 * VecU2)) * 0.5f;
	return(UVArea);
}



// return next vertex round this triangle if it's on an edge

RKVertex* RKFace::GetNextTrailing(RKVertex* pVertIn, float* pLength)
{
	float DX, DY, DZ;


	if(pVert1 == pVertIn && pVert2->m_OnEdge == true && this->pEdge1Cut == true)
	{
		DX = pVert2->X - pVertIn->X;  DY = pVert2->Y - pVertIn->Y;  DZ = pVert2->Z - pVertIn->Z;
		float len = sqrt(DX * DX + DY * DY + DZ * DZ);
		*pLength += len;
		return(pVert2);
	}

	if(pVert2 == pVertIn && pVert3->m_OnEdge == true && this->pEdge2Cut == true)
	{
		DX = pVert3->X - pVertIn->X;  DY = pVert3->Y - pVertIn->Y;  DZ = pVert3->Z - pVertIn->Z;
		float len = sqrt(DX * DX + DY * DY + DZ * DZ);
		*pLength += len;
		return(pVert3);
	}

	if(pVert3 == pVertIn && pVert1->m_OnEdge == true && this->pEdge3Cut == true)
	{
		DX = pVert1->X - pVertIn->X;  DY = pVert1->Y - pVertIn->Y;  DZ = pVert1->Z - pVertIn->Z;
		float len = sqrt(DX * DX + DY * DY + DZ * DZ);
		*pLength += len;
		return(pVert1);
	}

	return(NULL);
}



RKVertex* RKFace::GetPrevTrailing(RKVertex* pVertIn, float* pLength)
{
	float DX, DY, DZ;


	if(pVert1 == pVertIn && pVert3->m_OnEdge == true && this->pEdge3Cut == true)
	{
		DX = pVert3->X - pVertIn->X;  DY = pVert3->Y - pVertIn->Y;  DZ = pVert3->Z - pVertIn->Z;
		float len = sqrt(DX * DX + DY * DY + DZ * DZ);
		*pLength += len;
		return(pVert3);
	}

	if(pVert2 == pVertIn && pVert1->m_OnEdge == true && this->pEdge1Cut == true)
	{
		DX = pVert1->X - pVertIn->X;  DY = pVert1->Y - pVertIn->Y;  DZ = pVert1->Z - pVertIn->Z;
		float len = sqrt(DX * DX + DY * DY + DZ * DZ);
		*pLength += len;
		return(pVert1);
	}

	if(pVert3 == pVertIn && pVert2->m_OnEdge == true && this->pEdge2Cut == true)
	{
		DX = pVert2->X - pVertIn->X;  DY = pVert2->Y - pVertIn->Y;  DZ = pVert2->Z - pVertIn->Z;
		float len = sqrt(DX * DX + DY * DY + DZ * DZ);
		*pLength += len;
		return(pVert2);
	}

	return(NULL);
}



float RKFace::GetUVLengthsConnected(RKVertex* pThisVert)
{
	float Sum = 0.0f;
	float DeltaU, DeltaV;

	if(pThisVert == this->pVert1)
	{
		DeltaU = pVert3->U - pVert1->U;
		DeltaV = pVert3->V - pVert1->V;
		Sum += sqrt(DeltaU * DeltaU + DeltaV * DeltaV);
		DeltaU = pVert2->U - pVert1->U;
		DeltaV = pVert2->V - pVert1->V;
		Sum = sqrt(DeltaU * DeltaU + DeltaV * DeltaV);
		return(Sum);
	}

	if(pThisVert == this->pVert2)
	{
		DeltaU = pVert1->U - pVert2->U;
		DeltaV = pVert1->V - pVert2->V;
		Sum += sqrt(DeltaU * DeltaU + DeltaV * DeltaV);
		DeltaU = pVert3->U - pVert2->U;
		DeltaV = pVert3->V - pVert2->V;
		Sum = sqrt(DeltaU * DeltaU + DeltaV * DeltaV);
		return(Sum);
	}

	if(pThisVert == this->pVert3)
	{
		DeltaU = pVert2->U - pVert3->U;
		DeltaV = pVert2->V - pVert3->V;
		Sum += sqrt(DeltaU * DeltaU + DeltaV * DeltaV);
		DeltaU = pVert1->U - pVert3->U;
		DeltaV = pVert1->V - pVert3->V;
		Sum = sqrt(DeltaU * DeltaU + DeltaV * DeltaV);
		return(Sum);
	}


	return(Sum);
}





bool RKFace::IntersectUVs(float U1In, float V1In, float U2In, float V2In, float& Len3D)
{
	float Cross1, Cross2;
	RKVertex* Vert1, *Vert2, *Vert3;


	if(NinetyAngle == 2)
	{
		Vert1 = pVert1;
		Vert2 = pVert2;
		Vert3 = pVert3;
	}
	
	if(NinetyAngle == 3)
	{
		Vert1 = pVert2;
		Vert2 = pVert3;
		Vert3 = pVert1;
	}

	if(NinetyAngle == 1)
	{
		Vert1 = pVert3;
		Vert2 = pVert1;
		Vert3 = pVert2;
	}


	float Vec1U = U1In - Vert1->U;
	float Vec1V = V1In - Vert1->V;
	float Vec2U = U2In - Vert1->U;
	float Vec2V = V2In - Vert1->V;
	Cross1 = (Vec1U * Vec2V) - (Vec1V * Vec2U);

	Vec1U = U1In - Vert2->U;
	Vec1V = V1In - Vert2->V;
	Vec2U = U2In - Vert2->U;
	Vec2V = V2In - Vert2->V;
	Cross2 = (Vec1U * Vec2V) - (Vec1V * Vec2U);

	if((Cross2 < 0 && Cross1 > 0) || (Cross2 > 0 && Cross1 < 0))
	{
		float VecX = Vert2->X - Vert1->X;
		float VecY = Vert2->Y - Vert1->Y;
		float VecZ = Vert2->Z - Vert1->Z;
		Len3D = sqrt(VecX * VecX + VecY * VecY + VecZ * VecZ);
		m_WhichEdge = 0;
		return(true);
	}


	Vec1U = U1In - Vert2->U;
	Vec1V = V1In - Vert2->V;
	Vec2U = U2In - Vert2->U;
	Vec2V = V2In - Vert2->V;
	Cross1 = (Vec1U * Vec2V) - (Vec1V * Vec2U);

	Vec1U = U1In - Vert3->U;
	Vec1V = V1In - Vert3->V;
	Vec2U = U2In - Vert3->U;
	Vec2V = V2In - Vert3->V;
	Cross2 = (Vec1U * Vec2V) - (Vec1V * Vec2U);

	if((Cross2 < 0 && Cross1 > 0) || (Cross2 > 0 && Cross1 < 0))
	{
		float VecX = Vert3->X - Vert2->X;
		float VecY = Vert3->Y - Vert2->Y;
		float VecZ = Vert3->Z - Vert2->Z;
		Len3D = sqrt(VecX * VecX + VecY * VecY + VecZ * VecZ);
		m_WhichEdge = 1;
		return(true);
	}

	return(false);
}



void RKFace::SetEdgeLen(float EdgeLen)
{
	if(m_WhichEdge == 0)
	{
		m_Edge1Len = EdgeLen;
		return;
	}

	m_Edge2Len = EdgeLen;
}


void RKFace::NinteyAngles()
{
	float pVert1[3];
	float pVert2[3];
	float pVert3[3];

	pVert1[0] = 0.0f;
	pVert1[1] = 0.0f;
	pVert1[2] = 0.0f;

	pVert2[0] = 0.0f;
	pVert2[1] = m_Edge2Len;
	pVert2[2] = 0.0f;

	pVert3[0] = m_Edge1Len;
	pVert3[1] = 0.0f;
	pVert3[2] = 0.0f;

	float Angle1 = VecAngle(pVert3, pVert1, pVert2);
	float Angle2 = VecAngle(pVert1, pVert2, pVert3);
	float Angle3 = (float)M_PI - Angle2 - Angle1;

	if(NinetyAngle == 2)
	{
		m_ABFAngle1 = Angle3;
		m_ABFAngle2 = Angle1;
		m_ABFAngle3 = Angle2;
	}
	
	if(NinetyAngle == 3)
	{
		m_ABFAngle1 = Angle2;
		m_ABFAngle2 = Angle3;
		m_ABFAngle3 = Angle1;
	}

	if(NinetyAngle == 1)
	{
		m_ABFAngle1 = Angle1;
		m_ABFAngle2 = Angle2;
		m_ABFAngle3 = Angle3;
	}
}






// tolerance?
bool RKFace::Compare(RKFace* pFaceIn)
{
	double Diff1, Diff2, Diff3;

	Diff1 = fabs(m_InteriorAngle1 - pFaceIn->m_InteriorAngle1);
	Diff2 = fabs(m_InteriorAngle2 - pFaceIn->m_InteriorAngle2);
	Diff3 = fabs(m_InteriorAngle3 - pFaceIn->m_InteriorAngle3);

	if(Diff1 < 0.04 && Diff2 < 0.04 && Diff3 < 0.04)
	{
		if(pEdge1Cut == pFaceIn->pEdge1Cut && pEdge2Cut == pFaceIn->pEdge2Cut && pEdge3Cut == pFaceIn->pEdge3Cut) 
		{
			return(true);
		}
	}


	Diff1 = fabs(m_InteriorAngle1 - pFaceIn->m_InteriorAngle2);
	Diff2 = fabs(m_InteriorAngle2 - pFaceIn->m_InteriorAngle3);
	Diff3 = fabs(m_InteriorAngle3 - pFaceIn->m_InteriorAngle1);

	if(Diff1 < 0.04 && Diff2 < 0.04 && Diff3 < 0.04)
	{
		if(pEdge1Cut == pFaceIn->pEdge2Cut && pEdge2Cut == pFaceIn->pEdge3Cut && pEdge3Cut == pFaceIn->pEdge1Cut) 
		{
			return(true);
		}
	}


	Diff1 = fabs(m_InteriorAngle1 - pFaceIn->m_InteriorAngle3);
	Diff2 = fabs(m_InteriorAngle2 - pFaceIn->m_InteriorAngle1);
	Diff3 = fabs(m_InteriorAngle3 - pFaceIn->m_InteriorAngle2);

	if(Diff1 < 0.04 && Diff2 < 0.04 && Diff3 < 0.04)
	{
		if(pEdge1Cut == pFaceIn->pEdge3Cut && pEdge2Cut == pFaceIn->pEdge1Cut && pEdge3Cut == pFaceIn->pEdge2Cut) 
		{
			return(true);
		}
	}

	return(false);
}



// tolerance?
bool RKPolygon::Compare(RKPolygon* pPolyIn)
{
	if(m_numberVertices != pPolyIn->m_numberVertices) return(false);		// different vertex count

	for(int RotStart = 0; RotStart < m_numberVertices; RotStart++)
	{
		double Diff;
		double MaxDiff = 0.0;
		double LenDiff;
		double MaxLenDiff = 0.0;
		int DstIndex = RotStart;

		for(int SrcIndex = 0; SrcIndex < m_numberVertices; SrcIndex++)
		{
			Diff = fabs(m_InteriorAngles[SrcIndex] - pPolyIn->m_InteriorAngles[DstIndex]);
			if(Diff > MaxDiff) MaxDiff = Diff;

			LenDiff = abs((m_pEdges[SrcIndex]->m_length / pPolyIn->m_pEdges[DstIndex]->m_length) - 1.0);
			if(LenDiff > MaxLenDiff) MaxLenDiff = LenDiff;

			pPolyIn->m_pVertices[DstIndex]->m_pMatches = m_pVertices[SrcIndex];

			if(Diff < gMatchAngleTolerence)
			{
				RKPolygon* OtherFaceA = m_pEdges[SrcIndex]->GetOtherFace(this);
				RKPolygon* OtherFaceB = pPolyIn->m_pEdges[DstIndex]->GetOtherFace(pPolyIn);
				if(OtherFaceA == NULL && OtherFaceB != NULL) MaxDiff = 1.0f;
				if(OtherFaceA != NULL && OtherFaceB == NULL) MaxDiff = 1.0f;
			}

			DstIndex++;
			if(DstIndex == m_numberVertices) DstIndex = 0;
		}

		if(gEdgeLen && MaxDiff < gMatchAngleTolerence && MaxLenDiff < gMatchEdgeTolerence)	return(true);
		if(gEdgeLen == false && MaxDiff < gMatchAngleTolerence) return(true);
	}

	return(false);
}



bool RKPolygon::CompareReverse(RKPolygon* pPolyIn)
{
	if(m_numberVertices != pPolyIn->m_numberVertices) return(false);		// different vertex count

	for(int RotStart = 0; RotStart < m_numberVertices; RotStart++)
	{
		double Diff;
		double MaxDiff = 0.0;
		double LenDiff;
		double MaxLenDiff = 0.0;

		int DstIndex = RotStart;

		for(int SrcIndex = 0; SrcIndex < m_numberVertices; SrcIndex++)
		{
			Diff = fabs(m_InteriorAngles[SrcIndex] - pPolyIn->m_InteriorAngles[DstIndex]);
			if(Diff > MaxDiff) MaxDiff = Diff;

			pPolyIn->m_pVertices[DstIndex]->m_pMatches = m_pVertices[SrcIndex];

			DstIndex--;
			if(DstIndex < 0) DstIndex = m_numberVertices-1;

			if(Diff < gMatchAngleTolerence)
			{
				RKPolygon* OtherFaceA = m_pEdges[SrcIndex]->GetOtherFace(this);
				RKPolygon* OtherFaceB = pPolyIn->m_pEdges[DstIndex]->GetOtherFace(pPolyIn);
				if(OtherFaceA == NULL && OtherFaceB != NULL) MaxDiff = 1.0f;
				if(OtherFaceA != NULL && OtherFaceB == NULL) MaxDiff = 1.0f;
			}

			LenDiff = abs((m_pEdges[SrcIndex]->m_length / pPolyIn->m_pEdges[DstIndex]->m_length) - 1.0);
			if(LenDiff > MaxLenDiff) MaxLenDiff = LenDiff;
		}

		if(gEdgeLen && MaxDiff < gMatchAngleTolerence && MaxLenDiff < gMatchEdgeTolerence)	return(true);
		if(gEdgeLen == false && MaxDiff < gMatchAngleTolerence) return(true);
	}

	return(false);
}



bool RKPolygon::Compare(RKPolygon* pPolyIn, int StartPointA, int StartPointB)
{
	if(m_numberVertices != pPolyIn->m_numberVertices) return(false);		// different vertex count

	int SrcIndex = StartPointA;
	int DstIndex = StartPointB;
	double Diff;
	double MaxDiff = 0.0;
	double LenDiff;
	double MaxLenDiff = 0.0;

	for(int VertCounter = 0; VertCounter < m_numberVertices; VertCounter++)
	{
		Diff = fabs(m_InteriorAngles[SrcIndex] - pPolyIn->m_InteriorAngles[DstIndex]);
		if(Diff > MaxDiff) MaxDiff = Diff;

		LenDiff = abs((m_pEdges[SrcIndex]->m_length / pPolyIn->m_pEdges[DstIndex]->m_length) - 1.0);
		if(LenDiff > MaxLenDiff) MaxLenDiff = LenDiff;

		pPolyIn->m_pVertices[DstIndex]->m_pMatches = m_pVertices[SrcIndex];
		
		SrcIndex++;
		DstIndex++;
		if(DstIndex == m_numberVertices) DstIndex = 0;
		if(SrcIndex == m_numberVertices) SrcIndex = 0;
	}
	if(gEdgeLen && MaxDiff < gMatchAngleTolerence && MaxLenDiff < gMatchEdgeTolerence)	return(true);
	if(gEdgeLen == false && MaxDiff < gMatchAngleTolerence) return(true);

	return(false);
}
	



bool RKPolygon::CompareReverse(RKPolygon* pPolyIn, int StartPointA, int StartPointB)
{
	if(m_numberVertices != pPolyIn->m_numberVertices) return(false);		// different vertex count

	int SrcIndex = StartPointA;
	int DstIndex = StartPointB;
	double Diff;
	double MaxDiff = 0.0;
	double LenDiff;
	double MaxLenDiff = 0.0;

	for(int VertCounter = 0; VertCounter < m_numberVertices; VertCounter++)
	{
		Diff = fabs(m_InteriorAngles[SrcIndex] - pPolyIn->m_InteriorAngles[DstIndex]);
		if(Diff > MaxDiff) MaxDiff = Diff;

		pPolyIn->m_pVertices[DstIndex]->m_pMatches = m_pVertices[SrcIndex];
		
		DstIndex--;
		if(DstIndex < 0) DstIndex = m_numberVertices-1;

		LenDiff = abs((m_pEdges[SrcIndex]->m_length / pPolyIn->m_pEdges[DstIndex]->m_length) - 1.0);
		if(LenDiff > MaxLenDiff) MaxLenDiff = LenDiff;

		SrcIndex++;
		if(SrcIndex == m_numberVertices) SrcIndex = 0;
	}
	if(gEdgeLen && MaxDiff < gMatchAngleTolerence && MaxLenDiff < gMatchEdgeTolerence)	return(true);
	if(gEdgeLen == false && MaxDiff < gMatchAngleTolerence) return(true);

	return(false);
}
