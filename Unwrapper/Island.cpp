//#include "RKMesh.h"
//#include "Face.h"


#include <maya/MGlobal.h>
#include <maya/MStringArray.h>
#include <maya/MObject.h>
#include <maya/MDagPath.h>

#include <maya/MFloatPointArray.h>
#include <maya/MPointArray.h>
#include <maya/MFnMesh.h>

#include <maya/MItSelectionList.h>
#include <maya/MItMeshPolygon.h>

#include "Vector4.h"
#include "Island.h"
#include "FillHole.h"
#include "Unwrap.h"

#include <math.h>
#include <assert.h>
#include <algorithm>
#include <list>


#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/Geometry> 


extern int gScale;
extern double gMatchAngleTolerence;

extern int gMatch;
extern int gSymMatch;



using namespace Eigen;




RKEdge* EdgeHash::AddEdge(int VertIndex1, int VertIndex2, float EdgeLen, RKPolygon* pPoly, vector<RKEdge*> &m_Edges)
{
	int Size = EdgeList.size();

	for(int Index = 0; Index < this->EdgeList.size(); Index++)
	{
		if(EdgeList[Index]->m_Vert2 == VertIndex2)
		{
			EdgeList[Index]->m_pFaceList.push_back(pPoly);
			return(EdgeList[Index]);
		}
	}

	RKEdge* pNewEdge = new RKEdge();
	pNewEdge->m_Vert1 = VertIndex1;
	pNewEdge->m_Vert2 = VertIndex2;
	pNewEdge->m_length = EdgeLen;
	pNewEdge->m_pFaceList.push_back(pPoly);
	m_Edges.push_back(pNewEdge);
	EdgeList.push_back(pNewEdge);
	return(pNewEdge);
}



bool Island::Straighten()
{
	if(PinFinder.m_Pin1 == -1) PinFinder.PinMeshExtremes(m_Vertices);

	SquareGrid();
	m_LSCM.Unwrap(m_Vertices, m_Faces, PinFinder.m_Pin1, PinFinder.m_Pin2, false);
	if(CheckStraight() == true)
	{
		AverageEdges();
		m_LSCM.Unwrap(m_Vertices, m_Faces, PinFinder.m_Pin1, PinFinder.m_Pin2, false);
	}

	return(true);
}



bool Island::Unwrap(bool UseCPMS)
{
	FillHoles();
	if(PinFinder.m_Pin1 == -1) PinFinder.PinMeshExtremes(m_Vertices);


	bool UseLSCM = true;


//	if(UseCPMS) UseLSCM = m_LABF.SolveLABF(m_Vertices, m_Faces);

	if(UseCPMS) UseLSCM = m_CPMS.SolveCPMS(m_Vertices, m_Faces);

	if(UseLSCM)
	{
		m_LSCM.Unwrap(m_Vertices, m_Faces, PinFinder.m_Pin1, PinFinder.m_Pin2, !UseCPMS);
	}
	else
	{
		m_LSCM.Unwrap(m_Vertices, m_Faces, PinFinder.m_Pin1, PinFinder.m_Pin2, true);
	}

	return(true);
}



bool Island::Minimise(int IslandNumber, int IslandCount)
{
	FillHoles();

	m_MINI.Minimise(m_Vertices);

	return(true);
}


void Island::InitEdgeHash(int MaxVertCount)
{
	if(m_pEdgeHash != NULL) delete m_pEdgeHash;
	m_pEdgeHash = new EdgeHash[MaxVertCount];
}


RKEdge* Island::AddEdge(int Vert1, int Vert2, float EdgeLen, RKPolygon* pFace)
{
	if(Vert2 < Vert1)
	{
		int Temp = Vert1;
		Vert1 = Vert2;
		Vert2 = Temp;
	}

	RKEdge* pRetEdge = m_pEdgeHash[Vert1].AddEdge(Vert1, Vert2, EdgeLen, pFace, m_Edges);
	return(pRetEdge);
}



void Island::AddPolygon(int PolyIndex)
{
	RKPolygon* pPoly = Unwrap::m_Polygons[PolyIndex];

	if(pPoly->m_visited == true) return;

	m_Polygons.push_back(pPoly);
	pPoly->m_visited = true;

	for(int Index = 0; Index < pPoly->m_numberVertices; Index++)
	{
		int prevIndex = Index-1;
		int nextIndex = Index+1;
		if(prevIndex < 0) prevIndex = pPoly->m_numberVertices-1;
		if(nextIndex == pPoly->m_numberVertices) nextIndex = 0;

		int VertIndex1 = pPoly->m_pVertices[Index]->m_UVIndex;			//m_OriginalVertIndex;
		int VertIndex2 = pPoly->m_pVertices[nextIndex]->m_UVIndex;		//m_OriginalVertIndex;

		pPoly->m_pVertices[Index]->m_LinkedPolygons.push_back(pPoly);

		RKVertex* pVert1 = pPoly->m_pVertices[prevIndex];
		RKVertex* pVert2 = pPoly->m_pVertices[Index];
		RKVertex* pVert3 = pPoly->m_pVertices[nextIndex];

		double Angle1, Angle2, Angle3;
		RKFace::TriangleAngles(pVert1, pVert2, pVert3, &Angle1, &Angle2, &Angle3);

		float DeltaX = pVert2->X - pVert3->X;
		float DeltaY = pVert2->Y - pVert3->Y;
		float DeltaZ = pVert2->Z - pVert3->Z;
		float EdgeLen = sqrt(DeltaX * DeltaX + DeltaY * DeltaY + DeltaZ * DeltaZ);

		RKEdge* pEdge = AddEdge(VertIndex1, VertIndex2, EdgeLen, pPoly);
		pPoly->m_pEdges.push_back(pEdge);
		pPoly->m_InteriorAngles.push_back(Angle2);
	}
}




void Island::AddTriangle(RKFace* pStartFace)
{
	m_Faces.push_back(pStartFace);
	pStartFace->m_Visited = true;

	if(pStartFace->pVert1->m_Added == false)
	{
		m_Vertices.push_back(pStartFace->pVert1);
		pStartFace->pVert1->m_Added = true;
	}

	if(pStartFace->pVert2->m_Added == false)
	{
		m_Vertices.push_back(pStartFace->pVert2);
		pStartFace->pVert2->m_Added = true;
	}

	if(pStartFace->pVert3->m_Added == false)
	{
		m_Vertices.push_back(pStartFace->pVert3);
		pStartFace->pVert3->m_Added = true;
	}

	AddPolygon(pStartFace->m_pOrigPolygonindex);

//	m_Vertices.insert(pStartFace->pVert1);
//	m_Vertices.insert(pStartFace->pVert2);
///	m_Vertices.insert(pStartFace->pVert3);
}



void Island::AddTriangles(RKFace* pStartFace)
{
	AddTriangle(pStartFace);

	list<RKFace*> pLinked1;
	list<RKFace*> pLinked2;
	list<RKFace*> pLinked3;
	list<RKFace*>::iterator FaceIterator;
	RKFace* pFace;

	int WorkIndex = 0;

	while(WorkIndex != m_Faces.size())
	{
		pFace = m_Faces[WorkIndex];

		if(pFace->pVert1->m_Visited == false)
		{
			pFace->pVert1->m_Visited = true;
			pLinked1 = pFace->pVert1->m_LinkedFaces;
			for(FaceIterator = pLinked1.begin(); FaceIterator != pLinked1.end(); FaceIterator++)
			{
				RKFace* pFace2 = *FaceIterator;
				if(pFace2->m_Visited == false)
				{
					AddTriangle(pFace2);
				}
			}
		}

		if(pFace->pVert2->m_Visited == false)
		{
			pFace->pVert2->m_Visited = true;
			pLinked2 = pFace->pVert2->m_LinkedFaces;
			for(FaceIterator = pLinked2.begin(); FaceIterator != pLinked2.end(); FaceIterator++)
			{
				RKFace* pFace2 = *FaceIterator;
				if(pFace2->m_Visited == false)
				{
					AddTriangle(pFace2);
				}
			}
		}

		if(pFace->pVert3->m_Visited == false)
		{
			pFace->pVert3->m_Visited = true;
			pLinked3 = pFace->pVert3->m_LinkedFaces;
			for(FaceIterator = pLinked3.begin(); FaceIterator != pLinked3.end(); FaceIterator++)
			{
				RKFace* pFace2 = *FaceIterator;
				if(pFace2->m_Visited == false)
				{
					AddTriangle(pFace2);
				}
			}
		}

		WorkIndex++;
	}
}




void Island::SetUpMesh()
{
	m_Vertices.clear();
	m_Faces.clear();

}


void Island::FinaliseMesh()
{
	vector<RKVertex*>::iterator iter = m_Vertices.begin();
	int VertexIndex = 0;

	while(iter != m_Vertices.end())
	{
		(*iter)->UnwrapIndex = VertexIndex++;
		(*iter)->m_Visited = false;
		iter++;
	}


	PinFinder.m_Pin1 = -1;				// just for now
	PinFinder.m_Pin2 = -1;
	PinFinder.FindBoundaries(m_Vertices);
}



void Island::FillHoles()
{
	if(PinFinder.m_Boundaries.size() == 0) return;				// no boundary
	if(PinFinder.m_Boundaries.size() == 1)
	{
		PinFinder.GetPins();
		return;				// no holes
	}

	for(int HoleIndex = 1; HoleIndex < PinFinder.m_Boundaries.size(); HoleIndex++)
	{
		iFillHole(PinFinder.m_Boundaries[HoleIndex]->pBoundary);
	}

	PinFinder.GetPins();
}




extern int gIndices[1000];				// temp!!
extern double* gProjVertices;			// temp

void Island::iFillHole(vector<RKVertex*> boundaryVerts)
{
	FillHole FillIt;
	FillIt.HoleFiller(boundaryVerts, m_Faces);

//	return;


//	MStringArray result;
//	MString MelString;
//	char number[40];

	int Type = gIndices[0];
	int Index = 1;
	int VertIndex1 = 0;
	int VertIndex2 = 0;
	int VertIndex3 = 0;
	VertIndex1 = gIndices[Index];				// in case it's a fan


//	MelString += "polyCreateFacet";
/*
	int VertIndex = 0;
	vector<RKVertex*>::iterator VertIterator;
	for(VertIterator = boundaryVerts.begin(); VertIterator != boundaryVerts.end(); VertIterator++)
	{
		double* vert1 = &gProjVertices[VertIndex*2];
//		RKVertex* pVertex = *VertIterator;
		MelString += " -p ";

		sprintf_s(number, "%f", vert1[0]);
		MelString += number;
		MelString += " ";

		sprintf_s(number, "%f", vert1[1]);
		MelString += number;
		MelString += " ";

		sprintf_s(number, "%f", 0.0f);
		MelString += number;

		VertIndex++;
	}

	MelString += " -n RKHole";
	MGlobal::executeCommand(MelString, result, false, false);

	return;
*/

	while(Type != -4)
	{
		if(Type == -1)
		{
			VertIndex2 = gIndices[Index+1];
			VertIndex3 = gIndices[Index+2];
			Index++;
		}

		if(Type == -2)
		{
			VertIndex1 = gIndices[Index];		// flip on the odd/even?
			VertIndex2 = gIndices[Index+1];
			VertIndex3 = gIndices[Index+2];
			Index++;
		}

		if(Type == -3)
		{
			VertIndex1 = gIndices[Index];
			VertIndex2 = gIndices[Index+1];
			VertIndex3 = gIndices[Index+2];
			Index += 3;
		}


// draw the 2D triangles
/*
		double* vert1 = &gProjVertices[VertIndex1*2];
		double* vert2 = &gProjVertices[VertIndex2*2];
		double* vert3 = &gProjVertices[VertIndex3*2];

		MelString = "polyCreateFacet";
		MelString += " -p ";
		sprintf_s(number, "%f", vert1[0]);  MelString += number;  MelString += " ";
		sprintf_s(number, "%f", vert1[1]);  MelString += number;  MelString += " ";
		sprintf_s(number, "%f", 0.0f);  MelString += number;  MelString += " -p ";

		sprintf_s(number, "%f", vert2[0]);  MelString += number;  MelString += " ";
		sprintf_s(number, "%f", vert2[1]);  MelString += number;  MelString += " ";
		sprintf_s(number, "%f", 0.0f);  MelString += number;  MelString += " -p ";

		sprintf_s(number, "%f", vert3[0]);  MelString += number;  MelString += " ";
		sprintf_s(number, "%f", vert3[1]);  MelString += number;  MelString += " ";
		sprintf_s(number, "%f", 0.0f);  MelString += number;
*/

		RKVertex* pVertex1 = boundaryVerts[VertIndex1];
		RKVertex* pVertex2 = boundaryVerts[VertIndex3];
		RKVertex* pVertex3 = boundaryVerts[VertIndex2];


		RKFace* newFace = new RKFace();

		newFace->pVert1 = pVertex1;			//boundaryVerts[VertIndex1];		//origVertIndex1];
		newFace->pVert2 = pVertex2;			//boundaryVerts[VertIndex2];		//origVertIndex2];
		newFace->pVert3 = pVertex3;			//boundaryVerts[VertIndex3];		//origVertIndex3];
/*		newFace->pEdge1Cut = false;
		newFace->pEdge2Cut = false;
		newFace->pEdge3Cut = false;
		newFace->pVert1->m_OnEdge = false;
		newFace->pVert2->m_OnEdge = false;
		newFace->pVert3->m_OnEdge = false;
*/
		newFace->m_Filler = true;


		// I have a concern that some triangles may be rejected for being too small / long..  but the edgecut and OnEdge flags have been turned off
		// ToDo:: link these triangles in properly!!

		float area = newFace->GetArea();
		bool smallAngle = false;
		if(newFace->m_InteriorAngle1 < 0.05f) smallAngle = true;
		if(newFace->m_InteriorAngle2 < 0.05f) smallAngle = true;
		if(newFace->m_InteriorAngle3 < 0.05f) smallAngle = true;

		if(area == 0.0f || smallAngle == true)
		{
//			newFace->pVert1->m_OnEdge = true;
//			newFace->pVert2->m_OnEdge = true;
//			newFace->pVert3->m_OnEdge = true;
			delete newFace;
		}
		else
		{
			m_Faces.push_back(newFace);
			newFace->pVert1->m_LinkedFaces.push_back(newFace);
			newFace->pVert2->m_LinkedFaces.push_back(newFace);
			newFace->pVert3->m_LinkedFaces.push_back(newFace);
		}

///  Draw the 3D triangles
/*
		MelString = "polyCreateFacet";
		MelString += " -p ";
		sprintf_s(number, "%f", pVertex1->X);  MelString += number;  MelString += " ";
		sprintf_s(number, "%f", pVertex1->Y);  MelString += number;  MelString += " ";
		sprintf_s(number, "%f", pVertex1->Z);  MelString += number;  MelString += " -p ";

		sprintf_s(number, "%f", pVertex2->X);  MelString += number;  MelString += " ";
		sprintf_s(number, "%f", pVertex2->Y);  MelString += number;  MelString += " ";
		sprintf_s(number, "%f", pVertex2->Z);  MelString += number;  MelString += " -p ";

		sprintf_s(number, "%f", pVertex3->X);  MelString += number;  MelString += " ";
		sprintf_s(number, "%f", pVertex3->Y);  MelString += number;  MelString += " ";
		sprintf_s(number, "%f", pVertex3->Z);  MelString += number;

		MelString += " -n RKHole";
		MGlobal::executeCommand(MelString, result, false, false);
*/


		if(gIndices[Index] < 0)
		{
			Type = gIndices[Index];
			Index++;
		}

		if(gIndices[Index+2] < 0)
		{
			Type = gIndices[Index+2];
			Index+=3;
			VertIndex1 = gIndices[Index];				// in case it's a fan
		}
	}

	vector<RKVertex*>::iterator VertIterator;
	for(VertIterator = boundaryVerts.begin(); VertIterator != boundaryVerts.end(); VertIterator++)
	{
		RKVertex* pVertex = *VertIterator;
		pVertex->LinkUp();
	}

	return;

/*
	MStringArray result;
	MString MelString;
	char number[40];

	MelString += "polyCreateFacet";

	vector<RKVertex*>::iterator VertIterator;
	for(VertIterator = boundaryVerts.begin(); VertIterator != boundaryVerts.end(); VertIterator++)
	{
		RKVertex* pVertex = *VertIterator;
		MelString += " -p ";

		sprintf_s(number, "%f", pVertex->X);
		MelString += number;
		MelString += " ";

		sprintf_s(number, "%f", pVertex->Y);
		MelString += number;
		MelString += " ";

		sprintf_s(number, "%f", pVertex->Z);
		MelString += number;
	}

	MelString += " -n RKHole";
	MGlobal::executeCommand(MelString, result, false, false);



	MSelectionList selList;
	MGlobal::getActiveSelectionList( selList );
	MItSelectionList selListIter( selList );


	for( ; !selListIter.isDone(); selListIter.next() )
	{
		MDagPath		dagPath;
		MObject			component;
		selListIter.getDagPath( dagPath, component );

		MObject			mesh = dagPath.node();
		MFnMesh			meshFn( mesh );
		MItMeshPolygon	polyIter( mesh );


		for( ; !polyIter.isDone(); polyIter.next() )
		{
			bool goodTriangles = polyIter.hasValidTriangulation();

			if(goodTriangles)
			{
				MPointArray polyVerts;
				MIntArray polyIndices;

				int triCount;
				polyIter.numTriangles(triCount);
				polyIter.getTriangles(polyVerts, polyIndices, MSpace::kObject);

				MIntArray orderedIndices;
				polyIter.getVertices(orderedIndices);


				for(unsigned int Index = 0; Index < triCount; Index++)
				{
					int VertIndex1 = polyIndices[Index*3];
					int VertIndex2 = polyIndices[(Index*3)+1];
					int VertIndex3 = polyIndices[(Index*3)+2];

					RKFace* newFace = new RKFace();

					newFace->pVert1 = boundaryVerts[VertIndex1];		//origVertIndex1];
					newFace->pVert3 = boundaryVerts[VertIndex2];		//origVertIndex2];
					newFace->pVert2 = boundaryVerts[VertIndex3];		//origVertIndex3];
					newFace->pEdge1Cut = false;
					newFace->pEdge2Cut = false;
					newFace->pEdge3Cut = false;
					newFace->pVert1->m_OnEdge = false;
					newFace->pVert2->m_OnEdge = false;
					newFace->pVert3->m_OnEdge = false;
					newFace->m_Filler = true;

					float area = newFace->GetArea();
					if(area == 0.0f)
					{
						delete newFace;
					}
					else
					{
						m_Faces.push_back(newFace);
						newFace->pVert1->m_LinkedFaces.push_back(newFace);
						newFace->pVert2->m_LinkedFaces.push_back(newFace);
						newFace->pVert3->m_LinkedFaces.push_back(newFace);
					}
				}
			}
		}
	}

	MGlobal::executeCommand("delete");
*/
}



void Island::GetAreas(float& UVArea, float& PolygonArea)
{
	PolygonArea = UVArea = 0.0f;

	for(int FaceIndex = 0; FaceIndex != m_Faces.size(); FaceIndex++)
	{
		RKFace* pFace = m_Faces[FaceIndex];

		Vector4 CrossResult;
		Vector4 Vert1, Vert2, Vert3;
		//Vector4 DeltaX, DeltaY;

		Vector4 DeltaX(pFace->pVert2->X - pFace->pVert1->X, pFace->pVert2->Y - pFace->pVert1->Y, pFace->pVert2->Z - pFace->pVert1->Z, 1.0f);
		Vector4 DeltaY(pFace->pVert3->X - pFace->pVert1->X, pFace->pVert3->Y - pFace->pVert1->Y, pFace->pVert3->Z - pFace->pVert1->Z, 1.0f);

		float Mag = Vector4::CrossProduct(DeltaY, DeltaX);
//		float Mag = CrossResult.Magnitude();
		PolygonArea += Mag * 0.5f;


		DeltaX.Set(pFace->pVert2->U - pFace->pVert1->U, pFace->pVert2->V - pFace->pVert1->V, 0.0f, 1.0f);

//		DeltaX = Vector4(pFace->pVert2->U - pFace->pVert1->U, pFace->pVert2->V - pFace->pVert1->V, 0.0f, 1.0f);

		DeltaY.Set(pFace->pVert3->U - pFace->pVert1->U, pFace->pVert3->V - pFace->pVert1->V, 0.0f, 1.0f);
		//DeltaY = Vector4(pFace->pVert3->U - pFace->pVert1->U, pFace->pVert3->V - pFace->pVert1->V, 0.0f, 1.0f);

		Mag = Vector4::CrossProduct(DeltaY, DeltaX);
//		Mag = CrossResult.Magnitude();
		UVArea += Mag * 0.5f;


	}
}



void Island::GetUVBBox(float& MinU, float& MinV, float& MaxU, float& MaxV)
{
	MinU = MinV = 10000.0f;
	MaxU = MaxV = -10000.0f;

	vector<RKVertex*>::iterator iter = m_Vertices.begin();

	while(iter != m_Vertices.end())
	{
		RKVertex* pVert = *iter;
		if(pVert->U < MinU) MinU = pVert->U;
		if(pVert->U > MaxU) MaxU = pVert->U;
		if(pVert->V < MinV) MinV = pVert->V;
		if(pVert->V > MaxV) MaxV = pVert->V;
		iter++;
	}
}




void Island::UVTranslate(float UTrans, float VTrans)
{
	vector<RKVertex*>::iterator iter = m_Vertices.begin();

	while(iter != m_Vertices.end())
	{
		RKVertex* pVert = *iter;
		pVert->U += UTrans;
		pVert->V += VTrans;
		iter++;
	}
}



void Island::UVScale(float UVScale)
{
	vector<RKVertex*>::iterator iter = m_Vertices.begin();

	while(iter != m_Vertices.end())
	{
		RKVertex* pVert = *iter;
		pVert->U *= UVScale;
		pVert->V *= UVScale;
		iter++;
	}
}



void Island::Rotate(float Angle)
{
	float Sine = sin(Angle);
	float Cosine = cos(Angle);

	vector<RKVertex*>::iterator iter = m_Vertices.begin();

	while(iter != m_Vertices.end())
	{
		RKVertex* pVert = *iter;
		float OldU1 = pVert->U;
		float OldV1 = pVert->V;

		pVert->U = Cosine*OldU1 - Sine*OldV1;
		pVert->V = Sine*OldU1 + Cosine*OldV1;

		iter++;
	}
}




void Island::DeleteProfile()
{
	if(m_pProfilePixels != NULL) delete[] m_pProfilePixels;
//	if(m_pUpperProfile != NULL) delete[] m_pUpperProfile;
//	if(m_pLowerProfile != NULL) delete[] m_pLowerProfile;
//	if(m_pScaledUpper != NULL) delete[] m_pScaledUpper;
//	if(m_pScaledLower != NULL) delete[] m_pScaledLower;

	m_pProfilePixels = NULL;
//	m_pUpperProfile = NULL;
//	m_pLowerProfile = NULL;
//	m_pScaledUpper = NULL;
//	m_pScaledLower = NULL;
}



void Island::Profiles(int textureScale, float borderScale)			//, char* pFileName)
{
	// calc texture width / height from UMax and VMax

	float MinU, MinV, MaxU, MaxV;
	GetUVBBox(MinU, MinV, MaxU, MaxV);

	m_profileWidth = (int)(MaxU * (float)textureScale);
	m_profileHeight = (int)(MaxV * (float)textureScale);

	m_profileWidth += gScale*2;
	m_profileHeight += gScale*2;							// add on a gutter

	if(m_profileWidth <= 0) m_profileWidth = 1;
	if(m_profileHeight <= 0) m_profileHeight = 1;

//	m_scaledProfileWidth = (int)((MaxU * (float)textureScale) + borderScale);
//	m_scaledProfileHeight = (int)((MaxV * (float)textureScale) + borderScale);

	m_pProfilePixels = new unsigned char[m_profileWidth * m_profileHeight];
//	m_pUpperProfile = new unsigned int[m_profileWidth];
//	m_pLowerProfile = new unsigned int[m_profileWidth];

//	m_pScaledUpper = new unsigned int[m_scaledProfileWidth];
//	m_pScaledLower = new unsigned int[m_scaledProfileWidth];


	unsigned char Filler = 0;
	if(m_profileWidth == 1) Filler = 0xff;

	for(int Index = 0; Index < m_profileWidth*m_profileHeight; Index++)
	{
		m_pProfilePixels[Index] = Filler;
	}

	PolyRender.setTexture(m_pProfilePixels, m_profileWidth, m_profileHeight);

	for(int FaceIndex = 0; FaceIndex < m_Faces.size(); FaceIndex++)
	{
		float X1 = m_Faces[FaceIndex]->pVert1->U * textureScale;
		float Y1 = m_Faces[FaceIndex]->pVert1->V * textureScale;
		float X2 = m_Faces[FaceIndex]->pVert2->U * textureScale;
		float Y2 = m_Faces[FaceIndex]->pVert2->V * textureScale;
		float X3 = m_Faces[FaceIndex]->pVert3->U * textureScale;
		float Y3 = m_Faces[FaceIndex]->pVert3->V * textureScale;

		if(!m_Faces[FaceIndex]->m_Filler)
		{
			PolyRender.RenderPolygon(X1, Y1, X2, Y2, X3, Y3);
		}
	}

	CreateBitProfiles();

/*
	for(int XIndex = 0; XIndex < m_profileWidth; XIndex++)
	{
		bool foundBottom = false;
		bool foundTop = false;
		int YIndexBottom = 0;
		int YIndexTop = m_profileHeight-1;

		m_pLowerProfile[XIndex] = -1;		// if nothing was found
		m_pUpperProfile[XIndex] = -1;		// if nothing was found

		while(YIndexBottom < m_profileHeight/2 && (foundTop == false || foundBottom == false))
		{
			if(m_pProfilePixels[YIndexBottom * m_profileWidth + XIndex] != 0)
			{
				if(foundBottom == false)
				{
					m_pLowerProfile[XIndex] = YIndexBottom;
					foundBottom = true;
				}
				if(foundTop == false) m_pUpperProfile[XIndex] = YIndexBottom;
			}

			if(m_pProfilePixels[YIndexTop * m_profileWidth + XIndex] != 0)
			{
				if(foundTop == false)
				{
					m_pUpperProfile[XIndex] = YIndexTop;
					foundTop = true;
				}
				if(foundBottom == false) m_pLowerProfile[XIndex] = YIndexTop;
			}

			YIndexTop--;
			YIndexBottom++;
		}
	}


	int StickLength = m_scaledProfileWidth - m_profileWidth;			// the difference between the two borders
	int StickPoint = -(StickLength/2);									// the current point where the stick is in the lower border
	float borderWidth = m_scaledProfileHeight - this->m_profileHeight;

	for(int Index = 0; Index < m_scaledProfileWidth; Index++)
	{
		int Min, Max;
		FindMinMax(StickLength, StickPoint, Min, Max);
		if(Min == -1 && Max == -1)
		{
			m_pScaledUpper[Index] = -1;
			m_pScaledLower[Index] = -1;
		}
		else
		{
			m_pScaledUpper[Index] = (unsigned int)((float)Max + borderWidth);
			m_pScaledLower[Index] = (unsigned int)Min;
		}
		StickPoint++;
	}


	if(m_pProfilePixels != NULL) delete[] m_pProfilePixels;
	m_pProfilePixels = NULL;
*/
}

/*
void Island::FindMinMax(int StickLength, int StickPoint, int& Min, int& Max)
{
	Min = m_profileHeight;
	Max = 0;
	int CurrentPoint = StickPoint - (StickLength/2);
	bool found = false;

	for(int Index = 0; Index < StickLength; Index++)
	{
		if(CurrentPoint >= 0 && CurrentPoint < m_profileWidth)
		{
			int UpperHeight = m_pUpperProfile[CurrentPoint];
			int LowerHeight = m_pLowerProfile[CurrentPoint];

			if(UpperHeight != -1 && LowerHeight != -1)
			{
				if(UpperHeight > Max) Max = UpperHeight;
				if(LowerHeight < Min) Min = LowerHeight;
				found = true;
			}
		}

		CurrentPoint++;
	}

	if(found == false) { Min = -1; Max = -1; };
}
*/



void Island::CreateBitProfiles()
{
	int BitsWidth = m_profileWidth;
	while(BitsWidth & 0x3f) BitsWidth++;
	m_BitsWidth = BitsWidth >> 6;
	m_BitsWidth++;						// line of zeros to rotate into

	m_pProfileBits = new uint64_t[m_BitsWidth * m_profileHeight];

	unsigned int count = (m_BitsWidth * m_profileHeight);

	for(unsigned int Index = 0; Index < count; Index++)
	{
		m_pProfileBits[Index] = 0;
	}

	for(unsigned int YIndex = 0; YIndex < m_profileHeight; YIndex++)
	{
		for(unsigned int XIndex = 0; XIndex < m_profileWidth; XIndex++)
		{
			int DestAddress = (YIndex * m_BitsWidth) + (XIndex >> 6);

			if(m_pProfilePixels[YIndex * m_profileWidth + XIndex])
			{
				int RotVal = XIndex & 0x3f;
				uint64_t Pixel = 0x8000000000000000;

				m_pProfileBits[DestAddress] |= Pixel >> RotVal;
			}
		}
	}
}

/*
void Island::FindMinMax(int StickLength, int StickPoint, int& Min, int& Max)
{
	Min = m_profileHeight;
	Max = 0;
	int CurrentPoint = StickPoint - (StickLength/2);
	bool found = false;

	for(int Index = 0; Index < StickLength; Index++)
	{
		if(CurrentPoint >= 0 && CurrentPoint < m_profileWidth)
		{
			int UpperHeight = m_pUpperProfile[CurrentPoint];
			int LowerHeight = m_pLowerProfile[CurrentPoint];

			if(UpperHeight != -1 && LowerHeight != -1)
			{
				if(UpperHeight > Max) Max = UpperHeight;
				if(LowerHeight < Min) Min = LowerHeight;
				found = true;
			}
		}

		CurrentPoint++;
	}

	if(found == false) { Min = -1; Max = -1; };
}
*/






bool ComparePoints(const RKPoint* a, const RKPoint* b)
{
	if(a->x == b->x)
	{
		return(a->y < b->y);
	}

	return(a->x < b->x);
}



void Island::MinAreaRotate()
{
	vector <RKPoint*> sortedPoints;
	int numberPoints;
	RKPoint* pPointsIn = NULL;
	RKPoint* pPointsOut = NULL;

	if(PinFinder.m_Boundaries.size() == 0)
	{
		numberPoints = this->m_Vertices.size();
		pPointsIn = new RKPoint[numberPoints];
		pPointsOut = new RKPoint[numberPoints+1];

		for(int Index = 0; Index < numberPoints; Index++)
		{
			pPointsIn[Index].x = m_Vertices[Index]->U;
			pPointsIn[Index].y = m_Vertices[Index]->V;
			sortedPoints.push_back(&pPointsIn[Index]);
		}
	}
	else
	{
		numberPoints = PinFinder.m_Boundaries[0]->pBoundary.size();
		pPointsIn = new RKPoint[numberPoints];
		pPointsOut = new RKPoint[numberPoints+1];

		for(int Index = 0; Index < numberPoints; Index++)
		{
			pPointsIn[Index].x = PinFinder.m_Boundaries[0]->pBoundary[Index]->U;
			pPointsIn[Index].y = PinFinder.m_Boundaries[0]->pBoundary[Index]->V;
			sortedPoints.push_back(&pPointsIn[Index]);
		}
	}

	// sort pPointsIn
	std::sort(sortedPoints.begin(), sortedPoints.end(), ComparePoints);
	RKPoint* pSorted = new RKPoint[numberPoints];

	for(int Index = 0; Index < numberPoints; Index++)
	{
		pSorted[Index].x = sortedPoints[Index]->x;
		pSorted[Index].y = sortedPoints[Index]->y;
	}

	int hullSize = chainHull_2D(pSorted, numberPoints, pPointsOut);

	// chainHull_2D sometimes duplicates the first / last vertices
	if(pPointsOut[hullSize-1].x == pPointsOut[0].x && pPointsOut[hullSize-1].y == pPointsOut[0].y)
	{
		hullSize--;
	}


	int sourcePoint = hullSize-1;
	int numberOutPoints = 0;

	while(sourcePoint >= 0)
	{
//		if(numberOutPoints == 0)
//		{
//			if(pPointsOut[sourcePoint].x != pPointsOut[0].x && pPointsOut[sourcePoint].y != pPointsOut[0].y)
//			{
				pSorted[numberOutPoints].x = pPointsOut[sourcePoint].x;
				pSorted[numberOutPoints].y = pPointsOut[sourcePoint].y;
				numberOutPoints++;
//			}
//		}
//		else
//		{
//			if(pPointsOut[sourcePoint].x != pSorted[numberOutPoints-1].x && pPointsOut[sourcePoint].y != pSorted[numberOutPoints-1].y)
//			{
//				pSorted[numberOutPoints].x = pPointsOut[sourcePoint].x;
//				pSorted[numberOutPoints].y = pPointsOut[sourcePoint].y;
//				numberOutPoints++;
//			}
//		}

		sourcePoint--;
	}

//	float deltaX1 = pSorted[1].x - pSorted[0].x;
//	float deltaY1 = pSorted[1].y - pSorted[0].y;
//	float deltaX2 = pSorted[2].x - pSorted[0].x;
//	float deltaY2 = pSorted[2].y - pSorted[0].y;


//	float clockwise = (deltaX1 * deltaY2) - (deltaY1 * deltaX2);


	// find furthest two points
//	float Angle = Diameter(pSorted, numberOutPoints);
	float Angle = MinimumArea(pSorted, numberOutPoints);

	this->Rotate(Angle);

	float MinU, MinV, MaxU, MaxV;
	this->GetUVBBox(MinU, MinV, MaxU, MaxV);

	float LenU = MaxU - MinU;
	float LenV = MaxV - MinV;
	if(LenU > LenV)					// longest edge vertical
	{
		this->Rotate(M_PI/2);
	}


	sortedPoints.clear();
	delete[] pSorted;
	delete[] pPointsIn;
	delete[] pPointsOut;
}



float Island::Diameter(RKPoint* Hull, int NumberPoints)
{
	float Diameter = 0.0f;

	int BestVert1 = 0;
	int BestVert2 = 0;

	int Calip1Index = 0;						// starts at minimum x
	int Calip2Index = 0;						// starts at maximum x
	float Calip1Angle = 0.0f;						// starts pointing up
	float Calip2Angle = 0.0f;						// starts pointing down
	float* Angles = new float[NumberPoints];
	RKPoint* p1, *p2, *p3;

	int Calip1Start = 0;
	int Calip2Start = 0;

	float MinX = 1000000;
	float MaxX = -1000000;
	for(int Index = 0; Index < NumberPoints; Index++)
	{
		p1 = (Index == 0)? &Hull[NumberPoints-1]: &Hull[Index-1];
		p2 = &Hull[Index];
		p3 = (Index == NumberPoints-1)? &Hull[0]: &Hull[Index+1];

		Angles[Index] = (float)M_PI - vec2_angle(p1, p2, p3);

		if(Hull[Index].x < MinX)
		{
			MinX = Hull[Index].x;
			Calip1Start = Calip1Index = Index;
		}

		if(Hull[Index].x > MaxX)
		{
			MaxX = Hull[Index].x;
			Calip2Start = Calip2Index = Index;
		}
	}

	RKPoint Start;
	Start.x = Hull[Calip1Start].x;
	Start.y = Hull[Calip1Start].y + 1.0f;
	Calip1Angle = vec2_angle(&Hull[(Calip1Start+1)%NumberPoints], &Hull[Calip1Start], &Start);

	Start.x = Hull[Calip2Start].x;
	Start.y = Hull[Calip2Start].y - 1.0f;
	Calip2Angle = vec2_angle(&Hull[(Calip2Start+1)%NumberPoints], &Hull[Calip2Start], &Start);

	float rotated = 0.0f;

	while(rotated < M_PI * 2.0f)
	{
		float dx = Hull[Calip2Index].x - Hull[Calip1Index].x;
		float dy = Hull[Calip2Index].y - Hull[Calip1Index].y;
		float len = sqrt(dx * dx + dy * dy);

		if(len > Diameter)
		{
			Diameter = len;
			BestVert1 = Calip1Index;
			BestVert2 = Calip2Index;
		}

		float Angle1 = Angles[Calip1Index];
		float Angle2 = Angles[Calip2Index];

		if(Angle1 < Angle2)
		{
			Calip1Index = (Calip1Index+1) % NumberPoints;
			Calip2Angle -= Angle1;
			Calip1Angle = Angles[Calip1Index];
			rotated += Angle1;
		}
		else
		{
			Calip2Index = (Calip2Index+1) % NumberPoints;
			Calip1Angle -= Angle2;
			Calip2Angle = Angles[Calip2Index];
			rotated += Angle2;
		}
	}


	float VVec = Hull[BestVert2].y - Hull[BestVert1].y;
	float UVec = Hull[BestVert2].x - Hull[BestVert1].x;
	float len = sqrt(UVec * UVec + VVec * VVec);

	UVec /= len;
	VVec /= len;

	float RetAngle = atan2(VVec, UVec);
//	RetAngle += M_PI * 0.5f;

	return(-RetAngle);
}


float Island::MinimumArea(RKPoint* Hull, int NumberPoints)
{
//	float Diameter = 0.0f;
//	int BestVert1, BestVert2;

	float* Angles = new float[NumberPoints];
	RKPoint* p1, *p2, *p3, *p4, *p1n;

	int CalipIndices[4];
	float CalipAngles[4];

	float MinY = 1e10;
	float MaxY = -1e10;
	float MinX = 1e10;
	float MaxX = -1e10;
	for(int Index = 0; Index < NumberPoints; Index++)
	{
		p1 = (Index == 0)? &Hull[NumberPoints-1]: &Hull[Index-1];
		p2 = &Hull[Index];
		p3 = (Index == NumberPoints-1)? &Hull[0]: &Hull[Index+1];

		Angles[Index] = (float)M_PI - vec2_angle(p1, p2, p3);

		if(Hull[Index].x < MinX)
		{
			MinX = Hull[Index].x;
			CalipIndices[0] = Index;
		}

		if(Hull[Index].x > MaxX)
		{
			MaxX = Hull[Index].x;
			CalipIndices[2] = Index;
		}

		if(Hull[Index].y < MinY)
		{
			MinY = Hull[Index].y;
			CalipIndices[3] = Index;
		}

		if(Hull[Index].y > MaxY)
		{
			MaxY = Hull[Index].y;
			CalipIndices[1] = Index;
		}
	}

	RKPoint Start;
	int pointIndex = CalipIndices[0];
	Start.x = Hull[pointIndex].x;
	Start.y = Hull[pointIndex].y + 1.0f;
	CalipAngles[0] = vec2_angle(&Hull[(pointIndex+1)%NumberPoints], &Hull[pointIndex], &Start);

	pointIndex = CalipIndices[1];
	Start.x = Hull[pointIndex].x + 1.0f;
	Start.y = Hull[pointIndex].y;
	CalipAngles[1] = vec2_angle(&Hull[(pointIndex+1)%NumberPoints], &Hull[pointIndex], &Start);

	pointIndex = CalipIndices[2];
	Start.x = Hull[pointIndex].x;
	Start.y = Hull[pointIndex].y - 1.0f;
	CalipAngles[2] = vec2_angle(&Hull[(pointIndex+1)%NumberPoints], &Hull[pointIndex], &Start);

	pointIndex = CalipIndices[3];
	Start.x = Hull[pointIndex].x - 1.0f;
	Start.y = Hull[pointIndex].y;
	CalipAngles[3] = vec2_angle(&Hull[(pointIndex+1)%NumberPoints], &Hull[pointIndex], &Start);

	float rotated = 0.0f;
	float minarea = 1e10;
	float minangle = 0.0f;

	while(rotated < M_PI/2)
	{
		int minimumIndex = 0;
		float minimumAngle = 1e10;

		for(int Index = 0; Index < 4; Index++)
		{
			if(CalipAngles[Index] < minimumAngle)
			{
				minimumAngle = CalipAngles[Index];
				minimumIndex = Index;
			}
		}

		rotated += minimumAngle;
		int nextIndex = (CalipIndices[minimumIndex]+1)%NumberPoints;

		CalipAngles[minimumIndex] = Angles[nextIndex];
		CalipAngles[(minimumIndex+1)%4] = CalipAngles[(minimumIndex+1)%4] - minimumAngle;
		CalipAngles[(minimumIndex+2)%4] = CalipAngles[(minimumIndex+2)%4] - minimumAngle;
		CalipAngles[(minimumIndex+3)%4] = CalipAngles[(minimumIndex+3)%4] - minimumAngle;

		p1 = &Hull[CalipIndices[minimumIndex]];
		p1n = &Hull[nextIndex];
		p2 = &Hull[CalipIndices[(minimumIndex+1)%4]];
		p3 = &Hull[CalipIndices[(minimumIndex+2)%4]];
		p4 = &Hull[CalipIndices[(minimumIndex+3)%4]];

		float VVec = p1n->y - p1->y;
		float UVec = p1n->x - p1->x;
		float len = sqrt(UVec * UVec + VVec * VVec);

		if (len > 0.0f)
		{
			len = 1.0f/len;
			Start.x = (p1n->x - p1->x);		// * len;
			Start.y = (p1n->y - p1->y);		// * len;

			float area = rectangleArea(p1, &Start, p2, p3, p4);

			if (area < minarea)				// smallest area
			{
				minarea = area;
				minangle = rotated;
			}
		}

		CalipIndices[minimumIndex] = nextIndex;
	}

	if (minangle > (float)M_PI/4)
		minangle -= (float)M_PI/2;

	delete[] Angles;
	return(minangle);

/*
	float VVec = Hull[BestVert2].y - Hull[BestVert1].y;
	float UVec = Hull[BestVert2].x - Hull[BestVert1].x;
	float len = sqrt(UVec * UVec + VVec * VVec);

	UVec /= len;
	VVec /= len;

	float RetAngle = atan2(VVec, UVec);
//	RetAngle += M_PI * 0.5f;
*/
//	return(-RetAngle);
}



float Island::rectangleArea(RKPoint* p1, RKPoint* dir, RKPoint* p2, RKPoint* p3, RKPoint* p4)
{
	/* given 4 points on the rectangle edges and the direction of on edge,
	   compute the area of the rectangle */

//	float orthodir[2], corner1[2], corner2[2], corner3[2];
	RKPoint orthodir, corner1, corner2, corner3;

	orthodir.x = dir->y;
	orthodir.y = -dir->x;

//	if (!intersectLine2D(p1, dir, p2, &orthodir, &corner1))
	if(!intersect(p2, p1, &orthodir, &corner1))
		return 1e10;

//	if (!intersectLine2D(p1, dir, p4, &orthodir, &corner2))
	if(!intersect(p4, p1, &orthodir, &corner2))
		return 1e10;

//	if (!intersectLine2D(p3, dir, p4, &orthodir, &corner3))
	if(!intersect(p3, &corner2, dir, &corner3))
		return 1e10;

	return VecLen(&corner1, &corner2)*VecLen(&corner2, &corner3);
}


bool Island::intersect(RKPoint* p1, RKPoint* pointOnLine, RKPoint* dir, RKPoint* isect)
{
	float XKJ = p1->x -	pointOnLine->x;		//Point1X - EndPointX;
	float YKJ = p1->y -	pointOnLine->y;		//Point1Y - EndPointY;
	float XLK = dir->x;						//Point2X - Point1X;
	float YLK = dir->y;						//Point2Y - Point1Y;
	float Denom = XLK * XLK + YLK * YLK;

	if(Denom < 0.000001)
	{
		return(false);				// line has hardly any length
	}

	float T = -(XKJ*XLK + YKJ*YLK)/Denom;
	isect->x = pointOnLine->x + XKJ + T*XLK;
	isect->y = pointOnLine->y + YKJ + T*YLK;

	return(true);
}


bool Island::intersectLine2D(RKPoint *v1, RKPoint *dir1, RKPoint *v2, RKPoint *dir2, RKPoint *isect)
{
	float lmbda, div;

	div = dir2->x * dir1->y - dir2->y * dir1->x;

	if(div == 0.0f) return(false);

	lmbda= ((v1->y - v2->y) * dir1->x - (v1->x - v2->x) * dir1->y) / div;
	isect->x = v1->x + lmbda * dir2->x;
	isect->y = v1->y + lmbda * dir2->y;

	return(true);
}



float Island::VecLen(RKPoint* v1, RKPoint* v2)
{
	float x, y;

	x = v1->x - v2->x;
	y = v1->y - v2->y;
	return (float)sqrt(x*x+y*y);
}


float Island::vec2_angle(RKPoint* v1, RKPoint* v2, RKPoint* v3)
{
	float u1[2], u2[2], u3[2];
	float d1[2], d2[2];
	float dot;

	u1[0] = v1->x; u1[1] = v1->y;
	u2[0] = v2->x; u2[1] = v2->y;
	u3[0] = v3->x; u3[1] = v3->y;

	d1[0] = u1[0] - u2[0];
	d1[1] = u1[1] - u2[1];
	d2[0] = u3[0] - u2[0];
	d2[1] = u3[1] - u2[1];

	float len = sqrt((d1[0] * d1[0]) + (d1[1] * d1[1]));
	float len2 = sqrt((d2[0] * d2[0]) + (d2[1] * d2[1]));

	d1[0] /= len;   d1[1] /= len;
	d2[0] /= len2;   d2[1] /= len2;

	dot = d1[0]*d2[0] + d1[1]*d2[1];

	if (dot <= -1.0f)
		return (float)M_PI;
	else if (dot >= 1.0f)
		return 0.0f;
	else
		return (float)acos(dot);
}



// Copyright 2001, softSurfer (www.softsurfer.com)
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from its use.
// Users of this code must verify correctness for their application.


// Assume that a class is already given for the object:
//    Point with coordinates {float x, y;}
//===================================================================


// isLeft(): tests if a point is Left|On|Right of an infinite line.
//    Input:  three points P0, P1, and P2
//    Return: >0 for P2 left of the line through P0 and P1
//            =0 for P2 on the line
//            <0 for P2 right of the line
//    See: the January 2001 Algorithm on Area of Triangles
inline float isLeft( RKPoint P0, RKPoint P1, RKPoint P2 )
{
    return (P1.x - P0.x)*(P2.y - P0.y) - (P2.x - P0.x)*(P1.y - P0.y);
}
//===================================================================


// chainHull_2D(): Andrew's monotone chain 2D convex hull algorithm
//     Input:  P[] = an array of 2D points
//                   presorted by increasing x- and y-coordinates
//             n = the number of points in P[]
//     Output: H[] = an array of the convex hull vertices (max is n)
//     Return: the number of points in H[]

int Island::chainHull_2D( RKPoint* P, int n, RKPoint* H )
{
    // the output array H[] will be used as the stack
    int    bot=0, top=(-1);											// indices for bottom and top of the stack
    int    i;														// array scan index

    // Get the indices of points with min x-coord and min|max y-coord
    int minmin = 0, minmax;
    float xmin = P[0].x;

    for (i=1; i<n; i++)
        if (P[i].x != xmin) break;

	minmax = i-1;
    if (minmax == n-1)
	{       // degenerate case: all x-coords == xmin
        H[++top] = P[minmin];
        if (P[minmax].y != P[minmin].y) // a nontrivial segment
            H[++top] = P[minmax];

		H[++top] = P[minmin];           // add polygon endpoint
        return top+1;
    }

    // Get the indices of points with max x-coord and min|max y-coord
    int maxmin, maxmax = n-1;
    float xmax = P[n-1].x;
    for (i=n-2; i>=0; i--)
        if (P[i].x != xmax) break;
    maxmin = i+1;

    // Compute the lower hull on the stack H
    H[++top] = P[minmin];      // push minmin point onto stack
    i = minmax;
    while (++i <= maxmin)
    {
        // the lower line joins P[minmin] with P[maxmin]
        if (isLeft( P[minmin], P[maxmin], P[i]) >= 0 && i < maxmin)
            continue;          // ignore P[i] above or on the lower line

        while (top > 0)        // there are at least 2 points on the stack
        {
            // test if P[i] is left of the line at the stack top
            if (isLeft( H[top-1], H[top], P[i]) > 0)
                break;         // P[i] is a new hull vertex
            else
                top--;         // pop top point off stack
        }
        H[++top] = P[i];       // push P[i] onto stack
    }

    // Next, compute the upper hull on the stack H above the bottom hull
    if (maxmax != maxmin)      // if distinct xmax points
        H[++top] = P[maxmax];  // push maxmax point onto stack
    bot = top;                 // the bottom point of the upper hull stack
    i = maxmin;
    while (--i >= minmax)
    {
        // the upper line joins P[maxmax] with P[minmax]
        if (isLeft( P[maxmax], P[minmax], P[i]) >= 0 && i > minmax)
            continue;          // ignore P[i] below or on the upper line

        while (top > bot)    // at least 2 points on the upper stack
        {
            // test if P[i] is left of the line at the stack top
            if (isLeft( H[top-1], H[top], P[i]) > 0)
                break;         // P[i] is a new hull vertex
            else
                top--;         // pop top point off stack
        }
        H[++top] = P[i];       // push P[i] onto stack
    }
    if (minmax != minmin)
        H[++top] = P[minmin];  // push joining endpoint onto stack

    return top+1;
}





bool CompareBlocks(const WalkBlock* a, const WalkBlock* b)
{
	return a->m_DegreeOfDiff < b->m_DegreeOfDiff;
}

/****************************************************
*				Island comparison					*
****************************************************/

bool Island::Compare(Island* pIsland)
{
	if(pIsland->m_Polygons.size() != m_Polygons.size()) return(false);		// not the same if face count is different

	m_TestSym = true;
	//	bool Symmetrical = TestInteriorSymmetry();
	//	if(Symmetrical)	return(true);

	if(gSymMatch)
	{
		bool SymmetricalWith = TestSymmetrical(pIsland);					// symmetrical with another island
		if(SymmetricalWith) return(true);
	}

	m_TestSym = false;

	
	if(gMatch == false) return(false);


	vector<WalkBlock*> StartFaces = FindStartFaces(pIsland);			// delete StartFaces walkblocks
	if(StartFaces.size() == 0) return(false);
	if(StartFaces.size() > 1)
	{
		std::sort(StartFaces.begin(), StartFaces.end(), CompareBlocks);
	}

	// from start face, walk out

	for(int StartIndex = 0; StartIndex < StartFaces.size(); StartIndex++)
	{
		bool worked = WalkOutMatch(pIsland, StartFaces[StartIndex]);
		if(worked == true) 
		{
			pIsland->m_pSameAs = this;
			return(true);
		}
	}

	return(false);
}



// from starting matching face, walk out through vertices

bool Island::WalkOutMatch(Island* pIsland, WalkBlock* pStartBlock)
{
	list<WalkBlock*> toDo;

	// clear all m_Visited flags please
	for(int FaceIndex = 0; FaceIndex < m_Polygons.size(); FaceIndex++)
	{
		m_Polygons[FaceIndex]->m_visited = false;
		pIsland->m_Polygons[FaceIndex]->m_visited2 = false;
	}

	toDo.push_back(pStartBlock);
	pStartBlock->m_pFace1->m_visited = true;
	pStartBlock->m_pFace2->m_visited2 = true;

	while(toDo.size() != 0)
	{
		WalkBlock* pThisBlock = toDo.front();
		toDo.pop_front();
		bool matched = GetConnectedFaces(pThisBlock, toDo);
		delete pThisBlock;

		if(matched == false)
		{
			while(toDo.size() != 0)							// empty the todo list
			{
				WalkBlock* pThisBlock = toDo.front();
				toDo.pop_front();
				delete pThisBlock;
			}
			return(false);
		}
	}

	return(true);
}



bool Island::GetConnectedFaces(WalkBlock* pThisBlock, list<WalkBlock*> &toDo)
{
	bool match = false;

	pThisBlock->m_pFace1->m_visited = true;
	pThisBlock->m_pFace2->m_visited2 = true;

	int Edge1Index = pThisBlock->m_Edge1Index;
	int Edge2Index = pThisBlock->m_Edge2Index;
	RKPolygon* pPolyA = pThisBlock->m_pFace1;
	RKPolygon* pPolyB = pThisBlock->m_pFace2;

	if(m_TestSym)
	{
		match = pPolyA->CompareReverse(pPolyB, Edge1Index, Edge2Index);
		Edge2Index--;
		if(Edge2Index < 0) Edge2Index = pPolyB->m_numberVertices-1;
	}
	else
	{
		match = pPolyA->Compare(pPolyB, Edge1Index, Edge2Index);
	}
	if(match == false) return(false);

	int NumberOfEdges = pPolyA->m_numberVertices;

	for(int EdgeCount = 0; EdgeCount < NumberOfEdges; EdgeCount++)
	{
		RKPolygon* pOtherFaceA = pPolyA->m_pEdges[Edge1Index]->GetOtherFace(pPolyA);
		RKPolygon* pOtherFaceB = pPolyB->m_pEdges[Edge2Index]->GetOtherFace(pPolyB);
		
		if(pOtherFaceA != NULL && pOtherFaceB == NULL) return(false);
		if(pOtherFaceA == NULL && pOtherFaceB != NULL) return(false);

		if(pOtherFaceA && pOtherFaceB)			// may be null
		{
			if(pOtherFaceA->m_visited == false && pOtherFaceB->m_visited2 == false)
			{
				int OtherFaceAEdgeIndex;
				int OtherFaceBEdgeIndex;

				RKEdge* pFace1Edge = pPolyA->m_pEdges[Edge1Index];
				RKEdge* pFace2Edge = pPolyB->m_pEdges[Edge2Index];
				int MatchedVertCount = pOtherFaceA->m_numberVertices;

				for(int Index = 0; Index < MatchedVertCount; Index++)
				{
					if(pOtherFaceA->m_pEdges[Index] == pFace1Edge) OtherFaceAEdgeIndex = Index;
					if(pOtherFaceB->m_pEdges[Index] == pFace2Edge) OtherFaceBEdgeIndex = Index;
				}

				if(m_TestSym)
				{
					OtherFaceBEdgeIndex++;
					if(OtherFaceBEdgeIndex == MatchedVertCount) OtherFaceBEdgeIndex = 0;
					match = pOtherFaceA->CompareReverse(pOtherFaceB, OtherFaceAEdgeIndex, OtherFaceBEdgeIndex);
				}
				else
				{
					match = pOtherFaceA->Compare(pOtherFaceB, OtherFaceAEdgeIndex, OtherFaceBEdgeIndex);
				}

				if(match == true)
				{
					WalkBlock* pNewBlock = new WalkBlock;
					pNewBlock->m_pFace1 = pOtherFaceA;
					pNewBlock->m_pFace2 = pOtherFaceB;
					pNewBlock->m_Edge1Index = OtherFaceAEdgeIndex;
					pNewBlock->m_Edge2Index = OtherFaceBEdgeIndex;

					pOtherFaceA->m_visited = true;
					pOtherFaceB->m_visited2 = true;
					toDo.push_back(pNewBlock);
				}
				else
				{
					return(false);
				}
			}
		}

		if(m_TestSym)
		{
			Edge2Index--;
			if(Edge2Index < 0) Edge2Index = pPolyB->m_numberVertices-1;
		}
		else
		{
			Edge2Index++;
			if(Edge2Index == NumberOfEdges) Edge2Index = 0;
		}

		Edge1Index++;
		if(Edge1Index == NumberOfEdges) Edge1Index = 0;
	}

	return(true);
}




// OK, we're looking for two faces, with the same interior angles.
// and the same edge selections.  That should sort, cubes and spheres
// If we find none, then the meshes don't match.
// If we find more than one match, make a list, keep going.

vector<WalkBlock*> Island::FindStartFaces(Island* pIsland)
{
	// equilateral triangles might get added thrice, isocolese, twice
	vector<RKPolygon*> MinimumSet;
	vector<RKPolygon*> ThisSet;
	vector<WalkBlock*> RetList;
	int MinimumSize = 10000;

	MinimumSet.clear();

	for(int SrcIndex = 0; SrcIndex < m_Polygons.size(); SrcIndex++)
	{
		ThisSet.clear();
		RKPolygon* pFace1 = m_Polygons[SrcIndex];

		for(int DestIndex = 0; DestIndex < pIsland->m_Polygons.size(); DestIndex++)
		{
			RKPolygon* pFace2 = pIsland->m_Polygons[DestIndex];
			
			if(m_TestSym)
			{
				if(pFace1->CompareReverse(pFace2))
				{
					ThisSet.push_back(pFace1);
					ThisSet.push_back(pFace2);
				}
			}
			else
			{
				if(pFace1->Compare(pFace2))
				{
					ThisSet.push_back(pFace1);
					ThisSet.push_back(pFace2);
				}
			}
		}

		if(ThisSet.size() < 3) 
		{
			RetList = FirstMatches(ThisSet);
			ThisSet.clear();
			MinimumSet.clear();
			return(RetList);
		}

		if(MinimumSet.size() == 0 || ThisSet.size() < MinimumSize)				//MinimumSet.size())
		{
			MinimumSet.clear();
			MinimumSet = ThisSet;
			MinimumSize = MinimumSet.size();
		}
	}

	RetList = FirstMatches(MinimumSet);			//ThisSet);
	ThisSet.clear();
	MinimumSet.clear();
	return(RetList);
}



vector<WalkBlock*> Island::FindStartFacesSym()
{
	// equilateral triangles might get added thrice, isocolese, twice
	vector<RKPolygon*> MinimumSet;
	vector<RKPolygon*> ThisSet;
	vector<WalkBlock*> RetList;
	int MinimumSize = 10000;

	MinimumSet.clear();

	for(int SrcIndex = 0; SrcIndex < m_Polygons.size(); SrcIndex++)
	{
		ThisSet.clear();
		RKPolygon* pFace1 = m_Polygons[SrcIndex];

		for(int DestIndex = 0; DestIndex < m_Polygons.size(); DestIndex++)
		{
			if(SrcIndex != DestIndex)
			{
				RKPolygon* pFace2 = m_Polygons[DestIndex];
				
				if(pFace1->CompareReverse(pFace2))
				{
					ThisSet.push_back(pFace1);
					ThisSet.push_back(pFace2);
				}
			}
		}

		// only two polys matching in island
		if(ThisSet.size() == 2)				//< 3) 
		{
			RetList = FirstMatches(ThisSet);
			ThisSet.clear();
			MinimumSet.clear();
			return(RetList);
		}

		// 0 here could be polygon only matching with itself, but that could happen in plane of reflection and still be symmetrical
		if(ThisSet.size() != 0 && ThisSet.size() < MinimumSize)				// < MinimumSet.size())
		{
			MinimumSet.clear();
			MinimumSet = ThisSet;
			MinimumSize = MinimumSet.size();
		}
	}

	RetList = FirstMatches(MinimumSet);			//ThisSet);
	ThisSet.clear();
	MinimumSet.clear();
	return(RetList);
}



vector<WalkBlock*> Island::FirstMatches(vector<RKPolygon*> matchedFaces)
{
	vector<WalkBlock*> matchedBlocks;


	for(int StartIndex = 0; StartIndex < matchedFaces.size(); StartIndex+=2)
	{
		RKPolygon* pFaceA = matchedFaces[StartIndex];
		RKPolygon* pFaceB = matchedFaces[StartIndex+1];


		for(int RotStart = 0; RotStart < pFaceA->m_numberVertices; RotStart++)
		{
			double Diff;
			double MaxDiff = 0.0;
			int DstIndex = RotStart;

			for(int SrcIndex = 0; SrcIndex < pFaceA->m_numberVertices; SrcIndex++)
			{
				Diff = fabs(pFaceA->m_InteriorAngles[SrcIndex] - pFaceB->m_InteriorAngles[DstIndex]);
				if(Diff > MaxDiff) MaxDiff = Diff;
				
				if(m_TestSym)
				{
					DstIndex--;
					if(DstIndex < 0) DstIndex = pFaceA->m_numberVertices-1;
				}
				else
				{
					DstIndex++;
					if(DstIndex == pFaceA->m_numberVertices) DstIndex = 0;
				}
			}

			if(Diff < gMatchAngleTolerence)
			{
				WalkBlock* pNewBlock = new WalkBlock();
				pNewBlock->m_pFace1 = pFaceA;
				pNewBlock->m_pFace2 = pFaceB;
				pNewBlock->m_Edge1Index = 0;
				pNewBlock->m_Edge2Index = RotStart;
				MatchAngles(pNewBlock);
				matchedBlocks.push_back(pNewBlock);
			}
		}
	}

	return(matchedBlocks);
}




void Island::MatchAngles(WalkBlock* pBlock)
{
	RKPolygon* pPoly1 = pBlock->m_pFace1;
	RKPolygon* pPoly2 = pBlock->m_pFace2;
	int VertFace2 = pBlock->m_Edge2Index - pBlock->m_Edge1Index;
	int BestVertFace1 = 0;
	int BestVertFace2 = VertFace2;
	double BestAngle = 1.57;

	int NumberOfVertices = pBlock->m_pFace1->m_numberVertices;
	if(NumberOfVertices > 3)			// find the angle closest to 90 degrees
	{
		for(int Index = 0; Index < NumberOfVertices; Index++)
		{
			double angle = pBlock->m_pFace1->m_InteriorAngles[Index];
			double test = abs(1.5707 - angle);
			if(test < BestAngle)
			{
				BestAngle = test;
				BestVertFace1 = Index;
				BestVertFace2 = VertFace2;
			}

			if(m_TestSym)
			{
				VertFace2--;
				if(VertFace2 < 0) VertFace2 = NumberOfVertices-1;			
			}
			else
			{
				VertFace2++;
				if(VertFace2 == NumberOfVertices) VertFace2 = 0;
			}
		}
	}


	float LookAtFace1[3][3];
	float LookAtFace2[3][3];
	float LookAtFace3[3][3];
	float LookAtFace4[3][3];
	float UpY[3] = {0, -1.0f, 0};
	float UpZ[3] = {0, 0, 1.0f};

	bool Face1Valid = LookAtPoly(pPoly1, BestVertFace1, LookAtFace1, UpY, 0);
	bool Face2Valid = LookAtPoly(pPoly2, BestVertFace2, LookAtFace2, UpY, m_TestSym);
	bool Face3Valid = LookAtPoly(pPoly1, BestVertFace1, LookAtFace3, UpZ, 0);
	bool Face4Valid = LookAtPoly(pPoly2, BestVertFace2, LookAtFace4, UpZ, m_TestSym);

	float* pProjFace1 = new float[NumberOfVertices*3];
	float* pProjFace2 = new float[NumberOfVertices*3];
	float* pProjFace3 = new float[NumberOfVertices*3];
	float* pProjFace4 = new float[NumberOfVertices*3];
	VertFace2 = pBlock->m_Edge2Index - pBlock->m_Edge1Index;

	for(int Index = 0; Index < NumberOfVertices; Index++)
	{
		pProjFace1[Index*3+0] = LookAtFace1[0][0] * pPoly1->m_pVertices[Index]->X + LookAtFace1[1][0] * pPoly1->m_pVertices[Index]->Y + LookAtFace1[2][0] * pPoly1->m_pVertices[Index]->Z;
		pProjFace1[Index*3+1] = LookAtFace1[0][1] * pPoly1->m_pVertices[Index]->X + LookAtFace1[1][1] * pPoly1->m_pVertices[Index]->Y + LookAtFace1[2][1] * pPoly1->m_pVertices[Index]->Z;
		pProjFace1[Index*3+2] = LookAtFace1[0][2] * pPoly1->m_pVertices[Index]->X + LookAtFace1[1][2] * pPoly1->m_pVertices[Index]->Y + LookAtFace1[2][2] * pPoly1->m_pVertices[Index]->Z;

		pProjFace2[Index*3+0] = LookAtFace2[0][0] * pPoly2->m_pVertices[VertFace2]->X + LookAtFace2[1][0] * pPoly2->m_pVertices[VertFace2]->Y + LookAtFace2[2][0] * pPoly2->m_pVertices[VertFace2]->Z;
		pProjFace2[Index*3+1] = LookAtFace2[0][1] * pPoly2->m_pVertices[VertFace2]->X + LookAtFace2[1][1] * pPoly2->m_pVertices[VertFace2]->Y + LookAtFace2[2][1] * pPoly2->m_pVertices[VertFace2]->Z;
		pProjFace2[Index*3+2] = LookAtFace2[0][2] * pPoly2->m_pVertices[VertFace2]->X + LookAtFace2[1][2] * pPoly2->m_pVertices[VertFace2]->Y + LookAtFace2[2][2] * pPoly2->m_pVertices[VertFace2]->Z;
	
		pProjFace3[Index*3+0] = LookAtFace3[0][0] * pPoly1->m_pVertices[Index]->X + LookAtFace3[1][0] * pPoly1->m_pVertices[Index]->Y + LookAtFace3[2][0] * pPoly1->m_pVertices[Index]->Z;
		pProjFace3[Index*3+1] = LookAtFace3[0][1] * pPoly1->m_pVertices[Index]->X + LookAtFace3[1][1] * pPoly1->m_pVertices[Index]->Y + LookAtFace3[2][1] * pPoly1->m_pVertices[Index]->Z;
		pProjFace3[Index*3+2] = LookAtFace3[0][2] * pPoly1->m_pVertices[Index]->X + LookAtFace3[1][2] * pPoly1->m_pVertices[Index]->Y + LookAtFace3[2][2] * pPoly1->m_pVertices[Index]->Z;

		pProjFace4[Index*3+0] = LookAtFace4[0][0] * pPoly2->m_pVertices[VertFace2]->X + LookAtFace4[1][0] * pPoly2->m_pVertices[VertFace2]->Y + LookAtFace4[2][0] * pPoly2->m_pVertices[VertFace2]->Z;
		pProjFace4[Index*3+1] = LookAtFace4[0][1] * pPoly2->m_pVertices[VertFace2]->X + LookAtFace4[1][1] * pPoly2->m_pVertices[VertFace2]->Y + LookAtFace4[2][1] * pPoly2->m_pVertices[VertFace2]->Z;
		pProjFace4[Index*3+2] = LookAtFace4[0][2] * pPoly2->m_pVertices[VertFace2]->X + LookAtFace4[1][2] * pPoly2->m_pVertices[VertFace2]->Y + LookAtFace4[2][2] * pPoly2->m_pVertices[VertFace2]->Z;
	
		if(m_TestSym)
		{
			VertFace2--;
			if(VertFace2 < 0) VertFace2 = NumberOfVertices-1;
		}
		else
		{
			VertFace2++;
			if(VertFace2 == NumberOfVertices) VertFace2 = 0;
		}
	}

	double LowestDiff = 10000.0;
	double Difference;

	if(Face1Valid && Face2Valid)
	{
		Difference = TestDifference(NumberOfVertices, pProjFace1, pProjFace2);
		if(Difference < LowestDiff) LowestDiff = Difference;
	}

	if(Face1Valid && Face4Valid)
	{
		Difference = TestDifference(NumberOfVertices, pProjFace1, pProjFace4);
		if(Difference < LowestDiff) LowestDiff = Difference;
	}

	if(Face3Valid && Face2Valid)
	{
		Difference = TestDifference(NumberOfVertices, pProjFace3, pProjFace2);
		if(Difference < LowestDiff) LowestDiff = Difference;
	}

	if(Face3Valid && Face4Valid)
	{
		Difference = TestDifference(NumberOfVertices, pProjFace3, pProjFace4);
		if(Difference < LowestDiff) LowestDiff = Difference;
	}

	pBlock->m_DegreeOfDiff = LowestDiff;

	delete[] pProjFace1;
	delete[] pProjFace2;
	delete[] pProjFace3;
	delete[] pProjFace4;
}


double Island::TestDifference(int NumberOfVertices, float* pProjFace1, float* pProjFace2)
{
	double Difference = 0.0;

	int Vert2 = NumberOfVertices-1;
	for(int Index = 0; Index < NumberOfVertices; Index++)
	{
		float Vec1X = pProjFace1[Vert2*3] - pProjFace1[Index*3];
		float Vec1Y = pProjFace1[Vert2*3+1] - pProjFace1[Index*3+1];
		float Vec2X = pProjFace2[Vert2*3] - pProjFace2[Index*3];
		float Vec2Y = pProjFace2[Vert2*3+1] - pProjFace2[Index*3+1];
		float Len1 = sqrt(Vec1X * Vec1X + Vec1Y * Vec1Y);
		float Len2 = sqrt(Vec2X * Vec2X + Vec2Y * Vec2Y);
		Vec1X /= Len1;	Vec1Y /= Len1;
		Vec2X /= Len2;	Vec2Y /= Len2;

		float Dot = (Vec1X * Vec2X) + (Vec1Y * Vec2Y);
		if(Dot > 1.0f) Dot = 1.0f;
		if(Dot < -1.0f) Dot = -1.0f;
		Difference += acos(Dot);
		Vert2 = Index;
	}

	return(Difference);
}



// returns false if angle it too close to up

bool Island::LookAtPoly(RKPolygon* pPoly, int BestVert, float LookAt[3][3], float Up[3], bool SymTest)
{
	int PrevVert, ThisVert, NextVert;

	if(SymTest)
	{
		NextVert = BestVert - 1;
		if(NextVert < 0) NextVert = pPoly->m_numberVertices-1;
		ThisVert = BestVert;
		PrevVert = BestVert + 1;
		if(PrevVert == pPoly->m_numberVertices) PrevVert = 0;
	}
	else
	{
		PrevVert = BestVert - 1;
		if(PrevVert < 0) PrevVert = pPoly->m_numberVertices-1;
		ThisVert = BestVert;
		NextVert = BestVert + 1;
		if(NextVert == pPoly->m_numberVertices) NextVert = 0;
	}


	float Vec1X = pPoly->m_pVertices[NextVert]->X - pPoly->m_pVertices[ThisVert]->X;
	float Vec1Y = pPoly->m_pVertices[NextVert]->Y - pPoly->m_pVertices[ThisVert]->Y;
	float Vec1Z = pPoly->m_pVertices[NextVert]->Z - pPoly->m_pVertices[ThisVert]->Z;
	float Vec2X = pPoly->m_pVertices[PrevVert]->X - pPoly->m_pVertices[ThisVert]->X;
	float Vec2Y = pPoly->m_pVertices[PrevVert]->Y - pPoly->m_pVertices[ThisVert]->Y;
	float Vec2Z = pPoly->m_pVertices[PrevVert]->Z - pPoly->m_pVertices[ThisVert]->Z;
	float Face1NormX = (Vec1Y * Vec2Z) - (Vec1Z * Vec2Y);
	float Face1NormY = (Vec1Z * Vec2X) - (Vec1X * Vec2Z);
	float Face1NormZ = (Vec1X * Vec2Y) - (Vec1Y * Vec2X);

	float Len = sqrt(Face1NormX * Face1NormX + Face1NormY * Face1NormY + Face1NormZ * Face1NormZ);
	Face1NormX /= Len;  Face1NormY /= Len;  Face1NormZ /= Len; 
	LookAt[0][2] = Face1NormX;
	LookAt[1][2] = Face1NormY;
	LookAt[2][2] = Face1NormZ;

	float Dot = Face1NormX * Up[0] + Face1NormY * Up[1] + Face1NormZ * Up[2];
	if(Dot > 0.996f) return(false);
	if(Dot < -0.996) return(false);


	float XVecX = (Up[1] * Face1NormZ) - (Up[2] * Face1NormY);
	float XVecY = (Up[2] * Face1NormX) - (Up[0] * Face1NormZ);
	float XVecZ = (Up[0] * Face1NormY) - (Up[1] * Face1NormX);

	Len = sqrt(XVecX * XVecX + XVecY * XVecY + XVecZ * XVecZ);
	XVecX /= Len;  XVecY /= Len;  XVecZ /= Len;
	LookAt[0][0] = XVecX;
	LookAt[1][0] = XVecY;
	LookAt[2][0] = XVecZ;

	float YVecX = (Face1NormY * XVecZ) - (Face1NormZ * XVecY);
	float YVecY = (Face1NormZ * XVecX) - (Face1NormX * XVecZ);
	float YVecZ = (Face1NormX * XVecY) - (Face1NormY * XVecX);

	Len = sqrt(YVecX * YVecX + YVecY * YVecY + YVecZ * YVecZ);
	YVecX /= Len;  YVecY /= Len;  YVecZ /= Len;
	LookAt[0][1] = YVecX;
	LookAt[1][1] = YVecY;
	LookAt[2][1] = YVecZ;

	return(true);
}




void Island::PropagateUVs()
{
	if(m_pSameAs == NULL && this->m_pSymmetricalWith == NULL) return;

	for(int VertIndex = 0; VertIndex < m_Vertices.size(); VertIndex++)
	{
		RKVertex* pSrcVert = m_Vertices[VertIndex]->m_pMatches;
		m_Vertices[VertIndex]->U = pSrcVert->U;
		m_Vertices[VertIndex]->V = pSrcVert->V;
	}
}




bool Island::TestSymmetrical(Island* pIsland)
{
	vector<WalkBlock*> StartFaces = FindStartFaces(pIsland);
	if(StartFaces.size() == 0) return(false);
	if(StartFaces.size() > 1)
	{
		std::sort(StartFaces.begin(), StartFaces.end(), CompareBlocks);
	}

	// from start face, walk out

	for(int StartIndex = 0; StartIndex < StartFaces.size(); StartIndex++)
	{
		bool worked = WalkOutMatch(pIsland, StartFaces[StartIndex]);
		if(worked == true) 
		{
//			m_pSymmetricalWith = pIsland;
			pIsland->m_pSymmetricalWith = this;
			return(true);
		}
	}

	return(false);
}



bool Island::TestInteriorSymmetry()
{
	DivisionPlane* pSplitPlane = FindSymmetryPlane();

	if(pSplitPlane)
	{
		m_Symmetrical = true;
		m_DivisionPlane.m_A = pSplitPlane->m_A;
		m_DivisionPlane.m_B = pSplitPlane->m_B;
		m_DivisionPlane.m_C = pSplitPlane->m_C;
		m_DivisionPlane.m_D = pSplitPlane->m_D;
		return(true);
	}

	return(false);
}




DivisionPlane* Island::FindSymmetryPlane()
{
	vector<DivisionPlane*> m_Planes;

	m_CentreX = m_CentreY = m_CentreZ = 0.0;
	float dummy;
	bool Symmetrical;

	int NumberOfVertices = m_Vertices.size();
	for(int Index = 0; Index < NumberOfVertices; Index++)
	{
		m_CentreX += m_Vertices[Index]->X;
		m_CentreY += m_Vertices[Index]->Y;
		m_CentreZ += m_Vertices[Index]->Z;
	}

	m_CentreX /= NumberOfVertices;
	m_CentreY /= NumberOfVertices;
	m_CentreZ /= NumberOfVertices;


	RKVertex* pStartVertex = PinFinder.m_Boundaries[0]->pBoundary[0];
	RKVertex* pThisVertex = pStartVertex;
	RKVertex* pNextVertex = pStartVertex->GetNextTrailing(&dummy);
	

	// only need to travel half the boundary

	while(pNextVertex != pStartVertex)						// walk around the boundary
	{
		// what to do if it's not a split vertex

		double VNormX, VNormY, VNormZ;				// vertex norm
		double FNormX, FNormY, FNormZ;				// face normal
		double CNormX, CNormY, CNormZ;				// normal to centre of object
		double ENormX, ENormY, ENormZ;				// Edge Normal
		double MidX, MidY, MidZ;					// Edge Midpoint

		// edge mid point
		MidX = (pThisVertex->X + pNextVertex->X) * 0.5f;
		MidY = (pThisVertex->Y + pNextVertex->Y) * 0.5f;
		MidZ = (pThisVertex->Z + pNextVertex->Z) * 0.5f;

		// Mid edge to centre of object normal, needs normalising
		CNormX = MidX - m_CentreX;
		CNormY = MidY - m_CentreY;
		CNormZ = MidZ - m_CentreZ;

		pThisVertex->getFaceNormal(FNormX, FNormY, FNormZ, pNextVertex);
		pThisVertex->getVertexNormal(VNormX, VNormY, VNormZ);

		// Face Normal and Normal to Centre plane test
		Symmetrical = PlaneTest(FNormX, FNormY, FNormZ, CNormX, CNormY, CNormZ, m_Planes);


		// Vertex to centre of object normal, needs normalising
		CNormX = pThisVertex->X - m_CentreX;
		CNormY = pThisVertex->Y - m_CentreY;
		CNormZ = pThisVertex->Z - m_CentreZ;

		// Edge normal, needs normalising
		ENormX = pThisVertex->X - pNextVertex->X;
		ENormY = pThisVertex->Y - pNextVertex->Y;
		ENormZ = pThisVertex->Z - pNextVertex->Z;

		// Vertex Norm and	vertex to centre Normal plane test
		Symmetrical = PlaneTest(VNormX, VNormY, VNormZ, CNormX, CNormY, CNormZ, m_Planes);

		// if the edge contains a split vertex, this is the only test required.
		// Edge Norm and vertex to centre Normal plane test
		Symmetrical = PlaneTest(ENormX, ENormY, ENormZ, CNormX, CNormY, CNormZ, m_Planes);


		pThisVertex = pNextVertex;
		pNextVertex = pNextVertex->GetNextTrailing(&dummy);
	}

	if(m_Planes.size() == 0) return(NULL);
	if(m_Planes.size() == 1)
	{
		DivisionPlane* pThePlane = m_Planes.front();
		m_Planes.clear();
		return(pThePlane);
	}

	
	// TODO
	// get plane that is closest to one of the major axis'
	DivisionPlane* pThePlane = m_Planes.front();
	m_Planes.clear();
	return(pThePlane);
}



bool Island::PlaneTest(double Norm1X, double Norm1Y, double Norm1Z, double Norm2X, double Norm2Y, double Norm2Z, vector<DivisionPlane*> &SymPlanes)
{
	double PlaneA = (Norm1Y * Norm2Z) - (Norm2Y * Norm1Z);
	double PlaneB = (Norm1Z * Norm2X) - (Norm2Z * Norm1X);
	double PlaneC = (Norm1X * Norm2Y) - (Norm2X * Norm1Y);
	double maxDist = 0.0;
	int FurthestLeft = -1;

	double Length = (PlaneA * PlaneA) + (PlaneB * PlaneB) + (PlaneC * PlaneC);
	Length = sqrt(Length);
		
	PlaneA /= Length;
	PlaneB /= Length;
	PlaneC /= Length;
	double PlaneD = -((PlaneA * (m_CentreX)) + (PlaneB * (m_CentreY)) + (PlaneC * (m_CentreZ)));

	// do some extra shizzle for when we end up with a line instead of a plane..
	// use the normal of the polygon attached to the edge to make plane with

	int LeftCount = 0;
	int RightCount = 0;
	int NumberOfVertices = m_Vertices.size();

	for(int Index = 0; Index < NumberOfVertices; Index++)
	{
		double Result = PlaneA * m_Vertices[Index]->X + PlaneB * m_Vertices[Index]->Y + PlaneC * m_Vertices[Index]->Z + PlaneD;
		double absResult = fabs(Result);

		if(absResult > maxDist && Result > 0.0)
		{
			maxDist = absResult;
			FurthestLeft = Index;
		}

		m_Vertices[Index]->m_distFromPlane = Result;

		if(absResult > 0.00001)
		{
			if(Result < 0) 
				LeftCount++;
			else
				RightCount++;
		}
	}

	if(LeftCount == RightCount)
	{
		if(FurthestLeft == -1) return(false);
		bool result = MirrorTest(PlaneA, PlaneB, PlaneC, PlaneD, FurthestLeft, maxDist);
		if(result == true)
		{
			DivisionPlane* Plane = new DivisionPlane;
			Plane->m_A = PlaneA;
			Plane->m_B = PlaneB;
			Plane->m_C = PlaneC;
			Plane->m_D = PlaneD;

			SymPlanes.push_back(Plane);
			return(true);
		}
	}

	return(false);
}





bool Island::MirrorTest(double PlaneA, double PlaneB, double PlaneC, double PlaneD, int FurthestLeft, double maxDist)
{
	int NumberOfVertices = m_Vertices.size();


	for(int VertIndex = 0; VertIndex < NumberOfVertices; VertIndex++)
	{
		double RightDist = m_Vertices[VertIndex]->m_distFromPlane;

		if(RightDist < 0.0)
		{
			double testDist = fabs(fabs(RightDist) - maxDist);
			if(testDist < 0.01)												// tolerence
			{
				bool result = MirrorOverlap(FurthestLeft, VertIndex);
				if(result == true) return(true);
			}
		}
	}

	return(false);
}



bool Island::MirrorOverlap(int LeftIndex, int RightIndex)
{
	int MinimumSize = 10000;
	m_TestSym = true;
	WalkBlock testBlock;

	vector<WalkBlock*> StartFaces;
	vector<RKPolygon*> ThisSet;
	vector<RKPolygon*> MinimumSet;

	list<RKPolygon*>::iterator LeftFacesIterator;
	list<RKPolygon*>::iterator RightFacesIterator;

	for(LeftFacesIterator = m_Vertices[LeftIndex]->m_LinkedPolygons.begin(); LeftFacesIterator != m_Vertices[LeftIndex]->m_LinkedPolygons.end(); LeftFacesIterator++)
	{
		ThisSet.clear();
		RKPolygon* pFace1 = *LeftFacesIterator;

		for(RightFacesIterator = m_Vertices[RightIndex]->m_LinkedPolygons.begin(); RightFacesIterator != m_Vertices[RightIndex]->m_LinkedPolygons.end(); RightFacesIterator++)
		{
			RKPolygon* pFace2 = *RightFacesIterator;
				
			if(pFace1->CompareReverse(pFace2))
			{
				ThisSet.push_back(pFace1);
				ThisSet.push_back(pFace2);
			}
		}

		// only two polys matching in island
		if(ThisSet.size() == 2)
		{
			StartFaces = FirstMatches(ThisSet);
			ThisSet.clear();
			MinimumSet.clear();
			goto Test;
		}


		// 0 here could be polygon only matching with itself, but that could happen in plane of reflection and still be symmetrical
		if(ThisSet.size() != 0 && ThisSet.size() < MinimumSize)				// < MinimumSet.size())
		{
			MinimumSet.clear();
			MinimumSet = ThisSet;
			MinimumSize = MinimumSet.size();
		}
	}

	StartFaces = FirstMatches(MinimumSet);


Test:

	ThisSet.clear();
	MinimumSet.clear();


	if(StartFaces.size() == 0) return(false);

	// from start face, walk out

	for(int StartIndex = 0; StartIndex < StartFaces.size(); StartIndex++)
	{
		bool worked = WalkOutMatch(this, StartFaces[StartIndex]);
		if(worked == true) 
		{
			// delete unused StartFaces
			return(true);
		}
	}


	return(false);
}




// find the outer UV points that cross the plane of symmetry

void Island::SymUVPoints(float& U1, float& V1, float& U2, float& V2)
{
	vector<float> CrossPoints;
	float dummy;

	RKVertex* pStartVertex = PinFinder.m_Boundaries[0]->pBoundary[0];
	RKVertex* pThisVertex = pStartVertex;
	RKVertex* pNextVertex = pStartVertex->GetNextTrailing(&dummy);


	while(pNextVertex != pStartVertex)						// walk around the boundary
	{
		double  LeftResult = m_DivisionPlane.m_A * pThisVertex->X + m_DivisionPlane.m_B * pThisVertex->Y + m_DivisionPlane.m_C * pThisVertex->Z + m_DivisionPlane.m_D;
		double RightResult = m_DivisionPlane.m_A * pNextVertex->X + m_DivisionPlane.m_B * pNextVertex->Y + m_DivisionPlane.m_C * pNextVertex->Z + m_DivisionPlane.m_D;

		if(LeftResult <= 0 && RightResult > 0 || LeftResult >= 0 && RightResult < 0 )
		{
			double Cross = (RightResult - LeftResult);
			Cross = LeftResult / Cross;

			float UPoint = pNextVertex->U - pThisVertex->U;
			float VPoint = pNextVertex->V - pThisVertex->V;
			UPoint = (UPoint*Cross) + pNextVertex->U;
			VPoint = (VPoint*Cross) + pNextVertex->V;

			CrossPoints.push_back(UPoint);
			CrossPoints.push_back(VPoint);
		}

		pThisVertex = pNextVertex;
		pNextVertex = pNextVertex->GetNextTrailing(&dummy);
	}


	if(CrossPoints.size() == 4)
	{
		// need to do something cunning with the uprighter to get the 'top' and 'bottom' vertices
		U1 = CrossPoints[0];
		V1 = CrossPoints[1];
		U2 = CrossPoints[2];
		V2 = CrossPoints[3];
	}


	// find the two furthest apart points
	// don't need to do this!, any two points form a line
}




// Two symmetrical islands, which one is to the left of the plane?
// find the two nearest and furthest points in 3D space
// island to the left of the plane should have furthest U to the right of nearest U
// so if the UV island is upside down it'll flip.   Shit..

bool Island::LeftSymTest()
{
	float MinDist = 1e18;
	float MaxDist = -1;
	RKVertex* pNearest = NULL;
	RKVertex* pFurthest = NULL;


	for(int Index = 0; Index < m_Vertices.size(); Index++)
	{
		RKVertex* pVert = m_Vertices[Index];

		if(pVert->m_pMatches != NULL)
		{

			// does m_pMatches always point to something?
			float DeltaX = pVert->X - pVert->m_pMatches->X;
			float DeltaY = pVert->Y - pVert->m_pMatches->Y;
			float DeltaZ = pVert->Z - pVert->m_pMatches->Z;
	}
	}

	return(true);
}





void Island::SquareGrid()
{
	float Angle1, Angle2, Angle3;

	for(int FaceIndex = 0; FaceIndex < m_Faces.size(); FaceIndex++)
	{
		if(m_Faces[FaceIndex]->NinetyAngle == 2) {Angle2 =  1.57079633;  Angle1 = 0.785398163;  Angle3 = 0.785398163; }
		if(m_Faces[FaceIndex]->NinetyAngle == 3) {Angle3 =  1.57079633;  Angle1 = 0.785398163;  Angle2 = 0.785398163; }
		if(m_Faces[FaceIndex]->NinetyAngle == 1) {Angle1 =  1.57079633;  Angle2 = 0.785398163;  Angle3 = 0.785398163; }

		m_Faces[FaceIndex]->m_ABFAngle1 = Angle1;
		m_Faces[FaceIndex]->m_ABFAngle2 = Angle2;
		m_Faces[FaceIndex]->m_ABFAngle3 = Angle3;
	}
}



bool Island::CheckStraight()
{
	for(int FaceIndex = 0; FaceIndex < m_Faces.size(); FaceIndex++)
	{
		float Angle1, Angle2, Angle3;

		m_Faces[FaceIndex]->GetUVAngles(Angle1, Angle2, Angle3);

		float Diff1 = abs(Angle1 - 1.57079633f);
		float Diff2 = abs(Angle2 - 1.57079633f);
		float Diff3 = abs(Angle3 - 1.57079633f);

		if(!(Diff1 < 0.034906585 || Diff2 < 0.034906585 || Diff3 < 0.034906585)) return(false);		// more than 2 degrees out of whack
	}

	return(true);
}




void Island::AverageEdges()
{
	vector<RKVertex*> Boundary = PinFinder.m_Boundaries[0]->pBoundary;
	RKVertex* pVert = Boundary[0];
	RKVertex* pVertNext;
	float Dummy;
	float EdgeLen3D;
	float AverageLen;

	vector<int> HitFacesList;

	int BoundarySize = Boundary.size();		///2 + 2;

	for(int Index = 0; Index < BoundarySize; Index++)
	{
		HitFacesList.clear();
		AverageLen = 0.0f;

		pVertNext = pVert->GetNextTrailing(&Dummy);

		// midpoint of the edge and a point 90 degrees away from the edge

		float EdgeU1 = (pVert->U + pVertNext->U) * 0.5f;
		float EdgeV1 = (pVert->V + pVertNext->V) * 0.5f;
		float EdgeU2 = (pVertNext->V - pVert->V) + EdgeU1;
		float EdgeV2 = -(pVertNext->U - pVert->U) + EdgeV1;


		for(int FaceIndex = 0; FaceIndex < m_Faces.size(); FaceIndex++)
		{
			bool Hit = m_Faces[FaceIndex]->IntersectUVs(EdgeU1, EdgeV1, EdgeU2, EdgeV2, EdgeLen3D);
			if(Hit) 
			{
				HitFacesList.push_back(FaceIndex);
				AverageLen += EdgeLen3D;
			}
		}

		int NumHitFaces = HitFacesList.size();

		if(NumHitFaces != 0)
		{
			AverageLen /= NumHitFaces;

			for(int FaceIndex  = 0; FaceIndex < NumHitFaces; FaceIndex++)
			{
				int ThisFace = HitFacesList[FaceIndex];
				m_Faces[ThisFace]->SetEdgeLen(AverageLen);
			}
		}

		pVert = pVertNext;
	}



	// work out interior angles

	for(int FaceIndex = 0; FaceIndex < m_Faces.size(); FaceIndex++)
	{
		m_Faces[FaceIndex]->NinteyAngles();
	}
}




void Island::MaintainUVs()
{
	float CentreAX = 0.0f;
	float CentreAY = 0.0f;
	float CentreBX = 0.0f;
	float CentreBY = 0.0f;

	unsigned int NumberOfVertices = m_Vertices.size();


	for (unsigned int Index = 0; Index < NumberOfVertices; Index++)
	{
		CentreAX += m_Vertices[Index]->m_OldU;
		CentreAY += m_Vertices[Index]->m_OldV;
		CentreBX += m_Vertices[Index]->m_NewU;
		CentreBY += m_Vertices[Index]->m_NewV;
	}

	CentreAX /= NumberOfVertices;
	CentreAY /= NumberOfVertices;
	CentreBX /= NumberOfVertices;
	CentreBY /= NumberOfVertices;


	double R00, R01, R10, R11;
	R00 = R01 = R10 = R11 = 0.0f;

	double f_sd2_tar = 0, f_sd2 = 0;

	for (unsigned int Index = 0; Index < m_Vertices.size(); Index++)
	{
		float V1X = m_Vertices[Index]->m_OldU - CentreAX;
		float V1Y = m_Vertices[Index]->m_OldV - CentreAY;
		float V2X = m_Vertices[Index]->m_NewU - CentreBX;
		float V2Y = m_Vertices[Index]->m_NewV - CentreBY;

		R00 += V1X * V2X;
		R01 += V1X * V2Y;
		R10 += V1Y * V2X;
		R11 += V1Y * V2Y;

		f_sd2 += (V1X * V1X + V1Y * V1Y);
		f_sd2_tar += (V2X * V2X + V2Y * V2Y);
	}


	Matrix2d AMat(2, 2);

	AMat(0, 0) = R00;
	AMat(1, 0) = R01;
	AMat(0, 1) = R10;
	AMat(1, 1) = R11;

	JacobiSVD<Matrix2d> svdOfA(AMat, ComputeFullU | ComputeFullV);

	Eigen::Matrix2d R = svdOfA.matrixV() * svdOfA.matrixU().transpose();

	// calculate determinant of V*U^T to disambiguate rotation sign
	double f_det = R.determinant();
	Eigen::Vector2d e(1, (f_det < 0) ? -1 : 1);

	// recompute the rotation part if the determinant was negative
	if (f_det < 0)
		R.noalias() = svdOfA.matrixV() * e.asDiagonal() * svdOfA.matrixU().transpose();

	// renormalize the rotation (not needed but gives slightly more orthogonal transformations)
	//R = Eigen::Quaterniond(R).normalized().toRotationMatrix();

	// calculate the scale

	Eigen::Vector2d s = svdOfA.singularValues();

	double f_scale = s.dot(e) / f_sd2_tar;
	double f_inv_scale = s.dot(e) / f_sd2;

	R *= f_scale;

	R00 = R(0, 0);
	R01 = R(1, 0);
	R10 = R(0, 1);
	R11 = R(1, 1);


	for (unsigned int Index = 0; Index < m_Vertices.size(); Index++)
	{
		float NewU = m_Vertices[Index]->m_NewU - CentreBX;
		float NewV = m_Vertices[Index]->m_NewV - CentreBY;

		float FittedU = (NewU * R00) + (NewV * R10);
		float FittedV = (NewU * R01) + (NewV * R11);

		m_Vertices[Index]->U = FittedU + CentreAX;
		m_Vertices[Index]->V = FittedV + CentreAY;
	}
}