#include <algorithm>

#include <maya/MGlobal.h>			// temp
#include <math.h>

#include "RKProgress.h"
#include "Unwrap.h"


int gFunction;
int gScale;
int gUpright;
int gMaintainUVs = false;
int gMatch = false;
int gEdgeLen = false;
int gSymMatch = false;

double gMatchAngleTolerence = 4.0;
double gMatchEdgeTolerence = 0.1;

int gOuterBorderSwitch;
int gOuterBorderPixels;
int gOuterBorderTMapSize;


vector<RKPolygon*> Unwrap::m_Polygons;


/****************************************************
*		Create UV islands and unwrap them			*
****************************************************/

void Unwrap::DoUnwrap(bool UseCPMS)
{
    RKProgress::Get().SetText("Creating Islands");      //CreatingIslandsString);
	RKProgress::Get().SetProgress(0);

	// Scale the meshes to be in a 1.0 sphere, because packer is broken
	ScaleMeshes();

	CreateIslands();

//	return;										// temp!!
    
    
    //if((gFunction == 2 || gFunction == 3) && gMaintainUVs) return;

	bool Failed = false;
	float Scale = (float)gScale / 100.0f;
	Scale += 1.0f; 

	if(gOuterBorderSwitch)
	{
		float OuterRatio = (float)gOuterBorderPixels / (float)gOuterBorderTMapSize;
		if(OuterRatio > 0.45)
		{
			gOuterBorderPixels = (int)((float)gOuterBorderTMapSize * 0.45);
		}

		double Extra = (double)gOuterBorderPixels / (double)gOuterBorderTMapSize;
		Scale += Extra / 2;
	}

//	MGlobal::displayError( "unwrap islands" );
	if(gFunction == 0 || gFunction == 1 || gFunction == 7)
	{
		Failed = UnwrapIslands(gFunction);				// this Failed is not working.  Fix!!
	}

	if(gFunction == 4)
	{
		Failed = MinimiseIslands();
		return;
	}


	RKProgress::Get().SetText("Packing Islands");
	RKProgress::Get().SetProgress(0);
	if(Failed == false) 
	{
		if((bool)gMatch == true || (bool)gSymMatch == true)
		{
			vector<Island*> matchedIslands;

			IslandMatcher.MatchIslands(m_listOfIslands, matchedIslands);
			IslandPacker.PackIslands(matchedIslands, Scale);
			IslandMatcher.Propagate(m_listOfIslands);
		}
		else
		{
			IslandPacker.PackIslands(m_listOfIslands, Scale);
		}
	}


	if (gMaintainUVs && (gFunction == 0 || gFunction == 1 || gFunction == 7))
	{
		FixZeroArea();
		return;
	}
		


	// final check and scale of UVs
	float MinU = 1.0f;
	float MinV = 1.0f;
	float MaxU = 0.0f;
	float MaxV = 0.0f;
	int numVertices = (int)m_Vertices.size();

	for(int Index = 0; Index < numVertices; Index++)
	{
		RKVertex* pVert = m_Vertices[Index];
		if(pVert->U < MinU) MinU = pVert->U;
		if(pVert->V < MinV) MinV = pVert->V;
		if(pVert->U > MaxU) MaxU = pVert->U;
		if(pVert->V > MaxV) MaxV = pVert->V;
	}

	float OffsetU = 0.0f;
	float OffsetV = 0.0f;
	float ScaleUV = 1.0f;

	if(MinU < 0.0f || MinV < 0.0f)
	{
		OffsetU = -MinU;
		OffsetV = -MinV;
	}

	float OverU = MaxU-MinU;
	float OverV = MaxV-MinV;

//	if(OverU > 1.0f || OverV > 1.0f)
//	{
		if(OverU > OverV) ScaleUV = 1.0f/OverU;
		else
			ScaleUV = 1.0f/OverV;
//	}

	for(int Index = 0; Index < numVertices; Index++)
	{
		RKVertex* pVert = m_Vertices[Index];
		pVert->U += OffsetU;
		pVert->V += OffsetV;
		pVert->U *= ScaleUV;
		pVert->V *= ScaleUV;
	}

	if(gOuterBorderSwitch)
	{
		double shrink = (double)(gOuterBorderTMapSize - gOuterBorderPixels*2) / (double)gOuterBorderTMapSize;

		for(int Index = 0; Index < numVertices; Index++)
		{
			RKVertex* pVert = m_Vertices[Index];
			pVert->U -= 0.5f;
			pVert->V -= 0.5f;
			pVert->U *= shrink;
			pVert->V *= shrink;
			pVert->U += 0.5f;
			pVert->V += 0.5f;
		}
	}


	FixZeroArea();
}




void Unwrap::ScaleMeshes()
{
	float CentreX = 0.0f;
	float CentreY = 0.0f;
	float CentreZ = 0.0f;
	float Radius = 0.0f;


	unsigned NumberOfVertices = (int)m_Vertices.size();

	for(unsigned int VertIndex = 0; VertIndex < NumberOfVertices; VertIndex++)
	{
		CentreX += m_Vertices[VertIndex]->X;
		CentreY += m_Vertices[VertIndex]->Y;
		CentreZ += m_Vertices[VertIndex]->Z;
	}

	CentreX /= (float)NumberOfVertices;
	CentreY /= (float)NumberOfVertices;
	CentreZ /= (float)NumberOfVertices;


	for(unsigned int VertIndex = 0; VertIndex < NumberOfVertices; VertIndex++)
	{
		float DistX, DistY, DistZ;
		float Result;

		DistX = CentreX - m_Vertices[VertIndex]->X;
		DistY = CentreY - m_Vertices[VertIndex]->Y;
		DistZ = CentreZ - m_Vertices[VertIndex]->Z;

		Result = DistX * DistX + DistY * DistY + DistZ * DistZ;
		Result = sqrt(Result);

		if(Result > Radius) Radius = Result;
	}


	for(unsigned int VertIndex = 0; VertIndex < NumberOfVertices; VertIndex++)
	{
		m_Vertices[VertIndex]->X -= CentreX;				// centralise
		m_Vertices[VertIndex]->Y -= CentreY;
		m_Vertices[VertIndex]->Z -= CentreZ;

		m_Vertices[VertIndex]->X /= Radius;					// normalise
		m_Vertices[VertIndex]->Y /= Radius;
		m_Vertices[VertIndex]->Z /= Radius;
	}
}


/****************************************************
*		Unwrap all islands							*
****************************************************/

bool Unwrap::UnwrapIslands(bool UseCPMS)
{

	int IslandCount = (int)m_listOfIslands.size();
	int IslandIndex = 1;
	char Message[255];

	vector<Island*>::iterator IslandIterator;


	for(IslandIterator = m_listOfIslands.begin(); IslandIterator != m_listOfIslands.end(); IslandIterator++)
	{
		sprintf(Message, "Unwrapping Island: %d/%d", IslandIndex, IslandCount);
		RKProgress::Get().SetText(Message);
		RKProgress::Get().SetProgress(0);

		if(gFunction == 7)
		{
			(*IslandIterator)->Straighten();
		}
		else
		{
			(*IslandIterator)->Unwrap(UseCPMS);
		}

		//bool failed = (*IslandIterator)->Unwrap(UseCPMS);
		IslandIndex++;
	}

	return(false);
}



bool Unwrap::MinimiseIslands()
{
	int IslandCount = (int)m_listOfIslands.size();
	int IslandIndex = 1;
	char Message[255];

	vector<Island*>::iterator IslandIterator;

	for(IslandIterator = m_listOfIslands.begin(); IslandIterator != m_listOfIslands.end(); IslandIterator++)
	{
		sprintf(Message, "Minimising Island: %d/%d", IslandIndex, IslandCount);
		RKProgress::Get().SetText(Message);
		RKProgress::Get().SetProgress(0);

		//bool failed =
        (*IslandIterator)->Minimise(IslandIndex, IslandCount);
		IslandIndex++;
	}

	return(false);
}



/****************************************************
*		Set the number of UVs in mesh				*
****************************************************/

int Unwrap::SetVertexCount(int Size)
{
	int Base = (int)m_Vertices.size();
	m_Vertices.resize(Base + Size);


	for(int Index = 0; Index < Size; Index++)
	{
		m_Vertices[Base + Index] = new RKVertex();
	}

	return(Base);
}



/****************************************************
*		Add a vertex to the list					*
****************************************************/

void Unwrap::AddVertex(int UVIndex, int OriginalVertIndex, float X, float Y, float Z, float U, float V)
{
	m_Vertices[UVIndex]->X = X;
	m_Vertices[UVIndex]->Y = Y;
	m_Vertices[UVIndex]->Z = Z;
	m_Vertices[UVIndex]->U = U;
	m_Vertices[UVIndex]->V = V;
	m_Vertices[UVIndex]->m_OldU = U;
	m_Vertices[UVIndex]->m_OldV = V;
	m_Vertices[UVIndex]->m_Used = true;
	m_Vertices[UVIndex]->m_UVIndex = UVIndex;
	m_Vertices[UVIndex]->m_OriginalVertIndex = OriginalVertIndex;
}



/****************************************************
*		Add a face, also add face to each vertex	*
****************************************************/

bool Unwrap::AddFace(int VertIndex1, int VertIndex2, int VertIndex3, int PolygonIndex, bool ZeroArea, int NinetyIndex)
{
	RKFace* newFace = new RKFace();

	newFace->NinetyAngle = NinetyIndex;

	newFace->pVert1 = m_Vertices[VertIndex1];
	newFace->pVert2 = m_Vertices[VertIndex2];
	newFace->pVert3 = m_Vertices[VertIndex3];
	newFace->m_pOrigPolygonindex = PolygonIndex;

	float area = newFace->GetArea();
	if((area == 0.0f || ZeroArea == true) && (gFunction == 0 || gFunction == 1))
	{
		// zero area face, set all uv's to 0.0f.
		// If after unwrap one of this triangles verts is 0.0
		// Place UV along UV edge, proportional to point on 3D line

		newFace->pVert1->U = 0.0f;
		newFace->pVert1->V = 0.0f;
		newFace->pVert2->U = 0.0f;
		newFace->pVert2->V = 0.0f;
		newFace->pVert3->U = 0.0f;
		newFace->pVert3->V = 0.0f;

		m_ZeroArea.push_back(newFace);
		//delete newFace;
		return(true);
	}

	m_Vertices[VertIndex1]->m_LinkedFaces.push_back(newFace);
	m_Vertices[VertIndex2]->m_LinkedFaces.push_back(newFace);
	m_Vertices[VertIndex3]->m_LinkedFaces.push_back(newFace);

	m_Triangles.push_back(newFace);

	return(false);
}


int Unwrap::AddPolygon(vector<int> VertIndices)
{
	RKPolygon* newPoly = new RKPolygon();

	newPoly->m_numberVertices = (int)VertIndices.size();

	for(int Index = 0; Index < newPoly->m_numberVertices; Index++)
	{
		newPoly->m_pVertices.push_back(m_Vertices[VertIndices[Index]]);
	}

	int retVal = (int)m_Polygons.size();

	m_Polygons.push_back(newPoly);

	return(retVal);
}



/****************************************************
*		Delete and clear island list				*
****************************************************/

void Unwrap::DeleteIslands()
{
	for(int Index = 0; Index < m_listOfIslands.size(); Index++) 
	{
		delete(m_listOfIslands[Index]);
	}

	m_listOfIslands.clear();
}


void Unwrap::DeleteVertsAndFaces()
{
	DeleteIslands();

	for(int Index = 0; Index < m_ZeroArea.size(); Index++)
	{
		delete m_ZeroArea[Index];
	}

	m_ZeroArea.clear();

/*	for(int Index = 0; Index < m_Triangles.size(); Index++)
	{
		delete m_Triangles[Index];
	}
*/	m_Triangles.clear();


	for(int Index = 0; Index < m_Polygons.size(); Index++)
	{
		delete m_Polygons[Index];
	}
	m_Polygons.clear();


	for(int Index = 0; Index < m_Vertices.size(); Index++)
	{
		delete m_Vertices[Index];
	}
	m_Vertices.clear();
}




/****************************************************
*		Create a list of islands					*
****************************************************/

void Unwrap::CreateIslands()
{
	m_listOfIslands.clear();

	for(int Index = 0; Index < m_Triangles.size(); Index++) 
	{
		m_Triangles[Index]->m_Visited = false;
	}

	LinkEdges();

	int StartFace = FindStartFace();
	int TempRegionCount = 0;

	while(StartFace != -1)
	{
		Island* pIsland = new Island();
		pIsland->SetUpMesh();
		pIsland->InitEdgeHash((int)m_Vertices.size());

		pIsland->AddTriangles(m_Triangles[StartFace]);
		m_listOfIslands.push_back(pIsland);
		TempRegionCount++;

		StartFace = FindStartFace();
	}

	for(int Index = 0; Index < m_listOfIslands.size(); Index++)
	{
		m_listOfIslands[Index]->FinaliseMesh();
	}
}


/********************************************************
*		Find an unvisited face to start a new island	*
********************************************************/

int Unwrap::FindStartFace()
{
	for(int Index = 0; Index < m_Triangles.size(); Index++)
	{
		if(m_Triangles[Index]->m_Visited == false) return(Index);
	}

	return(-1);
}




/********************************************************
*		Link connected triangles via edges				*
********************************************************/

void Unwrap::LinkEdges()
{
	for(int Index = 0; Index < m_Vertices.size(); Index++)
	{
		m_Vertices[Index]->LinkUp();
	}
}



void Unwrap::FixZeroArea()
{
	if(m_ZeroArea.size() == 0) return;


	for(int Limit = 0; Limit < 1000; Limit++)
	{
		bool NewUV = false;

		for(int Index = 0; Index < m_ZeroArea.size(); Index++)
		{
			float DeltaX, DeltaY, DeltaZ, DeltaU, DeltaV;
			float Dist1, Dist2;
			bool UV1Moved, UV2Moved, UV3Moved;
			UV1Moved = UV2Moved = UV3Moved = false;
			RKVertex* pVert1 = m_ZeroArea[Index]->pVert1;
			RKVertex* pVert2 = m_ZeroArea[Index]->pVert2;
			RKVertex* pVert3 = m_ZeroArea[Index]->pVert3;

			if(pVert1->U != 0.0f || pVert1->V != 0.0f) UV1Moved = true;
			if(pVert2->U != 0.0f || pVert2->V != 0.0f) UV2Moved = true;
			if(pVert3->U != 0.0f || pVert3->V != 0.0f) UV3Moved = true;

			if(UV1Moved && UV2Moved && UV3Moved == false)
			{
				// delta UV3
				DeltaX = pVert2->X - pVert1->X;   DeltaY = pVert2->Y - pVert1->Y;  DeltaZ = pVert2->Z - pVert1->Z;
				Dist1 = sqrt(DeltaX * DeltaX + DeltaY * DeltaY + DeltaZ * DeltaZ);
				DeltaX = pVert3->X - pVert1->X;   DeltaY = pVert3->Y - pVert1->Y;  DeltaZ = pVert3->Z - pVert1->Z;
				Dist2 = sqrt(DeltaX * DeltaX + DeltaY * DeltaY + DeltaZ * DeltaZ);

				if(Dist1 != 0.0f)
				{
					float Scalar = Dist2 / Dist1;
					DeltaU = pVert2->U - pVert1->U;
					DeltaV = pVert2->V - pVert1->V;
					pVert3->U = pVert1->U + (DeltaU * Scalar);
					pVert3->V = pVert1->V + (DeltaV * Scalar);
					NewUV = true;
				}
			}

			if(UV1Moved == false && UV2Moved && UV3Moved)
			{
				// delta UV1
				DeltaX = pVert3->X - pVert2->X;   DeltaY = pVert3->Y - pVert2->Y;  DeltaZ = pVert3->Z - pVert2->Z;
				Dist1 = sqrt(DeltaX * DeltaX + DeltaY * DeltaY + DeltaZ * DeltaZ);
				DeltaX = pVert1->X - pVert2->X;   DeltaY = pVert1->Y - pVert2->Y;  DeltaZ = pVert1->Z - pVert2->Z;
				Dist2 = sqrt(DeltaX * DeltaX + DeltaY * DeltaY + DeltaZ * DeltaZ);

				if(Dist1 != 0.0f)
				{
					float Scalar = Dist2 / Dist1;
					DeltaU = pVert3->U - pVert2->U;
					DeltaV = pVert3->V - pVert2->V;
					pVert1->U = pVert2->U + (DeltaU * Scalar);
					pVert1->V = pVert2->V + (DeltaV * Scalar);
					NewUV = true;
				}
			}

			if(UV1Moved && UV2Moved == false && UV3Moved)
			{
				// delta UV2
				DeltaX = pVert1->X - pVert3->X;   DeltaY = pVert1->Y - pVert3->Y;  DeltaZ = pVert1->Z - pVert3->Z;
				Dist1 = sqrt(DeltaX * DeltaX + DeltaY * DeltaY + DeltaZ * DeltaZ);
				DeltaX = pVert2->X - pVert3->X;   DeltaY = pVert2->Y - pVert3->Y;  DeltaZ = pVert2->Z - pVert3->Z;
				Dist2 = sqrt(DeltaX * DeltaX + DeltaY * DeltaY + DeltaZ * DeltaZ);

				if(Dist1 != 0.0f)
				{
					float Scalar = Dist2 / Dist1;
					DeltaU = pVert1->U - pVert3->U;
					DeltaV = pVert1->V - pVert3->V;
					pVert2->U = pVert3->U + (DeltaU * Scalar);
					pVert2->V = pVert3->V + (DeltaV * Scalar);
					NewUV = true;
				}
			}
		}

		if(NewUV == false) return;
	}
}
