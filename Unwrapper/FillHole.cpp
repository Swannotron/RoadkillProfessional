#include <maya/M3dView.h>

#ifdef MAC_PLUGIN
#include <OpenGL/glu.h>
#else
#include <GL/glu.h>
#endif

#include "FillHole.h"
#include "bestfit.h"

#include <float.h>
#include <math.h>




bool gComplex;

int gVertIndex;
double* gVertices;
int gIndices[10000];							// don't like this one bit..  remove!

double* gProjVertices;
double* gTripleProjVerts;



void FillHole::HoleFiller(vector<RKVertex*> &Boundary, vector<RKFace*> &Faces)
{
	// check to see if it works without splitting

	gIndices[0] = -4;				// just in case it goes all funny!
//	m_BoundaryLength = Boundary.size();

	if(Boundary.size() > 500) return;		// too big, fuck off!

	HoleChunk* pWholeBoundary = new HoleChunk();
	for(int VertexIndex = 0; VertexIndex < Boundary.size(); VertexIndex++)
	{
		m_pVertices.push_back(Boundary[VertexIndex]);
		pWholeBoundary->m_pVertIndices.push_back(VertexIndex);
		//pWholeBoundary->m_pVertices.push_back(Boundary[VertexIndex]);
	}

	CheckChunk(pWholeBoundary, false);

	if(pWholeBoundary->m_Complex == true)
	{
		BreakChunk(pWholeBoundary);
	}


	ReOrder(pWholeBoundary);

	gIndices[0] = -3;
	int VertIndex = 1;

	for(int Index = 0; Index < pWholeBoundary->m_Triangles.size(); Index++)
	{
		gIndices[VertIndex++] = pWholeBoundary->m_Triangles[Index];
	}

	gIndices[VertIndex++] = -4;

	if(pWholeBoundary != NULL) delete pWholeBoundary;
	m_pVertices.clear();
}



void FillHole::ReOrder(HoleChunk* pWholeBoundary)
{
	vector<int> CwTris;											// clockwise triangles

	int NumberOfTris = pWholeBoundary->m_Triangles.size() / 3;
	int finalJoin = m_pVertices.size() - 1;


	for(int TriIndex = 0; TriIndex < NumberOfTris; TriIndex++)
	{
		int VertIndex1 = pWholeBoundary->m_Triangles[TriIndex*3];
		int VertIndex2 = pWholeBoundary->m_Triangles[TriIndex*3+1];
		int VertIndex3 = pWholeBoundary->m_Triangles[TriIndex*3+2];
		int Delta1 = VertIndex2 - VertIndex1;
		int Delta2 = VertIndex3 - VertIndex2;
		int Delta3 = VertIndex1 - VertIndex3;

		// triangle anti clockwise
		if(Delta1 == -1 || Delta2 == -1 || Delta3 == -1 || Delta1 == finalJoin || Delta2 == finalJoin || Delta3 == finalJoin)
		{
			CwTris.push_back(VertIndex1);
			CwTris.push_back(VertIndex3);
			CwTris.push_back(VertIndex2);
			pWholeBoundary->m_Triangles[TriIndex*3] = -1;
			pWholeBoundary->m_Triangles[TriIndex*3+1] = -1;
			pWholeBoundary->m_Triangles[TriIndex*3+2] = -1;
		}

		if(Delta1 == 1 || Delta2 == 1 || Delta3 == 1 || Delta1 == -finalJoin || Delta2 == -finalJoin || Delta3 == -finalJoin)
		{
			CwTris.push_back(VertIndex1);
			CwTris.push_back(VertIndex2);
			CwTris.push_back(VertIndex3);
			pWholeBoundary->m_Triangles[TriIndex*3] = -1;
			pWholeBoundary->m_Triangles[TriIndex*3+1] = -1;
			pWholeBoundary->m_Triangles[TriIndex*3+2] = -1;
		}
	}

	bool MoreToFlip = true;

	while(MoreToFlip == true)
	{
		MoreToFlip = false;

		for(int TriIndex = 0; TriIndex < NumberOfTris; TriIndex++)
		{
			int SrcVertIndex1 = pWholeBoundary->m_Triangles[TriIndex*3];

			if(SrcVertIndex1 != -1)
			{
				int SrcVertIndex2 = pWholeBoundary->m_Triangles[TriIndex*3+1];
				int SrcVertIndex3 = pWholeBoundary->m_Triangles[TriIndex*3+2];

				int NumberDestTris = CwTris.size() / 3;
				int action = 0;				// 0 do nothing,  1 copy tri,  2 flip and copy


				for(int ListIndex = 0; ListIndex < NumberDestTris; ListIndex++)
				{
					int ListVertIndex1 = CwTris[ListIndex*3];
					int ListVertIndex2 = CwTris[ListIndex*3+1];
					int ListVertIndex3 = CwTris[ListIndex*3+2];

					// source edge 1
					if(SrcVertIndex1 == ListVertIndex2 && SrcVertIndex2 == ListVertIndex1) action = 1;
					if(SrcVertIndex1 == ListVertIndex1 && SrcVertIndex2 == ListVertIndex2) action = 2;
					if(SrcVertIndex1 == ListVertIndex3 && SrcVertIndex2 == ListVertIndex2) action = 1;
					if(SrcVertIndex1 == ListVertIndex2 && SrcVertIndex2 == ListVertIndex3) action = 2;
					if(SrcVertIndex1 == ListVertIndex1 && SrcVertIndex2 == ListVertIndex3) action = 1;
					if(SrcVertIndex1 == ListVertIndex3 && SrcVertIndex2 == ListVertIndex1) action = 2;

					// source edge 2
					if(SrcVertIndex2 == ListVertIndex2 && SrcVertIndex3 == ListVertIndex1) action = 1;
					if(SrcVertIndex2 == ListVertIndex1 && SrcVertIndex3 == ListVertIndex2) action = 2;
					if(SrcVertIndex2 == ListVertIndex3 && SrcVertIndex3 == ListVertIndex2) action = 1;
					if(SrcVertIndex2 == ListVertIndex2 && SrcVertIndex3 == ListVertIndex3) action = 2;
					if(SrcVertIndex2 == ListVertIndex1 && SrcVertIndex3 == ListVertIndex3) action = 1;
					if(SrcVertIndex2 == ListVertIndex3 && SrcVertIndex3 == ListVertIndex1) action = 2;

					// source edge 3
					if(SrcVertIndex3 == ListVertIndex2 && SrcVertIndex1 == ListVertIndex1) action = 1;
					if(SrcVertIndex3 == ListVertIndex1 && SrcVertIndex1 == ListVertIndex2) action = 2;
					if(SrcVertIndex3 == ListVertIndex3 && SrcVertIndex1 == ListVertIndex2) action = 1;
					if(SrcVertIndex3 == ListVertIndex2 && SrcVertIndex1 == ListVertIndex3) action = 2;
					if(SrcVertIndex3 == ListVertIndex1 && SrcVertIndex1 == ListVertIndex3) action = 1;
					if(SrcVertIndex3 == ListVertIndex3 && SrcVertIndex1 == ListVertIndex1) action = 2;
				}

				if(action == 0) MoreToFlip = true;

				if(action == 1)
				{
					CwTris.push_back(SrcVertIndex1);
					CwTris.push_back(SrcVertIndex2);
					CwTris.push_back(SrcVertIndex3);
					pWholeBoundary->m_Triangles[TriIndex*3] = -1;
					pWholeBoundary->m_Triangles[TriIndex*3+1] = -1;
					pWholeBoundary->m_Triangles[TriIndex*3+2] = -1;
				}

				if(action == 2)
				{
					CwTris.push_back(SrcVertIndex1);
					CwTris.push_back(SrcVertIndex3);
					CwTris.push_back(SrcVertIndex2);
					pWholeBoundary->m_Triangles[TriIndex*3] = -1;
					pWholeBoundary->m_Triangles[TriIndex*3+1] = -1;
					pWholeBoundary->m_Triangles[TriIndex*3+2] = -1;
				}
			}
		}
	}

	pWholeBoundary->m_Triangles.clear();
	pWholeBoundary->m_Triangles.insert(pWholeBoundary->m_Triangles.begin(), CwTris.begin(), CwTris.end());

	CwTris.clear();
}





void FillHole::BreakChunk(HoleChunk* pWholeBoundary)
{
	vector<HoleChunk*> Chunks;

	Chunks.push_back(pWholeBoundary);


	while(Chunks[0]->m_Complex == true)
	{
		HoleChunk* pNewChunk = new HoleChunk();
		Chunks.push_back(pNewChunk);
		SplitChunks(Chunks[0], pNewChunk);
		CheckChunk(Chunks[0], false);
	}

	int size = Chunks.size();

	for(int Index = 1; Index < size; Index++)
	{
		pWholeBoundary->m_Triangles.insert(pWholeBoundary->m_Triangles.end(), Chunks[Index]->m_Triangles.begin(), Chunks[Index]->m_Triangles.end());
		delete Chunks[Index];
	}

	Chunks.clear();

	// find three verts in boundary that the tightest interior angle
	// walk away from that point, prev and next, check chunk for complexity over 4 verts..
	// when we hit complex, roll back one vertex, try adding verts from opposite end
	// add convex list to chunks

	// repeat on leftovers.
	// delaunay the lot.
	// add together
}




// find three verts in boundary that the tightest interior angle
// walk away from that point, prev and next, check chunk for complexity over 4 verts..
// when we hit complex, roll back one vertex, try adding verts from opposite end
// add convex list to chunks

// take out and triangulate planar regions first?

void FillHole::SplitChunks(HoleChunk* pSourceChunk, HoleChunk* pDestChunk)
{
	float vert1[3], vert2[3], vert3[3];
	double BestAngle = 1000.0f;
	list<int>::iterator BestVertIter;
	list<int>::iterator PrevVertIter;
	list<int>::iterator NextVertIter;


	list<int>::iterator VertIterator;
	for(VertIterator = pSourceChunk->m_pVertIndices.begin(); VertIterator != pSourceChunk->m_pVertIndices.end(); VertIterator++)
	{
		list<int>::iterator Vert2Iterator = GetNext(pSourceChunk->m_pVertIndices, VertIterator);
		list<int>::iterator Vert3Iterator = GetNext(pSourceChunk->m_pVertIndices, Vert2Iterator);

		RKVertex* pVertex1 = m_pVertices[*VertIterator];
		RKVertex* pVertex2 = m_pVertices[*Vert2Iterator];
		RKVertex* pVertex3 = m_pVertices[*Vert3Iterator];

		vert1[0] = pVertex1->X;   vert1[1] = pVertex1->Y;   vert1[2] = pVertex1->Z;
		vert2[0] = pVertex2->X;   vert2[1] = pVertex2->Y;   vert2[2] = pVertex2->Z;
		vert3[0] = pVertex3->X;   vert3[1] = pVertex3->Y;   vert3[2] = pVertex3->Z;

		double Angle = RKFace::VecAngle(vert1, vert2, vert3);
		if(Angle < BestAngle)
		{
			BestVertIter = Vert2Iterator;
			BestAngle = Angle;
		}
	}


	pDestChunk->m_Complex = false;
//	bool flipFlop = false;

	PrevVertIter = GetPrev(pSourceChunk->m_pVertIndices, BestVertIter);
	NextVertIter = GetNext(pSourceChunk->m_pVertIndices, BestVertIter);

	pDestChunk->m_pVertIndices.push_back(*BestVertIter);
	pDestChunk->m_pVertIndices.push_front(*PrevVertIter);
	pDestChunk->m_pVertIndices.push_back(*NextVertIter);
	pDestChunk->m_Triangles.push_back(*BestVertIter);
	pDestChunk->m_Triangles.push_back(*PrevVertIter);
	pDestChunk->m_Triangles.push_back(*NextVertIter);

	pSourceChunk->m_pVertIndices.remove(*BestVertIter);						// remove Best from source
	PrevVertIter = GetPrev(pSourceChunk->m_pVertIndices, PrevVertIter);
	NextVertIter = GetNext(pSourceChunk->m_pVertIndices, NextVertIter);		// get next


	while(pDestChunk->m_Complex == false)
	{
		HoleChunk LeftChunk;
		HoleChunk RightChunk;
		list<int>::iterator LeftSuccess;
		list<int>::iterator RightSuccess;
		int LeftStart, RightStart;

		LeftChunk.Copy(pDestChunk);
		RightChunk.Copy(pDestChunk);

		LeftSuccess = GetPrev(pSourceChunk->m_pVertIndices, NextVertIter);
		LeftChunk.m_pVertIndices.push_back(*NextVertIter);
		LeftStart = *NextVertIter;
		//NextVertIter = GetNext(pSourceChunk->m_pVertIndices, NextVertIter);

		RightSuccess = GetNext(pSourceChunk->m_pVertIndices, PrevVertIter);
		RightChunk.m_pVertIndices.push_front(*PrevVertIter);
		RightStart = *PrevVertIter;
		//PrevVertIter = GetPrev(pSourceChunk->m_pVertIndices, PrevVertIter);

		CheckChunk(&LeftChunk, true);
		CheckChunk(&RightChunk, true);

		bool WhichChunk = 0;			// 0 = left   1 == right

		if(LeftChunk.m_Complex == true && RightChunk.m_Complex == true) return;		// no more tris can go into this chunk
		if(LeftChunk.m_Complex == true && RightChunk.m_Complex == false) WhichChunk = 1;
		if(LeftChunk.m_Complex == false && RightChunk.m_Complex == false)
		{
			int TriIndex1, TriIndex2, TriIndex3;
			float SmallestInteriorLeft, SmallestInteriorRight;

			GetTri(&LeftChunk, *LeftSuccess, LeftStart, TriIndex1, TriIndex2, TriIndex3);
			SmallestInteriorLeft = GetSmallestIterior(TriIndex1, TriIndex2, TriIndex3);

			GetTri(&RightChunk, *RightSuccess, RightStart, TriIndex1, TriIndex2, TriIndex3);
			SmallestInteriorRight = GetSmallestIterior(TriIndex1, TriIndex2, TriIndex3);

			if(SmallestInteriorRight > SmallestInteriorLeft) WhichChunk = 1;

			// find left triangle, get smallest interior angle
			// find right triangle, get smallest interior angle

			// add direction that makes biggest triangle
		}


		if(WhichChunk == 0)
		{
			pDestChunk->Copy(&LeftChunk);
			pSourceChunk->m_pVertIndices.remove(*LeftSuccess);
			NextVertIter = GetNext(pSourceChunk->m_pVertIndices, NextVertIter);
		}
		else
		{
			pDestChunk->Copy(&RightChunk);
			pSourceChunk->m_pVertIndices.remove(*RightSuccess);
			PrevVertIter = GetPrev(pSourceChunk->m_pVertIndices, PrevVertIter);
		}

/*		flipFlop ^= 1;
		list<int>::iterator RemoveIfSuccess;
		list<int>::iterator RemoveIfFail;

		if(flipFlop)
		{
			RemoveIfSuccess = GetPrev(pSourceChunk->m_pVertIndices, NextVertIter);
			pDestChunk->m_pVertIndices.push_back(*NextVertIter);
			RemoveIfFail = NextVertIter;
			NextVertIter = GetNext(pSourceChunk->m_pVertIndices, NextVertIter);
		}
		else
		{
			RemoveIfSuccess = GetNext(pSourceChunk->m_pVertIndices, PrevVertIter);
			pDestChunk->m_pVertIndices.push_front(*PrevVertIter);
			RemoveIfFail = PrevVertIter;
			PrevVertIter = GetPrev(pSourceChunk->m_pVertIndices, PrevVertIter);
		}

		CheckChunk(pDestChunk);

		if(pDestChunk->m_Complex == true)
		{
			pDestChunk->m_pVertIndices.remove(*RemoveIfFail);				// remove previous from dest
			CheckChunk(pDestChunk);
			return;
		}
		else
		{
			pSourceChunk->m_pVertIndices.remove(*RemoveIfSuccess);			// remove RemoveIfSuccess from source
		}
*/

		if(pSourceChunk->m_pVertIndices.size() < 3) break;					// taken all the triangles
	}
}


void FillHole::GetTri(HoleChunk* pChunk, int EdgeIndex1, int EdgeIndex2, int& TriIndex1, int& TriIndex2, int& TriIndex3)
{
	int NumberOfTris = pChunk->m_Triangles.size() / 3;

	for(int TriIndex = 0; TriIndex < NumberOfTris; TriIndex++)
	{
		TriIndex1 = pChunk->m_Triangles[TriIndex*3];
		TriIndex2 = pChunk->m_Triangles[TriIndex*3+1];
		TriIndex3 = pChunk->m_Triangles[TriIndex*3+2];

		if(TriIndex1 == EdgeIndex1 && TriIndex2 == EdgeIndex2 || TriIndex1 == EdgeIndex2 && TriIndex2 == EdgeIndex1) return;
		if(TriIndex2 == EdgeIndex1 && TriIndex3 == EdgeIndex2 || TriIndex2 == EdgeIndex2 && TriIndex3 == EdgeIndex1) return;
		if(TriIndex3 == EdgeIndex1 && TriIndex1 == EdgeIndex2 || TriIndex3 == EdgeIndex2 && TriIndex1 == EdgeIndex1) return;
	}

	TriIndex1 = -1;
	TriIndex2 = -1;
	TriIndex3 = -1;
}



float FillHole::GetSmallestIterior(int TriIndex1, int TriIndex2, int TriIndex3)
{
	float Vert1[3], Vert2[3], Vert3[3];
	Vert1[0] = m_pVertices[TriIndex1]->X;	Vert1[1] = m_pVertices[TriIndex1]->Y; Vert1[2] = m_pVertices[TriIndex1]->Z;
	Vert2[0] = m_pVertices[TriIndex2]->X;	Vert2[1] = m_pVertices[TriIndex2]->Y; Vert2[2] = m_pVertices[TriIndex2]->Z;
	Vert3[0] = m_pVertices[TriIndex3]->X;	Vert3[1] = m_pVertices[TriIndex3]->Y; Vert3[2] = m_pVertices[TriIndex3]->Z;
	double Angle1 = RKFace::VecAngle(Vert3, Vert1, Vert2);
	double Angle2 = RKFace::VecAngle(Vert1, Vert2, Vert3);
	double Angle3 = (float)M_PI - Angle2 - Angle1;

	if(Angle1 < Angle2 && Angle1 < Angle3) return(Angle1);
	if(Angle2 < Angle1 && Angle2 < Angle3) return(Angle2);
	return(Angle3);
}



list<int>::iterator FillHole::GetNext(list<int> &Vertices, list<int>::iterator VertIterator)
{
	list<int>::iterator NextIter = VertIterator;
	NextIter++;

	if(NextIter == Vertices.end())
	{
		NextIter = Vertices.begin();
	}

	return(NextIter);
}



list<int>::iterator FillHole::GetPrev(list<int> &Vertices, list<int>::iterator VertIterator)
{
	list<int>::iterator PrevIter = VertIterator;

	if(PrevIter == Vertices.begin())
	{
		PrevIter = Vertices.end();					// --???
		PrevIter--;
	}
	else
	{
		PrevIter--;
	}

	return(PrevIter);
}



bool FillHole::CheckChunk(HoleChunk* pChunk, bool AngleCheck)
{
	double BestPlane[4];
	int numberVertices = pChunk->m_pVertIndices.size();
	double* weights = new double[numberVertices];
	int* ReMap = new int[numberVertices];
	gVertices = new double[numberVertices*3];
	gProjVertices = new double[numberVertices*3];
	gTripleProjVerts = new double[numberVertices*3];

	int VertIndex = 0;
	list<int>::iterator VertIterator;
	for(VertIterator = pChunk->m_pVertIndices.begin(); VertIterator != pChunk->m_pVertIndices.end(); VertIterator++)
	{
		int OrigVert = *VertIterator;
//		RKVertex* pVertex = *VertIterator;
		gVertices[VertIndex*3] = m_pVertices[OrigVert]->X;			//pVertex->X;
		gVertices[VertIndex*3+1] = m_pVertices[OrigVert]->Y;		//pVertex->Y;
		gVertices[VertIndex*3+2] = m_pVertices[OrigVert]->Z;		//pVertex->Z;
		weights[VertIndex] = 1.0f;
		ReMap[VertIndex] = OrigVert;
		VertIndex++;
	}

	BestFitPlane::getBestFitPlane(numberVertices, gVertices, 3*8, weights, 1*8, BestPlane);
	// rotate onto best plane

	projectPolygon( BestPlane, numberVertices );

	for(int Index = 0; Index < numberVertices; Index++)
	{
		gTripleProjVerts[Index*3] = gProjVertices[Index*2];
		gTripleProjVerts[Index*3+1] = gProjVertices[Index*2+1];
		gTripleProjVerts[Index*3+2] = 0.0;
	}

	bool ComplexBounds = CheckComplex(numberVertices, gTripleProjVerts);					// is the polygon intersecting itself on the best fit project plane
	pChunk->m_Complex = ComplexBounds;

//	if(numberVertices == 4) pChunk->m_Complex = true;			//  bodge to check the splitter works

	if(ComplexBounds == true)
	{
		delete[] gVertices;
		delete[] gProjVertices;
		delete[] gTripleProjVerts;
		delete[] ReMap;
		delete[] weights;
		return(false);			// going to have to split this up
	}


	Delaunay(numberVertices);
	pChunk->m_Triangles.clear();

	int Type = gIndices[0];
	int Index = 1;
	int VertIndex1, VertIndex2, VertIndex3;

	while(Type != -4)
	{
		VertIndex1 = gIndices[Index];
		VertIndex2 = gIndices[Index+1];
		VertIndex3 = gIndices[Index+2];
		Index += 3;

		int RealIndex1 = ReMap[VertIndex1];			//*V1Iter;		//pChunk->m_pVertIndices[VertIndex1];
		int RealIndex2 = ReMap[VertIndex2];			//*V2Iter;		//pChunk->m_pVertIndices[VertIndex2];
		int RealIndex3 = ReMap[VertIndex3];			//*V3Iter;		//pChunk->m_pVertIndices[VertIndex3];

		if(AngleCheck)
		{
			float Vert1[3], Vert2[3], Vert3[3];
			Vert1[0] = gProjVertices[VertIndex1*2];	Vert1[1] = gProjVertices[VertIndex1*2+1]; Vert1[2] = 0.0f;
			Vert2[0] = gProjVertices[VertIndex2*2]; Vert2[1] = gProjVertices[VertIndex2*2+1]; Vert2[2] = 0.0f;
			Vert3[0] = gProjVertices[VertIndex3*2]; Vert3[1] = gProjVertices[VertIndex3*2+1]; Vert3[2] = 0.0f;
			double Angle1_2D = RKFace::VecAngle(Vert3, Vert1, Vert2);
			double Angle2_2D = RKFace::VecAngle(Vert1, Vert2, Vert3);
			double Angle3_2D = (float)M_PI - Angle2_2D - Angle1_2D;


			Vert1[0] = m_pVertices[RealIndex1]->X;	Vert1[1] = m_pVertices[RealIndex1]->Y; Vert1[2] = m_pVertices[RealIndex1]->Z;
			Vert2[0] = m_pVertices[RealIndex2]->X;	Vert2[1] = m_pVertices[RealIndex2]->Y; Vert2[2] = m_pVertices[RealIndex2]->Z;
			Vert3[0] = m_pVertices[RealIndex3]->X;	Vert3[1] = m_pVertices[RealIndex3]->Y; Vert3[2] = m_pVertices[RealIndex3]->Z;
			double Angle1_3D = RKFace::VecAngle(Vert3, Vert1, Vert2);
			double Angle2_3D = RKFace::VecAngle(Vert1, Vert2, Vert3);
			double Angle3_3D = (float)M_PI - Angle2_3D - Angle1_3D;

			if(Angle1_2D < 0.1f) pChunk->m_Complex = true;			// very elongated polygon, time to rethink
			if(Angle2_2D < 0.1f) pChunk->m_Complex = true;
			if(Angle3_2D < 0.1f) pChunk->m_Complex = true;
			if(Angle1_3D < 0.1f) pChunk->m_Complex = true;
			if(Angle2_3D < 0.1f) pChunk->m_Complex = true;
			if(Angle3_3D < 0.1f) pChunk->m_Complex = true;
		}

		pChunk->m_Triangles.push_back(RealIndex1);		//VertIndex1);
		pChunk->m_Triangles.push_back(RealIndex2);		//VertIndex2);
		pChunk->m_Triangles.push_back(RealIndex3);		//VertIndex3);
		Type = gIndices[Index];
	}


	delete[] ReMap;
	delete[] weights;
	delete[] gVertices;
	delete[] gProjVertices;
	delete[] gTripleProjVerts;

	// pass plane into GLUTess
	// see if we have a complex polygon
	// see if we have any zero area triangles
	// see if all triangles wind the same way

	return(true);
}





#ifdef _WINDOWS
void ftglError( GLenum errCode, int* pIndices )			//FTMesh* mesh)
#else
void ftglError( GLenum errCode, int* pIndices )
#endif
{
	gComplex = true;
}


#ifdef _WINDOWS
void ftglVertex( void* data, int* pIndices )			//FTMesh* mesh)
#else
void ftglVertex( void* data, int* pIndices )
#endif
{
	size_t addr1 = (size_t)gVertices;
	size_t addr2 = (size_t)data;

	size_t index = (size_t)data - (size_t)gTripleProjVerts;
	index /= 24;

	pIndices[gVertIndex++] = (int)index;
}


#ifdef _WINDOWS
void ftglCombine( GLdouble coords[3], void* vertex_data[4], GLfloat weight[4], void** outData, int* pIndices )			//FTMesh* mesh)
#else
void ftglCombine( GLdouble coords[3], void* vertex_data[4], GLfloat weight[4], void** outData, int* pIndices )
#endif
{
	gComplex = true;
}


#ifdef _WINDOWS
void ftglBegin( GLenum type, int* pIndices )			//FTMesh* mesh)
#else
void ftglBegin( GLenum type, int* pIndices )
#endif
{
	if(type == GL_TRIANGLE_FAN)
	{
		pIndices[gVertIndex++] = -1;
	}

	if(type == GL_TRIANGLE_STRIP)
	{
		pIndices[gVertIndex++] = -2;
	}

	if(type == GL_TRIANGLES)
	{
		pIndices[gVertIndex++] = -3;
	}
}

#ifdef _WINDOWS
void ftglEnd( int* pIndices )			//FTMesh* mesh)
#else
void ftglEnd( int* pIndices )
#endif
{
	int a = 0;
	a++;
}


#ifdef _WINDOWS
void EdgeFlagDataCB(GLboolean flag, void* lpContext)
{
}
#else
void EdgeFlagDataCB(GLboolean flag, void* lpContext)
{
}
#endif



#if defined (_WIN32)
typedef GLvoid (*_GLUfuncptr)( );
#endif


bool FillHole::CheckComplex(int numberOfVertices, double* pVerts)
{
	gComplex = false;
	gVertIndex = 0;

    GLUtesselator* tobj = gluNewTess();

#ifdef MAC_PLUGIN
	gluTessCallback( tobj, GLU_TESS_BEGIN_DATA,     (GLvoid(*)()) ftglBegin);
    gluTessCallback( tobj, GLU_TESS_VERTEX_DATA,    (GLvoid(*)()) ftglVertex);
    gluTessCallback( tobj, GLU_TESS_COMBINE_DATA,   (GLvoid(*)()) ftglCombine);
    gluTessCallback( tobj, GLU_TESS_END_DATA,       (GLvoid(*)()) ftglEnd);
    gluTessCallback( tobj, GLU_TESS_ERROR_DATA,     (GLvoid(*)()) ftglError);
	gluTessCallback( tobj, GLU_TESS_EDGE_FLAG_DATA,	(GLvoid(*)()) EdgeFlagDataCB);
#else
    gluTessCallback( tobj, GLU_TESS_BEGIN_DATA,     (_GLUfuncptr) ftglBegin);
    gluTessCallback( tobj, GLU_TESS_VERTEX_DATA,    (_GLUfuncptr) ftglVertex);
    gluTessCallback( tobj, GLU_TESS_COMBINE_DATA,   (_GLUfuncptr) ftglCombine);
    gluTessCallback( tobj, GLU_TESS_END_DATA,       (_GLUfuncptr) ftglEnd);
    gluTessCallback( tobj, GLU_TESS_ERROR_DATA,     (_GLUfuncptr) ftglError);
	gluTessCallback( tobj, GLU_TESS_EDGE_FLAG_DATA,	(_GLUfuncptr) EdgeFlagDataCB);
#endif
	/*
    if( contourFlag & ft_outline_even_odd_fill) // ft_outline_reverse_fill
    {
        gluTessProperty( tobj, GLU_TESS_WINDING_RULE, GLU_TESS_WINDING_ODD);
    }
    else
    {
        gluTessProperty( tobj, GLU_TESS_WINDING_RULE, GLU_TESS_WINDING_NONZERO);
    }
*/


    gluTessProperty( tobj, GLU_TESS_TOLERANCE, 0);
 //   gluTessNormal( tobj, BestPlane[0], BestPlane[1], BestPlane[2]);
    gluTessBeginPolygon( tobj, gIndices );
	gluTessBeginContour( tobj);

	for(int vertIndex = 0; vertIndex < numberOfVertices; vertIndex++)
	{
		double* vert = &pVerts[vertIndex*3];
        gluTessVertex( tobj, vert, vert);
	}

	gluTessEndContour( tobj);
    gluTessEndPolygon( tobj);
    gluDeleteTess(tobj);

	gIndices[gVertIndex++] = -4;

	return(gComplex);
}



#define Dot(u,v)	(u[0]*v[0] + u[1]*v[1] + u[2]*v[2])

static void Normalize( double v[3] )
{
	double len = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];

//  assert( len > 0 );
	len = sqrt( len );
	v[0] /= len;
	v[1] /= len;
	v[2] /= len;
}


#undef	ABS
#define ABS(x)	((x) < 0 ? -(x) : (x))

static int LongAxis( double v[3] )
{
  int i = 0;

  if( ABS(v[1]) > ABS(v[0]) ) { i = 1; }
  if( ABS(v[2]) > ABS(v[i]) ) { i = 2; }
  return i;
}

//#define S_UNIT_X	0.50941539564955385	/* Pre-normalized */
//#define S_UNIT_Y	0.86052074622010633
#define S_UNIT_X	1.0
#define S_UNIT_Y	0.0

void FillHole::projectPolygon( double BestPlane[4], int NumberOfVertices )
{
//  GLUvertex *v, *vHead = &tess->mesh->vHead;
	double norm[3];
	double sUnit[3];
	double tUnit[3];
	int i;		// computedNormal = FALSE;

	norm[0] = BestPlane[0];
	norm[1] = BestPlane[1];
	norm[2] = BestPlane[2];

//	sUnit = tess->sUnit;
//	tUnit = tess->tUnit;
	i = LongAxis( norm );

//#if defined(FOR_TRITE_TEST_PROGRAM) || defined(TRUE_PROJECT)
	/* Choose the initial sUnit vector to be approximately perpendicular
	 * to the normal.
	 */
	Normalize( norm );			// does this normalise work?

	sUnit[i] = 0;
	sUnit[(i+1)%3] = S_UNIT_X;
	sUnit[(i+2)%3] = S_UNIT_Y;

	/* Now make it exactly perpendicular */
	double w = Dot( sUnit, norm );
	sUnit[0] -= w * norm[0];
	sUnit[1] -= w * norm[1];
	sUnit[2] -= w * norm[2];
	Normalize( sUnit );

	/* Choose tUnit so that (sUnit,tUnit,norm) form a right-handed frame */
	tUnit[0] = norm[1]*sUnit[2] - norm[2]*sUnit[1];
	tUnit[1] = norm[2]*sUnit[0] - norm[0]*sUnit[2];
	tUnit[2] = norm[0]*sUnit[1] - norm[1]*sUnit[0];
	Normalize( tUnit );
/*#else
	// Project perpendicular to a coordinate axis -- better numerically
	sUnit[i] = 0;
	sUnit[(i+1)%3] = S_UNIT_X;
	sUnit[(i+2)%3] = S_UNIT_Y;

	tUnit[i] = 0;
	tUnit[(i+1)%3] = (norm[i] > 0) ? -S_UNIT_Y : S_UNIT_Y;
	tUnit[(i+2)%3] = (norm[i] > 0) ? S_UNIT_X : -S_UNIT_X;
#endif
*/
	/* Project the vertices onto the sweep plane */
	double vert[3];
	for(int Index = 0; Index < NumberOfVertices; Index++)
	{
		vert[0] = gVertices[Index*3];
		vert[1] = gVertices[Index*3+1];
		vert[2] = gVertices[Index*3+2];

		gProjVertices[Index*2] = Dot( vert, sUnit );
		gProjVertices[Index*2+1] = Dot( vert, tUnit );
//		v->s = Dot( v->coords, sUnit );
//		v->t = Dot( v->coords, tUnit );
	}

//	if( computedNormal )
//	{
//		CheckOrientation( tess );
//	}
}






void FillHole::Delaunay(int NumberOfVertices)
{
	int NumberOfTriangles = NumberOfVertices-2;

	// rotate vertices to best plane

//	DelaunayVert* pVerts = new DelaunayVert[NumberOfVertices];
	DelaunayFace* pFaces = new DelaunayFace[NumberOfTriangles];
	vector<DelaunayEdge*> pEdges;

	for(int FaceIndex = 0; FaceIndex < NumberOfTriangles; FaceIndex++)
	{
		pFaces[FaceIndex].m_Edge1Index = EdgeIndex(pEdges, gIndices[FaceIndex*3+1], gIndices[FaceIndex*3+2], FaceIndex);
		pFaces[FaceIndex].m_Edge2Index = EdgeIndex(pEdges, gIndices[FaceIndex*3+2], gIndices[FaceIndex*3+3], FaceIndex);
		pFaces[FaceIndex].m_Edge3Index = EdgeIndex(pEdges, gIndices[FaceIndex*3+3], gIndices[FaceIndex*3+1], FaceIndex);
	}

	bool DidSwap = true;
	int bail = 100;

	while(DidSwap == true)
	{
		DidSwap = false;
		bail--;
		if(bail == 0) goto arse;

		for(int TriIndex = 0; TriIndex < NumberOfTriangles; TriIndex++)
		{
			int Edge1 = pFaces[TriIndex].m_Edge1Index;
			int Edge2 = pFaces[TriIndex].m_Edge2Index;
			int Edge3 = pFaces[TriIndex].m_Edge3Index;

			DidSwap |= DelaunayFlip(pFaces, pEdges, Edge1);
//			if(DidSwap == true) bail--;
//			if(bail == 0) goto arse;

			DidSwap |= DelaunayFlip(pFaces, pEdges, Edge2);
//			if(DidSwap == true) bail--;
//			if(bail == 0) goto arse;

			DidSwap |= DelaunayFlip(pFaces, pEdges, Edge3);
//			if(DidSwap == true) bail--;
//			if(bail == 0) goto arse;
		}
	}

arse:
	for(int FaceIndex = 0; FaceIndex < NumberOfTriangles; FaceIndex++)
	{
		int VertIndex1 = pEdges[pFaces[FaceIndex].m_Edge1Index]->m_Vert1Index;
		int VertIndex2 = pEdges[pFaces[FaceIndex].m_Edge1Index]->m_Vert2Index;
		int VertIndex3 = -1;

		if(pEdges[pFaces[FaceIndex].m_Edge2Index]->m_Vert1Index != VertIndex1 && pEdges[pFaces[FaceIndex].m_Edge2Index]->m_Vert1Index != VertIndex2)
		{
			VertIndex3 = pEdges[pFaces[FaceIndex].m_Edge2Index]->m_Vert1Index;
		}
		else
		{
			VertIndex3 = pEdges[pFaces[FaceIndex].m_Edge2Index]->m_Vert2Index;
		}

		gIndices[FaceIndex*3+1] = VertIndex1;			//pEdges[pFaces[FaceIndex].m_Edge1Index]->m_Vert1Index;
		gIndices[FaceIndex*3+2] = VertIndex2;			//pEdges[pFaces[FaceIndex].m_Edge2Index]->m_Vert1Index;
		gIndices[FaceIndex*3+3] = VertIndex3;			//pEdges[pFaces[FaceIndex].m_Edge3Index]->m_Vert1Index;

//		pFaces[FaceIndex].m_Edge1Index = EdgeIndex(pEdges, gIndices[FaceIndex*3], gIndices[FaceIndex*3+1], FaceIndex);
//		pFaces[FaceIndex].m_Edge2Index = EdgeIndex(pEdges, gIndices[FaceIndex*3+1], gIndices[FaceIndex*3+2], FaceIndex);
//		pFaces[FaceIndex].m_Edge3Index = EdgeIndex(pEdges, gIndices[FaceIndex*3+2], gIndices[FaceIndex*3], FaceIndex);
	}

	delete[] pFaces;
	for(int Index = 0; Index < pEdges.size(); Index++)
	{
		delete pEdges[Index];
	}

	pEdges.clear();
}



bool FillHole::DelaunayFlip(DelaunayFace* pFaces, vector<DelaunayEdge*> &pEdges, int FlipEdgeIndex)
{
	int VertIndex1 = 0;
	int VertIndex2 = 0;
	int VertIndex3 = 0;
	int VertIndex4 = 0;			// clockwise
	int EdgeIndex1 = 0;
	int EdgeIndex2 = 0;
	int EdgeIndex3 = 0;
	int EdgeIndex4 = 0;
	int FaceIndex1 = 0;
	int FaceIndex2 = 0;

	FaceIndex1 = pEdges[FlipEdgeIndex]->m_Face1Index;
	FaceIndex2 = pEdges[FlipEdgeIndex]->m_Face2Index;

	if(FaceIndex2 == -1) return(false);				// can't flip boundary edge

	VertIndex1 = pEdges[FlipEdgeIndex]->m_Vert1Index;
	VertIndex3 = pEdges[FlipEdgeIndex]->m_Vert2Index;
//	int OtherEdgeIndex = -1;

	if(pFaces[FaceIndex1].m_Edge1Index == FlipEdgeIndex) { EdgeIndex1 = pFaces[FaceIndex1].m_Edge2Index; EdgeIndex2 = pFaces[FaceIndex1].m_Edge3Index; }
	if(pFaces[FaceIndex1].m_Edge2Index == FlipEdgeIndex) { EdgeIndex1 = pFaces[FaceIndex1].m_Edge3Index; EdgeIndex2 = pFaces[FaceIndex1].m_Edge1Index; }
	if(pFaces[FaceIndex1].m_Edge3Index == FlipEdgeIndex) { EdgeIndex1 = pFaces[FaceIndex1].m_Edge1Index; EdgeIndex2 = pFaces[FaceIndex1].m_Edge2Index; }
	if(pEdges[EdgeIndex1]->m_Vert1Index != VertIndex1 && pEdges[EdgeIndex1]->m_Vert1Index != VertIndex3) VertIndex2 = pEdges[EdgeIndex1]->m_Vert1Index;
	if(pEdges[EdgeIndex1]->m_Vert2Index != VertIndex1 && pEdges[EdgeIndex1]->m_Vert2Index != VertIndex3) VertIndex2 = pEdges[EdgeIndex1]->m_Vert2Index;

	if(pFaces[FaceIndex2].m_Edge1Index == FlipEdgeIndex) { EdgeIndex3 = pFaces[FaceIndex2].m_Edge2Index; EdgeIndex4 = pFaces[FaceIndex2].m_Edge3Index; }
	if(pFaces[FaceIndex2].m_Edge2Index == FlipEdgeIndex) { EdgeIndex3 = pFaces[FaceIndex2].m_Edge3Index; EdgeIndex4 = pFaces[FaceIndex2].m_Edge1Index; }
	if(pFaces[FaceIndex2].m_Edge3Index == FlipEdgeIndex) { EdgeIndex3 = pFaces[FaceIndex2].m_Edge1Index; EdgeIndex4 = pFaces[FaceIndex2].m_Edge2Index; }
	if(pEdges[EdgeIndex3]->m_Vert1Index != VertIndex1 && pEdges[EdgeIndex3]->m_Vert1Index != VertIndex3) VertIndex4 = pEdges[EdgeIndex3]->m_Vert1Index;
	if(pEdges[EdgeIndex3]->m_Vert2Index != VertIndex1 && pEdges[EdgeIndex3]->m_Vert2Index != VertIndex3) VertIndex4 = pEdges[EdgeIndex3]->m_Vert2Index;

	// test for concaveness, return if concave
	bool Concave = TestConcave(VertIndex1, VertIndex2, VertIndex3, VertIndex4);
	if(Concave == true) return(false);

	bool PointsInCircle = false;

	PointsInCircle |= InCircle(VertIndex4, VertIndex1, VertIndex2, VertIndex3);
	PointsInCircle |= InCircle(VertIndex2, VertIndex1, VertIndex3, VertIndex4);

	if(PointsInCircle == true)
	{
		pEdges[FlipEdgeIndex]->m_Vert1Index = VertIndex2;
		pEdges[FlipEdgeIndex]->m_Vert2Index = VertIndex4;

		pFaces[FaceIndex1].m_Edge1Index = FlipEdgeIndex;
		pFaces[FaceIndex1].m_Edge2Index = EdgeIndex4;
		pFaces[FaceIndex1].m_Edge3Index = EdgeIndex1;

		pFaces[FaceIndex2].m_Edge1Index = FlipEdgeIndex;
		pFaces[FaceIndex2].m_Edge2Index = EdgeIndex2;
		pFaces[FaceIndex2].m_Edge3Index = EdgeIndex3;

		SetOpposites(pEdges, pFaces, FaceIndex1, FaceIndex2);

		pEdges[FlipEdgeIndex]->m_Face1Index = FaceIndex1;
		pEdges[FlipEdgeIndex]->m_Face2Index = FaceIndex2;

		return(true);					// we've flipped!
	}

	return(false);
}




void FillHole::SetOpposites(vector<DelaunayEdge*> &pEdges, DelaunayFace* pFaces, int FaceIndex1, int FaceIndex2)
{
	int F1E1 = pFaces[FaceIndex1].m_Edge1Index;
	int F1E2 = pFaces[FaceIndex1].m_Edge2Index;
	int F1E3 = pFaces[FaceIndex1].m_Edge3Index;

	int F2E1 = pFaces[FaceIndex2].m_Edge1Index;
	int F2E2 = pFaces[FaceIndex2].m_Edge2Index;
	int F2E3 = pFaces[FaceIndex2].m_Edge3Index;

	if(pEdges[F1E1]->m_Face1Index == FaceIndex1 || pEdges[F1E1]->m_Face1Index == FaceIndex2) pEdges[F1E1]->m_Face1Index = FaceIndex1;
	if(pEdges[F1E1]->m_Face2Index == FaceIndex1 || pEdges[F1E1]->m_Face2Index == FaceIndex2) pEdges[F1E1]->m_Face2Index = FaceIndex1;
	if(pEdges[F1E2]->m_Face1Index == FaceIndex1 || pEdges[F1E2]->m_Face1Index == FaceIndex2) pEdges[F1E2]->m_Face1Index = FaceIndex1;
	if(pEdges[F1E2]->m_Face2Index == FaceIndex1 || pEdges[F1E2]->m_Face2Index == FaceIndex2) pEdges[F1E2]->m_Face2Index = FaceIndex1;
	if(pEdges[F1E3]->m_Face1Index == FaceIndex1 || pEdges[F1E3]->m_Face1Index == FaceIndex2) pEdges[F1E3]->m_Face1Index = FaceIndex1;
	if(pEdges[F1E3]->m_Face2Index == FaceIndex1 || pEdges[F1E3]->m_Face2Index == FaceIndex2) pEdges[F1E3]->m_Face2Index = FaceIndex1;

	if(pEdges[F2E1]->m_Face1Index == FaceIndex1 || pEdges[F2E1]->m_Face1Index == FaceIndex2) pEdges[F2E1]->m_Face1Index = FaceIndex2;
	if(pEdges[F2E1]->m_Face2Index == FaceIndex1 || pEdges[F2E1]->m_Face2Index == FaceIndex2) pEdges[F2E1]->m_Face2Index = FaceIndex2;
	if(pEdges[F2E2]->m_Face1Index == FaceIndex1 || pEdges[F2E2]->m_Face1Index == FaceIndex2) pEdges[F2E2]->m_Face1Index = FaceIndex2;
	if(pEdges[F2E2]->m_Face2Index == FaceIndex1 || pEdges[F2E2]->m_Face2Index == FaceIndex2) pEdges[F2E2]->m_Face2Index = FaceIndex2;
	if(pEdges[F2E3]->m_Face1Index == FaceIndex1 || pEdges[F2E3]->m_Face1Index == FaceIndex2) pEdges[F2E3]->m_Face1Index = FaceIndex2;
	if(pEdges[F2E3]->m_Face2Index == FaceIndex1 || pEdges[F2E3]->m_Face2Index == FaceIndex2) pEdges[F2E3]->m_Face2Index = FaceIndex2;
}




// I'm sure this could be written a load better, unfortunatley I'm past giving a shit

bool FillHole::TestConcave(int VertIndex1, int VertIndex2, int VertIndex3, int VertIndex4)
{
	int NegCount = 0;
	int PosCount = 0;
	float VecX1, VecY1, VecX2, VecY2;
	float Cross;

	VecX1 = gProjVertices[VertIndex1*2] - gProjVertices[VertIndex2*2];
	VecY1 = gProjVertices[VertIndex1*2+1] - gProjVertices[VertIndex2*2+1];
	VecX2 = gProjVertices[VertIndex3*2] - gProjVertices[VertIndex2*2];
	VecY2 = gProjVertices[VertIndex3*2+1] - gProjVertices[VertIndex2*2+1];

	Cross = (VecX1 * VecY2) - (VecY1 * VecX2);
	if(Cross < 0) NegCount++;
	else
		PosCount++;


	VecX1 = gProjVertices[VertIndex2*2] - gProjVertices[VertIndex3*2];
	VecY1 = gProjVertices[VertIndex2*2+1] - gProjVertices[VertIndex3*2+1];
	VecX2 = gProjVertices[VertIndex4*2] - gProjVertices[VertIndex3*2];
	VecY2 = gProjVertices[VertIndex4*2+1] - gProjVertices[VertIndex3*2+1];

	Cross = (VecX1 * VecY2) - (VecY1 * VecX2);
	if(Cross < 0) NegCount++;
	else
		PosCount++;


	VecX1 = gProjVertices[VertIndex3*2] - gProjVertices[VertIndex4*2];
	VecY1 = gProjVertices[VertIndex3*2+1] - gProjVertices[VertIndex4*2+1];
	VecX2 = gProjVertices[VertIndex1*2] - gProjVertices[VertIndex4*2];
	VecY2 = gProjVertices[VertIndex1*2+1] - gProjVertices[VertIndex4*2+1];

	Cross = (VecX1 * VecY2) - (VecY1 * VecX2);
	if(Cross < 0) NegCount++;
	else
		PosCount++;


	VecX1 = gProjVertices[VertIndex4*2] - gProjVertices[VertIndex1*2];
	VecY1 = gProjVertices[VertIndex4*2+1] - gProjVertices[VertIndex1*2+1];
	VecX2 = gProjVertices[VertIndex2*2] - gProjVertices[VertIndex1*2];
	VecY2 = gProjVertices[VertIndex2*2+1] - gProjVertices[VertIndex1*2+1];

	Cross = (VecX1 * VecY2) - (VecY1 * VecX2);
	if(Cross < 0) NegCount++;
	else
		PosCount++;


	if(NegCount == 0 || PosCount == 0) return(false);

	return(true);
}




int FillHole::EdgeIndex(vector<DelaunayEdge*> &pEdges, int VertIndex1, int VertIndex2, int FaceIndex)
{
	for(int Index = 0; Index < pEdges.size(); Index++)
	{
		if(pEdges[Index]->m_Vert1Index == VertIndex1 && pEdges[Index]->m_Vert2Index == VertIndex2
			|| pEdges[Index]->m_Vert1Index == VertIndex2 && pEdges[Index]->m_Vert2Index == VertIndex1)
		{
			pEdges[Index]->m_Face2Index = FaceIndex;
			return(Index);
		}
	}

	DelaunayEdge* pNewEdge = new DelaunayEdge();
	pNewEdge->m_Vert1Index = VertIndex1;
	pNewEdge->m_Vert2Index = VertIndex2;
	pNewEdge->m_Face1Index = FaceIndex;
	pNewEdge->m_Face2Index = -1;
	pEdges.push_back(pNewEdge);

	return(pEdges.size()-1);
}



//Return TRUE if the point (xp,yp) lies inside the circumcircle
//made up by points (x1,y1) (x2,y2) (x3,y3)
//NOTE: A point on the edge is inside the circumcircle

bool FillHole::InCircle(int TestPointIndex, int P1Index, int P2Index, int P3Index)
{
	float Epsilon = FLT_EPSILON;
	float pX = gProjVertices[TestPointIndex*2];
	float pY = gProjVertices[TestPointIndex*2+1];
	float p1X = gProjVertices[P1Index*2];
	float p1Y = gProjVertices[P1Index*2+1];
	float p2X = gProjVertices[P2Index*2];
	float p2Y = gProjVertices[P2Index*2+1];
	float p3X = gProjVertices[P3Index*2];
	float p3Y = gProjVertices[P3Index*2+1];

	if (fabs(p1Y - p2Y ) < Epsilon && fabs(p2Y - p3Y) < Epsilon)
	{
		//INCIRCUM - F - Points are coincident !!
		return false;
	}

	float m1, m2, mx1, mx2, my1, my2, xc, yc;

	if(fabs(p2Y - p1Y) < Epsilon)
	{
		m2 = -(p3X - p2X) / (p3Y - p2Y);
		mx2 = (p2X + p3X) * 0.5;
		my2 = (p2Y + p3Y) * 0.5;

		//Calculate CircumCircle center (xc,yc)
		xc = (p2X + p1X) * 0.5;
		yc = m2 * (xc - mx2) + my2;
	}
	else if (fabs(p3Y - p2Y) < Epsilon)
	{
		m1 = -(p2X - p1X) / (p2Y - p1Y);
		mx1 = (p1X + p2X) * 0.5;
		my1 = (p1Y + p2Y) * 0.5;

		//Calculate CircumCircle center (xc,yc)

		xc = (p3X + p2X) * 0.5;
		yc = m1 * (xc - mx1) + my1;
	}
	else
	{
		m1 = -(p2X - p1X) / (p2Y - p1Y);
		m2 = -(p3X - p2X) / (p3Y - p2Y);
		mx1 = (p1X + p2X) * 0.5;

		mx2 = (p2X + p3X) * 0.5;
		my1 = (p1Y + p2Y) * 0.5;
		my2 = (p2Y + p3Y) * 0.5;

		//Calculate CircumCircle center (xc,yc)
		xc = (m1 * mx1 - m2 * mx2 + my2 - my1) / (m1 - m2);
		yc = m1 * (xc - mx1) + my1;
	}

	float dx = p2X - xc;
	float dy = p2Y - yc;
	float rsqr = dx * dx + dy * dy;

	//double r = Math.Sqrt(rsqr); //Circumcircle radius
	dx = pX - xc;
	dy = pY - yc;

	float drsqr = dx * dx + dy * dy;
	return ( drsqr < rsqr );
}










/********************************
*	Alternative Triangulator	*
********************************/
/*
CTriangulator::CTriangulator()
{
    m_points.clear();
}

///
///     Default destructor
///
CTriangulator::~CTriangulator()
{
    m_points.clear();
}

///
///     Adds a point to the contour
///
void CTriangulator::add(SVec3D* p)
{
    m_points.push_back(p);			//append(p);
}

///
///     Triangulates the contour
///
void CTriangulator::triangulate(vector<int> &indices)
{
    _process(indices);
}

///
///     Processes the triangulation
///
void CTriangulator::_process(vector<int> &indices)
{
    const int n = m_points.size();

    if (n < 3) return;

    int* V = new int[n];

	if (0.0f < _area())
    {
        for (int v = 0; v < n; v++) V[v] = v;
    }
    else
    {
        for (int v = 0; v < n; v++) V[v] = (n - 1) - v;
    }


	int nv = n;
    int count = 2 * nv;

	for (int m = 0, v = nv - 1; nv > 2;)
    {
        if (0 >= (count--)) return;

        int u = v;
        if (nv <= u)
		{
            u = 0;
		}

        v = u + 1;

		if (nv <= v)
		{
            v = 0;
		}

        int w = v + 1;
        if (nv <= w)
		{
            w = 0;
		}

        if (_snip(u, v, w, nv, V))
        {
            int a, b, c, s, t;
            a = V[u];
            b = V[v];
            c = V[w];
            indices.push_back(a);			//append(a);
            indices.push_back(b);			//->append(b);
            indices.push_back(c);			//->append(c);
            m++;

            for (s = v, t = v + 1; t < nv; s++, t++)
			{
                V[s] = V[t];
			}
            nv--;
            count = 2 * nv;
        }
    }

    delete V;
}

///
///     Returns the area of the contour
///

float CTriangulator::_area()
{
    int n = m_points.size();			//->getSize();
    float A = 0.0f;

    for (int p = n - 1, q = 0; q < n; p = q++)
    {
        SVec3D* pval = m_points[p];		//.get(p);
        SVec3D* qval = m_points[q];		//.get(q);
        A += pval->x * qval->y - qval->x * pval->y;
    }
    return(A * 0.5f);
}


bool CTriangulator::_snip(int u, int v, int w, int n, int *V)
{
	int p;
    SVec3D* A = m_points[V[u]];			//.get(V[u]);
    SVec3D* B = m_points[V[v]];		//->get(V[v]);
    SVec3D* C = m_points[V[w]];		//->get(V[w]);

    if (FLT_EPSILON > (((B->x - A->x) * (C->y - A->y)) - ((B->y - A->y) * (C->x - A->x))) )
	{
        return false;
	}

    for (p = 0; p < n; p++)
    {
        if ((p == u) || (p == v) || (p == w))
		{
            continue;
		}

        SVec3D* P = m_points[V[p]];		//->get(V[p]);

		if (_insideTriangle(A, B, C, P))
		{
            return false;
		}
    }
    return true;
}


///
///     Tests if a point is inside the given triangle
///
bool CTriangulator::_insideTriangle(SVec3D* A, SVec3D* B, SVec3D* C, SVec3D* P)
{
    float ax, ay, bx, by, cx, cy, apx, apy, bpx, bpy, cpx, cpy;
    float cCROSSap, bCROSScp, aCROSSbp;

    ax = C->x - B->x;  ay = C->y - B->y;
    bx = A->x - C->x;  by = A->y - C->y;
    cx = B->x - A->x;  cy = B->y - A->y;
    apx = P->x - A->x;  apy = P->y - A->y;
    bpx = P->x - B->x;  bpy = P->y - B->y;
    cpx = P->x - C->x;  cpy = P->y - C->y;

    aCROSSbp = ax * bpy - ay * bpx;
    cCROSSap = cx * apy - cy * apx;
    bCROSScp = bx * cpy - by * cpx;

    return ((aCROSSbp >= 0.0f) && (bCROSScp >= 0.0f) && (cCROSSap >= 0.0f));
}

*/
