#ifndef ISLAND_INCLUDED
#define ISLAND_INCLUDED

#include <vector>
#include <list>
#include <set>
#include <stdint.h>


#include "cpms.h"
#include "lscm.h"

#include "FindPins.h"
#include "PolyFiller.h"

#include "MinimiseStretch.h"
#include "Primitives.h"



using namespace std;



// for the convex hull code
class RKPoint
{
public:
	float x;
	float y;
};





class EdgeHash
{
public:
	EdgeHash() { EdgeList.clear(); }
	~EdgeHash() { EdgeList.clear(); }

	RKEdge* AddEdge(int VertIndex1, int VertIndex2, float EdgeLen, RKPolygon* pFace, vector<RKEdge*> &m_Edges);

	vector<RKEdge*> EdgeList;
};


class DivisionPlane
{
public:
	double m_A, m_B, m_C, m_D;
};



class WalkBlock
{
public:
	RKPolygon* m_pFace1;				// the face in the source Island
	RKPolygon* m_pFace2;				// the face in the dest Island
	int m_Edge1Index;					// The matching edge index for face 1
	int m_Edge2Index;					// The matching edge index for face 2
	float m_DegreeOfDiff;
};




class Island
{
public:
	Island()
	{
		UMin = VMin = UMax = VMax = 0.0f;
		m_PinnedCount = 0;

		m_pEdgeHash = NULL;
		m_pProfilePixels = NULL;
		m_pProfileBits = NULL;

//		m_pUpperProfile = NULL;
//		m_pLowerProfile = NULL;
//		m_pScaledUpper = NULL;
//		m_pScaledLower = NULL;
		m_pSameAs = NULL;
		m_pSymmetricalWith = NULL;
		m_Symmetrical = false;
		m_MidPoint = 0.0f;


		m_Vertices.clear();
		m_Faces.clear();
		m_Edges.clear();
	}

	~Island()
	{
		UMin = VMin = UMax = VMax = 0.0f;

		m_Vertices.clear();

		if(m_pEdgeHash != NULL) delete[] m_pEdgeHash;
		if(m_pProfilePixels != NULL) delete[] m_pProfilePixels;
		if(m_pProfileBits != NULL) delete[] m_pProfileBits;

//		if(m_pUpperProfile != NULL) delete[] m_pUpperProfile;
//		if(m_pLowerProfile != NULL) delete[] m_pLowerProfile;
//		if(m_pScaledUpper != NULL) delete[] m_pScaledUpper;
//		if(m_pScaledLower != NULL) delete[] m_pScaledLower;

		m_pProfilePixels = NULL;
		m_pProfileBits = NULL;

//		m_pUpperProfile = NULL;
//		m_pLowerProfile = NULL;
//		m_pScaledUpper = NULL;
//		m_pScaledLower = NULL;
		m_pSameAs = NULL;
		m_pSymmetricalWith = NULL;
		m_Symmetrical = false;
		m_MidPoint = 0.0f;

		for(int Index = 0; Index < m_Faces.size(); Index++)
		{
			delete m_Faces[Index];
		}
		m_Faces.clear();								// delete each face from here, to clear hole fillers

		for(int Index = 0; Index < m_Edges.size(); Index++)
		{
			delete m_Edges[Index];
		}
		m_Edges.clear();
	}


	Island* m_pSameAs;							// for the matching
	Island* m_pSymmetricalWith;
	bool Compare(Island* pIsland);
	bool TestInteriorSymmetry();


	bool WalkOutMatch(Island* pIsland, WalkBlock* pStartBlock);
	vector<WalkBlock*> FindStartFaces(Island* pIsland);
	vector<WalkBlock*> FindStartFacesSym();
	vector<WalkBlock*> FirstMatches(vector<RKPolygon*> matchedFaces);
	bool GetConnectedFaces(WalkBlock* pThisBlock, list<WalkBlock*> &toDo);
	void PropagateUVs();

	void InitEdgeHash(int MaxVertCount);


	int m_pixelX;										// for the profile packer
	int m_pixelY;

	float UMin, VMin, UMax, VMax;
	float m_Area;

	void SetUpMesh();
	void AddTriangle(RKFace* pStartFace);
	void AddTriangles(RKFace* pStartFace);
	void AddPolygon(int PolyIndex);

	void SetPins();
	void FinaliseMesh();

	void MaintainUVs();

	void GetAreas(float& UVArea, float& PolygonArea);
	void GetUVBBox(float& MinU, float& MinV, float& MaxU, float& MaxV);
	void UVTranslate(float UTrans, float VTrans);
	void UVScale(float UVScale);
	void Rotate(float Angle);

	void MinAreaRotate();

	bool Unwrap(bool Live);
	bool Straighten();
	bool Minimise(int IslandNumber, int IslandCount);


	void DeleteProfile();
	int m_PinnedCount;

	bool m_Symmetrical;
	float m_MidPoint;									// mid point
	vector<RKVertex*> m_Vertices;
	vector<RKFace*> m_Faces;
	vector<RKEdge*> m_Edges;
	vector<RKPolygon*> m_Polygons;

	FindPins PinFinder;

	void Profiles(int textureWidthHeight, float borderScale);		//, char* pFileName);
	void CreateBitProfiles();

//	int m_scaledProfileHeight;
//	int m_scaledProfileWidth;
	int m_profileWidth;
	int m_profileHeight;
	unsigned char* m_pProfilePixels;

	int m_BitsWidth;
	uint64_t* m_pProfileBits;

//	unsigned int* m_pUpperProfile;
//	unsigned int* m_pLowerProfile;

//	unsigned int* m_pScaledUpper;
//	unsigned int* m_pScaledLower;

	void FindMinMax(int StickLength, int StickPoint, int& Min, int& Max);
	void SymUVPoints(float& U1, float& V1, float& U2, float& V2);
	bool LeftSymTest();

private:

	RKFace* FindNotVisited(RKVertex* pVertex);


	float GetUpwards();

	void FillHoles();
	void iFillHole(vector<RKVertex*> boundaryVerts);

	int chainHull_2D( RKPoint* P, int n, RKPoint* H );
	float Diameter(RKPoint* Hull, int NumberPoints);
	float MinimumArea(RKPoint* Hull, int NumberPoints);

	float vec2_angle(RKPoint* v1, RKPoint* v2, RKPoint* v3);
	float rectangleArea(RKPoint* p1, RKPoint* dir, RKPoint* p2, RKPoint* p3, RKPoint* p4);
	bool intersectLine2D(RKPoint *v1, RKPoint *dir1, RKPoint *v2, RKPoint *dir2, RKPoint *isect);
	bool intersect(RKPoint* p1, RKPoint* pointOnLine, RKPoint* dir, RKPoint* isect);
	float VecLen(RKPoint* v1, RKPoint* v2);

	void SquareGrid();
	bool CheckStraight();
	void AverageEdges();


	void MatchAngles(WalkBlock* pBlock);
	bool LookAtPoly(RKPolygon* pPoly, int BestVert, float LookAt[3][3], float Up[3],  bool SymTest);
	double TestDifference(int NumberOfVertices, float* pProjFace1, float* pProjFace2);

	bool TestSymmetrical(Island* pIsland);
	DivisionPlane* FindSymmetryPlane();
	bool m_TestSym;
	bool PlaneTest(double Norm1X, double Norm1Y, double Norm1Z, double Norm2X, double Norm2Y, double Norm2Z, vector<DivisionPlane*> &SymPlanes);
	bool MirrorTest(double PlaneA, double PlaneB, double PlaneC, double PlaneD, int FurthestLeft, double maxDist);
	bool MirrorOverlap(int LeftIndex, int RightIndex);

	double m_CentreX;
	double m_CentreY;
	double m_CentreZ;			// the islands centre

	DivisionPlane m_DivisionPlane;

	RKEdge* AddEdge(int Vert1, int Vert2, float EdgeLen, RKPolygon* pFace);
	EdgeHash* m_pEdgeHash;


	PolyFiller PolyRender;

	CPMS m_CPMS;
	LSCMap m_LSCM;
	MinimiseStretch m_MINI;
};


#endif		// ISLAND_INCLUDED
