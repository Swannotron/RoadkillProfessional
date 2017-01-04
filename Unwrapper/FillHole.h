#ifndef FILLHOLE
#define FILLHOLE


#include "Primitives.h"

#include <list>
using namespace std;

class DelaunayEdge;		// forward reference


class DelaunayVert
{
public:
	float m_X;
	float m_Y;
	float m_Z;
};

class DelaunayFace
{
public:
	int m_Edge1Index;				// clockwise edges
	int m_Edge2Index;
	int m_Edge3Index;
};

class DelaunayEdge
{
public:
	DelaunayEdge()
	{
		m_Face1Index = -1;
		m_Face2Index = -1;
		m_Vert1Index = -1;
		m_Vert2Index = -1;
	}

	int m_Face1Index;
	int m_Face2Index;
	int m_Vert1Index;
	int m_Vert2Index;
};


class SVec3D
{
public:
	float x;
	float y;
	float z;
};






class HoleChunk
{
public:
	HoleChunk() { m_pVertIndices.clear(); m_Triangles.clear(); }
	~HoleChunk() { m_pVertIndices.clear(); m_Triangles.clear(); }

	list<int> m_pVertIndices;
	vector<int> m_Triangles;
	bool m_Complex;

	void Copy(HoleChunk* pChunkIn)
	{
		this->m_pVertIndices.clear();
		this->m_Triangles.clear();
		this->m_Complex = pChunkIn->m_Complex;
		this->m_pVertIndices.insert(this->m_pVertIndices.end(), pChunkIn->m_pVertIndices.begin(), pChunkIn->m_pVertIndices.end());
		this->m_Triangles.insert(this->m_Triangles.end(), pChunkIn->m_Triangles.begin(), pChunkIn->m_Triangles.end());
	}
};



class FillHole
{
public:

	FillHole() {};
	~FillHole() {};

	void HoleFiller(vector<RKVertex*> &Boundary, vector<RKFace*> &Faces);
	bool CheckChunk(HoleChunk* pChunk, bool AngleCheck);
	bool CheckComplex(int numberOfVertices, double* pVerts);		//, double BestPlane[4]);

	void BreakChunk(HoleChunk* pWholeBoundary);
	void SplitChunks(HoleChunk* pSourceChunk, HoleChunk* pDestChunk);
	void ReOrder(HoleChunk* pWholeBoundary);

	void GetTri(HoleChunk* pChunk, int EdgeIndex1, int EdgeIndex2, int& TriIndex1, int& TriIndex2, int& TriIndex3);
	float GetSmallestIterior(int TriIndex1, int TriIndex2, int TriIndex3);

	list<int>::iterator GetNext(list<int> &Vertices, list<int>::iterator VertIterator);
	list<int>::iterator GetPrev(list<int> &Vertices, list<int>::iterator VertIterator);

	void projectPolygon( double BestPlane[4], int NumberOfVertices );

	void Delaunay(int NumberOfVertices);
	int EdgeIndex(vector<DelaunayEdge*> &pEdges, int VertIndex1, int VertIndex2, int FaceIndex);

	bool DelaunayFlip(DelaunayFace* pFaces, vector<DelaunayEdge*> &pEdges, int EdgeIndex);
	bool TestConcave(int VertIndex1, int VertIndex2, int VertIndex3, int VertIndex4);
	bool InCircle(int TestPointIndex, int P1Index, int P2Index, int P3Index);
	void SetOpposites(vector<DelaunayEdge*> &pEdges, DelaunayFace* pFaces, int FaceIndex1, int FaceIndex2);


	vector<RKVertex*> m_pVertices;
};

/*
class CTriangulator
{
	vector <SVec3D*> m_points;
//    CVec3DArray	*m_points;              // points on contour

	bool _insideTriangle(SVec3D* A, SVec3D* B, SVec3D* C, SVec3D* P);

    float _area();
    bool _snip(int u, int v, int w, int n, int *V);
    void _process(vector<int> &indices);

public:

	CTriangulator();
    ~CTriangulator();

    void add(SVec3D* p);

    void triangulate(vector<int> &indices);			//CUIntArray *indices);

    ///     Returns the given point in the triangulator array
    inline SVec3D* get(const int id) { return(m_points[id]); }
};
*/


#endif		// FILLHOLE
