#ifndef NEWUNWRAP
#define NEWUNWRAP

#include <vector>

#include "Primitives.h"
#include "Island.h"
#include "IslandMatch.h"
#include "Pack.h"

using namespace std;



class Unwrap
{
public:

	Unwrap()
	{
		m_listOfIslands.clear();
		m_Vertices.clear();
		m_Triangles.clear();
		m_ZeroArea.clear();	
	}

	~Unwrap()
	{
		DeleteIslands();
		DeleteVertsAndFaces();
		m_Vertices.clear();
		m_Triangles.clear();
		m_ZeroArea.clear();	
	}

	void DoUnwrap(bool UseCPMS);

	int SetVertexCount(int Size);
	void AddVertex(int UVIndex, int OriginalVertIndex, float X, float Y, float Z, float U, float V);
	bool AddFace(int VertIndex1, int VertIndex2, int VertIndex3, int PolygonIndex, bool ZeroArea, int NinetyIndex);
	int Unwrap::AddPolygon(vector<int> VertIndices);
	void DeleteVertsAndFaces();

	vector<RKVertex*> m_Vertices;
	vector<RKFace*> m_Triangles;
	vector<RKFace*> m_ZeroArea;

	static vector<RKPolygon*> m_Polygons;


private:
	vector<Island*> m_listOfIslands;

	Pack IslandPacker;
	IslandMatch IslandMatcher;

	bool UnwrapIslands(bool UseCPMS);
	bool MinimiseIslands();

	void CreateIslands();
	void DeleteIslands();
	int FindStartFace();
	void LinkEdges();
	
	void FixZeroArea();

	void ScaleMeshes();
};



#endif
