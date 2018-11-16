#ifndef WORKINGMESH_INCLUDED
#define WORKINGMESH_INCLUDED

#include "Unwrapper/Unwrap.h"

#include <vector>
#include <set>

#include <maya/MFloatArray.h>
#include <maya/MIntArray.h>
#include <maya/MDagPath.h>
#include <maya/MString.h>
#include <maya/MFnMesh.h>
#include <maya/MFloatPointArray.h>

using namespace std;



struct LinkedEdge
{
	int EnteredEdgeIndex1;
	int PolygonIndex1;

	int EnteredEdgeIndex2;
	int PolygonIndex2;
};


class QuickFace
{
public:

	~QuickFace()
	{
		m_VertexIndices.setLength(0);
		m_UVIndices.setLength(0);
		m_EdgeIndices.setLength(0);
	}

	bool m_hasUVs;
	bool m_selected;
	int m_numberOfVerts;
	bool m_visited;
	MIntArray m_VertexIndices;
	MIntArray m_EdgeIndices;
	MIntArray m_UVIndices;
};


class QuickEdge
{
public:

	~QuickEdge()
	{
		m_ConnectedFaces.setLength(0);
	}

	int m_vertexIndex1;
	int m_vertexIndex2;
	bool m_Split;

	int m_thisIndex;
	int m_matchedWith;

	MIntArray m_ConnectedFaces;
};


class QuickVertex
{
public:
	~QuickVertex()
	{
		m_AttachedFaces.setLength(0);
		m_AttachedEdges.setLength(0);
	}

	MIntArray m_AttachedFaces;
	MIntArray m_AttachedEdges;
};


class Block
{
public:
	Block()
	{
		m_doneWith = false;
		m_selected = false;
		m_UVId = -1;
		m_outerEdge1 = -1;
		m_outerEdge2 = -1;
		m_outerFace1 = -1;
		m_outerFace2 = -1;
	}

	bool m_doneWith;
	bool m_selected;
	int m_UVId;
	int m_outerEdge1;
	int m_outerEdge2;
	int m_outerFace1;
	int m_outerFace2;
};








class WorkingMesh
{
public:
	WorkingMesh() 
	{
		m_QuickVerts = NULL;
		m_QuickEdges = NULL;
		m_QuickFaces = NULL;
		m_edgesToCut.clear();
		m_untextured = true;
	};

	~WorkingMesh() 
	{
		if(m_QuickVerts) delete[] m_QuickVerts;
		m_QuickVerts = NULL;

		if(m_QuickEdges) delete[] m_QuickEdges;
		m_QuickEdges = NULL;

		if(m_QuickFaces) delete[] m_QuickFaces;
		m_QuickFaces = NULL;

		m_edgesToCut.clear();
		m_UCoords.setLength(0);
		m_VCoords.setLength(0);
		m_SelectedFaces.setLength(0);
	};


	MDagPath m_MeshNode;
	MIntArray m_SelectedFaces;
	MIntArray m_SelectedEdges;
	MString	m_selUVSet;
	MString m_meshNodeName;

	void SplitSingletons();
	void issuePolyCut(bool hasHistory);
	void GetUVToFaceRatio();
	void ScaleUVs(double Scale);

	void SymmetryCut();

	void AddToUnwrapper(Unwrap* pUnwrapper);
	void SetNewUVs(Unwrap* pUnwrapper);

//	void DeselectEdges();
//	void SelectEdges(bool hasHistory, bool AddOrRemove);

	bool m_untextured;
	double m_UvToFaceRatio;

	MFloatArray m_UCoords;
	MFloatArray m_VCoords;



private:
	void splitWheel(int VertID);
	bool growBlock(Block* pBlock, int VertID, MIntArray& wheelFaces, MIntArray& wheelEdges);
	int getUVId(int faceIndex, int VertID);
	int getEdgeID(int faceIndex, int VertID, bool nextEdge, MIntArray& wheelEdges);
	bool CutZeroAreas(QuickFace* pTestFace, MFloatPointArray& vertices);

	int FindNinety(int VertIndex1, int VertIndex2, int VertIndex3);

	bool MatchUpEdges();
	int FindEdgeIndex(int PolyIndex, int EdgeIndex);

	void BodgyGetTriangles(MFnMesh& mesh, MIntArray& triangleCounts, MIntArray& triangleVertices);

	MString			ftransName;

	QuickVertex*	m_QuickVerts;
	QuickEdge*		m_QuickEdges;
	QuickFace*		m_QuickFaces;

	set<int>		m_edgesToCut;

	unsigned int	m_numUVs;
	int				m_UVBase;


};




#endif  // WORKINGMESH_INCLUDED

