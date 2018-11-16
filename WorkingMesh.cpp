#include "WorkingMesh.h"
#include "Unwrapper/RKProgress.h"


#include <algorithm>
#include <math.h>

#include <maya/MString.h>
#include <maya/MFnMesh.h>
#include <maya/MItMeshPolygon.h>
#include <maya/MItMeshVertex.h>
#include <maya/MItMeshEdge.h>

#include <maya/MPointArray.h>
#include <maya/MFloatPointArray.h>
#include <maya/MFnTransform.h>
#include <maya/MGlobal.h>
#include <maya/MFnDagNode.h>


#include <maya/MTransformationMatrix.h>
#include <maya/MMatrix.h>


// what to do with faces that have no UVs?

void WorkingMesh::SplitSingletons()
{
	MStatus			status;
	MObject			mesh = m_MeshNode.node();

	MObject			meshTrans = m_MeshNode.transform();
	MFnTransform	trans(meshTrans);
	ftransName =	trans.partialPathName();
//	MTransformationMatrix mat =  trans.transformation();
//	MMatrix			rotMat = mat.asMatrix();

	MFnMesh			meshfn( mesh );
	MItMeshVertex	vertIter( mesh );
	MItMeshPolygon	polyIter( mesh );
	MItMeshEdge		edgeIter( mesh );


//	MFnDependencyNode depNodeFn;
//	depNodeFn.setObject(mesh);
	m_meshNodeName = trans.fullPathName();			//depNodeFn.name();

	m_edgesToCut.clear();

	meshfn.getCurrentUVSetName(m_selUVSet);

	MFloatPointArray vertices;
	meshfn.getPoints(vertices, MSpace::kObject);
//	unsigned int numMeshVerts = vertices.length();
//	vertices.clear();


	unsigned int numMeshVerts = meshfn.numVertices();
	int numMeshEdges = meshfn.numEdges();
	int numMeshFaces = meshfn.numPolygons();

	m_QuickVerts = new QuickVertex[numMeshVerts+1];
	m_QuickEdges = new QuickEdge[numMeshEdges];
	m_QuickFaces = new QuickFace[numMeshFaces+1];

	for( ; !vertIter.isDone(); vertIter.next() )
	{
		int vertIndex = vertIter.index();
		vertIter.getConnectedFaces(m_QuickVerts[vertIndex].m_AttachedFaces);
		vertIter.getConnectedEdges(m_QuickVerts[vertIndex].m_AttachedEdges);
	}

	int numberOfTexturedFaces = 0;
	for( ; !polyIter.isDone(); polyIter.next() )
	{
		bool hasUVs = polyIter.hasUVs(m_selUVSet);
//		bool Holed = polyIter.isHoled();

		int polyIndex = polyIter.index();
		int numberOfVerts = polyIter.polygonVertexCount();

		m_QuickFaces[polyIndex].m_selected = false;
		m_QuickFaces[polyIndex].m_numberOfVerts = numberOfVerts;
		polyIter.getVertices(m_QuickFaces[polyIndex].m_VertexIndices);
		m_QuickFaces[polyIndex].m_UVIndices.setLength(numberOfVerts);

		if(hasUVs)					// && Holed == false)
		{
			m_QuickFaces[polyIndex].m_hasUVs = true;
			for(int vertIndex = 0; vertIndex < numberOfVerts; vertIndex++)
			{
				numberOfTexturedFaces++;
				polyIter.getUVIndex(vertIndex, m_QuickFaces[polyIndex].m_UVIndices[vertIndex], &m_selUVSet);
			}

			bool CutIt = CutZeroAreas(&m_QuickFaces[polyIndex], vertices);
			if(CutIt)
			{
				MIntArray EdgeList;
				polyIter.getEdges(EdgeList);
				for(int EdgeIndex = 0; EdgeIndex < EdgeList.length(); EdgeIndex++) 
				{
					m_edgesToCut.insert(EdgeList[EdgeIndex]);
				}
				EdgeList.clear();
			}
		}
		else
		{
			m_QuickFaces[polyIndex].m_hasUVs = false;
			for(int vertIndex = 0; vertIndex < numberOfVerts; vertIndex++)
			{
				m_QuickFaces[polyIndex].m_UVIndices[vertIndex] = -1;			// no uv's on this face
			}
		}
	}


	bool SelectedTextured = false;
	for(int Index = 0; Index < m_SelectedFaces.length(); Index++)
	{
		int SelectedFace = m_SelectedFaces[Index];
		m_QuickFaces[SelectedFace].m_selected = true;
		if(m_QuickFaces[SelectedFace].m_hasUVs == true) SelectedTextured = true;
	}

	if(m_SelectedFaces.length() == 0) SelectedTextured = true;


	if(numberOfTexturedFaces == 0 || SelectedTextured == false) 
	{
		vertices.clear();
		m_edgesToCut.clear();
		m_untextured = true;
		delete[] m_QuickVerts;
		delete[] m_QuickEdges;
		delete[] m_QuickFaces;
		m_QuickVerts = NULL;
		m_QuickEdges = NULL;
		m_QuickFaces = NULL;
		return;		// trying to map a set of unmapped faces
	}

	m_untextured = false;


	int edgeIndex = 0;
	for( ; !edgeIter.isDone(); edgeIter.next() )
	{
		int vertIndex1 = edgeIter.index(0);
		int vertIndex2 = edgeIter.index(1);
		m_QuickEdges[edgeIndex].m_vertexIndex1 = vertIndex1;
		m_QuickEdges[edgeIndex].m_vertexIndex2 = vertIndex2;
		edgeIter.getConnectedFaces(m_QuickEdges[edgeIndex].m_ConnectedFaces);

		if(m_QuickEdges[edgeIndex].m_ConnectedFaces.length() == 0)
		{
			int a = 0;
			a++;
		}
		

		m_QuickEdges[edgeIndex].m_Split = true;				// 1 face is a split /  >2 faces == split

		if(m_QuickEdges[edgeIndex].m_ConnectedFaces.length() == 2)
		{
			int face1 = m_QuickEdges[edgeIndex].m_ConnectedFaces[0];
			int face2 = m_QuickEdges[edgeIndex].m_ConnectedFaces[1];

			int Face1UV1 = getUVId(face1, vertIndex1);
			int Face1UV2 = getUVId(face2, vertIndex1);
			int Face2UV1 = getUVId(face1, vertIndex2);
			int Face2UV2 = getUVId(face2, vertIndex2);

			if(Face1UV1 == Face1UV2 && Face2UV1 == Face2UV2) 
			{
				m_QuickEdges[edgeIndex].m_Split = false;				// no UVs split
			}

			// one selected, one not
			if(m_QuickFaces[face1].m_selected ^ m_QuickFaces[face2].m_selected)
			{
				if(Face1UV1 == Face1UV2 && Face2UV1 == Face2UV2)			//m_edgesToCut.insert(edgeIndex);
				{
					m_QuickEdges[edgeIndex].m_Split = true;		// no potential split
				}
			}
		}

		edgeIndex++;
	}


	for(int vertIndex = 0; vertIndex < numMeshVerts; vertIndex++)
	{
		int faceCount = m_QuickVerts[vertIndex].m_AttachedFaces.length();

		if(faceCount > 1)																// can't get splits with only one face attached
		{
			splitWheel(vertIndex);
		}
	}

	delete[] m_QuickVerts;
	delete[] m_QuickEdges;
	delete[] m_QuickFaces;
	m_QuickVerts = NULL;
	m_QuickEdges = NULL;
	m_QuickFaces = NULL;
}



bool WorkingMesh::CutZeroAreas(QuickFace* pTestFace, MFloatPointArray& vertices)
{
	int EdgeCount = pTestFace->m_UVIndices.length();

	for(int Index = 0; Index < EdgeCount-2; Index++)
	{
		int LocalIndex1 = Index-1;
		int LocalIndex2 = Index;
		int LocalIndex3 = Index+1;
		if(LocalIndex1 < 0) LocalIndex1 = EdgeCount-1;

		int VertIndex1 = pTestFace->m_VertexIndices[LocalIndex1];
		int VertIndex2 = pTestFace->m_VertexIndices[LocalIndex2];
		int VertIndex3 = pTestFace->m_VertexIndices[LocalIndex3];


		//bool ZeroArea = false;
		float v1[3], v2[3], v3[3];
		double Angle1, Angle2, Angle3;
		v1[0] = vertices[VertIndex1].x;  v1[1] = vertices[VertIndex1].y;  v1[2] = vertices[VertIndex1].z;
		v2[0] = vertices[VertIndex2].x;  v2[1] = vertices[VertIndex2].y;  v2[2] = vertices[VertIndex2].z;
		v3[0] = vertices[VertIndex3].x;  v3[1] = vertices[VertIndex3].y;  v3[2] = vertices[VertIndex3].z;

		RKFace::TriangleAngles(v1, v2, v3, &Angle1, &Angle2, &Angle3);

		if(Angle1 < 0.00001f) return(true);
		if(Angle2 < 0.00001f) return(true);
		if(Angle3 < 0.00001f) return(true);
	}

	return(false);
}



/*************************************************************************
* Work through each mesh vertex.  
* Choose edge, walk outwards to find block.
* Block is a fan of faces, all with the same UV Id
* After we've found a block, if there's another block with the same UVId
* Select trailing edges from both blocks
*************************************************************************/

void WorkingMesh::splitWheel(int VertID)
{
	vector<Block*> BlockList;
	MIntArray wheelEdges = m_QuickVerts[VertID].m_AttachedEdges;
	MIntArray wheelFaces = m_QuickVerts[VertID].m_AttachedFaces;

	if(wheelFaces.length() == 1) return;				// 1 face, no splits possible

	while(wheelFaces.length() != 0)
	{
		Block* pBlock = new Block();					// time for a new block
		int startFace = wheelFaces[0];
		wheelFaces.remove(0);
		BlockList.push_back(pBlock);

		pBlock->m_UVId = getUVId(startFace, VertID);
		pBlock->m_outerFace1 = startFace;
		pBlock->m_outerFace2 = startFace;

		bool expand = true;

		while(expand)
		{
			expand = growBlock(pBlock, VertID, wheelFaces, wheelEdges);
		}
	}

	if(BlockList.size() > 1)
	{
		for(int BlockIndex = 0; BlockIndex < BlockList.size()-1; BlockIndex++)
		{
			if(BlockList[BlockIndex]->m_doneWith == false)
			{
				int thisUVId = BlockList[BlockIndex]->m_UVId;

				for(int OtherBlock = BlockIndex+1; OtherBlock < BlockList.size(); OtherBlock++)
				{
					if(thisUVId == BlockList[OtherBlock]->m_UVId)				// separate blocks, same UVId, split them off
					{
						m_edgesToCut.insert(BlockList[BlockIndex]->m_outerEdge1);
						m_edgesToCut.insert(BlockList[BlockIndex]->m_outerEdge2);
						m_edgesToCut.insert(BlockList[OtherBlock]->m_outerEdge1);
						m_edgesToCut.insert(BlockList[OtherBlock]->m_outerEdge2);
						BlockList[OtherBlock]->m_doneWith = true;
					}
				}
			}
		}
	}

	// delete the blocklist
	for(int Index = 0; Index < BlockList.size(); Index++)
	{
		Block* pBlock = BlockList[Index];
		delete pBlock;
	}

	BlockList.clear();
}


/*********************************************************
*	try to expand the block left and right
*********************************************************/

bool WorkingMesh::growBlock(Block* pBlock, int VertID, MIntArray& wheelFaces, MIntArray& wheelEdges)
{
	int outerFace1Index, prevFaceIndex = -1;
	int outerFace2Index, nextFaceIndex = -1;
	QuickEdge* prevEdge;
	QuickEdge* nextEdge;

	pBlock->m_outerEdge1 = getEdgeID(pBlock->m_outerFace1, VertID, false, wheelEdges);		// get the previous edge
	
	if(pBlock->m_outerEdge1 != -1)
	{
		prevEdge = &m_QuickEdges[pBlock->m_outerEdge1];
		outerFace1Index = pBlock->m_outerFace1;
		prevFaceIndex = prevEdge->m_ConnectedFaces[0];
		if(prevFaceIndex == outerFace1Index) prevFaceIndex = prevEdge->m_ConnectedFaces[1];		// what to do with face count > 2
		if(prevEdge->m_Split) prevFaceIndex = -1;	
	}


	pBlock->m_outerEdge2 = getEdgeID(pBlock->m_outerFace2, VertID, true, wheelEdges);		// get the next edge
	
	if(pBlock->m_outerEdge2 != -1)
	{
		nextEdge = &m_QuickEdges[pBlock->m_outerEdge2];
		outerFace2Index = pBlock->m_outerFace2;
		nextFaceIndex = nextEdge->m_ConnectedFaces[0];
		if(nextFaceIndex == outerFace2Index) nextFaceIndex = nextEdge->m_ConnectedFaces[1];
		if(nextEdge->m_Split) nextFaceIndex = -1;								//  no connecting face
	}






	if(prevFaceIndex == -1 && nextFaceIndex == -1) return(false);			// reached the limits of this block

	bool stillExpanding = false;

	if(prevFaceIndex != -1)
	{
		int foundIndex = -1;

		for(int faceIndex = 0; faceIndex < wheelFaces.length(); faceIndex++)
		{
			if(wheelFaces[faceIndex] == prevFaceIndex) foundIndex = faceIndex;
		}

		if(foundIndex != -1)
		{
			wheelFaces.remove(foundIndex);
			pBlock->m_outerFace1 = prevFaceIndex;
			stillExpanding = true;
		}
	}

	if(nextFaceIndex != -1)
	{
		int foundIndex = -1;

		for(int faceIndex = 0; faceIndex < wheelFaces.length(); faceIndex++)
		{
			if(wheelFaces[faceIndex] == nextFaceIndex) foundIndex = faceIndex;
		}

		if(foundIndex != -1)
		{
			wheelFaces.remove(foundIndex);
			pBlock->m_outerFace2 = nextFaceIndex;
			stillExpanding = true;
		}
	}

	return(stillExpanding);
}




/**********************************************************
*	get the UV ID at a vertex on a face	
**********************************************************/

int WorkingMesh::getUVId(int faceIndex, int VertID)
{
	MIntArray verts = m_QuickFaces[faceIndex].m_VertexIndices;

	for(int Index = 0; Index < verts.length(); Index++)
	{
		if(verts[Index] == VertID) return(m_QuickFaces[faceIndex].m_UVIndices[Index]);
	}

	return(-1);
}


/***********************************************************
*	get the previous or next edge around a vertex on a face
***********************************************************/

int WorkingMesh::getEdgeID(int faceIndex, int VertID, bool nextEdge, MIntArray& wheelEdges)
{
	MIntArray verts = m_QuickFaces[faceIndex].m_VertexIndices;
	int vertsLength = verts.length();

	for(int Index = 0; Index < verts.length(); Index++)
	{
		if(verts[Index] == VertID)
		{
			int nextIndex = Index + 1;
			int prevIndex = Index - 1;
			if(nextIndex > vertsLength-1) nextIndex = 0;
			if(prevIndex < 0) prevIndex = vertsLength-1;

			int prevVertIndex = verts[prevIndex];
			int nextVertIndex = verts[nextIndex];

			for(int EdgeIndex = 0; EdgeIndex < wheelEdges.length(); EdgeIndex++)
			{
				int thisEdgeIndex = wheelEdges[EdgeIndex];
				int EdgeVertIndex1 = m_QuickEdges[thisEdgeIndex].m_vertexIndex1;
				int EdgeVertIndex2 = m_QuickEdges[thisEdgeIndex].m_vertexIndex2;

				if(nextEdge == true)
				{
					if(EdgeVertIndex1 == VertID && EdgeVertIndex2 == nextVertIndex) return(thisEdgeIndex);
					if(EdgeVertIndex2 == VertID && EdgeVertIndex1 == nextVertIndex) return(thisEdgeIndex);
				}
				else
				{
					if(EdgeVertIndex1 == VertID && EdgeVertIndex2 == prevVertIndex) return(thisEdgeIndex);
					if(EdgeVertIndex2 == VertID && EdgeVertIndex1 == prevVertIndex) return(thisEdgeIndex);
				}
			}
		}
	}

	return(-1);
}




void WorkingMesh::issuePolyCut(bool hasHistory)
{
	if(m_edgesToCut.size() == 0) return;


	MString polyMapCutCommand;
	polyMapCutCommand.clear();

	if(!hasHistory)
	{
		polyMapCutCommand = "polyMapCut -ch 0";
	}
	else
	{
		polyMapCutCommand = "polyMapCut -ch 1";
	}

	int start = -1;
	int end = 0;
	char Number[20];
	set<int>::iterator it;

	for(it = m_edgesToCut.begin(); it != m_edgesToCut.end(); it++)
	{
		int next = *it;

		if(start == -1) 
		{
			start = next;
			end = next;
		}
		else
		{
			if(next == end+1)
			{
				end = next;
			}
			else
			{
				polyMapCutCommand += " ";
				polyMapCutCommand += ftransName;

				if(start == end)
				{
					sprintf(Number, ".e[%d]", start);
					polyMapCutCommand += Number;
				}
				else
				{
					sprintf(Number, ".e[%d:%d]", start,end);
					polyMapCutCommand += Number;
				}

				start = next;
				end = next;
			}
		}
	}

	polyMapCutCommand += " ";
	polyMapCutCommand += ftransName;

	if(start == end)
	{
		sprintf(Number, ".e[%d]", start);
		polyMapCutCommand += Number;
	}
	else
	{
		sprintf(Number, ".e[%d:%d]", start,end);
		polyMapCutCommand += Number;
	}

	m_edgesToCut.clear();

	int result;
	MGlobal::executeCommand(polyMapCutCommand, result, false, true);
}




void WorkingMesh::SymmetryCut()
{
	if(m_SelectedEdges.length() != 1)
	{
		MGlobal::displayError( "Roadkill: Symmetry cut: mesh has more than one edge selected");
		return;
	}


	MStatus			status;
	MObject			mesh = m_MeshNode.node();

	MObject			meshTrans = m_MeshNode.transform();
	MFnTransform	trans(meshTrans);
	ftransName =	trans.partialPathName();

	MFnMesh			meshfn( mesh );
	MItMeshVertex	vertIter( mesh );
	MItMeshPolygon	polyIter( mesh );
	MItMeshEdge		edgeIter( mesh );


	MFnDependencyNode depNodeFn;
	depNodeFn.setObject(mesh);
	m_meshNodeName = depNodeFn.name();

	m_edgesToCut.clear();

	meshfn.getCurrentUVSetName(m_selUVSet);

	MFloatPointArray vertices;
	meshfn.getPoints(vertices, MSpace::kObject);

	unsigned int numMeshVerts = meshfn.numVertices();
	int numMeshEdges = meshfn.numEdges();
	int numMeshFaces = meshfn.numPolygons();

	m_QuickVerts = new QuickVertex[numMeshVerts+1];
	m_QuickEdges = new QuickEdge[numMeshEdges];
	m_QuickFaces = new QuickFace[numMeshFaces+1];

	for( ; !vertIter.isDone(); vertIter.next() )
	{
		int vertIndex = vertIter.index();
		vertIter.getConnectedFaces(m_QuickVerts[vertIndex].m_AttachedFaces);
		vertIter.getConnectedEdges(m_QuickVerts[vertIndex].m_AttachedEdges);
	}

	int numberOfTexturedFaces = 0;
	for( ; !polyIter.isDone(); polyIter.next() )
	{
		bool hasUVs = polyIter.hasUVs(m_selUVSet);

		int polyIndex = polyIter.index();
		int numberOfVerts = polyIter.polygonVertexCount();

		m_QuickFaces[polyIndex].m_visited = false;
		m_QuickFaces[polyIndex].m_selected = false;
		m_QuickFaces[polyIndex].m_numberOfVerts = numberOfVerts;
		polyIter.getVertices(m_QuickFaces[polyIndex].m_VertexIndices);
		m_QuickFaces[polyIndex].m_UVIndices.setLength(numberOfVerts);
		polyIter.getEdges(m_QuickFaces[polyIndex].m_EdgeIndices);

		if(hasUVs)
		{
			m_QuickFaces[polyIndex].m_hasUVs = true;
			for(int vertIndex = 0; vertIndex < numberOfVerts; vertIndex++)
			{
				numberOfTexturedFaces++;
				polyIter.getUVIndex(vertIndex, m_QuickFaces[polyIndex].m_UVIndices[vertIndex], &m_selUVSet);
			}

			bool CutIt = CutZeroAreas(&m_QuickFaces[polyIndex], vertices);
			if(CutIt)
			{
				MIntArray EdgeList;
				polyIter.getEdges(EdgeList);
				for(int EdgeIndex = 0; EdgeIndex < EdgeList.length(); EdgeIndex++) 
				{
					m_edgesToCut.insert(EdgeList[EdgeIndex]);
				}
				EdgeList.clear();
			}
		}
		else
		{
			m_QuickFaces[polyIndex].m_hasUVs = false;
			for(int vertIndex = 0; vertIndex < numberOfVerts; vertIndex++)
			{
				m_QuickFaces[polyIndex].m_UVIndices[vertIndex] = -1;			// no uv's on this face
			}
		}
	}



	int edgeIndex = 0;
	for( ; !edgeIter.isDone(); edgeIter.next() )
	{
		int vertIndex1 = edgeIter.index(0);
		int vertIndex2 = edgeIter.index(1);
		m_QuickEdges[edgeIndex].m_vertexIndex1 = vertIndex1;
		m_QuickEdges[edgeIndex].m_vertexIndex2 = vertIndex2;
		m_QuickEdges[edgeIndex].m_matchedWith = -1;
		m_QuickEdges[edgeIndex].m_thisIndex = edgeIndex;
		edgeIter.getConnectedFaces(m_QuickEdges[edgeIndex].m_ConnectedFaces);
		

		m_QuickEdges[edgeIndex].m_Split = true;				// 1 face is a split /  >2 faces == split

		if(m_QuickEdges[edgeIndex].m_ConnectedFaces.length() == 2)
		{
			int face1 = m_QuickEdges[edgeIndex].m_ConnectedFaces[0];
			int face2 = m_QuickEdges[edgeIndex].m_ConnectedFaces[1];

			int Face1UV1 = getUVId(face1, vertIndex1);
			int Face1UV2 = getUVId(face2, vertIndex1);
			int Face2UV1 = getUVId(face1, vertIndex2);
			int Face2UV2 = getUVId(face2, vertIndex2);

			if(Face1UV1 == Face1UV2 && Face2UV1 == Face2UV2) 
			{
				m_QuickEdges[edgeIndex].m_Split = false;				// no UVs split
			}

			// one selected, one not
			if(m_QuickFaces[face1].m_selected ^ m_QuickFaces[face2].m_selected)
			{
				if(Face1UV1 == Face1UV2 && Face2UV1 == Face2UV2)			//m_edgesToCut.insert(edgeIndex);
				{
					m_QuickEdges[edgeIndex].m_Split = true;		// no potential split
				}
			}
		}

		edgeIndex++;
	}

	bool result = MatchUpEdges();


	if(result == true)
	{
		for(int EdgeIndex = 0; EdgeIndex < numMeshEdges; EdgeIndex++)
		{
			if(m_QuickEdges[EdgeIndex].m_matchedWith != -1 && m_QuickEdges[EdgeIndex].m_Split == false)
			{
				int MatchedEdge = m_QuickEdges[EdgeIndex].m_matchedWith;
				if(m_QuickEdges[MatchedEdge].m_Split == true)
				{
					m_edgesToCut.insert(EdgeIndex);
				}
			}
		}

		issuePolyCut(true);
	}


	m_edgesToCut.clear();
	m_untextured = true;
	delete[] m_QuickVerts;
	delete[] m_QuickEdges;
	delete[] m_QuickFaces;
	m_QuickVerts = NULL;
	m_QuickEdges = NULL;
	m_QuickFaces = NULL;
}




bool WorkingMesh::MatchUpEdges()
{
	list<LinkedEdge> EdgeList;

	int StartEdge = m_SelectedEdges[0];
	QuickEdge* pThisEdge = &m_QuickEdges[StartEdge];

	if(pThisEdge->m_ConnectedFaces.length() != 2)
	{
		MGlobal::displayError( "Roadkill: Symmetry cut: Selected Edge only connected to one polygon");
		return(false);
	}


	LinkedEdge newEdge;
	newEdge.EnteredEdgeIndex1 = StartEdge;
	newEdge.EnteredEdgeIndex2 = StartEdge;
	newEdge.PolygonIndex1 = pThisEdge->m_ConnectedFaces[0];
	newEdge.PolygonIndex2 = pThisEdge->m_ConnectedFaces[1];
	m_QuickFaces[pThisEdge->m_ConnectedFaces[0]].m_visited = true;
	m_QuickFaces[pThisEdge->m_ConnectedFaces[1]].m_visited = true;
	m_QuickEdges[StartEdge].m_matchedWith = StartEdge;
	EdgeList.push_back(newEdge);


	while(EdgeList.size() != 0)
	{
		LinkedEdge thisEdge = EdgeList.front();
		EdgeList.pop_front();

		QuickEdge* pEdge1 = &m_QuickEdges[thisEdge.EnteredEdgeIndex1];
		QuickEdge* pEdge2 = &m_QuickEdges[thisEdge.EnteredEdgeIndex2];
		QuickFace* pFace1 = &m_QuickFaces[thisEdge.PolygonIndex1];
		QuickFace* pFace2 = &m_QuickFaces[thisEdge.PolygonIndex2];

		pEdge1->m_matchedWith = pEdge2->m_thisIndex;
		pEdge2->m_matchedWith = pEdge1->m_thisIndex;

		if(pFace1->m_numberOfVerts != pFace2->m_numberOfVerts)
		{
			EdgeList.clear();
			MGlobal::displayError( "Roadkill: Symmetry cut: Mesh is not symmetrical");
			return(false);				// vertex count of linked faces does not match
		}

		int Edge1Index = FindEdgeIndex(thisEdge.PolygonIndex1, thisEdge.EnteredEdgeIndex1);
		int Edge2Index = FindEdgeIndex(thisEdge.PolygonIndex2, thisEdge.EnteredEdgeIndex2);
		int EdgeCount = pFace1->m_numberOfVerts;

		for(int Edge = 0; Edge < EdgeCount; Edge++)
		{
			int EdgeFace1 = pFace1->m_EdgeIndices[Edge1Index];
			int EdgeFace2 = pFace2->m_EdgeIndices[Edge2Index];

			int OtherFace1, OtherFace2;
			OtherFace1 = OtherFace2 = -1;

			if(m_QuickEdges[EdgeFace1].m_matchedWith == -1)
			{
				if(m_QuickEdges[EdgeFace1].m_ConnectedFaces.length() == 2)
				{
					if(m_QuickEdges[EdgeFace1].m_ConnectedFaces[0] == thisEdge.PolygonIndex1) OtherFace1 = m_QuickEdges[EdgeFace1].m_ConnectedFaces[1];
					if(m_QuickEdges[EdgeFace1].m_ConnectedFaces[1] == thisEdge.PolygonIndex1) OtherFace1 = m_QuickEdges[EdgeFace1].m_ConnectedFaces[0];
				}

				if(m_QuickEdges[EdgeFace2].m_ConnectedFaces.length() == 2)
				{
					if(m_QuickEdges[EdgeFace2].m_ConnectedFaces[0] == thisEdge.PolygonIndex2) OtherFace2 = m_QuickEdges[EdgeFace2].m_ConnectedFaces[1];
					if(m_QuickEdges[EdgeFace2].m_ConnectedFaces[1] == thisEdge.PolygonIndex2) OtherFace2 = m_QuickEdges[EdgeFace2].m_ConnectedFaces[0];
				}
			}

			if((OtherFace1 == -1 && OtherFace2 != -1) ||  (OtherFace1 != -1 && OtherFace2 == -1))
			{
				EdgeList.clear();
				MGlobal::displayError( "Roadkill: Symmetry cut: Mesh is not symmetrical");
				return(false);				// Edges are connected differently
			}

			if(OtherFace1 != -1 && OtherFace2 != -1)
			{
				m_QuickEdges[EdgeFace1].m_matchedWith = m_QuickEdges[EdgeFace2].m_thisIndex;
				m_QuickEdges[EdgeFace2].m_matchedWith = m_QuickEdges[EdgeFace1].m_thisIndex;

				if(m_QuickFaces[OtherFace1].m_visited == false && m_QuickFaces[OtherFace2].m_visited == false)
				{
					m_QuickFaces[OtherFace1].m_visited = true;
					m_QuickFaces[OtherFace2].m_visited = true;

					LinkedEdge newEdge;
					newEdge.EnteredEdgeIndex1 = EdgeFace1;
					newEdge.EnteredEdgeIndex2 = EdgeFace2;
					newEdge.PolygonIndex1 = OtherFace1;
					newEdge.PolygonIndex2 = OtherFace2;
					EdgeList.push_back(newEdge);
				}
			}

			Edge1Index++;
			Edge2Index--;

			if(Edge1Index == EdgeCount) Edge1Index = 0;
			if(Edge2Index == -1) Edge2Index = EdgeCount-1;
		}
	}

	return(true);
}



int WorkingMesh::FindEdgeIndex(int PolyIndex, int EdgeIndex)
{
	int EdgeCount = m_QuickFaces[PolyIndex].m_numberOfVerts;

	for(int OutIndex = 0; OutIndex < EdgeCount; OutIndex++)
	{
		if(m_QuickFaces[PolyIndex].m_EdgeIndices[OutIndex] == EdgeIndex) return(OutIndex);
	}

	return(-1);
}




void WorkingMesh::GetUVToFaceRatio()
{
	MFloatPointArray vertices;

	MDagPath		dagPath = m_MeshNode;
//	MObject			mesh = dagPath.node();
	MFnMesh			meshFn( dagPath );			//mesh );mesh );
	MItMeshPolygon	polyIter( dagPath );			//mesh );mesh );

/*
	double			Scale[3];
	MObject			meshTrans = m_MeshNode.transform();
	MFnTransform	trans(meshTrans);
	MTransformationMatrix mat =  trans.transformation();
	//MMatrix			rotMat = mat.asMatrix();
	mat.getScale(Scale, MSpace::kObject);
*/

	RKProgress::Get().SetText("Getting Polygons");
	RKProgress::Get().SetProgress(0);

	meshFn.getPoints(vertices, MSpace::kWorld);		//kObject);kObject);
/*	int numVerts = vertices.length();
	for(int VertIndex = 0; VertIndex < numVerts; VertIndex++)
	{
//		float X = vertices[VertIndex].x;
//		float Y = vertices[VertIndex].y;
//		float Z = vertices[VertIndex].z;

//		float outX = (X * rotMat.matrix[0][0]) + (Y * rotMat.matrix[1][0]) + (Z * rotMat.matrix[2][0]);
//		float outY = (X * rotMat.matrix[0][1]) + (Y * rotMat.matrix[1][1]) + (Z * rotMat.matrix[2][1]);
//		float outZ = (X * rotMat.matrix[0][2]) + (Y * rotMat.matrix[1][2]) + (Z * rotMat.matrix[2][2]);

		vertices[VertIndex].x *= Scale[0];
		vertices[VertIndex].y *= Scale[1];
		vertices[VertIndex].z *= Scale[2];
	}
*/

	RKProgress::Get().SetProgress(50);

	m_UCoords.setLength(0);
	m_VCoords.setLength(0);
	meshFn.getUVs(m_UCoords, m_VCoords, &m_selUVSet);
	RKProgress::Get().SetProgress(100);


	vector<int> selectedFaces;					// this is crap!

	for(int Index = 0; Index < m_SelectedFaces.length(); Index++)
	{
		int SelectedFace = m_SelectedFaces[Index];
		selectedFaces.push_back(SelectedFace);
	}

	int FaceCount = 0;
	double m_MeshAverage = 0;

	for( ; !polyIter.isDone(); polyIter.next() )
	{
//		bool goodTriangles = polyIter.hasValidTriangulation();			// bad triangle?  Highlight and error.
		bool hasUVs = polyIter.hasUVs();

		if(hasUVs == true)
		{
			MPointArray polyVerts;
			MIntArray polyIndices;

			int triCount;
			polyIter.numTriangles(triCount);
			polyIter.getTriangles(polyVerts, polyIndices, MSpace::kObject);

			MIntArray orderedIndices;
			polyIter.getVertices(orderedIndices);

//			unsigned int numVerts = polyVerts.length();
			unsigned int numOrderedVerts = orderedIndices.length();
			
			int faceIndex = polyIter.index();
			bool selected = true;

			if(selectedFaces.size() != 0)
			{
				vector<int>::iterator inSelected = find(selectedFaces.begin(), selectedFaces.end(), faceIndex);
				if(inSelected == selectedFaces.end())
				{
					selected = false;
				}
			}

			for(unsigned int Index = 0; Index < triCount; Index++)
			{
				int VertIndex1 = polyIndices[Index*3];
				int VertIndex2 = polyIndices[(Index*3)+1];
				int VertIndex3 = polyIndices[(Index*3)+2];

				int origVertIndex1 = 0;
				int origVertIndex2 = 0;
				int origVertIndex3 = 0;

				for(int tIndex = 0; tIndex < numOrderedVerts; tIndex++)
				{
					if(orderedIndices[tIndex] == VertIndex1) origVertIndex1 = tIndex;
				}

				for(int tIndex = 0; tIndex < numOrderedVerts; tIndex++)
				{
					if(orderedIndices[tIndex] == VertIndex2) origVertIndex2 = tIndex;
				}

				for(int tIndex = 0; tIndex < numOrderedVerts; tIndex++)
				{
					if(orderedIndices[tIndex] == VertIndex3) origVertIndex3 = tIndex;
				}


				int UV1Index, UV2Index, UV3Index;

				polyIter.getUVIndex(origVertIndex1, UV1Index, &m_selUVSet);
				polyIter.getUVIndex(origVertIndex2, UV2Index, &m_selUVSet);
				polyIter.getUVIndex(origVertIndex3, UV3Index, &m_selUVSet);


				float deltaX1 = vertices[VertIndex2].x - vertices[VertIndex1].x;
				float deltaY1 = vertices[VertIndex2].y - vertices[VertIndex1].y;
				float deltaZ1 = vertices[VertIndex2].z - vertices[VertIndex1].z;
				float deltaX2 = vertices[VertIndex3].x - vertices[VertIndex1].x;
				float deltaY2 = vertices[VertIndex3].y - vertices[VertIndex1].y;
				float deltaZ2 = vertices[VertIndex3].z - vertices[VertIndex1].z;
				
				float crossX = deltaY1 * deltaZ2 - deltaZ1 * deltaY2;
				float crossY = deltaX1 * deltaZ2 - deltaZ1 * deltaX2;
				float crossZ = deltaX1 * deltaY2 - deltaY1 * deltaX2;
				float Mag = sqrt((crossX * crossX) + (crossY * crossY) + (crossZ * crossZ));
				float FaceArea = Mag * 0.5f;


				deltaX1 = m_UCoords[UV2Index] - m_UCoords[UV1Index];
				deltaY1 = m_VCoords[UV2Index] - m_VCoords[UV1Index];
				deltaX2 = m_UCoords[UV3Index] - m_UCoords[UV1Index];
				deltaY2 = m_VCoords[UV3Index] - m_VCoords[UV1Index];

				crossZ = deltaX1 * deltaY2 - deltaY1 * deltaX2;
				float UVArea = fabs(crossZ * 0.5f);


				if(FaceArea > 0.00000001f && UVArea > 0.00000001f)
				{
					float StretchVal = UVArea / FaceArea;
					m_MeshAverage += StretchVal;
					FaceCount++;
				}
			}
		}
	}


	selectedFaces.clear();

	if(FaceCount != 0)
	{
		m_UvToFaceRatio = m_MeshAverage / (double)FaceCount;
	}
	else
	{
		m_UvToFaceRatio = 0.0;
	}
}




void WorkingMesh::ScaleUVs(double Scale)
{
	int numUVs = m_UCoords.length();

	for(int Index = 0; Index < numUVs; Index++)
	{
		m_UCoords[Index] *= Scale;
		m_VCoords[Index] *= Scale;
	}
}



int WorkingMesh::FindNinety(int VertIndex1, int VertIndex2, int VertIndex3)
{
	if(VertIndex1 == 0 && VertIndex2 == 1 && VertIndex3 == 2) return(2);
	if(VertIndex1 == 1 && VertIndex2 == 2 && VertIndex3 == 0) return(1);
	if(VertIndex1 == 2 && VertIndex2 == 0 && VertIndex3 == 1) return(3);
	if(VertIndex1 == 2 && VertIndex2 == 1 && VertIndex3 == 0) return(2);
	if(VertIndex1 == 0 && VertIndex2 == 2 && VertIndex3 == 1) return(3);
	if(VertIndex1 == 1 && VertIndex2 == 0 && VertIndex3 == 2) return(1);

	if(VertIndex1 == 0 && VertIndex2 == 1 && VertIndex3 == 3) return(1);
	if(VertIndex1 == 1 && VertIndex2 == 3 && VertIndex3 == 0) return(3);
	if(VertIndex1 == 3 && VertIndex2 == 0 && VertIndex3 == 1) return(2);
	if(VertIndex1 == 3 && VertIndex2 == 1 && VertIndex3 == 0) return(3);
	if(VertIndex1 == 0 && VertIndex2 == 3 && VertIndex3 == 1) return(1);
	if(VertIndex1 == 1 && VertIndex2 == 0 && VertIndex3 == 3) return(2);

	if(VertIndex1 == 0 && VertIndex2 == 2 && VertIndex3 == 3) return(3);
	if(VertIndex1 == 2 && VertIndex2 == 3 && VertIndex3 == 0) return(2);
	if(VertIndex1 == 3 && VertIndex2 == 0 && VertIndex3 == 2) return(1);
	if(VertIndex1 == 3 && VertIndex2 == 2 && VertIndex3 == 0) return(1);
	if(VertIndex1 == 0 && VertIndex2 == 3 && VertIndex3 == 3) return(2);
	if(VertIndex1 == 2 && VertIndex2 == 0 && VertIndex3 == 3) return(3);

	if(VertIndex1 == 1 && VertIndex2 == 2 && VertIndex3 == 3) return(2);
	if(VertIndex1 == 2 && VertIndex2 == 3 && VertIndex3 == 1) return(1);
	if(VertIndex1 == 3 && VertIndex2 == 1 && VertIndex3 == 2) return(3);
	if(VertIndex1 == 3 && VertIndex2 == 2 && VertIndex3 == 1) return(2);
	if(VertIndex1 == 1 && VertIndex2 == 3 && VertIndex3 == 2) return(3);
	if(VertIndex1 == 2 && VertIndex2 == 1 && VertIndex3 == 3) return(1);

	return(1);
}



void WorkingMesh::BodgyGetTriangles(MFnMesh& mesh, MIntArray& triangleCounts, MIntArray& triangleVertices)
{
	MItMeshPolygon	polyIter( mesh.object() );

	int NumPolygons = mesh.numPolygons();
	triangleCounts.setLength(NumPolygons);
	triangleVertices.clear();				//setLength(NumPolygons*6);

	int polyIndex = 0;
	//int triIndex = 0;

	for( ; !polyIter.isDone(); polyIter.next() )
	{
		MPointArray polyVerts;
		MIntArray polyIndices;

		int triCount;
		polyIter.numTriangles(triCount);
		polyIter.getTriangles(polyVerts, polyIndices, MSpace::kObject);
		triangleCounts[polyIndex++] = triCount;

		for(int Index = 0; Index < triCount; Index++)
		{
			triangleVertices.append(polyIndices[Index*3]);
			triangleVertices.append(polyIndices[Index*3+1]);
			triangleVertices.append(polyIndices[Index*3+2]);
		}
	}
}




void WorkingMesh::AddToUnwrapper(Unwrap* pUnwrapper)
{
	MFloatPointArray vertices;

	MDagPath		dagPath = m_MeshNode;
	MObject			mesh = dagPath.node();
	MFnMesh			meshFn( dagPath );		//mesh );
	MItMeshPolygon	polyIter( dagPath );		//mesh );
/*
	double			Scale[3];
	MObject			meshTrans = m_MeshNode.transform();
	MFnTransform	trans(meshTrans);
	MTransformationMatrix mat =  trans.transformation();
//	MMatrix			rotMat = mat.asMatrix();
	mat.getScale(Scale, MSpace::kObject);
*/

	RKProgress::Get().SetText("Getting Polygons");
	RKProgress::Get().SetProgress(0);
	RKProgress::Get().SetProgress(50);

	m_UCoords.setLength(0);
	m_VCoords.setLength(0);
	meshFn.getUVs(m_UCoords, m_VCoords, &m_selUVSet);
	m_numUVs = m_UCoords.length();	

	m_UVBase = pUnwrapper->SetVertexCount(m_numUVs);			// very slow!
	RKProgress::Get().SetProgress(100);


	vector<int> selectedFaces;					// this is crap!

	for(int Index = 0; Index < m_SelectedFaces.length(); Index++)
	{
		int SelectedFace = m_SelectedFaces[Index];
		selectedFaces.push_back(SelectedFace);
	}


	MIntArray triangleCounts;
	MIntArray triangleVertices;
	MIntArray uvCounts;
	MIntArray uvIDs;
	MIntArray polyVertIndices;

	meshFn.getPoints(vertices, MSpace::kWorld);		//kObject);
//	

#ifdef MAYA7
	BodgyGetTriangles(meshFn, triangleCounts, triangleVertices);
#else
	meshFn.getTriangles(triangleCounts, triangleVertices);
#endif

	meshFn.getAssignedUVs(uvCounts, uvIDs, &m_selUVSet);

	int QuadTriVerts[6];					// only for use on quads

	int BaseUVIndex = 0;
	int TriangleIndex = 0;
	int NumberOfPolygons = triangleCounts.length();

	vector<int> UVIndices;

	for(int PolyIndex = 0; PolyIndex < NumberOfPolygons; PolyIndex++)
	{
		meshFn.getPolygonVertices(PolyIndex, polyVertIndices);

		int numberPolyVerts = polyVertIndices.length();
		int UVCount = uvCounts[PolyIndex];
		int numberPolyTris = triangleCounts[PolyIndex];

		bool selected = true;

		if(selectedFaces.size() != 0)
		{
			vector<int>::iterator inSelected = find(selectedFaces.begin(), selectedFaces.end(), PolyIndex);
			if(inSelected == selectedFaces.end())
			{
				selected = false;
			}
		}

		if(UVCount != 0)
		{
			if(numberPolyVerts == 4)
			{
				for(int Test = 0; Test < 6; Test++)
				{
					int VertIndex = triangleVertices[TriangleIndex*3 + Test];

					if(VertIndex == polyVertIndices[0]) QuadTriVerts[Test] = 0;
					if(VertIndex == polyVertIndices[1]) QuadTriVerts[Test] = 1;
					if(VertIndex == polyVertIndices[2]) QuadTriVerts[Test] = 2;
					if(VertIndex == polyVertIndices[3]) QuadTriVerts[Test] = 3;
				}
			}

			UVIndices.clear();
			//int vertCount = polyIter.polygonVertexCount();

			for(int VertIndex = 0; VertIndex < numberPolyVerts; VertIndex++)
			{
				int UVIndex = uvIDs[BaseUVIndex+VertIndex];
				UVIndices.push_back(m_UVBase + UVIndex);
			}

			int PolyIndex = pUnwrapper->AddPolygon(UVIndices);
			UVIndices.clear();



			for(int Index = 0; Index < numberPolyTris; Index++)
			{
				int VertIndex1 = triangleVertices[(TriangleIndex+Index)*3];
				int VertIndex2 = triangleVertices[(TriangleIndex+Index)*3+1];
				int VertIndex3 = triangleVertices[(TriangleIndex+Index)*3+2];

				int origVertIndex1 = 0;
				int origVertIndex2 = 0;
				int origVertIndex3 = 0;

				for(int tIndex = 0; tIndex < numberPolyVerts; tIndex++)
				{
					if(polyVertIndices[tIndex] == VertIndex1) origVertIndex1 = tIndex;
				}

				for(int tIndex = 0; tIndex < numberPolyVerts; tIndex++)
				{
					if(polyVertIndices[tIndex] == VertIndex2) origVertIndex2 = tIndex;
				}

				for(int tIndex = 0; tIndex < numberPolyVerts; tIndex++)
				{
					if(polyVertIndices[tIndex] == VertIndex3) origVertIndex3 = tIndex;
				}

				int UVIndex1 = uvIDs[BaseUVIndex+origVertIndex1];
				int UVIndex2 = uvIDs[BaseUVIndex+origVertIndex2];
				int UVIndex3 = uvIDs[BaseUVIndex+origVertIndex3];

				
				bool ZeroArea = false;
				float v1[3], v2[3], v3[3];
				double Angle1, Angle2, Angle3;
				v1[0] = vertices[VertIndex1].x;  v1[1] = vertices[VertIndex1].y;  v1[2] = vertices[VertIndex1].z;
				v2[0] = vertices[VertIndex2].x;  v2[1] = vertices[VertIndex2].y;  v2[2] = vertices[VertIndex2].z;
				v3[0] = vertices[VertIndex3].x;  v3[1] = vertices[VertIndex3].y;  v3[2] = vertices[VertIndex3].z;

				RKFace::TriangleAngles(v1, v2, v3, &Angle1, &Angle2, &Angle3);

				if(Angle1 < 0.0001f) ZeroArea = true;
				if(Angle2 < 0.0001f) ZeroArea = true;
				if(Angle3 < 0.0001f) ZeroArea = true;

				int NinetyIndex = 1;

				if(numberPolyVerts == 4)
				{
					if(Index == 0)
					{
						NinetyIndex = FindNinety(QuadTriVerts[0], QuadTriVerts[1], QuadTriVerts[2]);

						if(QuadTriVerts[0] == 0 && QuadTriVerts[1] == 1 && QuadTriVerts[2] == 2) NinetyIndex = 1;
						if(QuadTriVerts[0] == 0 && QuadTriVerts[1] == 1 && QuadTriVerts[2] == 2) NinetyIndex = 1;

						if(QuadTriVerts[0] == 0 && QuadTriVerts[1] == 2 && QuadTriVerts[2] == 3) NinetyIndex = 2;
					}
					else
					{
						NinetyIndex = FindNinety(QuadTriVerts[3], QuadTriVerts[4], QuadTriVerts[5]);

						if(QuadTriVerts[1] == 0 && QuadTriVerts[4] == 2 && QuadTriVerts[5] == 3) NinetyIndex = 1;
						if(QuadTriVerts[1] == 1 && QuadTriVerts[4] == 2 && QuadTriVerts[5] == 3) NinetyIndex = 2;
					}
				}

					
					pUnwrapper->AddVertex(m_UVBase + UVIndex1, VertIndex1, vertices[VertIndex1].x, vertices[VertIndex1].y, vertices[VertIndex1].z, m_UCoords[UVIndex1], m_VCoords[UVIndex1]);
					pUnwrapper->AddVertex(m_UVBase + UVIndex2, VertIndex2, vertices[VertIndex2].x, vertices[VertIndex2].y, vertices[VertIndex2].z, m_UCoords[UVIndex2], m_VCoords[UVIndex2]);
					pUnwrapper->AddVertex(m_UVBase + UVIndex3, VertIndex3, vertices[VertIndex3].x, vertices[VertIndex3].y, vertices[VertIndex3].z, m_UCoords[UVIndex3], m_VCoords[UVIndex3]);

				if(selected == true)
				{	
					ZeroArea = pUnwrapper->AddFace(m_UVBase + UVIndex1, m_UVBase + UVIndex2, m_UVBase + UVIndex3, PolyIndex, ZeroArea, NinetyIndex);
				}
			}
		}

		BaseUVIndex += uvCounts[PolyIndex];
		TriangleIndex += numberPolyTris;
	}


/*
	vector<int> UVIndices;


	for( ; !polyIter.isDone(); polyIter.next() )
	{
//		bool goodTriangles = polyIter.hasValidTriangulation();			// bad triangle?  Highlight and error.
		bool hasUVs = polyIter.hasUVs();

		if(hasUVs == true)
		{
			UVIndices.clear();
			int vertCount = polyIter.polygonVertexCount();

			for(int VertIndex = 0; VertIndex < vertCount; VertIndex++)
			{
				int UVIndex;
				polyIter.getUVIndex(VertIndex, UVIndex, &m_selUVSet);
				UVIndices.push_back(m_UVBase + UVIndex);
			}

			int PolyIndex = pUnwrapper->AddPolygon(UVIndices);
			UVIndices.clear();



			MPointArray polyVerts;
			MIntArray polyIndices;

			int triCount;
			polyIter.numTriangles(triCount);
			polyIter.getTriangles(polyVerts, polyIndices, MSpace::kObject);

			MIntArray orderedIndices;
			polyIter.getVertices(orderedIndices);

//			unsigned int numVerts = polyVerts.length();
			unsigned int numOrderedVerts = orderedIndices.length();
			
			int faceIndex = polyIter.index();
			bool selected = true;

			if(selectedFaces.size() != 0)
			{
				vector<int>::iterator inSelected = find(selectedFaces.begin(), selectedFaces.end(), faceIndex);
				if(inSelected == selectedFaces.end())
				{
					selected = false;
				}
			}

			for(unsigned int Index = 0; Index < triCount; Index++)
			{
				int VertIndex1 = polyIndices[Index*3];
				int VertIndex2 = polyIndices[(Index*3)+1];
				int VertIndex3 = polyIndices[(Index*3)+2];

				int origVertIndex1 = 0;
				int origVertIndex2 = 0;
				int origVertIndex3 = 0;

				for(int tIndex = 0; tIndex < numOrderedVerts; tIndex++)
				{
					if(orderedIndices[tIndex] == VertIndex1) origVertIndex1 = tIndex;
				}

				for(int tIndex = 0; tIndex < numOrderedVerts; tIndex++)
				{
					if(orderedIndices[tIndex] == VertIndex2) origVertIndex2 = tIndex;
				}

				for(int tIndex = 0; tIndex < numOrderedVerts; tIndex++)
				{
					if(orderedIndices[tIndex] == VertIndex3) origVertIndex3 = tIndex;
				}


				int UVIndex1, UVIndex2, UVIndex3;

				polyIter.getUVIndex(origVertIndex1, UVIndex1, &m_selUVSet);
				polyIter.getUVIndex(origVertIndex2, UVIndex2, &m_selUVSet);
				polyIter.getUVIndex(origVertIndex3, UVIndex3, &m_selUVSet);



				bool ZeroArea = false;
				float v1[3], v2[3], v3[3];
				double Angle1, Angle2, Angle3;
				v1[0] = vertices[VertIndex1].x;  v1[1] = vertices[VertIndex1].y;  v1[2] = vertices[VertIndex1].z;
				v2[0] = vertices[VertIndex2].x;  v2[1] = vertices[VertIndex2].y;  v2[2] = vertices[VertIndex2].z;
				v3[0] = vertices[VertIndex3].x;  v3[1] = vertices[VertIndex3].y;  v3[2] = vertices[VertIndex3].z;

				RKFace::TriangleAngles(v1, v2, v3, &Angle1, &Angle2, &Angle3);

				if(Angle1 < 0.0001f) ZeroArea = true;
				if(Angle2 < 0.0001f) ZeroArea = true;
				if(Angle3 < 0.0001f) ZeroArea = true;
					
				if(selected == true)			// && ZeroArea == false) 
				{	
					pUnwrapper->AddVertex(m_UVBase + UVIndex1, VertIndex1, vertices[VertIndex1].x, vertices[VertIndex1].y, vertices[VertIndex1].z, m_UCoords[UVIndex1], m_VCoords[UVIndex1]);
					pUnwrapper->AddVertex(m_UVBase + UVIndex2, VertIndex2, vertices[VertIndex2].x, vertices[VertIndex2].y, vertices[VertIndex2].z, m_UCoords[UVIndex2], m_VCoords[UVIndex2]);
					pUnwrapper->AddVertex(m_UVBase + UVIndex3, VertIndex3, vertices[VertIndex3].x, vertices[VertIndex3].y, vertices[VertIndex3].z, m_UCoords[UVIndex3], m_VCoords[UVIndex3]);
					ZeroArea = pUnwrapper->AddFace(m_UVBase + UVIndex1, m_UVBase + UVIndex2, m_UVBase + UVIndex3, PolyIndex, ZeroArea);
				}
			}
		}
	}
*/

	selectedFaces.clear();
}




void WorkingMesh::SetNewUVs(Unwrap* pUnwrapper)
{
	for(int Index = 0; Index < m_numUVs; Index++)
	{
		m_UCoords[Index] = pUnwrapper->m_Vertices[m_UVBase + Index]->U;
		m_VCoords[Index] = pUnwrapper->m_Vertices[m_UVBase + Index]->V;
	}
}
