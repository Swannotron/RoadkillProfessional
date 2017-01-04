#include "FindPins.h"

#include <algorithm>
#include <math.h>



bool CompareBoundaries(const BoundaryEdge* a, const BoundaryEdge* b)
{
	return a->m_Length > b->m_Length;
}


/********************************************************************************
*	It's all gone a bit wrong, pin the extremes of the mesh as a last resort	*
********************************************************************************/

void FindPins::PinMeshExtremes(vector<RKVertex*> &rVertices)
{
	float MinX, MinY, MinZ;
	float MaxX, MaxY, MaxZ;
	int MinIndices[3];
	int MaxIndices[3];

	MinX = MinY = MinZ = (float)1e20;
	MaxX = MaxY = MaxZ = (float)-1e20;

	vector<RKVertex*>::iterator iter = rVertices.begin();

	while(iter != rVertices.end())
	{
		RKVertex* pVert = *iter;
		if(pVert->X < MinX) { MinX  = pVert->X; MinIndices[0] = pVert->UnwrapIndex; }
		if(pVert->Y < MinY) { MinY  = pVert->Y; MinIndices[1] = pVert->UnwrapIndex; }
		if(pVert->Z < MinZ) { MinZ  = pVert->Z; MinIndices[2] = pVert->UnwrapIndex; }

		if(pVert->X > MaxX) { MaxX  = pVert->X; MaxIndices[0] = pVert->UnwrapIndex; }
		if(pVert->Y > MaxY) { MaxY  = pVert->Y; MaxIndices[1] = pVert->UnwrapIndex; }
		if(pVert->Z > MaxZ) { MaxZ  = pVert->Z; MaxIndices[2] = pVert->UnwrapIndex; }
		iter++;
	}

	float MaxDistance = (float)-1e20;
	int DirIndex;

	if(MaxX - MinX > MaxDistance)
	{
		MaxDistance = MaxX - MinX;
		DirIndex = 0;
	}

	if(MaxY - MinY > MaxDistance)
	{
		MaxDistance = MaxY - MinY;
		DirIndex = 1;
	}

	if(MaxZ - MinZ > MaxDistance)
	{
		MaxDistance = MaxZ - MinZ;
		DirIndex = 2;
	}

	m_Pin1 = MinIndices[DirIndex];
	m_Pin2 = MaxIndices[DirIndex];
}







/************************************************************
*	Attempt to find the mesh boundary.						*
*	Find an edge only connected to one face					*
*	Walk round boundary from that edge, flagging as we go	*
*	Measure length of all connected edges					*
*	Longest set of connected edges is the boundary			*
************************************************************/

bool FindPins::FindBoundaries(vector<RKVertex*> &rVertices)				//, RKEdge* pEdgeList)
{
	m_Pin1 = -1;
	m_Pin2 = -1;


	vector<RKVertex*>::iterator iter = rVertices.begin();

	while(iter != rVertices.end())
	{
		RKVertex* pThisVertex = (*iter);

		if(pThisVertex->m_OnEdge == true && pThisVertex->m_Visited == false)
		{
			BoundaryEdge* NewBoundary = new BoundaryEdge;
			float* pLength = &NewBoundary->m_Length;

			NewBoundary->pBoundary.push_back(pThisVertex);
			pThisVertex->m_Visited = true;

			RKVertex* pNextVertex = pThisVertex->GetNextTrailing(pLength);

			while(pNextVertex != pThisVertex)
			{
				if(pNextVertex == NULL)			// all gone wrong, try to get the length of all the edges involved
				{
					return(false);
				}

				pNextVertex->m_Visited = true;
				NewBoundary->pBoundary.push_back(pNextVertex);
				pNextVertex = pNextVertex->GetNextTrailing(pLength);
			}

			m_Boundaries.push_back(NewBoundary);
		}

		iter++;
	}

	if(m_Boundaries.size() > 1)
	{
		std::sort(m_Boundaries.begin(), m_Boundaries.end(), CompareBoundaries);
	}

	return(true);

/*	m_Pin1 = -1;
	m_Pin2 = -1;
	float MaxLength = 0.0f;
	RKEdge* pStartEdge = NULL;
	RKEdge* pOppositeEdge = NULL;
	m_BoundaryEdge = NULL;


	for(int EdgeIndex = 0; EdgeIndex < NumberOfEdges; EdgeIndex++) pEdgeList[EdgeIndex].m_Visited = false;

	for(int EdgeIndex = 0; EdgeIndex < NumberOfEdges; EdgeIndex++)
	{
		if(pEdgeList[EdgeIndex].m_Visited == false && pEdgeList[EdgeIndex].m_pFace2 == NULL)
		{
			RKEdge* pReturnedOpposite;

			float Length = WalkEdge(&pEdgeList[EdgeIndex], pReturnedOpposite);
			BoundaryEdge NewBoundary;
			NewBoundary.m_pStartEdge = &pEdgeList[EdgeIndex];
			NewBoundary.m_Length = Length;
			m_Boundaries.push_back(NewBoundary);

			if(Length > MaxLength)
			{
				MaxLength = Length;
				pStartEdge = &pEdgeList[EdgeIndex];
				pOppositeEdge = pReturnedOpposite;
				m_BoundaryEdge = pStartEdge;
			}
		}
	}

	m_OuterLength = MaxLength;
*/
}


struct ltstr
{
  bool operator()(const int s1, const int s2) const
  {
    return(s1 == s2);
  }
};



void FindPins::SplitVerts()
{
	map<const int, vector<RKVertex*> > splitList;
	// test if biggest boundary is fooked..
	vector<RKVertex*> pBoundary = m_Boundaries[0]->pBoundary;

	for(int Index = 0; Index < pBoundary.size(); Index++)
	{
		RKVertex* thisVert = pBoundary[Index];
		int Key = thisVert->m_OriginalVertIndex;

		map<const int, vector<RKVertex*> >::iterator iter = splitList.find(Key);

		if(iter == splitList.end())
		{
			vector<RKVertex*> newList;		// = new list<RKVertex*>;
			newList.push_back(thisVert);
			splitList[Key] = newList;
		}
		else
		{
			(*iter).second.push_back(thisVert);;
		}
	}

	map<const int, vector<RKVertex*> >::iterator iter = splitList.begin();

	while(iter != splitList.end())
	{
		int size = (*iter).second.size();

		if(size > 1)
		{
			for(int Index = 0; Index < size; Index++)
			{
				(*iter).second[Index]->m_splitVert = true;
			}
		}

		(*iter).second.clear();

		iter++;
	}
}





// edge lengths no longer calculated.
// must now walk from vertex to vertex

void FindPins::GetPins()
{
	m_Pin1 = -1;
	m_Pin2 = -1;

	if(m_Boundaries.size() == 0) return;

	SplitVerts();
	// check longest bondary is OK


	vector<RKVertex*> Boundary = m_Boundaries[0]->pBoundary;
	RKVertex* pClearVert = NULL;
	RKVertex* pStartLongest = NULL;
	RKVertex* pEndLongest = NULL;
	RKVertex* pPrevVert = NULL;
	RKVertex* pCurrVert = NULL;
	RKVertex* pStartCurrent = NULL;
	float currentLength = 0.0f;
	float longestLength = 0.0f;
	float LengthPrev = 0.0f;
	float LengthNext = 0.0f;
	float dist = 0.0f;

	bool hasSplits = false;

	// find unsplit vertex
	for(int Index = 0; Index < Boundary.size()-1; Index++)
	{
		if(Boundary[Index]->m_splitVert == true && Boundary[Index+1]->m_splitVert == true) hasSplits = true;
		if(Boundary[Index]->m_splitVert == false && Boundary[Index+1]->m_splitVert == false)
		{
			pClearVert = Boundary[Index];
			break;
		}
	}

	if(pClearVert == NULL)
	{
		if(hasSplits == true) Mixture();
		return;								// no split verts
	}

	pPrevVert = pClearVert;
	pCurrVert = pClearVert->GetNextTrailing(&dist);

	while(pCurrVert != pClearVert)
	{
		float nextDist = 0.0f;
		RKVertex* pNextVert = pCurrVert->GetNextTrailing(&nextDist);

		if(pCurrVert->m_splitVert || pPrevVert->m_splitVert & pNextVert->m_splitVert)
		{
			if(currentLength == 0)
			{
				pStartCurrent = pCurrVert;
			}
			currentLength += dist;
		}
		else
		{
			if(currentLength > longestLength)
			{
				pStartLongest = pStartCurrent;
				pEndLongest = pNextVert;
				longestLength  = currentLength;
				currentLength = 0.0f;
			}
		}

		dist = nextDist;
		pPrevVert = pCurrVert;
		pCurrVert = pNextVert;
	}

	if(currentLength > longestLength)
	{
		pStartLongest = pStartCurrent;
		pEndLongest = pCurrVert;
		longestLength  = currentLength;
	}

	if(pStartLongest == NULL)				// we don't have a set of vert splits
	{
		return;
	}

	if(longestLength < 0.5f * m_Boundaries[0]->m_Length)
	{
		Mixture();
		return;
	}


	while(pStartLongest != pEndLongest)
	{
		float EdgeLen = 0.0f;

		if(LengthNext < LengthPrev)
		{
			pStartLongest = pStartLongest->GetNextTrailing(&EdgeLen);
			LengthNext += EdgeLen;
		}
		else
		{
			pEndLongest = pEndLongest->GetPrevTrailing(&EdgeLen);
			LengthPrev += EdgeLen;
		}
	}

	m_Pin1 = pStartLongest->UnwrapIndex;



	LengthPrev = 0.0f;
	LengthNext = 0.0f;
	pStartLongest = pStartLongest->GetPrevTrailing(&LengthNext);


	while(pStartLongest != pEndLongest)
	{
		float EdgeLen = 0.0f;

		if(LengthNext < LengthPrev)
		{
			pStartLongest = pStartLongest->GetPrevTrailing(&EdgeLen);
			LengthNext += EdgeLen;
		}
		else
		{
			pEndLongest = pEndLongest->GetNextTrailing(&EdgeLen);
			LengthPrev += EdgeLen;
		}
	}

	m_Pin2 = pStartLongest->UnwrapIndex;


/*

//	FindBoundaries(NumberOfEdges, pEdgeList);
	if(m_BoundaryEdge == NULL) return;				// can't find boundary edge

	RKEdge* pClearEdge = NULL;
	RKEdge* pNextEdge = m_BoundaryEdge->GetNextTrailing();
	float currentLength = 0.0f;
	float longestLength = 0.0f;
	RKEdge* pStartLongest = NULL;
	RKEdge* pEndLongest = NULL;
	RKEdge* pStartCurrent = NULL;

	float LengthPrev = 0.0f;
	float LengthNext = 0.0f;


	while(pNextEdge != m_BoundaryEdge)
	{
		if(pNextEdge->m_pVertex1->m_splitVert == false)
		{
			pClearEdge = pNextEdge;
			break;
		}
		pNextEdge = pNextEdge->GetNextTrailing();
	}

	float totalLength = 0.0f;

	if(pClearEdge != NULL)
	{
		pNextEdge = pClearEdge->GetNextTrailing();


		while(pNextEdge != pClearEdge)
		{
			totalLength += pNextEdge->m_Length;
			if(pNextEdge->m_pVertex2->m_splitVert == true)
			{
				if(currentLength == 0)
				{
					pStartCurrent = pNextEdge;
				}
				currentLength += pNextEdge->m_Length;
			}
			else
			{
				if(currentLength > longestLength)
				{
					pStartLongest = pStartCurrent;
					pEndLongest = pNextEdge;
					longestLength  = currentLength;
					currentLength = 0.0f;
				}
			}

			pNextEdge = pNextEdge->GetNextTrailing();
		}
	}



	if(pStartLongest == NULL || longestLength < 0.5f * totalLength)				// we don't have a set of vert splits
	{
		return;
//		pStartLongest = m_BoundaryEdge;
//		pEndLongest = m_BoundaryEdge->GetPrevTrailing();
	}


	while(pStartLongest != pEndLongest)
	{
		if(LengthNext < LengthPrev)
		{
			LengthNext += pStartLongest->m_Length;
			pStartLongest = pStartLongest->GetNextTrailing();
		}
		else
		{
			LengthPrev += pEndLongest->m_Length;
			pEndLongest = pEndLongest->GetPrevTrailing();
		}


	}

	m_Pin1 = pStartLongest->m_pVertex1->UnwrapIndex;
	pStartLongest->m_pVertex1->m_Pinned = true;			// testing


	LengthPrev = 0.0f;
	LengthNext = 0.0f;
	pStartLongest = pStartLongest->GetPrevTrailing();


	while(pStartLongest != pEndLongest)
	{
		if(LengthNext < LengthPrev)
		{
			LengthNext += pStartLongest->m_Length;
			pStartLongest = pStartLongest->GetPrevTrailing();
		}
		else
		{
			LengthPrev += pEndLongest->m_Length;
			pEndLongest = pEndLongest->GetNextTrailing();
		}
	}

	m_Pin2 = pStartLongest->m_pVertex2->UnwrapIndex;

	pStartLongest->m_pVertex2->m_Pinned = true;			// testing
*/
}



/*
float FindPins::WalkEdge(RKEdge* pEdge, RKEdge* &pOppositeEdge)
{
	float LengthPrev = 0.0f;
	float LengthNext = 0.0f;
	RKEdge* pPrevEdge = pEdge->GetPrevTrailing();
	RKEdge* pNextEdge = pEdge;

	if(pPrevEdge == NULL || pNextEdge == NULL) return(0.0f);

	pNextEdge->m_Visited = true;
	pPrevEdge->m_Visited = true;

	while(pPrevEdge != pNextEdge)							// will Edge always come back on itself?
	{
		if(LengthNext < LengthPrev)
		{
			LengthNext += pNextEdge->m_Length;
			pNextEdge = pNextEdge->GetNextTrailing();
			if(pNextEdge == NULL) return(LengthPrev + LengthNext);
			if(pNextEdge->m_Visited == true) break;
			pNextEdge->m_Visited = true;
		}
		else
		{
			LengthPrev += pPrevEdge->m_Length;
			pPrevEdge = pPrevEdge->GetPrevTrailing();
			if(pPrevEdge == NULL) return(LengthPrev + LengthNext);
			if(pPrevEdge->m_Visited == true) break;
			pPrevEdge->m_Visited = true;
		}
	}

	LengthNext += pNextEdge->m_Length;									// temp test!!!!!
	pOppositeEdge = pNextEdge;
	return(LengthPrev + LengthNext);

}
*/






void FindPins::Mixture()
{
	vector<RKVertex*> Dimples;

	float EdgeLen = 0.0f;
	RKVertex* pLoopStart = NULL;
	RKVertex* pLoopEnd = NULL;

	// should we start in a split?  I think so!
	RKVertex* pStartVertex = m_Boundaries[0]->pBoundary[0];

	while(pStartVertex->m_splitVert == false)
	{
		pStartVertex = pStartVertex->GetNextTrailing(&EdgeLen);
	}

	RKVertex* pPrevVertex = pStartVertex;
	RKVertex* pThisVertex = pStartVertex->GetNextTrailing(&EdgeLen);
	RKVertex* pNextVertex = pThisVertex->GetNextTrailing(&EdgeLen);
	pPrevVertex->m_Visited = false;

	while(pThisVertex != pStartVertex)
	{
		pThisVertex->m_Visited = false;			// re-use as dimple pointer

		if(pPrevVertex->m_splitVert == true && pThisVertex->m_splitVert == false && pNextVertex->m_splitVert == true)
		{
			Dimples.push_back(pThisVertex);
			pThisVertex->m_Visited = true;
		}

		if(pPrevVertex->m_splitVert == true && pThisVertex->m_splitVert == false && pNextVertex->m_splitVert == false)
		{
			pLoopStart = pThisVertex;
		}

		if(pPrevVertex->m_splitVert == false && pThisVertex->m_splitVert == false && pNextVertex->m_splitVert == true)
		{
			if(pLoopStart != NULL)
			{
				if(pLoopStart->m_OriginalVertIndex != pThisVertex->m_OriginalVertIndex)
				{
					pLoopEnd = pThisVertex;
					float DistA = 0.0f;
					float DistB = 0.0f;

					// find midpoint of loop
					while(pLoopStart != pLoopEnd)
					{
						if(DistA < DistB)
						{
							pLoopStart = pLoopStart->GetNextTrailing(&DistA);
							//DistA += EdgeLen;
						}
						else
						{
							pLoopEnd = pLoopEnd->GetPrevTrailing(&DistB);
							//DistB += EdgeLen;
						}
					}

					Dimples.push_back(pLoopStart);
					pLoopStart->m_Visited = true;
				}
			}

			pLoopStart = NULL;
		}

		pPrevVertex = pThisVertex;
		pThisVertex = pNextVertex;
		pNextVertex = pThisVertex->GetNextTrailing(&EdgeLen);
	}

	RKVertex* pVert1 = NULL;
	RKVertex* pVert2 = NULL;
	float Longest = 0.0f;
	float BoundaryLength = m_Boundaries[0]->m_Length;


	if(Dimples.size() == 0) return;

	if(Dimples.size() == 1)				// find other UV opposite this one
	{
		float ThisLength = 0.0f;
		pVert1 = Dimples[0];
		pVert2 = pVert1->GetNextTrailing(&ThisLength);

		while(ThisLength < BoundaryLength * 0.5f)
		{
			pVert2 = pVert2->GetNextTrailing(&ThisLength);
			//ThisLength += EdgeLen;
		}

		m_Pin1 = pVert1->UnwrapIndex;
		m_Pin2 = pVert2->UnwrapIndex;
		return;
	}

	if(Dimples.size() == 2)
	{
		m_Pin1 = Dimples[0]->UnwrapIndex;
		m_Pin2 = Dimples[1]->UnwrapIndex;
		return;
	}


	for(int DimpleIndex = 0; DimpleIndex < Dimples.size(); DimpleIndex++)
	{
		float ThisLength = 0.0f;
		RKVertex* pDimple1 = Dimples[DimpleIndex];
		RKVertex* pThis = pDimple1->GetNextTrailing(&ThisLength);

		while(ThisLength < BoundaryLength * 0.5f)
		{
			if(pThis->m_Visited == true)
			{
				if(ThisLength > Longest)
				{
					Longest = ThisLength;
					pVert1 = pDimple1;
					pVert2 = pThis;
				}
			}

			pThis = pThis->GetNextTrailing(&ThisLength);
			//ThisLength += EdgeLen;
		}
	}

	if(pVert1 == NULL) return;

	m_Pin1 = pVert1->UnwrapIndex;
	m_Pin2 = pVert2->UnwrapIndex;
}
