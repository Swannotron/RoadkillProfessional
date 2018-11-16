#ifndef FIND_PINS
#define FIND_PINS

#include "Primitives.h"

#include <list>
#include <map>

using namespace std;


class BoundaryEdge
{
public:
	BoundaryEdge()
	{
		pBoundary.clear();
		m_Length = 0.0f;
	}

	~BoundaryEdge()
	{
		pBoundary.clear();
		m_Length = 0.0f;
	}


	vector<RKVertex*> pBoundary;
	float m_Length;
};




class FindPins
{
public:

	FindPins()
	{
		m_Pin1 = -1;
		m_Pin2 = -1;
	}

	~FindPins()
	{
		m_Pin1 = -1;
		m_Pin2 = -1;

		for(int Index = 0; Index < m_Boundaries.size(); Index++)
		{
			delete m_Boundaries[Index];
		}

		m_Boundaries.clear();
	}

	void GetPins();
	bool FindBoundaries(vector<RKVertex*> &rVertices);
	void PinMeshExtremes(vector<RKVertex*> &rVertices);

	int m_Pin1;						// index of first pinned vertex
	int m_Pin2;						// index of second pinned vertex

//	RKEdge* m_BoundaryEdge;
	float m_OuterLength;
	vector<BoundaryEdge*> m_Boundaries;
//	float VecAngle(RKVertex *pVert1, RKVertex *pVert2, RKVertex *pVert3);			// should be somewhere else


private:

	void SplitVerts();

	void Mixture();

//	void FindSymmetry();
//	float Dist3DPointLine(float PX, float PY, float PZ, float X0, float Y0, float Z0, float X1, float Y1, float Z1);

//	float WalkEdge(RKEdge* pEdge, RKEdge* &pOppositeEdge);



//	float Normalise(float *n);
//	float VecAngleCos(float *pVert1, float *pVert2, float *pVert3);



};





#endif //FIND_PINS
