#ifndef PRIMITIVES
#define PRIMITIVES


#include <vector>
#include <list>
#include <set>


using namespace std;


class RKVertex;
class RKFace;
class RKPolygon;


//#define M_PI		3.14159265358979323846


class RKEdge
{
public:
	RKEdge()
	{
		m_Vert1 = -1;
		m_Vert2 = -1;
		m_pFaceList.clear();
	}

	~RKEdge()
	{
		m_pFaceList.clear();
	}

	double m_length;
	int m_Vert1;
	int m_Vert2;

	vector<RKPolygon*> m_pFaceList;

	RKPolygon* GetOtherFace(RKPolygon* pFaceIn)
	{
		for(int Index = 0; Index < m_pFaceList.size(); Index++)
		{
			if(m_pFaceList[Index] != pFaceIn) return(m_pFaceList[Index]);
		}

		return(NULL);
	}
};


class RKFace
{
public:

	RKFace()
	{
		m_Visited = false;
		pVert1 = NULL;
		pVert2 = NULL;
		pVert3 = NULL;
		pEdge1Cut = true;			// welded together later
		pEdge2Cut = true;
		pEdge3Cut = true;
		m_Filler = false;
	}

	~RKFace()
	{
		m_Visited = false;
		pVert1 = NULL;
		pVert2 = NULL;
		pVert3 = NULL;
		pEdge1Cut = true;			// welded together later
		pEdge2Cut = true;
		pEdge3Cut = true;
		m_Filler = false;
	}


	RKVertex* pVert1;
	RKVertex* pVert2;
	RKVertex* pVert3;

	float NX, NY, NZ;									// triangle Normal

	bool m_Visited;
	bool m_Filler;

	double m_InteriorAngle1;
	double m_InteriorAngle2;
	double m_InteriorAngle3;

	double m_ABFAngle1;
	double m_ABFAngle2;
	double m_ABFAngle3;

	bool pEdge1Cut;
	bool pEdge2Cut;
	bool pEdge3Cut;

	int m_UnwrapIndex;

	int NinetyAngle;
	float m_Edge1Len;
	float m_Edge2Len;
	int m_WhichEdge;


	int m_pOrigPolygonindex;



	inline bool operator==(const RKFace &rhs) const
	{
		if(this == &rhs) return(true);
		return(false);
	}

	void GetFaceAngles();
	float GetArea();
	float GetUVArea();

	RKVertex* GetNextTrailing(RKVertex* pVertIn, float* pLength);
	RKVertex* GetPrevTrailing(RKVertex* pVertIn, float* pLength);

	float GetUVLengthsConnected(RKVertex* pThisVert);

	static void TriangleAngles(float *pVert1, float *pVert2, float *pVert3, double *Angle1, double *Angle2, double *Angle3);
	static void TriangleAngles(RKVertex *pVert1, RKVertex *pVert2, RKVertex *pVert3, double *Angle1, double *Angle2, double *Angle3);
	static double VecAngle(float *pVert1, float *pVert2, float *pVert3);

	void GetUVAngles(float& Angle1, float& Angle2, float& Angle3);
	bool IntersectUVs(float U1, float V1, float U2, float V2, float& Len3D);
	void SetEdgeLen(float EdgeLen);
	void NinteyAngles();

	bool Compare(RKFace* pFace);

private:

	static double VecAngleCos(float *pVert1, float *pVert2, float *pVert3);
	static double Normalise(double *n);
};




class RKVertex
{
public:

	RKVertex()
	{
		U = V = 0.0f;
		m_NewU = m_NewV = 0.0f;
		m_Pinned = false;
		m_splitVert = false;
		m_OnEdge = false;
		m_Visited = false;
		m_Added = false;
		m_Used = false;
		m_OriginalVertIndex = -1;
		m_InteriorIndex = -1;
		m_pMatches = NULL;
//		m_LinkedFaces.clear();
	}

	~RKVertex()
	{
		U = V = 0.0f;
		m_NewU = m_NewV = 0.0f;
		m_Pinned = false;
		m_splitVert = false;
		m_OnEdge = false;
		m_Visited = false;
		m_Added = false;
		m_Used = false;
		m_InteriorIndex = -1;
		m_pMatches = NULL;
		m_LinkedFaces.clear();
		if(m_LinkedEdges.size() != 0)
		{
			for(int Index = 0; Index < m_LinkedEdges.size(); Index++)
			{
				delete(m_LinkedEdges[Index]);
				m_LinkedEdges[Index] = NULL;
			}
		}
		m_LinkedEdges.clear();
	}


	float X, Y, Z;
	float U, V;											// calculated by the unwrapper

	float m_OldU, m_OldV;
	float m_NewU, m_NewV;

	list<RKPolygon*> m_LinkedPolygons;
	list<RKFace*> m_LinkedFaces;
	vector<RKEdge*> m_LinkedEdges;

	bool m_Used;
	bool m_OnEdge;
	bool m_splitVert;
	bool m_Pinned;
	bool m_Visited;
	bool m_Added;

	int UnwrapIndex;
	int m_OriginalVertIndex;
	int m_UVIndex;
	int m_InteriorIndex;
	RKVertex* m_pMatches;

	void LinkUp();

	RKEdge* getEdge(int Vert2, float Length);

	void getVertexNormal(double& VNX, double& VNY, double& VNZ);
	void getFaceNormal(double& VNX, double& VNY, double& VNZ, RKVertex* pNextVertex);

	double m_distFromPlane;


	RKVertex* GetNextTrailing(float* pLength)
	{
		RKVertex* pReturnVert = NULL;

		list<RKFace*>::iterator FaceIterator;
		for(FaceIterator = m_LinkedFaces.begin(); FaceIterator != m_LinkedFaces.end(); FaceIterator++)
		{
			RKVertex* pTrailing = (*FaceIterator)->GetNextTrailing(this, pLength);
			if(pTrailing != NULL)
			{
				if(pReturnVert != NULL) return(NULL);			// two trialing edges from one vertex, it's all gone wrong..
				pReturnVert = pTrailing;
			}
		}

		return(pReturnVert);
	}


	RKVertex* GetPrevTrailing(float* pLength)
	{
		RKVertex* pReturnVert = NULL;

		list<RKFace*>::iterator FaceIterator;
		for(FaceIterator = m_LinkedFaces.begin(); FaceIterator != m_LinkedFaces.end(); FaceIterator++)
		{
			RKVertex* pTrailing = (*FaceIterator)->GetPrevTrailing(this, pLength);
			if(pTrailing != NULL)
			{
				if(pReturnVert != NULL) return(NULL);			// two trialing edges from one vertex, it's all gone wrong..
				pReturnVert = pTrailing;
			}
		}

		return(pReturnVert);
	}


/*
	RKEdge* GetTrailingEdge(RKVertex* pVert2)
	{
		RKEdge* pReturnEdge = NULL;

		list<RKFace*>::iterator FaceIterator;
		for(FaceIterator = m_LinkedFaces.begin(); FaceIterator != m_LinkedFaces.end(); FaceIterator++)
		{
			pReturnEdge = (*FaceIterator)->FindTrailing(this, pVert2);
			if(pReturnEdge != NULL) return(pReturnEdge);
		}

		return(NULL);
	}
*/
};



class RKPolygon
{
public:
	RKPolygon()
	{
		m_visited = false;
		m_visited2 = false;
		m_numberVertices = -1;
		m_InteriorAngles.clear();
		m_pVertices.clear();
		m_pEdges.clear();
	}

	~RKPolygon()
	{
		m_visited = false;
		m_visited2 = false;
		m_numberVertices = -1;
		m_InteriorAngles.clear();
		m_pVertices.clear();
		m_pEdges.clear();
	}

	bool m_visited;
	bool m_visited2;					// for local symmetry finder
	int m_numberVertices;
	vector<double> m_InteriorAngles;
	vector<RKVertex*> m_pVertices;
	vector<RKEdge*> m_pEdges;

	bool Compare(RKPolygon* pPolyIn);
	bool CompareReverse(RKPolygon* pPolyIn);

	bool Compare(RKPolygon* pPolyIn, int StartPointA, int StartPointB);
	bool CompareReverse(RKPolygon* pPolyIn, int StartPointA, int StartPointB);
};




#endif

