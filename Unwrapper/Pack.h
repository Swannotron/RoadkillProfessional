#ifndef PACK
#define PACK

#include <vector>

#include "Primitives.h"

using namespace std;

class Island;


#define MAX2(x,y)		( (x)>(y) ? (x) : (y) )
#define MAX3(x,y,z)		MAX2( MAX2((x),(y)) , (z) )

#define SHIFT3(type, a, b, c) { type tmp; tmp = a; a = c; c = b; b = tmp; }


class Vec2D
{
public:
	float X, Y;

};


class Pack
{
public:
	Pack()
	{
		SortedListOfIslands.clear();
		m_pCurrentProfile = NULL;
	};

	~Pack()
	{
		SortedListOfIslands.clear();

		if(m_pCurrentProfile != NULL) delete[] m_pCurrentProfile;
		m_pCurrentProfile = NULL;
	};

	void PackIslands(vector<Island*>& listOfIslands, float borderScale);

private:

	void MaintainUVs(vector<Island*>& listOfIslands);

	void FinalBounds();

	float CalcAreas(vector<Island*>& listOfIslands, float borderScale);
	bool ProfilePack(int TextureSize);
//	bool PackTry(float side);
	float GetUpwards(Island* pIsland);
	float GetUpwards2(Island* pIsland);

	void AddFirstProfile(int TextureSize);
	void CopyFromTo(unsigned int* pProfileIn, int InWidth, int InHeight);
	void MergeProfile(Island* pIsland);
	void ExpandProfile(int NewSize);

	bool BestFit(int IslandIndex);
	bool GetStatsScaled(Island* pIsland, int ProfileX, int& AreaUnder, int& ExpansionRequired, int& ProfileYOut);
	void GetStats(Island* pIsland, int ProfileX, int& AreaUnder, int& ExpansionRequired, int& ProfileYOut);

	int FindLast();

	int TestFit(int TestX, int TestY, Island* pIsland, int& Yindex, int NumberTests);
	void BurnIn(int PosX, int PosY, int fits, Island* pIsland);

	vector<Island*> SortedListOfIslands;

	unsigned int m_TextureBitsWidth;
	unsigned int m_TextureWidthHeight;
//	unsigned char* m_pCurrentProfile;
	uint64_t* m_pCurrentProfile;

	float m_TotalUVArea;

//	void addToSet(float VertUpX, float VertUpY, float VertDownX, float VertDownY, int& FaceCount, int& FlipCount, float& VecX, float& VecY);
//	float GetAngle(float VecX, float VecY);

//	vector<float> m_outVecs;
/*
	int m_FrontCount;
	int m_BackCount;
	int m_LeftCount;
	int m_RightCount;
	int m_TopCount;
	int m_BottomCount;


	int m_FrontFlipCount;
	int m_BackFlipCount;
	int m_LeftFlipCount;
	int m_RightFlipCount;
	int m_TopFlipCount;
	int m_BottomFlipCount;

	float m_VecXFront;
	float m_VecYFront;
	float m_VecXBack;
	float m_VecYBack;

	float m_VecXLeft;
	float m_VecYLeft;
	float m_VecXRight;
	float m_VecYRight;

	float m_VecXTop;
	float m_VecYTop;
	float m_VecXBottom;
	float m_VecYBottom;
*/


	void GetUpwards3(Island* pIsland);
	void GetUVsForPlane(RKFace* pFace, float& MidPoint, float& U1Out, float& V1Out, float& U2Out, float& V2Out, int Plane);
	void DeltaCalc(RKVertex* pVert1, RKVertex* pVert2, float& MidPoint, float& U1Out, float& V1Out, float& U2Out, float& V2Out, bool& FoundFirst, int Plane);


	void FindMeanLine(vector<Vec2D> &VecsX, float& VecOutX, float& VecOutY);

	void FinalRotate(Island* pIsland);
	float FindBestRot(int TestDelta, RKVertex* pVertMin, RKVertex* pVertMax, int testWidthHeight);


};


#endif
