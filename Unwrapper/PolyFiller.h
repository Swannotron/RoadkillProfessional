#ifndef POLYFFILLER
#define POLYFFILLER


class EdgeCalc
{
public:
	EdgeCalc() {};
	~EdgeCalc() {};


	float m_U1, m_V1;
	float m_U2, m_V2;
//	float m_PosX1, m_PosY1, m_PosZ1;
//	float m_PosX2, m_PosY2, m_PosZ2;
//	float m_NormX1, m_NormY1, m_NormZ1;
//	float m_NormX2, m_NormY2, m_NormZ2;
	int m_TextureHeight;
//	bool m_OnEdge;

	int m_CurrentY;
	float m_CurrentU;
//	float m_CurrentPosX, m_CurrentPosY, m_CurrentPosZ;
//	float m_CurrentNormX, m_CurrentNormY, m_CurrentNormZ;

	bool ClipY();
	bool CalcDeltas();
	bool NextLine();
	void HighestVFirst();

private:

	int m_Height;

	float m_AddonU;
//	float m_AddonPosX, m_AddonPosY, m_AddonPosZ;
//	float m_AddonNormX, m_AddonNormY, m_AddonNormZ;
};



class PolyFiller
{
public:
	PolyFiller() {};
	~PolyFiller() {};

	void RenderPolygon(float X1, float Y1, float X2, float Y2, float X3, float Y3);

	void setTexture(unsigned char* PixMap, int PixMapWidth, int PixMapHeight)
	{
		m_pPixMap = PixMap;
		m_textureHeight = PixMapHeight;
		m_textureWidth = PixMapWidth;
	}


private:

	void RenderLine(float X1, float Y1, float X2, float Y2);

	float m_leftEdge[4096];
	float m_rightEdge[4096];

	unsigned char* m_pPixMap;
	int m_textureHeight;
	int m_textureWidth;

	void CalcEdge(EdgeCalc& WorkingEdge);
	void WritePixels();
	void WriteFatPixel(int XPos, int YPos);
};


#endif
