#include "PolyFiller.h"	
#include <math.h>

extern int gScale;


void PolyFiller::RenderPolygon(float X1, float Y1, float X2, float Y2, float X3, float Y3)
{
	for(int Index = 0; Index < m_textureHeight; Index++)
	{
		m_leftEdge[Index] = (float)-1e18;
		m_rightEdge[Index] = (float)-1e18;
	}

	EdgeCalc WorkingEdge;
	WorkingEdge.m_TextureHeight = this->m_textureHeight;

	WorkingEdge.m_U1 = X1;
	WorkingEdge.m_V1 = Y1;
	WorkingEdge.m_U2 = X2;
	WorkingEdge.m_V2 = Y2;
	CalcEdge(WorkingEdge);
	RenderLine(X1, Y1, X2, Y2);

	WorkingEdge.m_U1 = X2;
	WorkingEdge.m_V1 = Y2;
	WorkingEdge.m_U2 = X3;
	WorkingEdge.m_V2 = Y3;
	CalcEdge(WorkingEdge);
	RenderLine(X2, Y2, X3, Y3);

	WorkingEdge.m_U1 = X3;
	WorkingEdge.m_V1 = Y3;
	WorkingEdge.m_U2 = X1;
	WorkingEdge.m_V2 = Y1;
	CalcEdge(WorkingEdge);
	RenderLine(X3, Y3, X1, Y1);

	WritePixels();
}



void PolyFiller::RenderLine(float X1, float Y1, float X2, float Y2)
{
	if(Y2 < Y1)
	{
		float TempX, TempY;
		TempX = X1; TempY = Y1;
		X1 = X2;  Y1 = Y2;
		X2 = TempX; Y2 = TempY;
	}

	float DeltaX = X2 - X1;
	float SpanX = fabs(DeltaX);
	float DeltaY = Y2 - Y1;
	float XPos = X1;
	float YPos = Y1;

	if(DeltaY > SpanX)		// horzontal bias
	{
		int Count = (int)(Y2 - Y1);
		float XAdder = (X2 - X1) / Count;

		do
		{
			int XPixel = (int)XPos;
			int YPixel = (int)YPos;

//			if(YPixel >= 0 && YPixel < m_textureHeight && XPixel >= 0 && XPixel < m_textureWidth)
//			{
				WriteFatPixel(XPixel, YPixel);
//				m_pPixMap[(YPixel * m_textureWidth) + XPixel] = 0xff0000;
//			}
			YPos += 1.0f;
			XPos += XAdder;
			Count--;

		}while(Count > 0);

	}
	else
	{
		int Count = (int)SpanX;									// vertical bais
		float YAdder = (Y2 - Y1) / Count;
		int XPixel = (int)XPos;
		int XAdder = (int)((X2 - X1) / Count);

		do
		{
			int YPixel = (int)YPos;

//			if(YPixel >= 0 && YPixel < m_textureHeight && XPixel >= 0 && XPixel < m_textureWidth)
//			{
				WriteFatPixel(XPixel, YPixel);
//				m_pPixMap[(YPixel * m_textureWidth) + XPixel] = 0xff;
//			}
			XPixel += XAdder;
			YPos += YAdder;
			Count--;

		}while(Count > 0);
	}
}




void EdgeCalc::HighestVFirst()
{
	if(m_V2 < m_V1)
	{
		float TempU, TempV;
		TempU = m_U1; TempV = m_V1;
		m_U1 = m_U2; m_V1 = m_V2;
		m_U2 = TempU; m_V2 = TempV;
	}
}


bool EdgeCalc::CalcDeltas()
{
	m_AddonU = (m_U2 - m_U1) / (m_V2 - m_V1);

	int Y1TSpace = (int)(m_V1+0.5f);						// 0.5 is there because I want to hit the centre of the pixel
	int Y2TSpace = (int)(m_V2+0.5f);

	float PixelCentre = ((float)Y1TSpace)+0.5f;
	float DiffY1 = PixelCentre - m_V1;

	m_CurrentU = (DiffY1 * m_AddonU) + m_U1;

	m_CurrentY = Y1TSpace;
	m_Height = Y2TSpace;

	if(m_Height - m_CurrentY == 0) return(false);

	return(true);
}





bool EdgeCalc::NextLine()
{
	m_CurrentU += m_AddonU;

	m_CurrentY++;
	if(m_CurrentY >= m_Height) return(false);
	return(true);
}




bool EdgeCalc::ClipY()
{
	if(m_V2 < 0.0f) return(false);
	if(m_V1 > m_TextureHeight) return(false);

	if(m_V1 < 0.0f)
	{
		float DeltaU = m_U2 - m_U1;
		float DeltaV = m_V2 - m_V1;
		float PercentageClipped = -m_V1 / (DeltaV);
		m_U1 = (DeltaU * PercentageClipped) + m_U1;
		m_V1 = (DeltaV * PercentageClipped) + m_V1;
	}


	if(m_V2 > m_TextureHeight)
	{
		float DeltaU = m_U2 - m_U1;
		float DeltaV = m_V2 - m_V1;
		float PercentageClipped = (m_TextureHeight-m_V2) / DeltaV;
		m_U2 = (DeltaU * PercentageClipped) + m_U2;
		m_V2 = (DeltaV * PercentageClipped) + m_V2;
	}

	return(true);
}







void PolyFiller::CalcEdge(EdgeCalc& WorkingEdge)
{
	WorkingEdge.HighestVFirst();
	bool Result = WorkingEdge.ClipY();
	if(Result == false) return;

	Result = WorkingEdge.CalcDeltas();
	if(Result == false) return;


	do
	{
		float EdgeX = WorkingEdge.m_CurrentU;
		int YIndex = WorkingEdge.m_CurrentY;

		if(m_rightEdge[YIndex] <= EdgeX)					// new edge is less than the one in rightEdge
		{
			m_leftEdge[YIndex] = m_rightEdge[YIndex];			// move to rightEdge
			m_rightEdge[YIndex] = (float)-1e18;
		}

		if(EdgeX <= m_rightEdge[YIndex])
		{
			m_leftEdge[YIndex] = EdgeX;
		}
		else
		{
			m_rightEdge[YIndex] = EdgeX;

		}
	} while(WorkingEdge.NextLine());
}




void PolyFiller::WritePixels()
{
	for(int LineCount = 0; LineCount < m_textureHeight; LineCount++)
	{
		if(m_rightEdge[LineCount] > 0.0f && m_leftEdge[LineCount] < m_textureWidth)
		{
			int X1TSpace = (int)(m_leftEdge[LineCount]+0.5f);
			int X2TSpace = (int)(m_rightEdge[LineCount]+0.5f);

			int Width = X2TSpace - X1TSpace;
			if(Width > 0)
			{
				for(int XIndex = X1TSpace; XIndex < X2TSpace; XIndex++)
				{
					WriteFatPixel(XIndex, LineCount);
				}
			}
		}
	}
}




void PolyFiller::WriteFatPixel(int XPos, int YPos)
{
	int WriteX = XPos;
	int WriteY = YPos;
	int pixelDim = gScale*2+1;

	// no writing outside the buffer
//	if(XPos <= 0 || XPos > m_textureWidth-((gScale*2)+1) || YPos <= 0 || YPos > m_textureHeight-((gScale*2)+1)) return;

	// write a fat pixel
	for(int FatXIndex = 0; FatXIndex < pixelDim; FatXIndex++)
	{
		WriteX = XPos;

		for(int FatYIndex = 0; FatYIndex < pixelDim; FatYIndex++)
		{
			if(WriteX >= 0 && WriteX < m_textureWidth && WriteY >= 0 && WriteY < m_textureHeight)
			{
				m_pPixMap[WriteY * m_textureWidth + WriteX] = 0xff;
			}
			WriteX++;
		}

		WriteY++;
	}
}