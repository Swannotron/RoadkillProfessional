#include <algorithm>
#include <math.h>

#include <stdint.h>

#include <maya/MGlobal.h>

#include "Pack.h"
#include "Island.h"
#include "Vector4.h"
#include "RKProgress.h"

extern int gFunction;
extern int gUpright;

extern int gMaintainUVs;


bool CompareIslands(const Island* a, const Island* b)
{
	return a->m_Area > b->m_Area;
}




void Pack::PackIslands(vector<Island*>& listOfIslands, float borderScale)
{
	if (gMaintainUVs)
	{
		MaintainUVs(listOfIslands);
		return;
	}


	if(borderScale > 2.0f) borderScale = 2.0f;

	float MaxSide = CalcAreas(listOfIslands, borderScale);						// get the maximum edge length of all islands

	if(SortedListOfIslands.size() == 0) return;

	if (SortedListOfIslands.size() == 1)
	{
		FinalBounds();
		return;
	}

	float MaxSideScaled = MaxSide * (borderScale + (borderScale - 1.0f) * 2.0f);
	float MinSide = MAX2(sqrt(m_TotalUVArea), MaxSide);


	MaxSide = (( (int)sqrt( (float)(SortedListOfIslands.size()-1) ) + 1) ) * MaxSideScaled;

	if(MaxSide > 4.0f) MaxSide = 4.0f; 

	bool firstPack = true;

	if (MinSide < MaxSide)											// should always be true
	{
		for (int i = 0; i < 15; i++)								// PACK_SEARCH_DEPTH			// I reckon can be a lot less
		{
			int TextureSize = (int)(((MinSide + MaxSide)*0.5f + (float)1e-5) * 1024.0f);			//512.0f );

			if (ProfilePack( (int) (((MinSide + MaxSide)*0.5f + (float)1e-5) * 1024.0f)) )			//256.0f )) )
			{
				MaxSide = (MinSide + MaxSide)*0.5f;
				firstPack = false;

				int MaxTexSize = (int)(((MaxSide)+(float)1e-5) * 256.0f);
				int MinTexSize = (int)(((MinSide)+(float)1e-5) * 256.0f);

				if (MaxTexSize - MinTexSize <= 2) break;
			}
			else
			{
				if(firstPack == true)
				{
					MaxSide *= 2.0f;
				}
				else
				{
					MinSide = (MinSide + MaxSide)*0.5f;
				}
			}

			RKProgress::Get().SetProgress(i * 6);
		}
	}

	RKProgress::Get().SetProgress(90);

	float Side = MaxSide + (float)1e-5;
	ProfilePack( (int)(Side*1024.0f) );			//256.0f) );
	float thisUVScale = (1024.0f / (float)m_TextureWidthHeight);			//256.0f / (float)m_TextureWidthHeight;


	for(unsigned int Index = 0; Index < SortedListOfIslands.size(); Index++)
	{
		Island* pIsland = SortedListOfIslands[Index];

		float thisUTrans = (float)pIsland->m_pixelX / (float)m_TextureWidthHeight;
		float thisVTrans = (float)pIsland->m_pixelY / (float)m_TextureWidthHeight;

		pIsland->UVScale(thisUVScale);
		pIsland->UVTranslate(thisUTrans, thisVTrans);
		pIsland->DeleteProfile();
	}

	FinalBounds();
}



void Pack::FinalBounds()
{
	if(SortedListOfIslands.size() == 0) return;

	// final check and scale of UVs
	float MinU = 1.0f;
	float MinV = 1.0f;
	float MaxU = 0.0f;
	float MaxV = 0.0f;

	for(unsigned int Index = 0; Index < SortedListOfIslands.size(); Index++)
	{
		Island* pIsland = SortedListOfIslands[Index];
		float ThisMinU, ThisMinV, ThisMaxU, ThisMaxV;

		pIsland->GetUVBBox(ThisMinU, ThisMinV, ThisMaxU, ThisMaxV);

		if(ThisMinU < MinU) MinU = ThisMinU;
		if(ThisMinV < MinV) MinV = ThisMinV;
		if(ThisMaxU > MaxU) MaxU = ThisMaxU;
		if(ThisMaxV > MaxV) MaxV = ThisMaxV;
	}

/*	
	int numVertices = m_Vertices.size();

	for(int Index = 0; Index < numVertices; Index++)
	{
		RKVertex* pVert = m_Vertices[Index];
		if(pVert->U < MinU) MinU = pVert->U;
		if(pVert->V < MinV) MinV = pVert->V;
		if(pVert->U > MaxU) MaxU = pVert->U;
		if(pVert->V > MaxV) MaxV = pVert->V;
	}
*/
	float OffsetU = 0.0f;
	float OffsetV = 0.0f;
	float ScaleUV = 1.0f;

	if(MinU < 0.0f || MinV < 0.0f)
	{
		OffsetU = -MinU;
		OffsetV = -MinV;
	}

	float OverU = MaxU-MinU;
	float OverV = MaxV-MinV;

	if(OverU > OverV) ScaleUV = 1.0f/OverU;
	else
		ScaleUV = 1.0f/OverV;


	for(unsigned int Index = 0; Index < SortedListOfIslands.size(); Index++)
	{
		Island* pIsland = SortedListOfIslands[Index];

		pIsland->UVScale(ScaleUV);
		pIsland->UVTranslate(OffsetU, OffsetV);
	}
/*
	for(int Index = 0; Index < numVertices; Index++)
	{
		RKVertex* pVert = m_Vertices[Index];
		pVert->U += OffsetU;
		pVert->V += OffsetV;
		pVert->U *= ScaleUV;
		pVert->V *= ScaleUV;
	}
*/
}



float Pack::CalcAreas(vector<Island*>& listOfIslands, float borderScale)
{
	float UVArea, PolygonArea;
	float MinU, MinV, MaxU, MaxV;

	float TotalArea = 0.0f;
	float MaxSide = 0.0f;
//	float MinSide = 0.0f;

	m_TotalUVArea = 0.0f;

//	char output[255];

	SortedListOfIslands.clear();


	for(unsigned int Index = 0; Index < listOfIslands.size(); Index++)
	{
		float Rescale = 1.0f;
		Island* pIsland = listOfIslands[Index];

		pIsland->GetAreas(UVArea, PolygonArea);

		if(UVArea != 0)
		{
			if(gUpright == 2)
			{
				GetUpwards3(pIsland);
			}
			
			if(gUpright == 3)
			{
				pIsland->MinAreaRotate();
			}
		}


		pIsland->GetUVBBox(MinU, MinV, MaxU, MaxV);

		if(gFunction != 3)			// pack but no scale
		{
			Rescale = (UVArea > 0.0f)? sqrt(PolygonArea)/sqrt(UVArea): 0.0f;
		}


		MinU = -MinU;
		MinV = -MinV;
		pIsland->UVTranslate(MinU, MinV);
		pIsland->UVScale(Rescale);

		pIsland->GetUVBBox(MinU, MinV, MaxU, MaxV);
		float UDelta = MaxU - MinU;
		float VDelta = MaxV - MinV;
		float IslandArea = UDelta * VDelta;

		pIsland->m_Area = UDelta;
		if(VDelta > UDelta) pIsland->m_Area = VDelta;
		//pIsland->m_Area = IslandArea;
		TotalArea += IslandArea;

		SortedListOfIslands.push_back(pIsland);
	}

	if (SortedListOfIslands.size() == 0) return(0);					// nothing to pack, everything is pinned down, or it's all gone to shite

	TotalArea *= borderScale;									// add on the border to Area
	float islandScale = 1.0f / sqrt(TotalArea);


	std::sort(SortedListOfIslands.begin(), SortedListOfIslands.end(), CompareIslands);

	float maxEdge = 0.0f;

	for(unsigned int Index = 0; Index < SortedListOfIslands.size(); Index++)
	{
		Island* pIsland = SortedListOfIslands[Index];
		pIsland->UVScale(islandScale);
		pIsland->GetUVBBox(MinU, MinV, MaxU, MaxV);
		float EdgeU = MaxU - MinU;
		float EdgeV = MaxV - MinV;

		if(EdgeU > maxEdge) maxEdge = EdgeU;			// any islands still have edge of over 1.0?
		if(EdgeV > maxEdge) maxEdge = EdgeV;
	}

	if(maxEdge > 1.0f)									// if so, scale everything down so longest edge == 1.0
	{
		islandScale = 1.0f / maxEdge;

		for(unsigned int Index = 0; Index < SortedListOfIslands.size(); Index++)
		{
			Island* pIsland = SortedListOfIslands[Index];
			pIsland->UVScale(islandScale);

		}
	}

	for(unsigned int Index = 0; Index < SortedListOfIslands.size(); Index++)
	{
		Island* pIsland = SortedListOfIslands[Index];
		pIsland->GetAreas(UVArea, PolygonArea);
		m_TotalUVArea += UVArea;
	}

	// was 1024.0f
	float borderPixels = 512.0f * (borderScale - 1.0f);			//128.0f * (borderScale - 1.0f);		//m_TotalUVArea * 1024.0f * (borderScale - 1.0f);

//	if(borderPixels < 1.0f) borderPixels = 1.0f;


	for(unsigned int Index = 0; Index < SortedListOfIslands.size(); Index++)
	{
		Island* pIsland = SortedListOfIslands[Index];
		pIsland->Profiles(1024, borderPixels);				// 256
		pIsland->GetUVBBox(MinU, MinV, MaxU, MaxV);
		MaxSide = MAX3(MaxSide, MaxU, MaxV);
	}

	return(MaxSide);
}



bool Pack::ProfilePack(int TextureSize)
{
	if(m_pCurrentProfile != NULL) delete[] m_pCurrentProfile;
	m_pCurrentProfile = NULL;

	m_TextureWidthHeight = TextureSize;

	while(TextureSize & 63) TextureSize++;
	m_TextureBitsWidth = TextureSize >> 6;
	m_pCurrentProfile = new uint64_t[m_TextureBitsWidth * m_TextureWidthHeight];

	int NumPixels = (m_TextureBitsWidth * m_TextureWidthHeight);
	for(int Index = 0; Index < NumPixels; Index++)
	{
		m_pCurrentProfile[Index] = 0;
	}

	for(unsigned int IslandIndex  = 0; IslandIndex < SortedListOfIslands.size(); IslandIndex++)
	{
		bool doesFit = BestFit(IslandIndex);
		if(doesFit == false) return(false);
	}

	return(true);
}


bool Pack::BestFit(int IslandIndex)
{
	Island* pIsland = SortedListOfIslands[IslandIndex];

	if(pIsland->m_profileWidth > m_TextureWidthHeight) return(false);
	if(pIsland->m_profileHeight > m_TextureWidthHeight) return(false);

//	int PotentialXPositions = (int)(m_TextureWidthHeight+8 >> 3) - pIsland->m_BitsWidth;		//m_profileWidth;
	int PotentialXPositions = (int)m_TextureWidthHeight - pIsland->m_profileWidth;
	int PotentialYPositions = (int)m_TextureWidthHeight - pIsland->m_profileHeight;

	int LowestYPos = PotentialYPositions;
	int BestXPos = -1;
	int BestYPos = -1;
	int BestFits = -1;
	int YIndex = 0;
	int XIndex = 0;


	while(PotentialXPositions > 0)
	{
		int Count = 64;
		if(PotentialXPositions < 64)
		{
			Count = PotentialXPositions;
		}
		int fits = TestFit(XIndex>>6, LowestYPos, pIsland, YIndex, Count);
		if(YIndex < LowestYPos && fits != -1)
		{
			LowestYPos = YIndex;
			BestXPos = XIndex>>6;
			BestYPos = YIndex;
			BestFits = fits;
		}

		PotentialXPositions -= 64;
		XIndex += 64;
	}

	if(BestXPos != -1)
	{
		BurnIn(BestXPos, BestYPos, BestFits, pIsland);
		return(true);
	}


	return(false);
}




int Pack::TestFit(int TestX, int BestLowest, Island* pIsland, int& YIndexOut, int NumberTests)
{
	int IslandWidth = pIsland->m_BitsWidth;
	int IslandHeight = pIsland->m_profileHeight;
	int TextureWidth = m_TextureBitsWidth;
	int TestY;
	int ProfileStartY;
	bool PrfoileFits;
	uint64_t* pProfileBits = pIsland->m_pProfileBits;
	uint64_t BadBits;
	int BestFitRotation = -1;

	int TestUp, XIndex, YIndex, Count;


	for (unsigned int RotNumber = 0; RotNumber < NumberTests; RotNumber++)
	{
		PrfoileFits = false;
		ProfileStartY = 0;
	ReTest1:

		while (ProfileStartY < BestLowest && PrfoileFits == false)
		{
			PrfoileFits = true;
			TestY = ProfileStartY;

			for (unsigned int YIndex = 0; YIndex < IslandHeight - 1; YIndex++)
			{
				uint64_t* ProfilePtr = &m_pCurrentProfile[TestY * TextureWidth + TestX];
				uint64_t* BitsPtr = &pProfileBits[YIndex * IslandWidth];
				uint64_t PreviousBits = 0;

				for (int XIndex = 0; XIndex < IslandWidth; XIndex++)
				{
					uint64_t TheseBits = PreviousBits << (63 - RotNumber) | (*BitsPtr) >> RotNumber;
					PreviousBits = (*BitsPtr);

					if (TheseBits & *ProfilePtr)				// bit clash
					{
						PrfoileFits = false;
						Count = BestLowest - ProfileStartY;

						for (TestUp = 0; TestUp < Count; TestUp++)
						{
							if (!(TheseBits & *ProfilePtr))
							{
								ProfileStartY += TestUp;
								goto ReTest1;
							}

							ProfilePtr += TextureWidth;			// = uintTextureWidth;
						}

						ProfileStartY = BestLowest;
					}

					BitsPtr++;
					ProfilePtr++;
				}
				TestY++;
			}
		}

		if (ProfileStartY < BestLowest)
		{
			YIndexOut = ProfileStartY;
			BestLowest = ProfileStartY;
			BestFitRotation = RotNumber;
		}
	}

	return(BestFitRotation);


/*
	int IslandWidth = pIsland->m_BitsWidth >> 2;
	int IslandHeight = pIsland->m_profileHeight;
	int TextureWidth = m_TextureBitsWidth;
	int uintTextureWidth = m_TextureBitsWidth >> 2;
	int TestY;

	unsigned int* pProfileBits1 = pIsland->m_pProfileBits1;
	unsigned int* pProfileBits2 = pIsland->m_pProfileBits2;
	unsigned int* pProfileBits3 = pIsland->m_pProfileBits3;
	unsigned int* pProfileBits4 = pIsland->m_pProfileBits4;
	unsigned int* pProfileBits5 = pIsland->m_pProfileBits5;
	unsigned int* pProfileBits6 = pIsland->m_pProfileBits6;
	unsigned int* pProfileBits7 = pIsland->m_pProfileBits7;
	unsigned int* pProfileBits8 = pIsland->m_pProfileBits8;

	int Profile1StartY, Profile2StartY, Profile3StartY, Profile4StartY, Profile5StartY, Profile6StartY, Profile7StartY, Profile8StartY;
	Profile1StartY = Profile2StartY = Profile3StartY = Profile4StartY = Profile5StartY = Profile6StartY = Profile7StartY = Profile8StartY = 0;

	bool P1Fits, P2Fits, P3Fits, P4Fits, P5Fits, P6Fits, P7Fits, P8Fits;
	P1Fits = P2Fits = P3Fits = P4Fits = P5Fits = P6Fits = P7Fits = P8Fits = false;

	unsigned int* BitsPtr;
	unsigned int* ProfilePtr;
	int TestUp, XIndex, YIndex, Count;

	unsigned int BadBits;


	if(NumberTests == 0) goto End;
ReTest1:
	while(Profile1StartY < BestLowest && P1Fits == false)
	{
		P1Fits = true;
		TestY = Profile1StartY;

		for(YIndex = 0; YIndex < IslandHeight-1; YIndex++)
		{
			ProfilePtr = (unsigned int*) (&m_pCurrentProfile[TestY * TextureWidth + TestX]);
			BitsPtr = &pProfileBits1[YIndex * IslandWidth];

			for(XIndex = 0; XIndex < IslandWidth; XIndex++)
			{
				if(*BitsPtr++ & *ProfilePtr++)				// bit clash
				{
					P1Fits = false;
					BadBits = BitsPtr[-1];
					ProfilePtr--;
					Count = BestLowest - Profile1StartY;

					for(TestUp = 0; TestUp < Count; TestUp++)
					{
						if(!(BadBits & *ProfilePtr))
						{
							Profile1StartY += TestUp;
							goto ReTest1;
						}

						ProfilePtr += uintTextureWidth;
					}

					Profile1StartY = BestLowest;
				}
			}
			TestY++;
		}
	}
	if(Profile1StartY < BestLowest) BestLowest = Profile1StartY;
	

	if(NumberTests == 1) goto End;
ReTest2:
	while(Profile2StartY < BestLowest && P2Fits == false)
	{
		P2Fits = true;
		TestY = Profile2StartY;

		for(YIndex = 0; YIndex < IslandHeight-1; YIndex++)
		{
			ProfilePtr = (unsigned int*) (&m_pCurrentProfile[TestY * TextureWidth + TestX]);
			BitsPtr = &pProfileBits2[YIndex * IslandWidth];

			for(XIndex = 0; XIndex < IslandWidth; XIndex++)
			{
				if(*BitsPtr++ & *ProfilePtr++)			// bit clash
				{
					P2Fits = false;
					BadBits = BitsPtr[-1];
					ProfilePtr--;
					Count = BestLowest - Profile2StartY;

					for(TestUp = 0; TestUp < Count; TestUp++)
					{
						if(!(BadBits & *ProfilePtr))
						{
							Profile2StartY += TestUp;
							goto ReTest2;
						}

						ProfilePtr += uintTextureWidth;
					}

					Profile2StartY = BestLowest;
				}
			}
			TestY++;
		}
	}
	if(Profile2StartY < BestLowest) BestLowest = Profile2StartY;


	if(NumberTests == 2) goto End;
ReTest3:
	while(Profile3StartY < BestLowest && P3Fits == false)
	{
		P3Fits = true;
		TestY = Profile3StartY;

		for(YIndex = 0; YIndex < IslandHeight-1; YIndex++)
		{
			ProfilePtr = (unsigned int*) (&m_pCurrentProfile[TestY * TextureWidth + TestX]);
			BitsPtr = &pProfileBits3[YIndex * IslandWidth];

			for(XIndex = 0; XIndex < IslandWidth; XIndex++)
			{
				if(*BitsPtr++ & *ProfilePtr++)			// bit clash
				{
					P3Fits = false;
					ProfilePtr--;
					Count = BestLowest - Profile3StartY;
					BadBits = BitsPtr[-1];

					for(TestUp = 0; TestUp < Count; TestUp++)
					{
						if(!(BadBits & *ProfilePtr)) 
						{
							Profile3StartY += TestUp;
							goto ReTest3;
						}

						ProfilePtr += uintTextureWidth;
					}

					Profile3StartY = BestLowest;
				}
			}
			TestY++;
		}
	}
	if(Profile3StartY < BestLowest) BestLowest = Profile3StartY;


	if(NumberTests == 3) goto End;
ReTest4:
	while(Profile4StartY < BestLowest && P4Fits == false)
	{
		P4Fits = true;
		TestY = Profile4StartY;

		for(YIndex = 0; YIndex < IslandHeight-1; YIndex++)
		{
			ProfilePtr = (unsigned int*) (&m_pCurrentProfile[TestY * TextureWidth + TestX]);
			BitsPtr = &pProfileBits4[YIndex * IslandWidth];

			for(XIndex = 0; XIndex < IslandWidth; XIndex++)
			{
				if(*BitsPtr++ & *ProfilePtr++)		// bit clash
				{
					P4Fits = false;
					BadBits = BitsPtr[-1];
					ProfilePtr--;
					Count = BestLowest - Profile4StartY;

					for(TestUp = 0; TestUp < Count; TestUp++)
					{
						if(!(BadBits & *ProfilePtr))
						{
							Profile4StartY += TestUp;
							goto ReTest4;
						}

						ProfilePtr += uintTextureWidth;
					}
					
					Profile4StartY = BestLowest;
				}
			}
			TestY++;
		}
	}
	if(Profile4StartY < BestLowest) BestLowest = Profile4StartY;


	if(NumberTests == 4) goto End;
ReTest5:
	while(Profile5StartY < BestLowest && P5Fits == false)
	{
		P5Fits = true;
		TestY = Profile5StartY;

		for(YIndex = 0; YIndex < IslandHeight-1; YIndex++)
		{
			ProfilePtr = (unsigned int*) (&m_pCurrentProfile[TestY * TextureWidth + TestX]);
			BitsPtr = &pProfileBits5[YIndex * IslandWidth];

			for(XIndex = 0; XIndex < IslandWidth; XIndex++)
			{	
				if(*BitsPtr++ & *ProfilePtr++)			// bit clash
				{
					P5Fits = false;
					ProfilePtr--;
					BadBits = BitsPtr[-1];
					Count = BestLowest - Profile5StartY;
					
					for(TestUp = 0; TestUp < Count; TestUp++)
					{
						if(!(BadBits & *ProfilePtr)) 
						{
							Profile5StartY += TestUp;
							goto ReTest5;
						}

						ProfilePtr += uintTextureWidth;
					}

					Profile5StartY = BestLowest;
				}
			}
			TestY++;
		}
	}
	if(Profile5StartY < BestLowest) BestLowest = Profile5StartY;


	if(NumberTests == 5) goto End;
ReTest6:
	while(Profile6StartY < BestLowest && P6Fits == false)
	{
		P6Fits = true;
		TestY = Profile6StartY;

		for(YIndex = 0; YIndex < IslandHeight-1; YIndex++)
		{
			ProfilePtr = (unsigned int*) (&m_pCurrentProfile[TestY * TextureWidth + TestX]);
			BitsPtr = &pProfileBits6[YIndex * IslandWidth];

			for(XIndex = 0; XIndex < IslandWidth; XIndex++)
			{
				if(*BitsPtr++ & *ProfilePtr++)			// bit clash
				{
					P6Fits = false;
					ProfilePtr--;
					BadBits = BitsPtr[-1];
					Count = BestLowest - Profile6StartY;

					for(TestUp = 0; TestUp < Count; TestUp++)
					{
						if(!(BadBits & *ProfilePtr))
						{
							Profile6StartY += TestUp;
							goto ReTest6;
						}

						ProfilePtr += uintTextureWidth;
					}

					Profile6StartY = BestLowest;
				}
			}
			TestY++;
		}
	}
	if(Profile6StartY < BestLowest) BestLowest = Profile6StartY;


	if(NumberTests == 6) goto End;
ReTest7:
	while(Profile7StartY < BestLowest && P7Fits == false)
	{
		P7Fits = true;
		TestY = Profile7StartY;

		for(YIndex = 0; YIndex < IslandHeight-1; YIndex++)
		{
			ProfilePtr = (unsigned int*) (&m_pCurrentProfile[TestY * TextureWidth + TestX]);
			BitsPtr = &pProfileBits7[YIndex * IslandWidth];

			for(XIndex = 0; XIndex < IslandWidth; XIndex++)
			{
				if(*BitsPtr++ & *ProfilePtr++)			// bit clash
				{
					P7Fits = false;
					ProfilePtr--;
					BadBits = BitsPtr[-1];
					Count = BestLowest - Profile7StartY;

					for(TestUp = 0; TestUp < Count; TestUp++)
					{
						if(!(BadBits & *ProfilePtr))
						{
							Profile7StartY += TestUp;
							goto ReTest7;
						}

						ProfilePtr += uintTextureWidth;
					}

					Profile7StartY = BestLowest;
				}
			}
			TestY++;
		}
	}
	if(Profile7StartY < BestLowest) BestLowest = Profile7StartY;


	if(NumberTests == 7) goto End;
ReTest8:
	while(Profile8StartY < BestLowest && P8Fits == false)
	{
		P8Fits = true;
		TestY = Profile8StartY;

		for(YIndex = 0; YIndex < IslandHeight-1; YIndex++)
		{
			ProfilePtr = (unsigned int*) (&m_pCurrentProfile[TestY * TextureWidth + TestX]);
			BitsPtr = &pProfileBits8[YIndex * IslandWidth];

			for(XIndex = 0; XIndex < IslandWidth; XIndex++)
			{
				if(*BitsPtr++ & *ProfilePtr++)			// bit clash
				{
					P8Fits = false;
					ProfilePtr--;
					BadBits = BitsPtr[-1];
					Count = BestLowest - Profile8StartY;

					for(TestUp = 0; TestUp < Count; TestUp++)
					{
						if(!(BadBits & *ProfilePtr))
						{
							Profile8StartY += TestUp;
							goto ReTest8;
						}
						
						ProfilePtr += uintTextureWidth;
					}

					Profile8StartY = BestLowest;
				}
			}
			TestY++;
		}
	}
	if(Profile8StartY < BestLowest) BestLowest = Profile8StartY;

End:
	if(BestLowest == Profile1StartY && P1Fits) { YIndexOut = Profile1StartY;  return(0); }
	if(BestLowest == Profile2StartY && P2Fits) { YIndexOut = Profile2StartY;  return(1); }
	if(BestLowest == Profile3StartY && P3Fits) { YIndexOut = Profile3StartY;  return(2); }
	if(BestLowest == Profile4StartY && P4Fits) { YIndexOut = Profile4StartY;  return(3); }
	if(BestLowest == Profile5StartY && P5Fits) { YIndexOut = Profile5StartY;  return(4); }
	if(BestLowest == Profile6StartY && P6Fits) { YIndexOut = Profile6StartY;  return(5); }
	if(BestLowest == Profile7StartY && P7Fits) { YIndexOut = Profile7StartY;  return(6); }
	if(BestLowest == Profile8StartY && P8Fits) { YIndexOut = Profile8StartY;  return(7); }
*/

	return(-1);			// fits
}





void Pack::BurnIn(int PosX, int PosY, int fits, Island* pIsland)
{
	uint64_t* pProfilePixels;
	uint64_t* pCurrentProfile;

	int IslandWidth = pIsland->m_BitsWidth;					//m_profileWidth;
	pProfilePixels = pIsland->m_pProfileBits;

	pIsland->m_pixelX = PosX * 64 + fits;
	pIsland->m_pixelY = PosY;


	int IslandHeight = pIsland->m_profileHeight;

	for (int YIndex = 0; YIndex < IslandHeight - 1; YIndex++)
	{
		int SrcYAdder = PosY * m_TextureBitsWidth;
		//int DstYAdder = YIndex * (IslandWidth >> 2);

		pCurrentProfile = &(m_pCurrentProfile[PosX + SrcYAdder]);
		uint64_t PreviousBits = 0;

		for (int XIndex = 0; XIndex < IslandWidth; XIndex++)
		{
			uint64_t TheseBits = PreviousBits << (63 - fits) | (*pProfilePixels) >> fits;
			PreviousBits = *pProfilePixels;

			*pCurrentProfile |= TheseBits;
			pCurrentProfile++;
			pProfilePixels++;
		}

		PosY++;
	}
}










// Third time lucky?

void Pack::GetUpwards3(Island* pIsland)
{
	vector<Vec2D> VecsX;
	vector<Vec2D> VecsY;
	vector<Vec2D> VecsZ;

	float VecXLineX = 0.0f;
	float VecYLineX = 0.0f;
	float VecXLineY = 0.0f;
	float VecYLineY = 0.0f;
	float VecXLineZ = 0.0f;
	float VecYLineZ = 0.0f;

	if(pIsland->PinFinder.m_Boundaries.size() == 0)
	{
		pIsland->MinAreaRotate();
		return;
	}


	for(unsigned int FaceIndex = 0; FaceIndex < pIsland->m_Faces.size(); FaceIndex++)
	{
		RKFace* pFace = pIsland->m_Faces[FaceIndex];

		float MidX = (pFace->pVert1->X + pFace->pVert2->X  + pFace->pVert3->X) / 3.0f;
		float MidY = (pFace->pVert1->Y + pFace->pVert2->Y  + pFace->pVert3->Y) / 3.0f;
		float MidZ = (pFace->pVert1->Z + pFace->pVert2->Z  + pFace->pVert3->Z) / 3.0f;
		float XPlaneU1, XPlaneV1, XPlaneU2, XPlaneV2;
		float YPlaneU1, YPlaneV1, YPlaneU2, YPlaneV2;
		float ZPlaneU1, ZPlaneV1, ZPlaneU2, ZPlaneV2;

		GetUVsForPlane(pFace, MidX, XPlaneU1, XPlaneV1, XPlaneU2, XPlaneV2, 1);
		GetUVsForPlane(pFace, MidY, YPlaneU1, YPlaneV1, YPlaneU2, YPlaneV2, 2);
		GetUVsForPlane(pFace, MidZ, ZPlaneU1, ZPlaneV1, ZPlaneU2, ZPlaneV2, 3);

		// or maybe all these points need to go into 3 linear regressors?

		Vec2D XLine;
		Vec2D YLine;
		Vec2D ZLine;

		XLine.X = XPlaneU2 - XPlaneU1; XLine.Y = XPlaneV2 - XPlaneV1;
		YLine.X = YPlaneU2 - YPlaneU1; YLine.Y = YPlaneV2 - YPlaneV1;
		ZLine.X = ZPlaneU2 - ZPlaneU1; ZLine.Y = ZPlaneV2 - ZPlaneV1;

//		if(XLine.x < 0.0f) { XLine.x = -XLine.x;  XLine.y = -XLine.y; }
//		if(YLine.x < 0.0f) { YLine.x = -YLine.x;  YLine.y = -YLine.y; }
//		if(ZLine.x < 0.0f) { ZLine.x = -ZLine.x;  ZLine.y = -ZLine.y; }

		float lineXLen = sqrt((XLine.X * XLine.X)  + (XLine.Y * XLine.Y));
		float lineYLen = sqrt((YLine.X * YLine.X)  + (YLine.Y * YLine.Y));
		float lineZLen = sqrt((ZLine.X * ZLine.X)  + (ZLine.Y * ZLine.Y));

		if(lineXLen != 0.0f)
		{
			XLine.X /= lineXLen;				// normalise
			XLine.Y /= lineXLen;
//			VecXLineX += XLine.x;
//			VecYLineX += XLine.y;
//			LeastSquaresXPlane.add_sample_point(XLine);
			VecsX.push_back(XLine);
		}

		if(lineYLen != 0.0f)
		{
			YLine.X /= lineYLen;
			YLine.Y /= lineYLen;
//			VecXLineY += YLine.x;
//			VecYLineY += YLine.y;
//			LeastSquaresYPlane.add_sample_point(YLine);
			VecsY.push_back(YLine);
		}

		if(lineZLen != 0.0f)
		{
			ZLine.X /= lineZLen;
			ZLine.Y /= lineZLen;
//			VecXLineZ += ZLine.x;
//			VecYLineZ += ZLine.y;
//			LeastSquaresZPlane.add_sample_point(ZLine);
			VecsZ.push_back(ZLine);
		}
	}
/*
	float LenXLine = sqrt((VecXLineX * VecXLineX) + (VecYLineX * VecYLineX));
	float LenYLine = sqrt((VecXLineY * VecXLineY) + (VecYLineY * VecYLineY));
	float LenZLine = sqrt((VecXLineZ * VecXLineZ) + (VecYLineZ * VecYLineZ));

	VecXLineX /= LenXLine;  VecYLineX /= LenXLine;
	VecXLineY /= LenYLine;  VecYLineY /= LenYLine;
	VecXLineZ /= LenZLine;  VecYLineZ /= LenZLine;
*/

	FindMeanLine(VecsX, VecXLineX, VecYLineX);
	FindMeanLine(VecsY, VecXLineY, VecYLineY);
	FindMeanLine(VecsZ, VecXLineZ, VecYLineZ);

	float ErrorX = 0.0f;
	float ErrorY = 0.0f;
	float ErrorZ = 0.0f;

	// what to do if there's only one face in a plane?
	for(unsigned int VecXIndex = 0; VecXIndex < VecsX.size(); VecXIndex++)
	{
		float Dot = VecXLineX * VecsX[VecXIndex].X + VecYLineX * VecsX[VecXIndex].Y;
		if(Dot < 0) Dot = -Dot;
		if(Dot > 1.0f) Dot = 1.0f;
		ErrorX += acos(Dot);
	}
	if(VecsX.size() == 0) ErrorX = 1000000.0f;

	for(unsigned int VecYIndex = 0; VecYIndex < VecsY.size(); VecYIndex++)
	{
		float Dot = VecXLineY * VecsY[VecYIndex].X + VecYLineY * VecsY[VecYIndex].Y;
		if(Dot < 0) Dot = -Dot;
		if(Dot > 1.0f) Dot = 1.0f;
		ErrorY += acos(Dot);
	}
	if(VecsY.size() == 0) ErrorY = 1000000.0f;

	for(unsigned int VecZIndex = 0; VecZIndex < VecsZ.size(); VecZIndex++)
	{
		float Dot = VecXLineZ * VecsZ[VecZIndex].X + VecYLineZ * VecsZ[VecZIndex].Y;
		if(Dot < 0) Dot = -Dot;
		if(Dot > 1.0f) Dot = 1.0f;
		ErrorZ += acos(Dot);
	}
	if(VecsZ.size() == 0) ErrorZ = 1000000.0f;

	float RetAngle = 0.0f;
	if(ErrorX <= ErrorY && ErrorX <= ErrorZ)
	{
//		RetAngle = acos(VecYLineX);			// - M_PI * 0.5f;
		RetAngle = -atan2(VecYLineX, VecXLineX);
		RetAngle -= (float)M_PI;				// * 0.5f;
	}

	if(ErrorY <= ErrorX && ErrorY <= ErrorZ)
	{
		//RetAngle = acos(VecYLineY);
		RetAngle = -atan2(VecYLineY, VecXLineY);
		RetAngle -= (float)M_PI;			// * 0.5f;
	}

	if(ErrorZ <= ErrorX && ErrorZ <= ErrorY)
	{
		//RetAngle = acos(VecYLineZ);
		RetAngle = -atan2(VecYLineZ, VecXLineZ);
		RetAngle -= (float)M_PI;		// * 0.5f;				// good for cube
	}

	// Plane with least angle of error is the winner, rotate UVs so they're parrallel with that plane
	// going to have to take the number of faces into account.  Because one face can lead to that being the upright plane.

	pIsland->Rotate(RetAngle - (float)M_PI * 0.5f);

	FinalRotate(pIsland);

	return;
}


void Pack::FinalRotate(Island* pIsland)
{
//	float BestRotate = 0.0f;
	float MinX = 1000000.0f;
	float MinY = 1000000.0f;
	float MinZ = 1000000.0f;
	float MaxX = -1000000.0f;
	float MaxY = -1000000.0f;
	float MaxZ = -1000000.0f;
	RKVertex *VertMinX = NULL;
	RKVertex *VertMinY = NULL;
	RKVertex *VertMinZ = NULL;
	RKVertex *VertMaxX = NULL;
	RKVertex *VertMaxY = NULL;
	RKVertex *VertMaxZ = NULL;
//	RKVertex* pVertUV1, pVertUV2;

	vector<RKVertex*>::iterator iter = pIsland->m_Vertices.begin();

	while(iter != pIsland->m_Vertices.end())
	{
		RKVertex* pVert = *iter;
		iter++;

		if(pVert->X < MinX) { MinX = pVert->X; VertMinX = pVert; }
		if(pVert->X > MaxX) { MaxX = pVert->X; VertMaxX = pVert; }
		if(pVert->Y < MinY) { MinY = pVert->Y; VertMinY = pVert; }
		if(pVert->Y > MaxY) { MaxY = pVert->Y; VertMaxY = pVert; }
		if(pVert->Z < MinZ) { MinZ = pVert->Z; VertMinZ = pVert; }
		if(pVert->Z > MaxZ) { MaxZ = pVert->Z; VertMaxZ = pVert; }
	}

	float DeltaX = MaxX - MinX;
	float DeltaY = MaxY - MinY;
	float DeltaZ = MaxZ - MinZ;
	float RotAngle = 0.0f;

	if(DeltaY > (DeltaX * 0.25) && DeltaY > (DeltaZ * 0.1f))
	{
		RotAngle = FindBestRot(1, VertMaxY, VertMinY, 2);			// top to bottom, find rot which maximises V Delta
	}
	else
	{
		if((DeltaX * 0.25) > (DeltaZ * 0.1f))
		{
			RotAngle = FindBestRot(2, VertMinX, VertMaxX, 1);		// left to right, find rot which maximises U Delta
		}
		else
		{
			RotAngle = FindBestRot(3, VertMinZ, VertMaxZ, 2);		// near to far, find rot which maximises V Delta
		}
	}

	pIsland->Rotate(RotAngle);
}



float Pack::FindBestRot(int TestDelta, RKVertex* pVertMin, RKVertex* pVertMax, int testWidthHeight)		// 1 = U, 2 = V
{
	float BestRotate = 0.0f;
	float BestDelta = 0.0f;
	float ThisAngle = 0.0f;
	float WorldDelta = 0.0f;

	if(TestDelta == 2) WorldDelta = pVertMax->X - pVertMin->X;
	if(TestDelta == 1) WorldDelta = pVertMax->Y - pVertMin->Y;
	if(TestDelta == 3) WorldDelta = pVertMax->Z - pVertMin->Z;


	for(int RotIndex = 0; RotIndex < 4; RotIndex++)
	{
		float Sine = sin(ThisAngle);
		float Cosine = cos(ThisAngle);
		float testUTop, testVTop, testUBot, testVBot;

		testUTop = Cosine*pVertMin->U - Sine*pVertMin->V;
		testVTop = Sine*pVertMin->U + Cosine*pVertMin->V;
		testUBot = Cosine*pVertMax->U - Sine*pVertMax->V;
		testVBot = Sine*pVertMax->U + Cosine*pVertMax->V;

		float ThisDeltaUV = 0.0f;

		if(testWidthHeight == 1)
		{
			ThisDeltaUV = testUBot - testUTop;
		}
		else
		{
			ThisDeltaUV = testVBot - testVTop;
		}

		if(ThisDeltaUV * WorldDelta > BestDelta)
		{
			BestDelta = ThisDeltaUV * WorldDelta;
			BestRotate = ThisAngle;
		}

		ThisAngle += (float)M_PI * 0.5f;
	}

	return(BestRotate);
}





void Pack::FindMeanLine(vector<Vec2D> &Vecs, float& VecOutX, float& VecOutY)
{
	if(Vecs.size() == 0) return;

	double BestFitX = Vecs[0].X;		// first best fit is the first Normal
	double BestFitY = Vecs[0].Y;


	for(int Iteration = 0; Iteration < 5; Iteration++)
	{
		double VecX = 0.0;
		double VecY = 0.0;
		unsigned int FlipCount = 0;

		for(unsigned int VecIndex = 0; VecIndex < Vecs.size(); VecIndex++)
		{
			double NormVecX = Vecs[VecIndex].X;
			double NormVecY = Vecs[VecIndex].Y;

			if(VecX == 0.0f && VecY == 0.0f)
			{
				VecX = NormVecX;
				VecY = NormVecY;
			}
			else
			{
				float DotResult = (float)((BestFitX * NormVecX) + (BestFitY * NormVecY));
				if(DotResult < 0.0f)
				{
					NormVecX = -NormVecX;
					NormVecY = -NormVecY;
					FlipCount++;
				}

				VecX += NormVecX;
				VecY += NormVecY;
			}
		}

		if(FlipCount > Vecs.size() / 2)
		{
			VecX = -VecX;
			VecY = -VecY;
		}

		double VecLen = sqrt(VecX * VecX + VecY * VecY);
		if(VecLen != 0.0)
		{
			VecX /= VecLen;
			VecY /= VecLen;
		}

		double Dot = VecX * BestFitX + VecY * BestFitY;
//		double error = acos(Dot);

		BestFitX = VecX;
		BestFitY = VecY;
		VecOutX = (float)VecX;
		VecOutY = (float)VecY;

		if(Dot > 0.999) return;			// less than 2 degrees difference from last vector
	}
}



void Pack::GetUVsForPlane(RKFace* pFace, float& MidPoint, float& U1Out, float& V1Out, float& U2Out, float& V2Out, int Plane)
{
	U1Out = V1Out = U2Out = V2Out = 0.0f;
	bool firstFound = false;

	if(Plane == 1)
	{
		float X1 = pFace->pVert1->X;  float X2 = pFace->pVert2->X;  float X3 = pFace->pVert3->X;
		if(MidPoint > X1 && MidPoint <= X2 || MidPoint > X2 && MidPoint <= X1) DeltaCalc(pFace->pVert1, pFace->pVert2, MidPoint, U1Out, V1Out, U2Out, V2Out, firstFound, 1);
		if(MidPoint > X2 && MidPoint <= X3 || MidPoint > X3 && MidPoint <= X2) DeltaCalc(pFace->pVert2, pFace->pVert3, MidPoint, U1Out, V1Out, U2Out, V2Out, firstFound, 1);
		if(MidPoint > X3 && MidPoint <= X1 || MidPoint > X1 && MidPoint <= X3) DeltaCalc(pFace->pVert3, pFace->pVert1, MidPoint, U1Out, V1Out, U2Out, V2Out, firstFound, 1);
	}

	if(Plane == 2)
	{
		float Y1 = pFace->pVert1->Y;  float Y2 = pFace->pVert2->Y;  float Y3 = pFace->pVert3->Y;
		if(MidPoint > Y1 && MidPoint <= Y2 || MidPoint > Y2 && MidPoint <= Y1) DeltaCalc(pFace->pVert1, pFace->pVert2, MidPoint, U1Out, V1Out, U2Out, V2Out, firstFound, 2);
		if(MidPoint > Y2 && MidPoint <= Y3 || MidPoint > Y3 && MidPoint <= Y2) DeltaCalc(pFace->pVert2, pFace->pVert3, MidPoint, U1Out, V1Out, U2Out, V2Out, firstFound, 2);
		if(MidPoint > Y3 && MidPoint <= Y1 || MidPoint > Y1 && MidPoint <= Y3) DeltaCalc(pFace->pVert3, pFace->pVert1, MidPoint, U1Out, V1Out, U2Out, V2Out, firstFound, 2);
	}

	if(Plane == 3)
	{
		float Z1 = pFace->pVert1->Z;  float Z2 = pFace->pVert2->Z;  float Z3 = pFace->pVert3->Z;
		if(MidPoint > Z1 && MidPoint <= Z2 || MidPoint > Z2 && MidPoint <= Z1) DeltaCalc(pFace->pVert1, pFace->pVert2, MidPoint, U1Out, V1Out, U2Out, V2Out, firstFound, 3);
		if(MidPoint > Z2 && MidPoint <= Z3 || MidPoint > Z3 && MidPoint <= Z2) DeltaCalc(pFace->pVert2, pFace->pVert3, MidPoint, U1Out, V1Out, U2Out, V2Out, firstFound, 3);
		if(MidPoint > Z3 && MidPoint <= Z1 || MidPoint > Z1 && MidPoint <= Z3) DeltaCalc(pFace->pVert3, pFace->pVert1, MidPoint, U1Out, V1Out, U2Out, V2Out, firstFound, 3);
	}
}



void Pack::DeltaCalc(RKVertex* pVert1, RKVertex* pVert2, float& MidPoint, float& U1Out, float& V1Out, float& U2Out, float& V2Out, bool& FoundFirst, int Plane)
{
	float Start = 0.0f;
	float End = 0.0f;

	if(Plane == 1) { Start = pVert1->X; End = pVert2->X; }
	if(Plane == 2) { Start = pVert1->Y; End = pVert2->Y; }
	if(Plane == 3) { Start = pVert1->Z; End = pVert2->Z; }

	float ratio = ((MidPoint - Start) / (End - Start));

	if(FoundFirst == false)
	{
		U1Out = ((pVert2->U - pVert1->U) * ratio) + pVert1->U;
		V1Out = ((pVert2->V - pVert1->V) * ratio) + pVert1->V;
		FoundFirst = true;
	}
	else
	{
		U2Out = ((pVert2->U - pVert1->U) * ratio) + pVert1->U;
		V2Out = ((pVert2->V - pVert1->V) * ratio) + pVert1->V;
	}
}


void Pack::MaintainUVs(vector<Island*>& listOfIslands)
{
	for (unsigned int Index = 0; Index < listOfIslands.size(); Index++)
	{
		Island* pIsland = listOfIslands[Index];

		pIsland->MaintainUVs();
	}
}