#include "IslandMatch.h"



void IslandMatch::MatchIslands(vector<Island*> listOfIslands, vector<Island*> &matchedIslands)
{
	matchedIslands.push_back(listOfIslands[0]);

	for(int SrcIslandIndex = 0; SrcIslandIndex < matchedIslands.size(); SrcIslandIndex++)
	{
		listOfIslands[SrcIslandIndex]->TestInteriorSymmetry();
	}

	for(int DestIslandIndex = 1; DestIslandIndex < listOfIslands.size(); DestIslandIndex++)
	{
		bool Matched = false;

		for(int SrcIslandIndex = 0; SrcIslandIndex < matchedIslands.size(); SrcIslandIndex++)
		{
			Matched = matchedIslands[SrcIslandIndex]->Compare(listOfIslands[DestIslandIndex]);
			if(Matched) break;
		}

		if(Matched == false)
		{
			matchedIslands.push_back(listOfIslands[DestIslandIndex]);
		}
	}
}



void IslandMatch::Propagate(vector<Island*> listOfIslands)
{
	for(int Index = 0; Index < listOfIslands.size(); Index++)
	{
		listOfIslands[Index]->PropagateUVs();
	}
}



