#ifndef ISLANDMATCH_INCLUDED
#define ISLANDMATCH_INCLUDED

#include <vector>
#include "Island.h"


class IslandMatch
{
public:
	IslandMatch() {};
	~IslandMatch() {};

	void MatchIslands(vector<Island*> listOfIslands, vector<Island*> &matchedIslands);
	void Propagate(vector<Island*> listOfIslands);


private:


};



#endif			// endif ISLANDMATCH