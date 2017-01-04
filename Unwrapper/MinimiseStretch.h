#ifndef MINIMISE_INCLUDED
#define MINIMISE_INCLUDED


#include "Primitives.h"


class MinimiseStretch
{
public:
	MinimiseStretch() {};
	~MinimiseStretch() {};

	void Minimise(vector<RKVertex*> &rVertices);

private:

	void StretchIteration(vector<RKVertex*> &rVertices);
	float GetFaceStretch(RKFace* pFace);
	float GetStretchVertex(RKVertex* pVert);
};



#endif		// MINIMISE_INCLUDED
