#include <math.h>

#include "RKProgress.h"

#include "MinimiseStretch.h"


extern int gScale;						// used as the iteration counter

#define P_STRETCH_ITER 20




void MinimiseStretch::Minimise(vector<RKVertex*> &rVertices)
{
	float Progress = 0.0f;
	float ProgressAdded = 100.0f / (float)gScale;

	for(int Index = 0; Index < gScale; Index++)
	{
		StretchIteration(rVertices);
		RKProgress::Get().SetProgress((int)Progress);
		Progress += ProgressAdded;
	}
}



float MinimiseStretch::GetFaceStretch(RKFace* pFace)
{
	float T, w;		//, tmp[3];
	float Ps[3], Pt[3];
	float a, c, area, delta;
//	PEdge *e1 = f->edge, *e2 = e1->next, *e3 = e2->next;
//	PVert *v1 = e1->vert, *v2 = e2->vert, *v3 = e3->vert;

	area = pFace->GetUVArea();					//p_face_uv_area_signed(f);

	if (area <= 0.0f)							// flipped face -> infinite stretch
	{
		return 1e10f;
	}

	w = 1.0f /(2.0f*area);

	/* compute derivatives */
//	VecCopyf(Ps, v1->co);
//	VecMulf(Ps, (v2->uv[1] - v3->uv[1]));
	delta = pFace->pVert2->V - pFace->pVert3->V;
	Ps[0] = pFace->pVert1->X * delta;  Ps[1] = pFace->pVert1->Y * delta;  Ps[2] = pFace->pVert1->Z * delta;

//	VecCopyf(tmp, v2->co);
//	VecMulf(tmp, (v3->uv[1] - v1->uv[1]));
//	VecAddf(Ps, Ps, tmp);
	delta = pFace->pVert3->V - pFace->pVert1->V;
	Ps[0] += pFace->pVert2->X * delta;  Ps[1] += pFace->pVert2->Y * delta;  Ps[2] += pFace->pVert2->Z * delta;

//	VecCopyf(tmp, v3->co);
//	VecMulf(tmp, (v1->uv[1] - v2->uv[1]));
//	VecAddf(Ps, Ps, tmp);
	delta = pFace->pVert1->V - pFace->pVert2->V;
	Ps[0] += pFace->pVert3->X * delta;  Ps[1] += pFace->pVert3->Y * delta;  Ps[2] += pFace->pVert3->Z * delta;

//	VecMulf(Ps, w);
	Ps[0] *= w;  Ps[1] *= w;  Ps[2] *= w;


//	VecCopyf(Pt, v1->co);
//	VecMulf(Pt, (v3->uv[0] - v2->uv[0]));
	delta = pFace->pVert3->U - pFace->pVert2->U;
	Pt[0] = pFace->pVert1->X * delta;  Pt[1] = pFace->pVert1->Y * delta;  Pt[2] = pFace->pVert1->Z * delta;

//	VecCopyf(tmp, v2->co);
//	VecMulf(tmp, (v1->uv[0] - v3->uv[0]));
//	VecAddf(Pt, Pt, tmp);
	delta = pFace->pVert1->U - pFace->pVert3->U;
	Pt[0] += pFace->pVert2->X * delta;  Pt[1] += pFace->pVert2->Y * delta;  Pt[2] += pFace->pVert2->Z * delta;

//	VecCopyf(tmp, v3->co);
//	VecMulf(tmp, (v2->uv[0] - v1->uv[0]));
//	VecAddf(Pt, Pt, tmp);
	delta = pFace->pVert2->U - pFace->pVert1->U;
	Pt[0] += pFace->pVert3->X * delta;  Pt[1] += pFace->pVert3->Y * delta;  Pt[2] += pFace->pVert3->Z * delta;

//	VecMulf(Pt, w);
	Pt[0] *= w;  Pt[1] *= w;  Pt[2] *= w;

	/* Sander Tensor */
//	a = Inpf(Ps, Ps);
	a = Ps[0]*Ps[0] + Ps[1]*Ps[1] + Ps[2]*Ps[2];
//	c = Inpf(Pt, Pt);
	c = Pt[0]*Pt[0] + Pt[1]*Pt[1] + Pt[2]*Pt[2];

	T =  sqrt(0.5f*(a + c));

//	if (f->flag & PFACE_FILLED)
	if(pFace->m_Filler)
	{
		T *= 0.2f;
	}

	return T;
}



float MinimiseStretch::GetStretchVertex(RKVertex* pVert)
{
	float Sum = 0.0f;

	list<RKFace*>::iterator FaceIterator;
	for(FaceIterator = pVert->m_LinkedFaces.begin(); FaceIterator != pVert->m_LinkedFaces.end(); FaceIterator++)
//	for(int faceIndex = 0; faceIndex < pVert->m_LinkedFaces.size(); faceIndex++)
	{
		Sum += GetFaceStretch(*FaceIterator);
	}

	return(Sum);

/*
	PEdge *e = v->edge;
	float sum = 0.0f;

	do {
		sum += p_face_stretch(e->face);
		e = p_wheel_edge_next(e);
	} while (e && e != (v->edge));

	return sum;
*/
}






void MinimiseStretch::StretchIteration(vector<RKVertex*> &rVertices)
{
//	PVert *v;
//	PEdge *e;
	int j, nedges;
	float orig_stretch, low, stretch_low, high, stretch_high, mid, stretch;
	float orig_uv[2], dir[2], random_angle, trusted_radius;

//	int pinnedCount = 0;

//	for(v=chart->verts; v; v=v->nextlink)
	for(int vertIndex = 0; vertIndex < rVertices.size(); vertIndex++)
	{
		RKVertex* pThisVert = rVertices[vertIndex];

//		if((v->flag & PVERT_PIN)) pinnedCount++;

		//if((v->flag & PVERT_PIN))
		//	continue;

		if(pThisVert->m_OnEdge) continue;					// do not stretch outer boundary

		orig_stretch = GetStretchVertex(pThisVert);			//p_stretch_compute_vertex(v);
		orig_uv[0] = pThisVert->U;			//v->uv[0];
		orig_uv[1] = pThisVert->V;			//v->uv[1];

		/* move vertex in a random direction */
		trusted_radius = 0.0f;
		nedges = 0;
		//e = v->edge;


		list<RKFace*>::iterator FaceIterator;
		for(FaceIterator = pThisVert->m_LinkedFaces.begin(); FaceIterator != pThisVert->m_LinkedFaces.end(); FaceIterator++)
//		for(int faceIndex = 0; faceIndex < pThisVert->m_LinkedFaces.size(); faceIndex++)
		{
			RKFace* pFace = *FaceIterator;			//pThisVert->m_LinkedFaces[faceIndex];
			trusted_radius += pFace->GetUVLengthsConnected(pThisVert);
			nedges += 2;		//++;			// do this better, list of edges?
		}
/*
		do
		{
			trusted_radius += p_edge_uv_length(e);
			nedges++;

			e = p_wheel_edge_next(e);
		} while (e && e != (v->edge));
*/
		trusted_radius /= 2 * nedges;

		random_angle = (float) rand()/RAND_MAX;
		random_angle = random_angle * 2.0f * (float)M_PI;

		//random_angle = rng_getFloat(rng) * 2.0f * (float)M_PI;
		dir[0] = trusted_radius * cos(random_angle);
		dir[1] = trusted_radius * sin(random_angle);

		/* calculate old and new stretch */
		low = 0;
		stretch_low = orig_stretch;

		pThisVert->U += dir[0];
		pThisVert->V += dir[1];
//		Vec2Addf(v->uv, orig_uv, dir);
		high = 1;
		stretch = stretch_high = GetStretchVertex(pThisVert);			//v);

		/* binary search for lowest stretch position */
		for (j = 0; j < P_STRETCH_ITER; j++)
		{
			mid = 0.5f * (low + high);
			pThisVert->U = orig_uv[0] + mid*dir[0];
			pThisVert->V = orig_uv[1] + mid*dir[1];
			stretch = GetStretchVertex(pThisVert);		//v);

			if (stretch_low < stretch_high)
			{
				high = mid;
				stretch_high = stretch;
			}
			else
			{
				low = mid;
				stretch_low = stretch;
			}
		}

		/* no luck, stretch has increased, reset to old values */
		if(stretch >= orig_stretch)
		{
			pThisVert->U = orig_uv[0];
			pThisVert->V = orig_uv[1];
			//Vec2Copyf(v->uv, orig_uv);
		}
	}
}
