#ifndef VECTOR4_INCLUDED
#define VECTOR4_INCLUDED

#include <vector>

using namespace std;

//#include "RKEdge.h"


class Vector4
{
public:

	float	m_X, m_Y, m_Z, m_W;

	Vector4()
	{
		m_X = 0.0f;
		m_Y = 0.0f;
		m_Z = 0.0f;
//		m_Count = 0;
	};
	Vector4(float XIn, float YIn, float ZIn, float WIn)
	{
		m_X = XIn; m_Y = YIn; m_Z = ZIn; m_W = WIn;
	}

	Vector4(Vector4* pVectorIn)
	{
		m_X = pVectorIn->m_X; m_Y = pVectorIn->m_Y;  m_Z = pVectorIn->m_Z;  m_W = pVectorIn->m_W;
	}

	Vector4(Vector4& rVectorIn)
	{
		m_X = rVectorIn.m_X; m_Y = rVectorIn.m_Y;  m_Z = rVectorIn.m_Z;  m_W = rVectorIn.m_W;
	}

	~Vector4()
	{
//		m_LinkedEdges.clear();
	};

/*
	vector <int> m_LinkedEdges;
	float m_Dijikstra;
	bool m_Visited;
	int m_Count;
*/

	void Set(float XIn, float YIn, float ZIn, float WIn)
	{
		m_X = XIn; m_Y = YIn; m_Z = ZIn; m_W = WIn;
	}

	static float CrossProduct(Vector4 &Vec1, Vector4& Vec2);
	float DotProduct(Vector4 &In);
	float DotProduct(float XIn, float YIn, float ZIn);
	float Magnitude();
	void Normalise();



	Vector4 &operator *= (float Mult)
	{
		m_X *= Mult;
		m_Y *= Mult;
		m_Z *= Mult;
		m_W *= Mult;
		return *this;
	}

	Vector4 &operator += (Vector4& VIn)
	{
		m_X	+= VIn.m_X;
		m_Y	+= VIn.m_Y;
		m_Z	+= VIn.m_Z;
		m_W += VIn.m_W;
		return *this;
	}

	Vector4 &operator = (Vector4& VIn)
	{
		m_X = VIn.m_X;
		m_Y = VIn.m_Y;
		m_Z = VIn.m_Z;
		m_W = VIn.m_W;
		return *this;
	}
/*
	Vector4 operator - (Vector4& VIn)
	{
//		Vector4 Temp;

//		Temp.m_X = m_X - VIn.m_X;
//		Temp.m_Y = m_Y - VIn.m_Y;
//		Temp.m_Z = m_Z - VIn.m_Z;
//		Temp.m_Y = m_W - VIn.m_W;
		return Vector4(m_X - VIn.m_X, m_Y - VIn.m_Y, m_Z - VIn.m_Z, m_W - VIn.m_W);
	}
*/

	bool operator == (Vector4& VIn)
	{
		if(VIn.m_X != m_X) return false;
		if(VIn.m_Y != m_Y) return false;
		if(VIn.m_Z != m_Z) return false;
		return true;
	}
};


#endif //VECTOR4_INCLUDED
