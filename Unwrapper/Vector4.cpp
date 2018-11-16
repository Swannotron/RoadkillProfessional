#include <math.h>
#include "Vector4.h"


/* *************************************************************************************************************************************
	Vector type

	Description:
		Point in 3d space + manipulation functions

************************************************************************************************************************************** */



/********************************
*	Cross product				*
********************************/

float Vector4::CrossProduct(Vector4 &Vec1, Vector4& Vec2)
{
	Vector4 Temp;

//	Out.X	= Y * In.Z - In.Y * Z;
//	Out.Y	= X * In.Z - In.X * Z;
//	Out.Z	= X * In.Y - In.X * Y;

//	Temp.m_X = (Vec1.m_Y * Vec2.m_Z) - (Vec1.m_Z * Vec2.m_Y);
//	Temp.m_Y = (Vec1.m_X * Vec2.m_Z) - (Vec1.m_Z * Vec2.m_X);
//	Temp.m_Z = (Vec1.m_X * Vec2.m_Y) - (Vec1.m_Y * Vec2.m_X);

	Temp.m_X = (Vec1.m_Y * Vec2.m_Z) - (Vec1.m_Z * Vec2.m_Y);
	Temp.m_Y = (Vec1.m_Z * Vec2.m_X) - (Vec1.m_X * Vec2.m_Z);
	Temp.m_Z = (Vec1.m_X * Vec2.m_Y) - (Vec1.m_Y * Vec2.m_X);
	Temp.m_W = 1.0f;

	float Length = (Temp.m_X * Temp.m_X) + (Temp.m_Y * Temp.m_Y) + (Temp.m_Z * Temp.m_Z);			//  Only work on 3d displacement
	Length = (float)sqrt(Length);

	return(Length);
}



/********************************
*	Dot product					*
********************************/

float Vector4::DotProduct(Vector4 &In)
{
	float Result;
	Result = (m_X * In.m_X) + (m_Y * In.m_Y) + (m_Z * In.m_Z);
	return(Result);
}



/********************************
*	Dot product					*
********************************/

float Vector4::DotProduct(float XIn, float YIn, float ZIn)
{
	float Result;
	Result = (m_X * XIn) + (m_Y * YIn) + (m_Z * ZIn);
	return(Result);
}


/********************************
*	Get the Vectors Magnitude	*
********************************/

float Vector4::Magnitude()
{
	float Length;

	Length = (m_X * m_X) + (m_Y * m_Y) + (m_Z * m_Z);			//  Only work on 3d displacement
	Length = (float)sqrt(Length);

	return(Length);
}


/****************************
*	Normalise the Vector	*
****************************/

void Vector4::Normalise()
{
	float Length = Magnitude();

	m_X /= Length;
	m_Y /= Length;
	m_Z /= Length;
	m_W  = Length;
}
