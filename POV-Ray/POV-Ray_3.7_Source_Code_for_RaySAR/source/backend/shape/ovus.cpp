/*******************************************************************************
 * ovus.cpp
 *
 * This module implements the ovus primitive.
 *
 * ---------------------------------------------------------------------------
 * Persistence of Vision Ray Tracer ('POV-Ray') version 3.7.
 * Copyright 1991-2013 Persistence of Vision Raytracer Pty. Ltd.
 *
 * POV-Ray is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * POV-Ray is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ---------------------------------------------------------------------------
 * POV-Ray is based on the popular DKB raytracer version 2.12.
 * DKBTrace was originally written by David K. Buck.
 * DKBTrace Ver 2.0-2.12 were written by David K. Buck & Aaron A. Collins.
 * ---------------------------------------------------------------------------
 * $File: //depot/public/povray/3.x/source/backend/shape/ovus.cpp $
 * $Revision: #1 $
 * $Change: 6069 $
 * $DateTime: 2013/11/06 11:59:40 $
 * $Author: chrisc $
 *******************************************************************************/

/****************************************************************************
*
*  Explanation:
*
*    -
*
*  Syntax:
*
*  ovus
*  {
*     bottom_radius,top_radius
*  }
*
*  The so long awaited 'Egg' forms.
*
*  Normally, the bottom_radius is bigger than the top_radius
*  the center of the top sphere is at the zenith of the bottom sphere
*  and the bottom sphere is centered at 0.0
*  The radius of the connecting surface is the double of the biggest radius
*   (yes, the biggest diameter is used as the curvature of the connection)
*
*  The hard part was just to find where the connection starts and ends.
*
*  The name is ovus because:
*    - it's rather like torus (similar interface pushs for similar ending)
*    - the 2D curve is called 'ove' (latin for egg)
*
*****************************************************************************/


// frame.h must always be the first POV file included (pulls in platform config)
#include "backend/frame.h"
#include "backend/math/vector.h"
#include "backend/math/polysolv.h"
#include "backend/math/matrices.h"
#include "backend/scene/objects.h"
#include "backend/shape/ovus.h"
#include "backend/scene/threaddata.h"

// this must be the last file included
#include "base/povdebug.h"

namespace pov
{

/*****************************************************************************
* Local preprocessor defines
******************************************************************************/

// Minimal depth for a valid intersection.
// TODO FIXME - can we use EPSILON or a similar more generic constant instead?
const DBL DEPTH_TOLERANCE = 1.0e-4;

// Tolerance used for order reduction during root finding.
// TODO FIXME - can we use EPSILON or a similar more generic constant instead?
const DBL ROOT_TOLERANCE = 1.0e-4;



void Ovus::Intersect_Ovus_Spheres(const VECTOR &P, const VECTOR &D,
                                  DBL * Depth1, DBL * Depth2, DBL * Depth3,
                                  DBL * Depth4, DBL * Depth5, DBL * Depth6,
                                  SceneThreadData *Thread) const
{
	DBL OCSquared, t_Closest_Approach, Half_Chord, t_Half_Chord_Squared;
	VECTOR Padj, Ipoint;
	VECTOR IPoint;
	DBL R2, r2, Py2, Dy2, PDy2, k1, k2, horizontal, vertical;
	DBL Rad1, Rad2;
	int n;
	int lcount=0;
	VECTOR Second_Center;
	DBL c[5], r[4];

	*Depth1 = *Depth2 = *Depth3 = *Depth4 = *Depth5 = *Depth6 = -100; // TODO FIXME - magic value
	// no hit unless...

	VScale(Padj, P, -1);
	Rad1 = Sqr(BottomRadius);
	Rad2 = Sqr(TopRadius);
	VDot(OCSquared, Padj, Padj);

	VDot(t_Closest_Approach, Padj, D);

	if ((OCSquared < Rad1) || (t_Closest_Approach > EPSILON))
	{

		t_Half_Chord_Squared = Rad1 - OCSquared + Sqr(t_Closest_Approach);

		if (t_Half_Chord_Squared > EPSILON)
		{
			Half_Chord = sqrt(t_Half_Chord_Squared);

			*Depth1 = t_Closest_Approach - Half_Chord;
			*Depth2 = t_Closest_Approach + Half_Chord;
			VEvaluateRay(IPoint, P, *Depth1, D);
			if (IPoint[Y] < BottomVertical)
			{
				lcount++;
			}
			VEvaluateRay(IPoint, P, *Depth2, D);
			if (IPoint[Y] < BottomVertical)
			{
				lcount++;
			}
		}
	}
	if (lcount > 1) return;
	Make_Vector(Second_Center, 0, BottomRadius, 0);
	VSub(Padj, Second_Center, P);

	VDot(OCSquared, Padj, Padj);

	VDot(t_Closest_Approach, Padj, D);

	if ((OCSquared < Rad2) || (t_Closest_Approach > EPSILON))
	{

		t_Half_Chord_Squared = Rad2 - OCSquared + Sqr(t_Closest_Approach);

		if (t_Half_Chord_Squared > EPSILON)
		{
			Half_Chord = sqrt(t_Half_Chord_Squared);

			*Depth3 = t_Closest_Approach - Half_Chord;
			*Depth4 = t_Closest_Approach + Half_Chord;
			VEvaluateRay(IPoint, P, *Depth3, D);
			if (IPoint[Y] > TopVertical)
			{
				lcount++;
			}
			VEvaluateRay(IPoint, P, *Depth4, D);
			if (IPoint[Y] > TopVertical)
			{
				lcount++;
			}

		}
	}
	if (lcount > 1) return;
	Make_Vector(Second_Center, 0, VerticalPosition, 0);
	VSub(Padj, P, Second_Center);
	R2 = Sqr(HorizontalPosition);
	r2 = Sqr(ConnectingRadius);
	// Notice : ConnectingRadius > HorizontalPosition here !

	Py2 = Padj[Y] * Padj[Y];
	Dy2 = D[Y] * D[Y];
	PDy2 = Padj[Y] * D[Y];

	k1 = Padj[X] * Padj[X] + Padj[Z] * Padj[Z] + Py2 - R2 - r2;
	k2 = Padj[X] * D[X] + Padj[Z] * D[Z] + PDy2;
	// this is just like a big torus
	c[0] = 1.0;

	c[1] = 4.0 * k2;

	c[2] = 2.0 * (k1 + 2.0 * (k2 * k2 + R2 * Dy2));

	c[3] = 4.0 * (k2 * k1 + 2.0 * R2 * PDy2);

	c[4] = k1 * k1 + 4.0 * R2 * (Py2 - r2);

	n = Solve_Polynomial(4, c, r, Test_Flag(this, STURM_FLAG), ROOT_TOLERANCE, Thread);
	while (n--)
	{
		// here we only keep the 'lemon' inside the torus
		// and dismiss the 'apple'
		// If you find a solution to resolve the rotation of
		//   (x + r)^2 + y^2 = R^2 around y (so replacing x by sqrt(x^2+z^2))
		// with something which is faster than a 4th degree polynome,
		// please feel welcome to update and share...

		VEvaluateRay(Ipoint, P, r[n], D);

		vertical = Ipoint[Y];
		if ((vertical > BottomVertical) && (vertical < TopVertical))
		{
			horizontal = sqrt(Sqr(Ipoint[X]) + Sqr(Ipoint[Z]));
			OCSquared = Sqr((horizontal + HorizontalPosition)) + Sqr((vertical - VerticalPosition));
			if (fabs(OCSquared - Sqr(ConnectingRadius)) < ROOT_TOLERANCE)
			{
				if (*Depth5 < 0)
				{
					*Depth5 = r[n];
				}
				else
				{
					*Depth6 = r[n];
				}
			}
		}
	}
}
/*****************************************************************************
*
* FUNCTION
*
*   All_Ovus_Intersections
*
* INPUT
*
*   Object      - Object
*   Ray         - Ray
*   Depth_Stack - Intersection stack
*
* OUTPUT
*
*   Depth_Stack
*
* RETURNS
*
*   int - true, if an intersection was found
*
* AUTHOR
*
*   Jerome Grimbert
*
* DESCRIPTION
*
*   Determine ray/ovus intersection and clip intersection found.
*
* CHANGES
*
*   Jun 1994 : Creation.
*
******************************************************************************/

bool Ovus::All_Intersections(const Ray& ray, IStack& Depth_Stack, SceneThreadData *Thread)
{
	bool Found = false;
	VECTOR Real_Normal, Real_Pt, INormal, IPoint;
	DBL Depth1, Depth2, Depth3, Depth4, Depth5, Depth6;
	DBL len, horizontal;
	VECTOR P,D;

	Thread->Stats()[Ray_Ovus_Tests]++;
	MInvTransPoint(P, ray.Origin, Trans);
	MInvTransDirection(D, ray.Direction, Trans);
	VLength(len, D);
	VInverseScaleEq(D, len);

	Intersect_Ovus_Spheres(P, D, &Depth1, &Depth2, &Depth3,
	                       &Depth4, &Depth5, &Depth6, Thread);
	if (Depth1 > EPSILON)
	{
		VEvaluateRay(IPoint, P, Depth1, D);
		if (IPoint[Y] < BottomVertical)
		{
			MTransPoint(Real_Pt, IPoint, Trans);
			if (Clip.empty()||(Point_In_Clip(Real_Pt, Clip, Thread)))
			{
				Assign_Vector(INormal, IPoint);
				VInverseScaleEq(INormal, BottomRadius);
				MTransNormal(Real_Normal, INormal, Trans);
				VNormalizeEq(Real_Normal);
				Depth_Stack->push(Intersection(Depth1/len, Real_Pt, Real_Normal, this));
				Found = true;
			}
		}
	}

	if (Depth2 > EPSILON)
	{
		VEvaluateRay(IPoint, P, Depth2, D);

		if (IPoint[Y] < BottomVertical)
		{
			MTransPoint(Real_Pt, IPoint, Trans);
			if (Clip.empty()||(Point_In_Clip(Real_Pt, Clip, Thread)))
			{
				Assign_Vector(INormal, IPoint);
				VInverseScaleEq(INormal, BottomRadius);
				MTransNormal(Real_Normal, INormal, Trans);
				VNormalizeEq(Real_Normal);
				Depth_Stack->push(Intersection(Depth2/len, Real_Pt, Real_Normal, this));
				Found = true;
			}
		}
	}

	if (Depth3 > EPSILON)
	{
		VEvaluateRay(IPoint, P, Depth3, D);

		if (IPoint[Y] > TopVertical)
		{
			MTransPoint(Real_Pt, IPoint, Trans);
			if (Clip.empty()||(Point_In_Clip(Real_Pt, Clip, Thread)))
			{
				Assign_Vector(INormal, IPoint);
				INormal[Y] -= BottomRadius;
				VInverseScaleEq(INormal, TopRadius);
				MTransNormal(Real_Normal, INormal, Trans);
				VNormalizeEq(Real_Normal);
				Depth_Stack->push(Intersection(Depth3/len, Real_Pt, Real_Normal, this));
				Found = true;
			}
		}
	}
	if (Depth4 > EPSILON)
	{
		VEvaluateRay(IPoint, P, Depth4, D);

		if (IPoint[Y] > TopVertical)
		{
			MTransPoint(Real_Pt, IPoint, Trans);
			if (Clip.empty()||(Point_In_Clip(Real_Pt, Clip, Thread)))
			{
				Assign_Vector(INormal, IPoint);
				INormal[Y] -= BottomRadius;
				VInverseScaleEq(INormal, TopRadius);
				MTransNormal(Real_Normal, INormal, Trans);
				VNormalizeEq(Real_Normal);
				Depth_Stack->push(Intersection(Depth4/len, Real_Pt, Real_Normal, this));
				Found = true;
			}
		}
	}

	if (Depth5 > EPSILON)
	{
		VEvaluateRay(IPoint, P, Depth5, D);
		MTransPoint(Real_Pt, IPoint, Trans);

		if (Clip.empty()||(Point_In_Clip(Real_Pt, Clip, Thread)))
		{
			Assign_Vector(INormal, IPoint);

			INormal[Y] -= VerticalPosition;
			horizontal = sqrt(Sqr(INormal[X]) + Sqr(INormal[Z]));
			INormal[X] += (INormal[X] * HorizontalPosition / horizontal);
			INormal[Z] += (INormal[Z] * HorizontalPosition / horizontal);
			VNormalizeEq(INormal);
			MTransNormal(Real_Normal, INormal, Trans);
			VNormalizeEq(Real_Normal);
			Depth_Stack->push(Intersection(Depth5/len, Real_Pt, Real_Normal, this));
			Found = true;
		}
	}
	if (Depth6 > EPSILON)
	{
		VEvaluateRay(IPoint, P, Depth6, D);
		MTransPoint(Real_Pt, IPoint, Trans);

		if (Clip.empty()||(Point_In_Clip(Real_Pt, Clip, Thread)))
		{
			Assign_Vector(INormal, IPoint);
			INormal[Y] -= VerticalPosition;
			horizontal = sqrt(Sqr(INormal[X]) + Sqr(INormal[Z]));
			INormal[X] += (INormal[X] * HorizontalPosition / horizontal);
			INormal[Z] += (INormal[Z] * HorizontalPosition / horizontal);
			VNormalizeEq(INormal);
			MTransNormal(Real_Normal, INormal, Trans);
			VNormalizeEq(Real_Normal);

			Depth_Stack->push(Intersection(Depth6/len, Real_Pt, Real_Normal, this));
			Found = true;
		}
	}
	if (Found)
	{
		Thread->Stats()[Ray_Ovus_Tests_Succeeded]++;
	}
	return (Found);
}



/*****************************************************************************
*
* FUNCTION
*
*   Inside_Ovus
*
* INPUT
*
*   IPoint - Intersection point
*
* OUTPUT
*
* RETURNS
*
*   int - true if inside
*
* AUTHOR
*
*   Jerome Grimbert
*
* DESCRIPTION
*
*   Test if a point lies inside the ovus.
*
* CHANGES
*
*
******************************************************************************/

bool Ovus::Inside(const VECTOR IPoint, TraceThreadData *Thread) const
{
	DBL OCSquared;
	DBL horizontal, vertical;
	bool INSide = false;
	VECTOR Origin, New_Point, Other;
	Make_Vector(Origin, 0, BottomRadius, 0);
	MInvTransPoint(New_Point, IPoint, Trans);
	VDot(OCSquared, New_Point, New_Point);
	if (OCSquared < Sqr(BottomRadius))
	{
		INSide = true;
	}
	VSub(Other, New_Point, Origin);
	VDot(OCSquared, Other, Other);
	if (OCSquared < Sqr(TopRadius))
	{
		INSide = true;
	}
	vertical = New_Point[Y];
	if ((vertical > BottomVertical) && (vertical < TopVertical))
	{
		horizontal = sqrt(Sqr(New_Point[X]) + Sqr(New_Point[Z]));
		OCSquared = Sqr(horizontal + HorizontalPosition) + Sqr((vertical - VerticalPosition));
		if (OCSquared < Sqr(ConnectingRadius))
		{
			INSide = true;
		}
	}
	if (Test_Flag(this, INVERTED_FLAG))
	{
		return !INSide;
	}
	else
	{
		return INSide;
	}
}



/*****************************************************************************
*
* FUNCTION
*
*   Ovus_Normal
*
* INPUT
*
*   Result - Normal vector
*   Object - Object
*   Inter  - Intersection found
*
* OUTPUT
*
*   Result
*
* RETURNS
*
* AUTHOR
*
*   Jerome Grimbert
*
* DESCRIPTION
*
*   Calculate the normal of the ovus in a given point.
*
* CHANGES
*
******************************************************************************/

void Ovus::Normal(VECTOR Result, Intersection *Inter, TraceThreadData *Thread) const
{
	Assign_Vector(Result, Inter->INormal);
}



/*****************************************************************************
*
* FUNCTION
*
*   Translate_Ovus
*
* INPUT
*
*   Object - Object
*   Vector - Translation vector
*
* OUTPUT
*
*   Object
*
* RETURNS
*
* AUTHOR
*
*   Jerome Grimbert
*
* DESCRIPTION
*
*   Translate an ovus.
*
* CHANGES
*
*
******************************************************************************/

void Ovus::Translate(const VECTOR, const TRANSFORM *tr)
{
	Transform(tr);
}



/*****************************************************************************
*
* FUNCTION
*
*   Rotate_Ovus
*
* INPUT
*
*   Object - Object
*   Vector - Rotation vector
*
* OUTPUT
*
*   Object
*
* RETURNS
*
* AUTHOR
*
*   Jerome Grimbert
*
* DESCRIPTION
*
*   Rotate an ovus.
*
* CHANGES
*
*
******************************************************************************/

void Ovus::Rotate(const VECTOR, const TRANSFORM *tr)
{
	Transform(tr);
}



/*****************************************************************************
*
* FUNCTION
*
*   Scale_Ovus
*
* INPUT
*
*   Object - Object
*   Vector - Scaling vector
*
* OUTPUT
*
*   Object
*
* RETURNS
*
* AUTHOR
*
*   Jerome Grimbert
*
* DESCRIPTION
*
*   Scale an ovus.
*
* CHANGES
*
*
******************************************************************************/

void Ovus::Scale(const VECTOR, const TRANSFORM *tr)
{
	Transform(tr);
}



/*****************************************************************************
*
* FUNCTION
*
*   Transform_Ovus
*
* INPUT
*
*   Object - Object
*   Trans  - Transformation to apply
*
* OUTPUT
*
*   Object
*
* RETURNS
*
* AUTHOR
*
*   Jerome Grimbert
*
* DESCRIPTION
*
*   Transform an ovus and recalculate its bounding box.
*
* CHANGES
*
*
******************************************************************************/

void Ovus::Transform(const TRANSFORM *tr)
{
	if(Trans == NULL)
		Trans = Create_Transform();

	Compose_Transforms(Trans, tr);

	Compute_BBox();
}



/*****************************************************************************
*
* FUNCTION
*
*   Invert_Ovus
*
* INPUT
*
*   Object - Object
*
* OUTPUT
*
*   Object
*
* RETURNS
*
* AUTHOR
*
*   Jerome Grimbert
*
* DESCRIPTION
*
*   Invert an ovus.
*
* CHANGES
*
*
******************************************************************************/

void Ovus::Invert()
{
	Invert_Flag(this, INVERTED_FLAG);
}



/*****************************************************************************
*
* FUNCTION
*
*   Create_Ovus
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
*   OVUS * - new ovus
*
* AUTHOR
*
*   Jerome Grimbert
*
* DESCRIPTION
*
*   Create a new ovus.
*
* CHANGES
*
*
******************************************************************************/

Ovus::Ovus() : ObjectBase(OVUS_OBJECT)
{
	Trans = Create_Transform();

	TopRadius = 0.0;
	BottomRadius = 0.0;
	HorizontalPosition = 0.0;
	VerticalPosition = 0.0;
	BottomVertical = 0.0;
	TopVertical = 0.0;
	ConnectingRadius = 0.0;
}



/*****************************************************************************
*
* FUNCTION
*
*   Copy_Ovus
*
* INPUT
*
*   Object - Object
*
* OUTPUT
*
* RETURNS
*
*   void * - New ovus
*
* AUTHOR
*
*   Jerome Grimbert
*
* DESCRIPTION
*
*   Copy an ovus.
*
* CHANGES
*
*
******************************************************************************/

ObjectPtr Ovus::Copy()
{
	Ovus *New = new Ovus();

	Destroy_Transform(New->Trans);
	*New = *this;
	New->Trans = Copy_Transform(Trans);

	return (New);
}



/*****************************************************************************
*
* FUNCTION
*
*   Destroy_Ovus
*
* INPUT
*
*   Object - Object
*
* OUTPUT
*
*   Object
*
* RETURNS
*
* AUTHOR
*
*   Jerome Grimbert
*
* DESCRIPTION
*
*   Destroy a ovus.
*
* CHANGES
*
*
******************************************************************************/

Ovus::~Ovus()
{
	Destroy_Transform(Trans);
}



/*****************************************************************************
*
* FUNCTION
*
*   Compute_Ovus_BBox
*
* INPUT
*
*   Ovus - Ovus
*
* OUTPUT
*
*   Ovus
*
* RETURNS
*
* AUTHOR
*
*   Jerome Grimbert
*
* DESCRIPTION
*
*   Calculate the bounding box of an ovus.
*
* CHANGES
*
*
******************************************************************************/

void Ovus::Compute_BBox()
{
	// Compute the biggest vertical cylinder radius
	DBL biggest;
	biggest = ConnectingRadius - HorizontalPosition;
	if (biggest < BottomRadius)
	{
		biggest = BottomRadius;
	}
	if (biggest < TopRadius)
	{
		biggest = TopRadius;
	}

	Make_BBox(BBox, -biggest, -BottomRadius, -biggest,
	          2.0 * biggest, 2.0 * BottomRadius + TopRadius, 2.0 * biggest);

	Recompute_BBox(&BBox, Trans);
}



/*****************************************************************************
*
* FUNCTION
*
*   Ovus_UVCoord
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Jerome Grimbert
*
* DESCRIPTION
*
*   -
*
* CHANGES
*
*   -
*
******************************************************************************/

void Ovus::UVCoord(UV_VECT Result, const Intersection *Inter, TraceThreadData *Thread) const
{
	CalcUV(Inter->IPoint, Result);
}


/*****************************************************************************
*
* FUNCTION
*
*   CalcUV
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Jerome Grimbert
*
* DESCRIPTION
*
*   Calculate the u/v coordinate of a point on an ovus
*
* CHANGES
*
******************************************************************************/

void Ovus::CalcUV(const VECTOR IPoint, UV_VECT Result) const
{
	DBL len, x, y, z;
	DBL phi, theta;
	VECTOR P;

	// Transform the ray into the ovus space.
	MInvTransPoint(P, IPoint, Trans);

	// the center of UV coordinate is the bottom center when top radius ->0
	// and it is the top center when top radius -> 2.0* bottom radius
	// when top radius == bottom radius, it is half-way between both center
	//
	// bottom center is <0,0,0>
	// top center is <0,BottomRadius,0>
	// TODO FIXME - comment doesn't seem to match the following code
	x = P[X];
//	y = P[Y] - BottomRadius*(TopRadius/(2.0*BottomRadius));
	y = P[Y] - (TopRadius/2.0);
	z = P[Z];

	// now assume it's just a sphere, for UV mapping/projection
	len = sqrt(x * x + y * y + z * z);

	if (len == 0.0)
		return;
	else
	{
		x /= len;
		y /= len;
		z /= len;
	}

	// Determine its angle from the x-z plane.
	phi = 0.5 + asin(y) / M_PI; // This will be from 0 to 1

	// Determine its angle from the point (1, 0, 0) in the x-z plane.
	len = x * x + z * z;

	if (len > EPSILON)
	{
		len = sqrt(len);
		if (z == 0.0)
		{
			if (x > 0)
				theta = 0.0;
			else
				theta = M_PI;
		}
		else
		{
			theta = acos(x / len);
			if (z < 0.0)
				theta = TWO_M_PI - theta;
		}

		theta /= TWO_M_PI; // This will be from 0 to 1
	}
	else
		// This point is at one of the poles. Any value of xcoord will be ok...
		theta = 0;

	Result[U] = theta;
	Result[V] = phi;

}

}
