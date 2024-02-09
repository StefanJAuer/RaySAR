/*******************************************************************************
 * torus.cpp
 *
 * This module implements functions that manipulate torii.
 *
 * This module was written by Dieter Bayer [DB].
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
 * $File: //depot/public/povray/3.x/source/backend/shape/torus.cpp $
 * $Revision: #1 $
 * $Change: 6069 $
 * $DateTime: 2013/11/06 11:59:40 $
 * $Author: chrisc $
 *******************************************************************************/

/****************************************************************************
*
*  Explanation:
*
*  ---
*
*  June 1994 : Creation. [DB]
*
*****************************************************************************/

// frame.h must always be the first POV file included (pulls in platform config)
#include "backend/frame.h"
#include "povray.h"
#include "backend/math/vector.h"
#include "backend/bounding/bbox.h"
#include "backend/math/polysolv.h"
#include "backend/math/matrices.h"
#include "backend/scene/objects.h"
#include "backend/shape/torus.h"
#include "backend/scene/threaddata.h"

// this must be the last file included
#include "base/povdebug.h"

namespace pov
{

/*****************************************************************************
* Local preprocessor defines
******************************************************************************/

/* Minimal depth for a valid intersection. */

const DBL DEPTH_TOLERANCE = 1.0e-4;

/* Tolerance used for order reduction during root finding. */

const DBL ROOT_TOLERANCE = 1.0e-4;



/*****************************************************************************
*
* FUNCTION
*
*   All_Torus_Intersections
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
*   Dieter Bayer
*   
* DESCRIPTION
*
*   Determine ray/torus intersection and clip intersection found.
*
* CHANGES
*
*   Jun 1994 : Creation.
*
******************************************************************************/

bool Torus::All_Intersections(const Ray& ray, IStack& Depth_Stack, SceneThreadData *Thread)
{
	int i, max_i, Found;
	DBL Depth[4];
	VECTOR IPoint;

	Found = false;

	if ((max_i = Intersect(ray, Depth, Thread)) > 0)
	{
		for (i = 0; i < max_i; i++)
		{
			if ((Depth[i] > DEPTH_TOLERANCE) && (Depth[i] < MAX_DISTANCE))
			{
				VEvaluateRay(IPoint, ray.Origin, Depth[i], ray.Direction);

				if (Clip.empty() || Point_In_Clip(IPoint, Clip, Thread))
				{
					Depth_Stack->push(Intersection(Depth[i], IPoint, this));

					Found = true;
				}
			}
		}
	}

	return(Found);
}



/*****************************************************************************
*
* FUNCTION
*
*   intersect_torus
*
* INPUT
*
*   Ray   - Ray
*   Torus - Torus
*   Depth - Intersections found
*   
* OUTPUT
*
*   Depth
*   
* RETURNS
*
*   int - Number of intersections found
*   
* AUTHOR
*
*   Dieter Bayer
*   
* DESCRIPTION
*
*   Determine ray/torus intersection.
*
*   Note that the torus is rotated about the y-axis!
*
* CHANGES
*
*   Jun 1994 : Creation.
*
******************************************************************************/

int Torus::Intersect(const Ray& ray, DBL *Depth, SceneThreadData *Thread) const
{
	int i, n;
	DBL len, R2, Py2, Dy2, PDy2, k1, k2;
	DBL y1, y2, r1, r2;
	DBL c[5];
	DBL r[4];
	VECTOR P, D;
	DBL DistanceP;            // Distance from P to torus center (origo).
	DBL BoundingSphereRadius; // Sphere fully (amply) enclosing torus.
	DBL Closer;               // P is moved Closer*D closer to torus.

	Thread->Stats()[Ray_Torus_Tests]++;

	/* Transform the ray into the torus space. */

	MInvTransPoint(P, ray.Origin, Trans);

	MInvTransDirection(D, ray.Direction, Trans);

	VLength(len, D);

	VInverseScaleEq(D, len);

	i = 0;

	y1 = -MinorRadius;
	y2 =  MinorRadius;
	r1 = Sqr(MajorRadius - MinorRadius);

	if ( MajorRadius < MinorRadius )
		r1 = 0;

	r2 = Sqr(MajorRadius + MinorRadius);

#ifdef TORUS_EXTRA_STATS
	Thread->Stats()[Torus_Bound_Tests]++;
#endif

	if (Test_Thick_Cylinder(P, D, y1, y2, r1, r2))
	{
#ifdef TORUS_EXTRA_STATS
		Thread->Stats()[Torus_Bound_Tests_Succeeded]++;
#endif

		// Move P close to bounding sphere to have more precise root calculation.
		// Bounding sphere radius is R + r, we add r once more to ensure
		// that P is safely outside sphere.
		BoundingSphereRadius = MajorRadius + MinorRadius + MinorRadius;
		DistanceP = VSumSqr(P); // Distance is currently squared.
		Closer = 0.0;
		if (DistanceP > Sqr(BoundingSphereRadius))
		{
			DistanceP = sqrt(DistanceP); // Now real distance.
			Closer = DistanceP - BoundingSphereRadius;
			VAddScaledEq(P, Closer, D);
		}

		R2   = Sqr(MajorRadius);
		r2   = Sqr(MinorRadius);

		Py2  = P[Y] * P[Y];
		Dy2  = D[Y] * D[Y];
		PDy2 = P[Y] * D[Y];

		k1   = P[X] * P[X] + P[Z] * P[Z] + Py2 - R2 - r2;
		k2   = P[X] * D[X] + P[Z] * D[Z] + PDy2;

		c[0] = 1.0;

		c[1] = 4.0 * k2;

		c[2] = 2.0 * (k1 + 2.0 * (k2 * k2 + R2 * Dy2));

		c[3] = 4.0 * (k2 * k1 + 2.0 * R2 * PDy2);

		c[4] = k1 * k1 + 4.0 * R2 * (Py2 - r2);

		n = Solve_Polynomial(4, c, r, Test_Flag(this, STURM_FLAG), ROOT_TOLERANCE, Thread);

		while(n--)
			Depth[i++] = (r[n] + Closer) / len;
	}

	if (i)
		Thread->Stats()[Ray_Torus_Tests_Succeeded]++;

	return(i);
}



/*****************************************************************************
*
* FUNCTION
*
*   Inside_Torus
*
* INPUT
*
*   IPoint - Intersection point
*   Object - Object
*   
* OUTPUT
*   
* RETURNS
*
*   int - true if inside
*   
* AUTHOR
*
*   Dieter Bayer
*   
* DESCRIPTION
*
*   Test if a point lies inside the torus.
*
* CHANGES
*
*   Jun 1994 : Creation.
*
******************************************************************************/

bool Torus::Inside(const VECTOR IPoint, TraceThreadData *Thread) const
{
	DBL r, r2;
	VECTOR P;

	/* Transform the point into the torus space. */

	MInvTransPoint(P, IPoint, Trans);

	r  = sqrt(Sqr(P[X]) + Sqr(P[Z]));

	r2 = Sqr(P[Y]) + Sqr(r - MajorRadius);

	if (r2 <= Sqr(MinorRadius))
	{
		return(!Test_Flag(this, INVERTED_FLAG));
	}
	else
	{
		return(Test_Flag(this, INVERTED_FLAG));
	}
}



/*****************************************************************************
*
* FUNCTION
*
*   Torus_Normal
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
*   Dieter Bayer
*   
* DESCRIPTION
*
*   Calculate the normal of the torus in a given point.
*
* CHANGES
*
*   Jun 1994 : Creation.
*
******************************************************************************/

void Torus::Normal(VECTOR Result, Intersection *Inter, TraceThreadData *Thread) const
{
	DBL dist;
	VECTOR P, N, M;

	/* Transform the point into the torus space. */

	MInvTransPoint(P, Inter->IPoint, Trans);

	/* Get normal from derivatives. */

	dist = sqrt(P[X] * P[X] + P[Z] * P[Z]);

	if (dist > EPSILON)
	{
		M[X] = MajorRadius * P[X] / dist;
		M[Y] = 0.0;
		M[Z] = MajorRadius * P[Z] / dist;
	}
	else
	{
		Make_Vector(M, 0.0, 0.0, 0.0);
	}

	VSub(N, P, M);

	/* Transform the normalt out of the torus space. */

	MTransNormal(Result, N, Trans);

	VNormalize(Result, Result);
}



/*****************************************************************************
*
* FUNCTION
*
*   Translate_Torus
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
*   Dieter Bayer
*
* DESCRIPTION
*
*   Translate a torus.
*
* CHANGES
*
*   Jun 1994 : Creation.
*
******************************************************************************/

void Torus::Translate(const VECTOR, const TRANSFORM *tr)
{
	Transform(tr);
}



/*****************************************************************************
*
* FUNCTION
*
*   Rotate_Torus
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
*   Dieter Bayer
*   
* DESCRIPTION
*
*   Rotate a torus.
*
* CHANGES
*
*   Jun 1994 : Creation.
*
******************************************************************************/

void Torus::Rotate(const VECTOR, const TRANSFORM *tr)
{
	Transform(tr);
}



/*****************************************************************************
*
* FUNCTION
*
*   Scale_Torus
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
*   Dieter Bayer
*   
* DESCRIPTION
*
*   Scale a torus.
*
* CHANGES
*
*   Jun 1994 : Creation.
*
******************************************************************************/

void Torus::Scale(const VECTOR, const TRANSFORM *tr)
{
	Transform(tr);
}



/*****************************************************************************
*
* FUNCTION
*
*   Transform_Torus
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
*   Dieter Bayer
*   
* DESCRIPTION
*
*   Transform a torus and recalculate its bounding box.
*
* CHANGES
*
*   Jun 1994 : Creation.
*
******************************************************************************/

void Torus::Transform(const TRANSFORM *tr)
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
*   Invert_Torus
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
*   Dieter Bayer
*   
* DESCRIPTION
*
*   Invert a torus.
*
* CHANGES
*
*   Jun 1994 : Creation.
*
******************************************************************************/

void Torus::Invert()
{
	Invert_Flag(this, INVERTED_FLAG);
}



/*****************************************************************************
*
* FUNCTION
*
*   Create_Torus
*
* INPUT
*   
* OUTPUT
*   
* RETURNS
*
*   TORUS * - new torus
*   
* AUTHOR
*
*   Dieter Bayer
*   
* DESCRIPTION
*
*   Create a new torus.
*
* CHANGES
*
*   Jun 1994 : Creation.
*
******************************************************************************/

Torus::Torus() : ObjectBase(TORUS_OBJECT)
{
	Trans = Create_Transform();

	MajorRadius = 0.0;
	MinorRadius = 0.0;
}



/*****************************************************************************
*
* FUNCTION
*
*   Copy_Torus
*
* INPUT
*
*   Object - Object
*   
* OUTPUT
*   
* RETURNS
*
*   void * - New torus
*   
* AUTHOR
*
*   Dieter Bayer
*   
* DESCRIPTION
*
*   Copy a torus.
*
* CHANGES
*
*   Jun 1994 : Creation.
*
*   Sep 1994 : fixed memory leakage [DB]
*
******************************************************************************/

ObjectPtr Torus::Copy()
{
	Torus *New = new Torus();

	Destroy_Transform(New->Trans);
	*New = *this;
	New->Trans = Copy_Transform(Trans);

	return (New);
}



/*****************************************************************************
*
* FUNCTION
*
*   Destroy_Torus
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
*   Dieter Bayer
*   
* DESCRIPTION
*
*   Destroy a torus.
*
* CHANGES
*
*   Jun 1994 : Creation.
*
******************************************************************************/

Torus::~Torus()
{
	Destroy_Transform(Trans);
}



/*****************************************************************************
*
* FUNCTION
*
*   Compute_Torus_BBox
*
* INPUT
*
*   Torus - Torus
*   
* OUTPUT
*
*   Torus
*   
* RETURNS
*   
* AUTHOR
*
*   Dieter Bayer
*   
* DESCRIPTION
*
*   Calculate the bounding box of a torus.
*
* CHANGES
*
*   Jun 1994 : Creation.
*
******************************************************************************/

void Torus::Compute_BBox()
{
	DBL r1, r2;

	r1 = MinorRadius;
	r2 = MajorRadius + MinorRadius;

	Make_BBox(BBox, -r2, -r1, -r2, 2.0 * r2, 2.0 * r1, 2.0 * r2);

	Recompute_BBox(&BBox, Trans);
}



/*****************************************************************************
*
* FUNCTION
*
*   Test_Thick_Cylinder
*
* INPUT
*
*   P  - Ray initial
*   D  - Ray direction
*   h1 - Height 1
*   h2 - Height 2
*   r1 - Square of inner radius
*   r2 - Square of outer radius
*   
* OUTPUT
*   
* RETURNS
*
*   int - true, if hit
*   
* AUTHOR
*
*   Dieter Bayer
*   
* DESCRIPTION
*
*   Test if a given ray defined in the lathe's coordinate system
*   intersects a "thick" cylinder (rotated about y-axis).
*
* CHANGES
*
*   Jun 1994 : Creation.
*
******************************************************************************/

bool Torus::Test_Thick_Cylinder(const VECTOR P, const VECTOR D, DBL h1, DBL h2, DBL r1, DBL r2) const
{
	DBL a, b, c, d;
	DBL u, v, k, r, h;

	if (fabs(D[Y]) < EPSILON)
	{
		if ((P[Y] < h1) || (P[Y] > h2))
		{
			return(false);
		}
	}
	else
	{
		/* Intersect ray with the cap-plane. */

		k = (h2 - P[Y]) / D[Y];

		u = P[X] + k * D[X];
		v = P[Z] + k * D[Z];

		if ((k > EPSILON) && (k < MAX_DISTANCE))
		{
			r = u * u + v * v;

			if ((r >= r1) && (r <= r2))
			{
				return(true);
			}
		}

		/* Intersectionersect ray with the base-plane. */

		k = (h1 - P[Y]) / D[Y];

		u = P[X] + k * D[X];
		v = P[Z] + k * D[Z];

		if ((k > EPSILON) && (k < MAX_DISTANCE))
		{
			r = u * u + v * v;

			if ((r >= r1) && (r <= r2))
			{
				return(true);
			}
		}
	}

	a = D[X] * D[X] + D[Z] * D[Z];

	if (a > EPSILON)
	{
		/* Intersect with outer cylinder. */

		b = P[X] * D[X] + P[Z] * D[Z];

		c = P[X] * P[X] + P[Z] * P[Z] - r2;

		d = b * b - a * c;

		if (d >= 0.0)
		{
			d = sqrt(d);

			k = (-b + d) / a;

			if ((k > EPSILON) && (k < MAX_DISTANCE))
			{
				h = P[Y] + k * D[Y];

				if ((h >= h1) && (h <= h2))
				{
					return(true);
				}
			}

			k = (-b - d) / a;

			if ((k > EPSILON) && (k < MAX_DISTANCE))
			{
				h = P[Y] + k * D[Y];

				if ((h >= h1) && (h <= h2))
				{
					return(true);
				}
			}
		}

		/* Intersect with inner cylinder. */

		c = P[X] * P[X] + P[Z] * P[Z] - r1;

		d = b * b - a * c;

		if (d >= 0.0)
		{
			d = sqrt(d);

			k = (-b + d) / a;

			if ((k > EPSILON) && (k < MAX_DISTANCE))
			{
				h = P[Y] + k * D[Y];

				if ((h >= h1) && (h <= h2))
				{
					return(true);
				}
			}

			k = (-b - d) / a;

			if ((k > EPSILON) && (k < MAX_DISTANCE))
			{
				h = P[Y] + k * D[Y];

				if ((h >= h1) && (h <= h2))
				{
					return(true);
				}
			}
		}
	}

	return(false);
}


/*****************************************************************************
*
* FUNCTION
*
*   Torus_UVCoord
*
* INPUT
*   
* OUTPUT
*   
* RETURNS
*   
* AUTHOR
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

void Torus::UVCoord(UV_VECT Result, const Intersection *Inter, TraceThreadData *Thread) const
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
* DESCRIPTION
*
*   Calculate the u/v coordinate of a point on a torus
*
* CHANGES
*
*   Fix with correct placing for intersection point. It have to be
*    untransformed before further calculations. [ABX]
*   Fix with correct space torus space. Change meaning of y and z. [ABX]
*
******************************************************************************/

void Torus::CalcUV(const VECTOR IPoint, UV_VECT Result) const
{
	DBL len, v, u, x, y, z;
	VECTOR P;

	// Transform the ray into the torus space.
	MInvTransPoint(P, IPoint, Trans);
	x = P[X];
	y = P[Y];
	z = P[Z];

	// Determine its angle from the y-axis.
	u = (1.0 - (atan2(z, x) + M_PI) / TWO_M_PI);

	len = sqrt(x * x + z * z);

	// Now rotate about the y-axis to get the point P into the x-z plane.
	x = len - MajorRadius;
	v = (atan2(y, x) + M_PI) / TWO_M_PI;

	Result[U] = u;
	Result[V] = v;
}

}
