/*******************************************************************************
 * discs.cpp
 *
 * This module implements the disc primitive.
 * This file was written by Alexander Enzmann.  He wrote the code for
 * discs and generously provided us these enhancements.
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
 * $File: //depot/public/povray/3.x/source/backend/shape/discs.cpp $
 * $Revision: #1 $
 * $Change: 6069 $
 * $DateTime: 2013/11/06 11:59:40 $
 * $Author: chrisc $
 *******************************************************************************/

// frame.h must always be the first POV file included (pulls in platform config)
#include "backend/frame.h"
#include "backend/math/vector.h"
#include "backend/bounding/bbox.h"
#include "backend/shape/discs.h"
#include "backend/math/matrices.h"
#include "backend/scene/objects.h"
#include "backend/scene/threaddata.h"

// this must be the last file included
#include "base/povdebug.h"

namespace pov
{

const DBL DEPTH_TOLERANCE = 1.0e-6;

/*****************************************************************************
*
* FUNCTION
*
*   All_Disc_Intersections
*
* INPUT
*   
* OUTPUT
*   
* RETURNS
*   
* AUTHOR
*
*   Alexander Enzmann
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

bool Disc::All_Intersections(const Ray& ray, IStack& Depth_Stack, TraceThreadData *Thread)
{
	int Intersection_Found;
	DBL Depth;
	VECTOR IPoint;

	Intersection_Found = false;

	Thread->Stats()[Ray_Disc_Tests]++;
	if (Intersect(ray, &Depth))
	{
		Thread->Stats()[Ray_Disc_Tests_Succeeded]++;
		VEvaluateRay(IPoint, ray.Origin, Depth, ray.Direction);

		if (Clip.empty() || Point_In_Clip (IPoint, Clip, Thread))
		{
			Depth_Stack->push(Intersection(Depth,IPoint,this));
			Intersection_Found = true;
		}
	}

	return (Intersection_Found);
}



/*****************************************************************************
*
* FUNCTION
*
*   Intersect_Disc
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Alexander Enzmann
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

bool Disc::Intersect(const Ray& ray, DBL *Depth) const
{
	DBL t, u, v, r2, len;
	VECTOR P, D;

	/* Transform the point into the discs space */

	MInvTransPoint(P, ray.Origin, Trans);
	MInvTransDirection(D, ray.Direction, Trans);

	VLength(len, D);
	VInverseScaleEq(D, len);

	if (fabs(D[Z]) > EPSILON)
	{
		t = -P[Z] / D[Z];

		if (t >= 0.0)
		{
			u = P[X] + t * D[X];
			v = P[Y] + t * D[Y];

			r2 = Sqr(u) + Sqr(v);

			if ((r2 >= iradius2) && (r2 <= oradius2))
			{
				*Depth = t / len;

				if ((*Depth > DEPTH_TOLERANCE) && (*Depth < MAX_DISTANCE))
					return (true);
			}
		}
	}

	return (false);
}



/*****************************************************************************
*
* FUNCTION
*
*   Inside_Disc
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Alexander Enzmann
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

bool Disc::Inside(const VECTOR IPoint, TraceThreadData *Thread) const
{
	VECTOR New_Point;

	/* Transform the point into the discs space */

	MInvTransPoint(New_Point, IPoint, Trans);

	if (New_Point[Z] >= 0.0)
	{
		/* We are outside. */

		return (Test_Flag(this, INVERTED_FLAG));
	}
	else
	{
		/* We are inside. */

		return (!Test_Flag(this, INVERTED_FLAG));
	}
}



/*****************************************************************************
*
* FUNCTION
*
*   Disc_Normal
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Alexander Enzmann
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

void Disc::Normal (VECTOR Result, Intersection *, TraceThreadData *) const
{
	Assign_Vector(Result, normal);
}



/*****************************************************************************
*
* FUNCTION
*
*   Translate_Disc
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Alexander Enzmann
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

void Disc::Translate(const VECTOR, const TRANSFORM *tr)
{
	Transform(tr);
}



/*****************************************************************************
*
* FUNCTION
*
*   Rotate_Disc
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Alexander Enzmann
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

void Disc::Rotate(const VECTOR, const TRANSFORM *tr)
{
	Transform(tr);
}



/*****************************************************************************
*
* FUNCTION
*
*   Scale_Disc
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Alexander Enzmann
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

void Disc::Scale(const VECTOR, const TRANSFORM *tr)
{
	Transform(tr);
}



/*****************************************************************************
*
* FUNCTION
*
*   Invert_Disc
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Alexander Enzmann
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

void Disc::Invert()
{
	Invert_Flag(this, INVERTED_FLAG);
}



/*****************************************************************************
*
* FUNCTION
*
*   Transform_Disc
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Alexander Enzmann
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

void Disc::Transform(const TRANSFORM *tr)
{
	MTransNormal(normal, normal, tr);

	VNormalize(normal, normal);

	Compose_Transforms(Trans, tr);

	/* Recalculate the bounds */

	Compute_BBox();
}



/*****************************************************************************
*
* FUNCTION
*
*   Create_Disc
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Alexander Enzmann
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

Disc::Disc() : ObjectBase(DISC_OBJECT)
{
	Make_Vector (center, 0.0, 0.0, 0.0);
	Make_Vector (normal, 0.0, 0.0, 1.0);

	iradius2 = 0.0;
	oradius2 = 1.0;

	d = 0.0;

	Trans = Create_Transform();

	/* Default bounds */

	Make_BBox(BBox, -1.0, -1.0, -SMALL_TOLERANCE, 2.0,  2.0, 2.0 * SMALL_TOLERANCE);
}



/*****************************************************************************
*
* FUNCTION
*
*   Copy_Disc
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Alexander Enzmann
*
* DESCRIPTION
*
*   -
*
* CHANGES
*
*   Sep 1994 : Fixed memory leakage [DB]
*
******************************************************************************/

ObjectPtr Disc::Copy()
{
	Disc *New = new Disc();
	Destroy_Transform(New->Trans);
	*New = *this;
	New->Trans = Copy_Transform(Trans);

	return (New);
}



/*****************************************************************************
*
* FUNCTION
*
*   Destroy_Disc
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Alexander Enzmann
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

Disc::~Disc()
{
	Destroy_Transform(Trans);
}



/*****************************************************************************
*
* FUNCTION
*
*   Compute_Disc
*
* INPUT
*
*   Disc - Disc
*
* OUTPUT
*
*   Disc
*
* RETURNS
*
* AUTHOR
*
*   Dieter Bayer
*
* DESCRIPTION
*
*   Calculate the transformation that scales, rotates, and translates
*   the disc to the desired location and orientation.
*
* CHANGES
*
*   Aug 1994 : Creation.
*
******************************************************************************/

void Disc::Compute_Disc()
{
	Compute_Coordinate_Transform(Trans, center, normal, 1.0, 1.0);

	Compute_BBox();
}



/*****************************************************************************
*
* FUNCTION
*
*   Compute_Disc_BBox
*
* INPUT
*
*   Disc - Disc
*
* OUTPUT
*
*   Disc
*
* RETURNS
*
* AUTHOR
*
*   Dieter Bayer
*
* DESCRIPTION
*
*   Calculate the bounding box of a disc.
*
* CHANGES
*
*   Aug 1994 : Creation.
*
******************************************************************************/

void Disc::Compute_BBox()
{
	DBL rad;

	rad = sqrt(oradius2);

	Make_BBox(BBox, -rad, -rad, -SMALL_TOLERANCE, 2.0*rad, 2.0*rad, 2.0*SMALL_TOLERANCE);

	Recompute_BBox(&BBox, Trans);
}

}
