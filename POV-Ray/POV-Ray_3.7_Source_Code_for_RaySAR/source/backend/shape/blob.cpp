/*******************************************************************************
 * blob.cpp
 *
 * This module implements functions that manipulate blobs.
 *
 * The original file was written by Alexander Enzmann.
 * He wrote the code for blobs and generously provided us these enhancements.
 *
 * Modifications and enhancements by Dieter Bayer [DB].
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
 * $File: //depot/public/povray/3.x/source/backend/shape/blob.cpp $
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
*    blob
*    {
*      threshold THRESHOLD_VALUE
*
*      component STRENGTH, RADIUS, <CENTER>
*
*      sphere { <CENTER>, RADIUS, [strength] STRENGTH
*        [ translate <VECTOR> ]
*        [ rotate <VECTOR> ]
*        [ scale <VECTOR> ]
*        [ finish { ... } ]
*        [ pigment { ... } ]
*        [ tnormal { ... } ]
*        [ texture { ... } ]
*      }
*
*      cylinder { <END1>, <END2>, RADIUS, [strength] STRENGTH
*        [ translate <VECTOR> ]
*        [ rotate <VECTOR> ]
*        [ scale <VECTOR> ]
*        [ finish { ... } ]
*        [ pigment { ... } ]
*        [ tnormal { ... } ]
*        [ texture { ... } ]
*      }
*
*      [ sturm ]
*      [ hierarchy FLAG ]
*    }
*
*  ---
*
*  Jul 1994 : Most functions rewritten, bounding hierarchy added. [DB]
*
*  Aug 1994 : Cylindrical blobs added. [DB]
*
*  Sep 1994 : Multi-texturing added (each component can have its own texture).
*             Translation, rotation and scaling of each component added. [DB]
*
*  Oct 1994 : Adopted the method for the bounding slab creation to build the
*             bounding sphere hierarchy of the blob to get a much better
*             hierarchy. Improved bounding sphere calculation for tighter
*             bounds. [DB]
*
*  Dec 1994 : Added code for dynamic blob queue allocation. [DB]
*
*  Feb 1995 : Moved bounding sphere stuff into a seperate file. [DB]
*
*****************************************************************************/

// frame.h must always be the first POV file included (pulls in platform config)
#include "backend/frame.h"
#include "backend/math/vector.h"
#include "backend/shape/blob.h"
#include "backend/bounding/bbox.h"
#include "backend/bounding/bsphere.h"
#include "backend/math/matrices.h"
#include "backend/scene/objects.h"
#include "backend/math/polysolv.h"
#include "backend/texture/texture.h"
#include "backend/scene/threaddata.h"
#include "base/pov_err.h"

#include <algorithm>

// this must be the last file included
#include "base/povdebug.h"

namespace pov
{

/*****************************************************************************
* Local preprocessor defines
******************************************************************************/

/* Minimal intersection depth for a valid intersection. */
const DBL DEPTH_TOLERANCE = 1.0e-2;

/* Tolerance for inside test. */
const DBL INSIDE_TOLERANCE = 1.0e-6;

/* Ray enters/exits a component. */
const int ENTERING = 0;
const int EXITING  = 1;


/*****************************************************************************
*
* FUNCTION
*
*   All_Blob_Intersections
*
* INPUT
*
*   Object      - Object
*   Ray         - Ray
*
* OUTPUT
*
*   Depth_Stack - Intersection stack
*
* RETURNS
*
*   int - true, if a intersection was found
*   
* AUTHOR
*
*   Alexander Enzmann
*   
* DESCRIPTION
*
*   Generate intervals of influence for each component. After these
*   are made, determine their aggregate effect on the ray. As the
*   individual intervals are checked, a quartic is generated
*   that represents the density at a particular point on the ray.
*
*   Explanation for spherical components:
*
*   After making the substitutions in MakeBlob, there is a formula
*   for each component that has the form:
*
*      c0 * r^4 + c1 * r^2 + c2.
*
*   In order to determine the influence on the ray of all of the
*   individual components, we start by determining the distance
*   from any point on the ray to the specified point.  This can
*   be found using the pythagorean theorem, using C as the center
*   of this component, P as the start of the ray, and D as the
*   direction of travel of the ray:
*
*      r^2 = (t * D + P - C) . (t * D + P - C)
*
*   we insert this equation for each appearance of r^2 in the
*   components' formula, giving:
*
*      r^2 = D.D t^2 + 2 t D . (P - C) + (P - C) . (P - C)
*
*   Since the direction vector has been normalized, D.D = 1.
*   Using the substitutions:
*
*      t0 = (P - C) . (P - C),
*      t1 = D . (P - C)
*
*   We can write the formula as:
*
*      r^2 = t0 + 2 t t1 + t^2
*
*   Taking r^2 and substituting into the formula for this component
*   of the Blob we get the formula:
*
*      density = c0 * (r^2)^2 + c1 * r^2 + c2,
*
*   or:
*
*      density = c0 * (t0 + 2 t t1 + t^2)^2 +
*                c1 * (t0 + 2 t t1 + t^2) +
*                c2
*
*   Expanding terms and collecting with respect to "t" gives:
*
*      t^4 * c0 +
*      t^3 * 4 c0 t1 +
*      t^2 * (c1 + 2 * c0 t0 + 4 c0 t1^2)
*      t   * 2 (c1 t1 + 2 c0 t0 t1) +
*            c2 + c1*t0 + c0*t0^2
*
*   This formula can now be solved for "t" by any of the quartic
*   root solvers that are available.
*
* CHANGES
*
*   Jul 1994 : Added code for cylindrical and ellipsoidical blobs. [DB]
*
*   Oct 1994 : Added code to convert polynomial into a bezier curve for
*              a quick test if an intersection exists in an interval. [DB]
*
*   Sep 1995 : Added code to avoid numerical problems with distant blobs. [DB]
*
******************************************************************************/

bool Blob::All_Intersections(const Ray& ray, IStack& Depth_Stack, TraceThreadData *Thread)
{
	int i, j, cnt;
	int root_count, in_flag;
	int Intersection_Found = false;
	DBL t0, t1, t2, c0, c1, c2, dist, len, start_dist;
	DBL *fcoeffs;
	DBL coeffs[5];
	DBL roots[4];
	VECTOR P, D, V1, PP, DD;
	VECTOR IPoint;
	const Blob_Element *Element;
	DBL newcoeffs[5], dk[5];
	DBL max_bound;
	DBL l, w;
	DBL depthTolerance = (ray.IsSubsurfaceRay()? 0 : DEPTH_TOLERANCE);

	Thread->Stats()[Ray_Blob_Tests]++;

	/* Transform the ray into blob space. */

	if (Trans != NULL)
	{
		MInvTransPoint(P, ray.Origin, Trans);
		MInvTransDirection(D, ray.Direction, Trans);

		VLength(len, D);
		VInverseScaleEq(D, len);
	}
	else
	{
		Assign_Vector(P, ray.Origin);
		Assign_Vector(D, ray.Direction);

		len = 1.0;
	}

	/* Get the intervals along the ray where each component has an effect. */
	Blob_Interval_Struct* intervals = Thread->Blob_Intervals;
	if ((cnt = determine_influences(P, D, depthTolerance, intervals, Thread)) == 0)
	{
		/* Ray doesn't hit any of the component elements. */
		return (false);
	}

	/* To avoid numerical problems we start at the first interval. */

	if ((start_dist = intervals[0].bound) < SMALL_TOLERANCE)
	{
		start_dist = 0.0;
	}

	for (i = 0; i < cnt; i++)
	{
		intervals[i].bound -= start_dist;
	}

	VAddScaledEq(P, start_dist, D);

	max_bound = intervals[0].bound;
	for (i = 0; i < cnt; i++)
	{
		if (intervals[i].bound > max_bound)
			max_bound = intervals[i].bound;
	}

	/* Get the new starting point. */

	if (max_bound != 0)
	{
		VScaleEq(D, max_bound);

		for (i = 0; i < cnt; i++)
			intervals[i].bound /= max_bound;
	}
	else
		max_bound = 1;

	/* Clear out the coefficients. */

	coeffs[0] =
	coeffs[1] =
	coeffs[2] =
	coeffs[3] = 0.0;
	coeffs[4] = -Data->Threshold;

	/*
	 * Step through the list of intersection points, adding the
	 * influence of each component as it appears. 
	 */

	fcoeffs = NULL;

	for (i = in_flag = 0; i < cnt; i++)
	{
		if ((intervals[i].type & 1) == ENTERING)
		{
			/*
			 * Something is just starting to influence the ray, so calculate
			 * its coefficients and add them into the pot. 
			 */

			in_flag++;

			Element = intervals[i].Element;
			fcoeffs = &Thread->Blob_Coefficients[Element->index * 5];

			switch (Element->Type)
			{
				case BLOB_SPHERE:

					VSub(V1, P, Element->O);

					VDot(t0, V1, V1);
					VDot(t1, V1, D);
					t2 = max_bound * max_bound;

					c0 = Element->c[0];
					c1 = Element->c[1];
					c2 = Element->c[2];

					fcoeffs[0] = c0 * t2 * t2;
					fcoeffs[1] = 4.0 * c0 * t1 * t2;
					fcoeffs[2] = 2.0 * c0 * (2.0 * t1 * t1 + t0 * t2) + c1 * t2;
					fcoeffs[3] = 2.0 * t1 * (2.0 * c0 * t0 + c1);
					fcoeffs[4] = t0 * (c0 * t0 + c1) + c2;

					break;

				case BLOB_ELLIPSOID:

					MInvTransPoint(PP, P, Element->Trans);
					MInvTransDirection(DD, D, Element->Trans);

					VSub(V1, PP, Element->O);

					VDot(t0, V1, V1);
					VDot(t1, V1, DD);
					VDot(t2, DD, DD);

					c0 = Element->c[0];
					c1 = Element->c[1];
					c2 = Element->c[2];

					fcoeffs[0] = c0 * t2 * t2;
					fcoeffs[1] = 4.0 * c0 * t1 * t2;
					fcoeffs[2] = 2.0 * c0 * (2.0 * t1 * t1 + t0 * t2) + c1 * t2;
					fcoeffs[3] = 2.0 * t1 * (2.0 * c0 * t0 + c1);
					fcoeffs[4] = t0 * (c0 * t0 + c1) + c2;

					break;

				case BLOB_BASE_HEMISPHERE:
				case BLOB_APEX_HEMISPHERE:

					MInvTransPoint(PP, P, Element->Trans);
					MInvTransDirection(DD, D, Element->Trans);

					if (Element->Type == BLOB_APEX_HEMISPHERE)
					{
						PP[Z] -= Element->len;
					}

					VDot(t0, PP, PP);
					VDot(t1, PP, DD);
					VDot(t2, DD, DD);

					c0 = Element->c[0];
					c1 = Element->c[1];
					c2 = Element->c[2];

					fcoeffs[0] = c0 * t2 * t2;
					fcoeffs[1] = 4.0 * c0 * t1 * t2;
					fcoeffs[2] = 2.0 * c0 * (2.0 * t1 * t1 + t0 * t2) + c1 * t2;
					fcoeffs[3] = 2.0 * t1 * (2.0 * c0 * t0 + c1);
					fcoeffs[4] = t0 * (c0 * t0 + c1) + c2;

					break;

				case BLOB_CYLINDER:

					/* Transform ray into cylinder space. */

					MInvTransPoint(PP, P, Element->Trans);
					MInvTransDirection(DD, D, Element->Trans);

					t0 = PP[X] * PP[X] + PP[Y] * PP[Y];
					t1 = PP[X] * DD[X] + PP[Y] * DD[Y];
					t2 = DD[X] * DD[X] + DD[Y] * DD[Y];

					c0 = Element->c[0];
					c1 = Element->c[1];
					c2 = Element->c[2];

					fcoeffs[0] = c0 * t2 * t2;
					fcoeffs[1] = 4.0 * c0 * t1 * t2;
					fcoeffs[2] = 2.0 * c0 * (2.0 * t1 * t1 + t0 * t2) + c1 * t2;
					fcoeffs[3] = 2.0 * t1 * (2.0 * c0 * t0 + c1);
					fcoeffs[4] = t0 * (c0 * t0 + c1) + c2;

					break;

				default:

					throw POV_EXCEPTION_STRING("Unknown blob component in All_Blob_Intersections().");
			}

			for (j = 0; j < 5; j++)
			{
				coeffs[j] += fcoeffs[j];
			}
		}
		else
		{
			/* 
			 * We are losing the influence of a component -->
			 * subtract off its coefficients. 
			 */
			fcoeffs = &Thread->Blob_Coefficients[intervals[i].Element->index * 5];

			for (j = 0; j < 5; j++)
			{
				coeffs[j] -= fcoeffs[j];
			}

			/* If no components are currently affecting the ray ---> skip ahead. */

			if (--in_flag == 0)
			{
				continue;
			}
		}

		/*
		 * If the following intersection lies close to the current intersection
		 * then first add/subtract next region before testing. [DB 7/94] 
		 */

		if ((i + 1 < cnt) && (fabs(intervals[i].bound - intervals[i + 1].bound) < EPSILON))
		{
			continue;
		}

		/*
		 * Transform polynomial in a way that the interval boundaries are moved
		 * to 0 and 1, i. e. the roots of interest are between 0 and 1. [DB 10/94]
		 */

		l = intervals[i].bound;
		w = intervals[i+1].bound - l;

		newcoeffs[0] = coeffs[0] * w * w * w * w;
		newcoeffs[1] = (coeffs[1] + 4.0 * coeffs[0] * l) * w * w * w;
		newcoeffs[2] = (3.0 * l * (2.0 * coeffs[0] * l + coeffs[1]) + coeffs[2]) * w * w;
		newcoeffs[3] = (2.0 * l * (2.0 * l * (coeffs[0] * l + 0.75 * coeffs[1]) + coeffs[2]) + coeffs[3]) * w;
		newcoeffs[4] = l * (l * (l * (coeffs[0] * l + coeffs[1]) + coeffs[2]) + coeffs[3]) + coeffs[4];

		/* Calculate coefficients of corresponding bezier curve. [DB 10/94] */

		dk[0] = newcoeffs[4];
		dk[1] = newcoeffs[4] + 0.25 * newcoeffs[3];
		dk[2] = newcoeffs[4] + 0.50 * (newcoeffs[3] + newcoeffs[2] / 3.0);
		dk[3] = newcoeffs[4] + 0.50 * (1.5 * newcoeffs[3] + newcoeffs[2] + 0.5 * newcoeffs[1]);
		dk[4] = newcoeffs[4] + newcoeffs[3] + newcoeffs[2] + newcoeffs[1] + newcoeffs[0];

		/*
		 * Skip this interval if the ray doesn't intersect the convex hull of the
		 * bezier curve, because no valid intersection will be found. [DB 10/94]
		 */

		if (((dk[0] >= 0.0) && (dk[1] >= 0.0) && (dk[2] >= 0.0) && (dk[3] >= 0.0) && (dk[4] >= 0.0)) ||
		    ((dk[0] <= 0.0) && (dk[1] <= 0.0) && (dk[2] <= 0.0) && (dk[3] <= 0.0) && (dk[4] <= 0.0)))
		{
			continue;
		}

		/*
		 * Now we could do bezier clipping to find the roots
		 * but I have no idea how this works. [DB 2/95]
		 */


		/* Solve polynomial. */

		root_count = Solve_Polynomial(4, coeffs, roots, Test_Flag(this, STURM_FLAG), 1.0e-11, Thread);

		/* See if any of the roots are valid. */

		for (j = 0; j < root_count; j++)
		{
			dist = roots[j];

			/*
			 * First see if the root is in the interval of influence of
			 * the currently active components.
			 */

			if ((dist >= intervals[i].bound) &&
			    (dist <= intervals[i+1].bound))
			{
				/* Correct distance. */

				dist = (dist * max_bound + start_dist) / len;

				if ((dist > depthTolerance) && (dist < MAX_DISTANCE))
				{
					VEvaluateRay(IPoint, ray.Origin, dist, ray.Direction);

					if (Clip.empty() || Point_In_Clip(IPoint, Clip, Thread))
					{
						Depth_Stack->push(Intersection(dist, IPoint, this));

						Intersection_Found = true;
					}
				}
			}
		}

		/*
		 * If the blob isn't used inside a CSG and we have found at least
		 * one intersection then we are ready, because all possible intersections
		 * will be further away (we have a sorted list!). [DB 7/94]
		 */

		if (!(Type & IS_CHILD_OBJECT) && (Intersection_Found))
		{
			break;
		}
	}

	if (Intersection_Found)
		Thread->Stats()[Ray_Blob_Tests_Succeeded]++;

	return (Intersection_Found);
}



/*****************************************************************************
*
* FUNCTION
*
*   insert_hit
*
* INPUT
*
*   Blob      - Pointer to blob structure
*   Element   - Element to insert
*   t0, t1    - Intersection depths
*
* OUTPUT
*
*   intervals - Pointer to sorted list of hits
*   cnt       - Number of hits in intervals
*
* RETURNS
*
* AUTHOR
*
*   Alexander Enzmann
*
* DESCRIPTION
*
*   Store the points of intersection. Keep track of: whether this is
*   the start or end point of the hit, which component was pierced
*   by the ray, and the point along the ray that the hit occured at.
*
* CHANGES
*
*   Oct 1994 : Modified to use memmove instead of loops for copying. [DB]
*   Sep 1995 : Changed to allow use of memcpy if memmove isn't available. [AED]
*   Jul 1996 : Changed to use POV_MEMMOVE, which can be memmove or pov_memmove.
*   Oct 1996 : Changed to avoid unnecessary compares. [DB]
*
******************************************************************************/

void Blob::insert_hit(const Blob_Element *Element, DBL t0, DBL t1, Blob_Interval_Struct *intervals, unsigned int *cnt)
{
	unsigned int k;

	/* We are entering the component. */

	intervals[*cnt].type    = Element->Type | ENTERING;
	intervals[*cnt].bound   = t0;
	intervals[*cnt].Element = Element;

	for (k = 0; t0 > intervals[k].bound; k++);

	if (k < *cnt)
	{
		/*
		 * This hit point is smaller than one that already exists -->
		 * bump the rest and insert it here.
		 */

		POV_MEMMOVE(&intervals[k+1], &intervals[k], (*cnt-k)*sizeof(Blob_Interval_Struct));

		/* We are entering the component. */

		intervals[k].type    = Element->Type | ENTERING;
		intervals[k].bound   = t0;
		intervals[k].Element = Element;

		(*cnt)++;

		/* We are exiting the component. */

		intervals[*cnt].type    = Element->Type | EXITING;
		intervals[*cnt].bound   = t1;
		intervals[*cnt].Element = Element;

		for (k = k + 1; t1 > intervals[k].bound; k++);

		if (k < *cnt)
		{
			POV_MEMMOVE(&intervals[k+1], &intervals[k], (*cnt-k)*sizeof(Blob_Interval_Struct));

			/* We are exiting the component. */

			intervals[k].type    = Element->Type | EXITING;
			intervals[k].bound   = t1;
			intervals[k].Element = Element;
		}

		(*cnt)++;
	}
	else
	{
		/* Just plop the start and end points at the end of the list */

		(*cnt)++;

		/* We are exiting the component. */

		intervals[*cnt].type    = Element->Type | EXITING;
		intervals[*cnt].bound   = t1;
		intervals[*cnt].Element = Element;

		(*cnt)++;
	}
}



/*****************************************************************************
*
* FUNCTION
*
*   intersect_cylinder
*
* INPUT
*
*   Element    - Pointer to element structure
*   P, D       - Ray = P + t * D
*   mindist    - Min. valid distance
*
* OUTPUT
*
*   tmin, tmax - Intersection depths found
*
* RETURNS
*
* AUTHOR
*
*   Dieter Bayer (with help from Alexander Enzmann)
*
* DESCRIPTION
*
*   -
*
* CHANGES
*
*   Jul 1994 : Creation.
*
******************************************************************************/

int Blob::intersect_cylinder(const Blob_Element *Element, const VECTOR P, const VECTOR D, DBL mindist, DBL *tmin, DBL *tmax)
{
	DBL a, b, c, d, t, u, v, w, len;
	VECTOR PP, DD;

	/* Transform ray into cylinder space. */

	MInvTransPoint(PP, P, Element->Trans);
	MInvTransDirection(DD, D, Element->Trans);

	VLength(len, DD);
	VInverseScaleEq(DD, len);

	/* Intersect ray with cylinder. */

	a = DD[X] * DD[X] + DD[Y] * DD[Y];

	if (a > EPSILON)
	{
		b = PP[X] * DD[X] + PP[Y] * DD[Y];
		c = PP[X] * PP[X] + PP[Y] * PP[Y] - Element->rad2;

		d = b * b - a * c;

		if (d > EPSILON)
		{
			d = sqrt(d);

			t = ( - b + d) / a;

			w = PP[Z] + t * DD[Z];

			if ((w >= 0.0) && (w <= Element->len))
			{
				if (t < *tmin) { *tmin = t; }
				if (t > *tmax) { *tmax = t; }
			}

			t = ( - b - d) / a;

			w = PP[Z] + t * DD[Z];

			if ((w >= 0.0) && (w <= Element->len))
			{
				if (t < *tmin) { *tmin = t; }
				if (t > *tmax) { *tmax = t; }
			}
		}
	}

	/* Intersect base/cap plane. */

	if (fabs(DD[Z]) > EPSILON)
	{
		/* Intersect base plane. */

		t = - PP[Z] / DD[Z];

		u = PP[X] + t * DD[X];
		v = PP[Y] + t * DD[Y];

		if ((u * u + v * v) <= Element->rad2)
		{
			if (t < *tmin) { *tmin = t; }
			if (t > *tmax) { *tmax = t; }
		}

		/* Intersect cap plane. */

		t = (Element->len - PP[Z]) / DD[Z];

		u = PP[X] + t * DD[X];
		v = PP[Y] + t * DD[Y];

		if ((u * u + v * v) <= Element->rad2)
		{
			if (t < *tmin) { *tmin = t; }
			if (t > *tmax) { *tmax = t; }
		}
	}

	/* Check if the intersections are valid. */

	*tmin /= len;
	*tmax /= len;

	if (*tmin < mindist) { *tmin = 0.0; }
	if (*tmax < mindist) { *tmax = 0.0; }

	if (*tmin >= *tmax)
	{
		return (false);
	}

	return (true);
}



/*****************************************************************************
*
* FUNCTION
*
*   intersect_ellipsoid
*
* INPUT
*
*   Element    - Pointer to element structure
*   P, D       - Ray = P + t * D
*   mindist    - Min. valid distance
*
* OUTPUT
*
*   tmin, tmax - Intersection depths found
*
* RETURNS
*
* AUTHOR
*
*   Dieter Bayer
*
* DESCRIPTION
*
*   -
*
* CHANGES
*
*   Sep 1994 : Creation.
*
******************************************************************************/

int Blob::intersect_ellipsoid(const Blob_Element *Element, const VECTOR P, const VECTOR D, DBL mindist, DBL *tmin, DBL *tmax)
{
	DBL b, d, t, len;
	VECTOR V1, PP, DD;

	MInvTransPoint(PP, P, Element->Trans);
	MInvTransDirection(DD, D, Element->Trans);

	VLength(len, DD);
	VInverseScaleEq(DD, len);

	VSub(V1, PP, Element->O);
	VDot(b, V1, DD);
	VDot(t, V1, V1);

	d = b * b - t + Element->rad2;

	if (d < EPSILON)
	{
		return (false);
	}

	d = sqrt(d);

	*tmax = ( - b + d) / len;  if (*tmax < mindist) { *tmax = 0.0; }
	*tmin = ( - b - d) / len;  if (*tmin < mindist) { *tmin = 0.0; }

	if (*tmax == *tmin)
	{
		return (false);
	}
	else
	{
		if (*tmax < *tmin)
		{
			d = *tmin;  *tmin = *tmax;  *tmax = d;
		}
	}

	return (true);
}



/*****************************************************************************
*
* FUNCTION
*
*   intersect_hemisphere
*
* INPUT
*
*   Element    - Pointer to element structure
*   P, D       - Ray = P + t * D
*   mindist    - Min. valid distance
*
* OUTPUT
*
*   tmin, tmax - Intersection depths found
*
* RETURNS
*
* AUTHOR
*
*   Dieter Bayer
*
* DESCRIPTION
*
*   -
*
* CHANGES
*
*   Jul 1994 : Creation (with help from Alexander Enzmann).
*
******************************************************************************/

int Blob::intersect_hemisphere(const Blob_Element *Element, const VECTOR P, const VECTOR D, DBL mindist, DBL *tmin, DBL *tmax)
{
	DBL b, d, t, z1, z2, len;
	VECTOR PP, DD;

	/* Transform ray into hemisphere space. */

	MInvTransPoint(PP, P, Element->Trans);
	MInvTransDirection(DD, D, Element->Trans);

	VLength(len, DD);
	VInverseScaleEq(DD, len);

	if (Element->Type == BLOB_BASE_HEMISPHERE)
	{
		VDot(b, PP, DD);
		VDot(t, PP, PP);

		d = b * b - t + Element->rad2;

		if (d < EPSILON)
		{
			return (false);
		}

		d = sqrt(d);

		*tmax = - b + d;
		*tmin = - b - d;

		if (*tmax < *tmin)
		{
			d = *tmin;  *tmin = *tmax;  *tmax = d;
		}

		/* Cut intersection at the plane. */

		z1 = PP[Z] + *tmin * DD[Z];
		z2 = PP[Z] + *tmax * DD[Z];

		/* If both points are inside --> no intersection */

		if ((z1 >= 0.0) && (z2 >= 0.0))
		{
			return (false);
		}

		/* If both points are outside --> intersections found */

		if ((z1 < 0.0) && (z2 < 0.0))
		{
			*tmin /= len;
			*tmax /= len;

			return (true);
		}

		/* Determine intersection with plane. */

		t = - PP[Z] / DD[Z];

		if (z1 >= 0.0)
		{
			/* Ray is crossing the plane from inside to outside. */

			*tmin = (t < mindist) ? 0.0 : t;
		}
		else
		{
			/* Ray is crossing the plane from outside to inside. */

			*tmax = (t < mindist) ? 0.0 : t;
		}

		*tmin /= len;
		*tmax /= len;

		return (true);
	}
	else
	{
		PP[Z] -= Element->len;

		VDot(b, PP, DD);
		VDot(t, PP, PP);

		d = b * b - t + Element->rad2;

		if (d < EPSILON)
		{
			return (false);
		}

		d = sqrt(d);

		*tmax = - b + d;
		*tmin = - b - d;

		if (*tmax < *tmin)
		{
			d = *tmin;  *tmin = *tmax;  *tmax = d;
		}

		/* Cut intersection at the plane. */

		z1 = PP[Z] + *tmin * DD[Z];
		z2 = PP[Z] + *tmax * DD[Z];

		/* If both points are inside --> no intersection */

		if ((z1 <= 0.0) && (z2 <= 0.0))
		{
			return (false);
		}

		/* If both points are outside --> intersections found */

		if ((z1 > 0.0) && (z2 > 0.0))
		{
			*tmin /= len;
			*tmax /= len;

			return (true);
		}

		/* Determine intersection with plane. */

		t = - PP[Z] / DD[Z];

		if (z1 <= 0.0)
		{
			/* Ray is crossing the plane from inside to outside. */

			*tmin = (t < mindist) ? 0.0 : t;
		}
		else
		{
			/* Ray is crossing the plane from outside to inside. */

			*tmax = (t < mindist) ? 0.0 : t;
		}

		*tmin /= len;
		*tmax /= len;

		return (true);
	}
}



/*****************************************************************************
*
* FUNCTION
*
*   intersect_sphere
*
* INPUT
*
*   Element    - Pointer to element structure
*   P, D       - Ray = P + t * D
*   mindist    - Min. valid distance
*
* OUTPUT
*
*   tmin, tmax - Intersection depths found
*
* RETURNS
*
* AUTHOR
*
*   Dieter Bayer
*
* DESCRIPTION
*
*   -
*
* CHANGES
*
*   Jul 1994 : Creation (with help from Alexander Enzmann).
*
******************************************************************************/

int Blob::intersect_sphere(const Blob_Element *Element, const VECTOR P, const VECTOR D, DBL mindist, DBL *tmin, DBL *tmax)
{
	DBL b, d, t;
	VECTOR V1;

	VSub(V1, P, Element->O);
	VDot(b, V1, D);
	VDot(t, V1, V1);

	d = b * b - t + Element->rad2;

	if (d < EPSILON)
	{
		return (false);
	}

	d = sqrt(d);

	*tmax = - b + d;  if (*tmax < mindist) { *tmax = 0.0; }
	*tmin = - b - d;  if (*tmin < mindist) { *tmin = 0.0; }

	if (*tmax == *tmin)
	{
		return (false);
	}
	else
	{
		if (*tmax < *tmin)
		{
			d = *tmin;  *tmin = *tmax;  *tmax = d;
		}
	}

	return (true);
}



/*****************************************************************************
*
* FUNCTION
*
*   intersect_element
*
* INPUT
*
*   P, D       - Ray = P + t * D
*   Element    - Pointer to element structure
*   mindist    - Min. valid distance
*
* OUTPUT
*
*   tmin, tmax - Intersection depths found
*
* RETURNS
*
* AUTHOR
*
*   Dieter Bayer
*
* DESCRIPTION
*
*   -
*
* CHANGES
*
*   Jul 1994 : Creation.
*
******************************************************************************/

int Blob::intersect_element(const VECTOR P, const VECTOR D, const Blob_Element *Element, DBL mindist, DBL *tmin, DBL *tmax, TraceThreadData *Thread)
{
#ifdef BLOB_EXTRA_STATS
	Thread->Stats()[Blob_Element_Tests]++;
#endif

	*tmin = BOUND_HUGE;
	*tmax = - BOUND_HUGE;

	switch (Element->Type)
	{
		case BLOB_SPHERE:

			if (!intersect_sphere(Element, P, D, mindist, tmin, tmax))
			{
				return (false);
			}

			break;

		case BLOB_ELLIPSOID:

			if (!intersect_ellipsoid(Element, P, D, mindist, tmin, tmax))
			{
				return (false);
			}

			break;

		case BLOB_BASE_HEMISPHERE:
		case BLOB_APEX_HEMISPHERE:

			if (!intersect_hemisphere(Element, P, D, mindist, tmin, tmax))
			{
				return (false);
			}

			break;

		case BLOB_CYLINDER:

			if (!intersect_cylinder(Element, P, D, mindist, tmin, tmax))
			{
				return (false);
			}

			break;
	}

#ifdef BLOB_EXTRA_STATS
	Thread->Stats()[Blob_Element_Tests_Succeeded]++;
#endif

	return (true);
}



/*****************************************************************************
*
* FUNCTION
*
*   determine_influences
*
* INPUT
*
*   P, D       - Ray = P + t * D
*   Blob       - Pointer to blob structure
*   mindist    - Min. valid distance
*
* OUTPUT
*
*   intervals  - Sorted list of intersections found
*
* RETURNS
*
*   int - Number of intersection found
*
* AUTHOR
*
*   Dieter Bayer
*
* DESCRIPTION
*
*   Make a sorted list of points along the ray at which the various blob
*   components start and stop adding their influence.
*
* CHANGES
*
*   Jul 1994 : Added code for bounding hierarchy traversal. [DB]
*
******************************************************************************/

int Blob::determine_influences(const VECTOR P, const VECTOR  D, DBL mindist, Blob_Interval_Struct *intervals, TraceThreadData *Thread) const
{
	int i;
	unsigned int cnt, size;
	DBL b, t, t0, t1;
	VECTOR V1;
	BSPHERE_TREE *Tree;
	BSPHERE_TREE **Queue = reinterpret_cast<BSPHERE_TREE **>(Thread->Blob_Queue);

	cnt = 0;

	if (Data->Tree == NULL)
	{
		/* There's no bounding hierarchy so just step through all elements. */

		for (i = 0; i < Data->Number_Of_Components; i++)
		{
			if (intersect_element(P, D, &Data->Entry[i], mindist, &t0, &t1, Thread))
			{
				insert_hit(&Data->Entry[i], t0, t1, intervals, &cnt);
			}
		}
	}
	else
	{
		/* Use blob's bounding hierarchy. */

		size = 0;

		Queue[size++] = Data->Tree;

		while (size > 0)
		{
			Tree = Queue[--size];

			/* Test if current node is a leaf. */

			if (Tree->Entries <= 0)
			{
				/* Test element. */

				if (intersect_element(P, D, reinterpret_cast<Blob_Element *>(Tree->Node), mindist, &t0, &t1, Thread))
				{
					insert_hit(reinterpret_cast<Blob_Element *>(Tree->Node), t0, t1, intervals, &cnt);
				}
			}
			else
			{
				/* Test all sub-nodes. */

				for (i = 0; i < (int)Tree->Entries; i++)
				{
#ifdef BLOB_EXTRA_STATS
					Thread->Stats()[Blob_Bound_Tests]++;
#endif

					VSub(V1, Tree->Node[i]->C, P);
					VDot(b, V1, D);
					VDot(t, V1, V1);

					if ((t - Sqr(b)) <= Tree->Node[i]->r2)
					{
#ifdef BLOB_EXTRA_STATS
						Thread->Stats()[Blob_Bound_Tests_Succeeded]++;
#endif

						if (insert_node(Tree->Node[i], &size, Thread))
							Queue = reinterpret_cast<BSPHERE_TREE **>(Thread->Blob_Queue);
					}
				}
			}
		}
	}

	return (cnt);
}



/*****************************************************************************
*
* FUNCTION
*
*   calculate_element_field
*
* INPUT
*
*   Element - Pointer to element structure
*   P       - Point whos field value is calculated
*
* OUTPUT
*
* RETURNS
*
*   DBL - Field value
*
* AUTHOR
*
*   Alexander Enzmann
*
* DESCRIPTION
*
*   Calculate the field value of a single element in a given point P
*   (which must already have been transformed into blob space).
*
* CHANGES
*
*   Jul 1994 : Added code for cylindrical and ellipsoidical blobs. [DB]
*
******************************************************************************/

DBL Blob::calculate_element_field(const Blob_Element *Element, const VECTOR P)
{
	DBL rad2, density;
	VECTOR V1, PP;

	density = 0.0;

	switch (Element->Type)
	{
		case BLOB_SPHERE:

			VSub(V1, P, Element->O);

			VDot(rad2, V1, V1);

			if (rad2 < Element->rad2)
			{
				density = rad2 * (rad2 * Element->c[0] + Element->c[1]) + Element->c[2];
			}

			break;

		case BLOB_ELLIPSOID:

			MInvTransPoint(PP, P, Element->Trans);

			VSub(V1, PP, Element->O);

			VDot(rad2, V1, V1);

			if (rad2 < Element->rad2)
			{
				density = rad2 * (rad2 * Element->c[0] + Element->c[1]) + Element->c[2];
			}

			break;

		case BLOB_BASE_HEMISPHERE:

			MInvTransPoint(PP, P, Element->Trans);

			if (PP[Z] <= 0.0)
			{
				VDot(rad2, PP, PP);

				if (rad2 <= Element->rad2)
				{
					density = rad2 * (rad2 * Element->c[0] + Element->c[1]) + Element->c[2];
				}
			}

			break;

		case BLOB_APEX_HEMISPHERE:

			MInvTransPoint(PP, P, Element->Trans);

			PP[Z] -= Element->len;

			if (PP[Z] >= 0.0)
			{
				VDot(rad2, PP, PP);

				if (rad2 <= Element->rad2)
				{
					density = rad2 * (rad2 * Element->c[0] + Element->c[1]) + Element->c[2];
				}
			}

			break;

		case BLOB_CYLINDER:

			MInvTransPoint(PP, P, Element->Trans);

			if ((PP[Z] >= 0.0) && (PP[Z] <= Element->len))
			{
				if ((rad2 = Sqr(PP[X]) + Sqr(PP[Y])) <= Element->rad2)
				{
					density = rad2 * (rad2 * Element->c[0] + Element->c[1]) + Element->c[2];
				}
			}

			break;
	}

	return (density);
}



/*****************************************************************************
*
* FUNCTION
*
*   calculate_field_value
*
* INPUT
*
*   Blob - Pointer to blob structure
*   P       - Point whos field value is calculated
*
* OUTPUT
*
* RETURNS
*
*   DBL - Field value
*
* AUTHOR
*
*   Dieter Bayer
*
* DESCRIPTION
*
*   Calculate the field value of a blob in a given point P
*   (which must already have been transformed into blob space).
*
* CHANGES
*
*   Jul 1994 : Added code for bounding hierarchy traversal. [DB]
*
******************************************************************************/

DBL Blob::calculate_field_value(const VECTOR P, TraceThreadData *Thread) const
{
	int i;
	unsigned int size;
	DBL density, rad2;
	VECTOR V1;
	BSPHERE_TREE *Tree;
	BSPHERE_TREE **Queue = reinterpret_cast<BSPHERE_TREE **>(Thread->Blob_Queue);

	density = 0.0;

	if (Data->Tree == NULL)
	{
		/* There's no tree --> step through all elements. */

		for (i = 0; i < Data->Number_Of_Components; i++)
		{
			density += calculate_element_field(&Data->Entry[i], P);
		}
	}
	else
	{
		/* A tree exists --> step through the tree. */

		size = 0;

		Queue[size++] = Data->Tree;

		while (size > 0)
		{
			Tree = Queue[--size];

			/* Test if current node is a leaf. */

			if (Tree->Entries <= 0)
			{
				density += calculate_element_field(reinterpret_cast<Blob_Element *>(Tree->Node), P);
			}
			else
			{
				/* Test all sub-nodes. */

				for (i = 0; i < (int)Tree->Entries; i++)
				{
					/* Insert sub-node if we are inside. */

					VSub(V1, P, Tree->Node[i]->C);

					VDot(rad2, V1, V1);

					if (rad2 <= Tree->Node[i]->r2)
						if (insert_node(Tree->Node[i], &size, Thread))
							Queue = reinterpret_cast<BSPHERE_TREE **>(Thread->Blob_Queue);
				}
			}
		}
	}

	return (density);
}



/*****************************************************************************
*
* FUNCTION
*
*   Inside_Blob
*
* INPUT
*
*   Test_Point - Point to test
*   Object     - Pointer to blob structure
*
* OUTPUT
*
* RETURNS
*
*   int - true if Test_Point is inside
*
* AUTHOR
*
*   Alexander Enzmann
*
* DESCRIPTION
*
*   Calculate the density at the given point and then compare to
*   the threshold to see if we are in or out of the blob.
*
* CHANGES
*
*   -
*
******************************************************************************/

bool Blob::Inside(const VECTOR Test_Point, TraceThreadData *Thread) const
{
	VECTOR New_Point;

	/* Transform the point into blob space. */

	if (Trans != NULL)
	{
		MInvTransPoint(New_Point, Test_Point, Trans);
	}
	else
	{
		Assign_Vector(New_Point, Test_Point);
	}

	if (calculate_field_value(New_Point, Thread) > Data->Threshold - INSIDE_TOLERANCE)
	{
		/* We are inside. */

		return (!Test_Flag(this, INVERTED_FLAG));
	}
	else
	{
		/* We are outside. */

		return (Test_Flag(this, INVERTED_FLAG));
	}
}



/*****************************************************************************
*
* FUNCTION
*
*   element_normal
*
* INPUT
*
*   P       - Surface point
*   Element - Pointer to element structure
*
* OUTPUT
*
*   Result  - Element's normal
*
* RETURNS
*
* AUTHOR
*
*   Dieter Bayer
*
* DESCRIPTION
*
*   Calculate the normal of a single element in the point P.
*
* CHANGES
*
*   Jul 1994 : Creation (with help from Alexander Enzmann).
*
******************************************************************************/

void Blob::element_normal(VECTOR Result, const VECTOR P, const Blob_Element *Element)
{
	DBL val, dist;
	VECTOR V1, PP;

	switch (Element->Type)
	{
		case BLOB_SPHERE:

			VSub(V1, P, Element->O);

			VDot(dist, V1, V1);

			if (dist <= Element->rad2)
			{
				val = -2.0 * Element->c[0] * dist - Element->c[1];

				VAddScaledEq(Result, val, V1);
			}

			break;

		case BLOB_ELLIPSOID:

			MInvTransPoint(PP, P, Element->Trans);

			VSub(V1, PP, Element->O);

			VDot(dist, V1, V1);

			if (dist <= Element->rad2)
			{
				val = -2.0 * Element->c[0] * dist - Element->c[1];

				MTransNormal(V1, V1, Element->Trans);

				VAddScaledEq(Result, val, V1);
			}

			break;

		case BLOB_BASE_HEMISPHERE:

			MInvTransPoint(PP, P, Element->Trans);

			if (PP[Z] <= 0.0)
			{
				VDot(dist, PP, PP);

				if (dist <= Element->rad2)
				{
					val = -2.0 * Element->c[0] * dist - Element->c[1];

					MTransNormal(PP, PP, Element->Trans);

					VAddScaledEq(Result, val, PP);
				}
			}

			break;

		case BLOB_APEX_HEMISPHERE:

			MInvTransPoint(PP, P, Element->Trans);

			PP[Z] -= Element->len;

			if (PP[Z] >= 0.0)
			{
				VDot(dist, PP, PP);

				if (dist <= Element->rad2)
				{
					val = -2.0 * Element->c[0] * dist - Element->c[1];

					MTransNormal(PP, PP, Element->Trans);

					VAddScaledEq(Result, val, PP);
				}
			}

			break;

		case BLOB_CYLINDER:

			MInvTransPoint(PP, P, Element->Trans);

			if ((PP[Z] >= 0.0) && (PP[Z] <= Element->len))
			{
				if ((dist = Sqr(PP[X]) + Sqr(PP[Y])) <= Element->rad2)
				{
					val = -2.0 * Element->c[0] * dist - Element->c[1];

					PP[Z] = 0.0;

					MTransNormal(PP, PP, Element->Trans);

					VAddScaledEq(Result, val, PP);
				}
			}

			break;
	}
}



/*****************************************************************************
*
* FUNCTION
*
*   Blob_Normal
*
* INPUT
*
*   Object  - Pointer to blob structure
*   Inter   - Pointer to intersection
*
* OUTPUT
*
*   Result  - Blob's normal
*
* RETURNS
*
* AUTHOR
*
*   Alexander Enzmann
*
* DESCRIPTION
*
*   Calculate the blob's surface normal in the intersection point.
*
* CHANGES
*
*   Jul 1994 : Added code for bounding hierarchy traversal. [DB]
*
******************************************************************************/

void Blob::Normal(VECTOR Result, Intersection *Inter, TraceThreadData *Thread) const
{
	int i;
	unsigned int size;
	DBL dist, val;
	VECTOR New_Point, V1;
	BSPHERE_TREE *Tree;
	BSPHERE_TREE **Queue = reinterpret_cast<BSPHERE_TREE **>(Thread->Blob_Queue);

	/* Transform the point into the blob space. */
	getLocalIPoint(New_Point, Inter);

	Make_Vector(Result, 0.0, 0.0, 0.0);

	/* For each component that contributes to this point, add its bit to the normal */

	if (Data->Tree == NULL)
	{
		/* There's no tree --> step through all elements. */

		for (i = 0; i < Data->Number_Of_Components; i++)
		{
			element_normal(Result, New_Point, &(Data->Entry[i]));
		}
	}
	else
	{
		/* A tree exists --> step through the tree. */

		size = 0;

		Queue[size++] = Data->Tree;

		while (size > 0)
		{
			Tree = Queue[--size];

			/* Test if current node is a leaf. */

			if (Tree->Entries <= 0)
			{
				element_normal(Result, New_Point, reinterpret_cast<Blob_Element *>(Tree->Node));
			}
			else
			{
				/* Test all sub-nodes. */

				for (i = 0; i < (int)Tree->Entries; i++)
				{
					/* Insert sub-node if we are inside. */

					VSub(V1, New_Point, Tree->Node[i]->C);

					VDot(dist, V1, V1);

					if (dist <= Tree->Node[i]->r2)
						if (insert_node(Tree->Node[i], &size, Thread))
							Queue = reinterpret_cast<BSPHERE_TREE **>(Thread->Blob_Queue);
				}
			}
		}
	}

	VDot(val, Result, Result);

	if (val == 0.0)
	{
		Make_Vector(Result, 1.0, 0.0, 0.0);
	}
	else
	{
		/* Normalize normal vector. */

		val = 1.0 / sqrt(val);

		VScaleEq(Result, val);
	}

	/* Transform back to world space. */

	if (Trans != NULL)
	{
		MTransNormal(Result, Result, Trans);

		VNormalize(Result, Result);
	}
}



/*****************************************************************************
*
* FUNCTION
*
*   Translate_Blob
*
* INPUT
*
*   Vector - Translation vector
*
* OUTPUT
*
*   Object - Pointer to blob structure
*
* RETURNS
*
* AUTHOR
*
*   Alexander Enzmann
*
* DESCRIPTION
*
*   Translate a blob.
*
* CHANGES
*
*   -
*
******************************************************************************/

void Blob::Translate(const VECTOR, const TRANSFORM *tr)
{
	Transform(tr);
}



/*****************************************************************************
*
* FUNCTION
*
*   Rotate_Blob
*
* INPUT
*
*   Vector - Rotation vector
*
* OUTPUT
*
*   Object - Pointer to blob structure
*
* RETURNS
*
* AUTHOR
*
*   Alexander Enzmann
*
* DESCRIPTION
*
*   Rotate a blob.
*
* CHANGES
*
*   -
*
******************************************************************************/

void Blob::Rotate(const VECTOR, const TRANSFORM *tr)
{
	Transform(tr);
}



/*****************************************************************************
*
* FUNCTION
*
*   Scale_Blob
*
* INPUT
*
*   Vector - Scaling vector
*
* OUTPUT
*
*   Object - Pointer to blob structure
*
* RETURNS
*
* AUTHOR
*
*   Alexander Enzmann
*
* DESCRIPTION
*
*   Scale a blob.
*
* CHANGES
*
*   -
*
******************************************************************************/

void Blob::Scale(const VECTOR, const TRANSFORM *tr)
{
	Transform(tr);
}



/*****************************************************************************
*
* FUNCTION
*
*   Transform_Blob
*
* INPUT
*
*   Trans  - Pointer to transformation
*
* OUTPUT
*
*   Object - Pointer to blob structure
*
* RETURNS
*   
* AUTHOR
*
*   Alexander Enzmann
*   
* DESCRIPTION
*
*   Transform a blob.
*
* CHANGES
*
*   -
*
******************************************************************************/

void Blob::Transform(const TRANSFORM *tr)
{
	int i;

	if(Trans == NULL)
		Trans = Create_Transform();

	Recompute_BBox(&BBox, tr);

	Compose_Transforms(Trans, tr);

	for(i = 0; i < Data->Number_Of_Components; i++)
		Transform_Textures(Element_Texture[i], tr);
}



/*****************************************************************************
*
* FUNCTION
*
*   Invert_Blob
*
* INPUT
*
*   Object - Pointer to blob structure
*
* OUTPUT
*
*   Object
*
* RETURNS
*
* AUTHOR
*
*   Alexander Enzmann
*
* DESCRIPTION
*
*   Invert a blob.
*
* CHANGES
*
*   -
*
******************************************************************************/

void Blob::Invert()
{
	Invert_Flag(this, INVERTED_FLAG);
}



/*****************************************************************************
*
* FUNCTION
*
*   Create_Blob
*
* INPUT
*
*   Object - Pointer to blob structure
*
* OUTPUT
*
*   Object
*
* RETURNS
*
* AUTHOR
*
*   Alexander Enzmann
*
* DESCRIPTION
*
*   Create a new blob.
*
* CHANGES
*
*   -
*
******************************************************************************/

Blob::Blob() : ObjectBase(BLOB_OBJECT)
{
	Set_Flag(this, HIERARCHY_FLAG);
	Trans = NULL;
	Element_Texture = NULL;
	Data = NULL ;
}

/*****************************************************************************
*
* FUNCTION
*
*   Copy_Blob
*
* INPUT
*
*   Object - Pointer to blob structure
*
* OUTPUT
*
*   Object
*
* RETURNS
*
* AUTHOR
*
*   Alexander Enzmann
*
* DESCRIPTION
*
*   Copy a blob.
*
*   NOTE: The components are not copied, only the number of references is
*         counted, so that Destroy_Blob() knows if they can be destroyed.
*
* CHANGES
*
*   Jul 1994 : Added code for blob data reference counting. [DB]
*
******************************************************************************/

ObjectPtr Blob::Copy()
{
	int i;
	Blob *New = new Blob();

	/* Copy blob. */

	Destroy_Transform(New->Trans);
	*New = *this;
	New->Trans = Copy_Transform(Trans);
	New->Data = Data->AcquireReference () ;
	New->Element_Texture = reinterpret_cast<TEXTURE **>(POV_MALLOC(New->Data->Number_Of_Components*sizeof(TEXTURE *), "blob texture list"));
	for (i = 0; i < New->Data->Number_Of_Components; i++)
		New->Element_Texture[i] = Copy_Textures(Element_Texture[i]);
	return (New);
}



/*****************************************************************************
*
* FUNCTION
*
*   Create_Blob_List_Element
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
*   Blob_List_Struct * - Pointer to blob element
*
* AUTHOR
*
*   Dieter Bayer
*
* DESCRIPTION
*
*   Create a new blob element in the component list used during parsing.
*
* CHANGES
*
*   Sep 1994 : Creation.
*
******************************************************************************/

Blob_List_Struct *Blob::Create_Blob_List_Element()
{
	return (new Blob_List_Struct);
}



/*****************************************************************************
*
* FUNCTION
*
*   Destroy_Blob
*
* INPUT
*
*   Object - Pointer to blob structure
*
* OUTPUT
*
*   Object
*
* RETURNS
*
* AUTHOR
*
*   Alexander Enzmann
*
* DESCRIPTION
*
*   Destroy a blob.
*
*   NOTE: The blob data is destroyed if they are no longer used by any copy.
*
* CHANGES
*
*   Jul 1994 : Added code for blob data reference counting. [DB]
*
*   Dec 1994 : Fixed memory leakage. [DB]
*
*   Aug 1995 : Fixed freeing of already freed memory. [DB]
*
******************************************************************************/

Blob::~Blob()
{
	Destroy_Transform(Trans);
	for (int i = 0; i < Data->Number_Of_Components; i++)
		Destroy_Textures(Element_Texture[i]);
	POV_FREE(Element_Texture);
	if (Data != NULL)
		Data->ReleaseReference () ;
}

/*****************************************************************************
*
* FUNCTION
*
*   Compute_Blob_BBox
*
* INPUT
*
*   Blob - Blob
*
* OUTPUT
*
*   Blob
*
* RETURNS
*
* AUTHOR
*
*   Dieter Bayer
*
* DESCRIPTION
*
*   Calculate the bounding box of a blob.
*
* CHANGES
*
*   Aug 1994 : Creation.
*
******************************************************************************/

void Blob::Compute_BBox()
{
	int i;
	DBL radius, radius2;
	VECTOR Center, Min, Max;

	Make_Vector(Min, BOUND_HUGE, BOUND_HUGE, BOUND_HUGE);
	Make_Vector(Max, - BOUND_HUGE, - BOUND_HUGE, - BOUND_HUGE);

	for (i = 0; i < Data->Number_Of_Components; i++)
	{
		if (Data->Entry[i].c[2] > 0.0)
		{
			get_element_bounding_sphere(&Data->Entry[i], Center, &radius2);

			radius = sqrt(radius2);

			Min[X] = min(Min[X], Center[X] - radius);
			Min[Y] = min(Min[Y], Center[Y] - radius);
			Min[Z] = min(Min[Z], Center[Z] - radius);
			Max[X] = max(Max[X], Center[X] + radius);
			Max[Y] = max(Max[Y], Center[Y] + radius);
			Max[Z] = max(Max[Z], Center[Z] + radius);
		}
	}

	Make_BBox_from_min_max(BBox, Min, Max);

	if (Trans != NULL)
	{
		Recompute_BBox(&BBox, Trans);
	}
}



/*****************************************************************************
*
* FUNCTION
*
*   get_element_bounding_sphere
*
* INPUT
*
*   Element - Pointer to element
*   Center  - Bounding sphere's center
*   Radius2 - Bounding sphere's squared radius
*
* OUTPUT
*
*   Center, Radius2
*
* RETURNS
*
* AUTHOR
*
*   Dieter Bayer
*
* DESCRIPTION
*
*   Calculate the bounding sphere of a blob element.
*
* CHANGES
*
*   Sep 1994 : Creation.
*
******************************************************************************/

void Blob::get_element_bounding_sphere(const Blob_Element *Element, VECTOR Center, DBL *Radius2)
{
	DBL r, r2 = 0.0;
	VECTOR C;
	BBOX local_BBox;

	switch (Element->Type)
	{
		case BLOB_SPHERE:
		case BLOB_ELLIPSOID:

			r2 = Element->rad2;

			Assign_Vector(C, Element->O);

			break;

		case BLOB_BASE_HEMISPHERE:

			r2 = Element->rad2;

			Make_Vector(C, 0.0, 0.0, 0.0);

			break;

		case BLOB_APEX_HEMISPHERE:

			r2 = Element->rad2;

			Make_Vector(C, 0.0, 0.0, Element->len);

			break;

		case BLOB_CYLINDER :

			Make_Vector(C, 0.0, 0.0, 0.5 * Element->len);

			r2 = Element->rad2 + Sqr(0.5 * Element->len);

			break;
	}

	/* Transform bounding sphere if necessary. */
	if (Element->Trans != NULL)
	{
		r = sqrt(r2);

		MTransPoint(C, C, Element->Trans);

		Make_BBox(local_BBox, 0, 0, 0, r, r, r);
		Recompute_BBox(&local_BBox, Element->Trans);
		r = max(max(fabs(local_BBox.Lengths[X]), fabs(local_BBox.Lengths[Y])),
		           fabs(local_BBox.Lengths[Z]));

		r2 = Sqr(r) + EPSILON;
	}

	Assign_Vector(Center, C);

	*Radius2 = r2;
}



/*****************************************************************************
*
* FUNCTION
*
*
* INPUT
*
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Dieter Bayer + others
*
* DESCRIPTION
*
*
* CHANGES
*
*   Sep 1994 : Creation.
*
******************************************************************************/

Blob_Element::Blob_Element (void)
{
	Type = 0;
	index = 0;
	len = 0.0;
	rad2 = 0.0;
	c[0] = 0.0;
	c[1] = 0.0;
	c[2] = 0.0;
	Make_Vector(O, 0.0, 0.0, 0.0);
	Texture = NULL;
	Trans = NULL;
}

Blob_Element::~Blob_Element ()
{
}

Blob_Data::Blob_Data (int Count)
{
	References = 1;
	Tree = NULL;
	Number_Of_Components = Count;
	Entry.resize (Count) ;
}

Blob_Data *Blob_Data::AcquireReference (void)
{
	References++;
	return (this) ;
}

void Blob_Data::ReleaseReference (void)
{
	if (--References == 0)
	{
		Destroy_Bounding_Sphere_Hierarchy(Tree);

		/*
		 * Make sure to destroy multiple references of a texture
		 * and/or transformation only once. Multiple references
		 * are only used with cylindrical blobs. Thus it's
		 * enough to ignore all cylinder caps.
		 */
		for (int i = 0; i < Number_Of_Components; i++)
			if ((Entry[i].Type == BLOB_SPHERE) || (Entry[i].Type == BLOB_ELLIPSOID) || (Entry[i].Type == BLOB_CYLINDER))
				Destroy_Transform(Entry[i].Trans);

		delete this ;
	}
}

Blob_Data::~Blob_Data ()
{
	// NOTE: ReleaseReference will call "delete this" if References == 1, so
	//       we need to ensure we don't recurse when that happens. Only call
	//       ReleaseReference() if References > 0.
	assert (References <= 1) ;
	if (References > 0)
		ReleaseReference () ;
}

/*****************************************************************************
*
* FUNCTION
*
*   Make_Blob
*
* INPUT
*
*   Blob       - Pointer to blob structure
*   threshold  - Blob's threshold
*   BlobList   - Pointer to elements
*   npoints    - Number of elements
*
* OUTPUT
*
*   Blob
*
* RETURNS
*
* AUTHOR
*
*   Alexander Enzmann
*
* DESCRIPTION
*
*   Create a blob after it was read from the scene file.
*
*   Starting with the density function: (1-r^2)^2, we have a field
*   that varies in strength from 1 at r = 0 to 0 at r = 1. By
*   substituting r/rad for r, we can adjust the range of influence
*   of a particular component. By multiplication by coeff, we can
*   adjust the amount of total contribution, giving the formula:
*
*     coeff * (1 - (r/rad)^2)^2
*
*   This varies in strength from coeff at r = 0, to 0 at r = rad.
*
* CHANGES
*
*   Jul 1994 : Added code for cylindrical and ellipsoidical blobs. [DB]
*
******************************************************************************/

int Blob::Make_Blob(DBL threshold, Blob_List_Struct *BlobList, int npoints, TraceThreadData *Thread)
{
	int i, count;
	DBL rad2, coeff;
	Blob_List_Struct *temp;

	if (npoints < 1)
		throw POV_EXCEPTION_STRING("Need at least one component in a blob.");

	/* Figure out how many components there will be. */

	for (i = 0, count = npoints, temp = BlobList; i < npoints; i++, temp = temp->next)
		if (temp->elem.Type & BLOB_CYLINDER)
			count += 2;

	/* Test for too many components. [DB 12/94] */
	if (count >= MAX_BLOB_COMPONENTS)
		throw POV_EXCEPTION_STRING("There are more than the maximum supported components in a blob.");

	/* Initialize the blob data. */

	Data->Threshold = threshold;
	Data->Number_Of_Components = count;
	Data->Entry.resize (count) ;

	vector<Blob_Element>::iterator it = Data->Entry.begin() ;
	for (i = 0; i < npoints; i++)
	{
		temp = BlobList;
		assert (it != Data->Entry.end()) ;
		Blob_Element* Entry = &*it++ ;

		if ((fabs(temp->elem.c[2]) < EPSILON) || (temp->elem.rad2 < EPSILON))
		{
;// TODO MESSAGE			Warning(0, "Degenerate Blob element");
		}

		/* Initialize component. */
		*Entry = temp->elem;

		/* We have a multi-texture blob. */
		if (Entry->Texture != NULL)
			Set_Flag(this, MULTITEXTURE_FLAG);

		/* Store blob specific information. */
		rad2 = temp->elem.rad2;
		coeff = temp->elem.c[2];
		Entry->c[0] = coeff / (rad2 * rad2);
		Entry->c[1] = -(2.0 * coeff) / rad2;
		Entry->c[2] = coeff;

		if (temp->elem.Type == BLOB_CYLINDER)
		{
			/* Create hemispherical component at the base. */
			Entry = &*it++ ;
			*Entry = temp->elem;
			Entry->Type = BLOB_BASE_HEMISPHERE;
			Entry->c[0] = coeff / (rad2 * rad2);
			Entry->c[1] = -(2.0 * coeff) / rad2;
			Entry->c[2] = coeff;

			/* Create hemispherical component at the apex. */
			Entry = &*it++ ;
			*Entry = temp->elem;
			Entry->Type = BLOB_APEX_HEMISPHERE;
			Entry->c[0] = coeff / (rad2 * rad2);
			Entry->c[1] = -(2.0 * coeff) / rad2;
			Entry->c[2] = coeff;
		}

		/* Get rid of texture non longer needed. */

		BlobList = BlobList->next;
		Destroy_Textures(temp->elem.Texture);
		delete temp ;
	}

	for (i = 0; i < count; i++)
		Data->Entry[i].index = i;

	/* Compute bounding box. */

	Compute_BBox();

	/* Create bounding sphere hierarchy. */

	if (Test_Flag(this, HIERARCHY_FLAG))
		build_bounding_hierarchy();

	if (count * 5 >= Thread->Blob_Coefficient_Count)
	{
		POV_FREE(Thread->Blob_Coefficients);
		Thread->Blob_Coefficient_Count = count * 7;
		Thread->Blob_Coefficients = reinterpret_cast<DBL *>(POV_MALLOC(sizeof(DBL) * Thread->Blob_Coefficient_Count, "Blob Coefficients"));
	}

	if (Data->Number_Of_Components * 2 >= Thread->Blob_Interval_Count)
	{
		delete[] Thread->Blob_Intervals;
		Thread->Blob_Interval_Count = Data->Number_Of_Components * 5 / 2;
		Thread->Blob_Intervals = new Blob_Interval_Struct [Thread->Blob_Interval_Count];
	}

	return (count) ;
}

/*****************************************************************************
*
* FUNCTION
*
*   Test_Blob_Opacity
*
* INPUT
*
*   Blob - Pointer to blob structure
*
* OUTPUT
*
*   Blob
*
* RETURNS
*
* AUTHOR
*
*   Dieter Bayer
*
* DESCRIPTION
*
*   Set the opacity flag of the blob according to the opacity
*   of the blob's texture(s).
*
* CHANGES
*
*   Apr 1996 : Creation.
*
******************************************************************************/

void Blob::Test_Blob_Opacity()
{
	int i;

	/* Initialize opacity flag to the opacity of the object's texture. */

	if ((Texture == NULL) || (Test_Opacity(Texture)))
	{
		Set_Flag(this, OPAQUE_FLAG);
	}

	if (Test_Flag(this, MULTITEXTURE_FLAG))
	{
		for (i = 0; i < Data->Number_Of_Components; i++)
		{
			if (Element_Texture[i] != NULL)
			{
				/* If component's texture isn't opaque the blob is neither. */

				if (!Test_Opacity(Element_Texture[i]))
				{
					Clear_Flag(this, OPAQUE_FLAG);
				}
			}
		}
	}
}



/*****************************************************************************
*
* FUNCTION
*
*   build_bounding_hierarchy
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Dieter Bayer
*
* DESCRIPTION
*
*   Create the bounding sphere hierarchy.
*
* CHANGES
*
*   Oct 1994 : Creation. (Derived from the bounding slab creation code)
*
******************************************************************************/

void Blob::build_bounding_hierarchy()
{
	int i, nElem, maxelements;
	BSPHERE_TREE **Elements;

	nElem = (int)Data->Number_Of_Components;

	maxelements = 2 * nElem;

	/*
	 * Now allocate an array to hold references to these elements.
	 */

	Elements = reinterpret_cast<BSPHERE_TREE **>(POV_MALLOC(maxelements*sizeof(BSPHERE_TREE *), "blob bounding hierarchy"));

	/* Init list with blob elements. */

	for (i = 0; i < nElem; i++)
	{
		Elements[i] = reinterpret_cast<BSPHERE_TREE *>(POV_MALLOC(sizeof(BSPHERE_TREE), "blob bounding hierarchy"));

		Elements[i]->Entries = 0;
		Elements[i]->Node    = reinterpret_cast<BSPHERE_TREE **>(&Data->Entry[i]);

		get_element_bounding_sphere(&Data->Entry[i], Elements[i]->C, &Elements[i]->r2);
	}

	Build_Bounding_Sphere_Hierarchy(&Data->Tree, nElem, &Elements);

	/* Get rid of the Elements array. */

	POV_FREE(Elements);
}



/*****************************************************************************
*
* FUNCTION
*
*   Determine_Blob_Textures
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Dieter Bayer
*
* DESCRIPTION
*
*   Determine the textures and weights of all components affecting
*   the given intersection point. The weights are calculated from
*   the field values and sum to 1.
*
* CHANGES
*
*   Sep 1994 : Creation.
*
*   Mar 1996 : Make the call to resize the textures/weights list just once
*              at the beginning instead of doing it for every element. [DB]
*
******************************************************************************/

void Blob::Determine_Textures(Intersection *isect, bool hitinside, WeightedTextureVector& textures, TraceThreadData *Thread)
{
	int i;
	unsigned int size;
	DBL rad2;
	VECTOR V1, P;
	Blob_Element *Element;
	BSPHERE_TREE *Tree;
	size_t firstinserted = textures.size();
	BSPHERE_TREE **Queue = reinterpret_cast<BSPHERE_TREE **>(Thread->Blob_Queue);

	/* Transform the point into the blob space. */
	getLocalIPoint(P, isect);

	if (Data->Tree == NULL)
	{
		/* There's no tree --> step through all elements. */

		for (i = 0; i < Data->Number_Of_Components; i++)
		{
			Element = &Data->Entry[i];
			determine_element_texture(Element, Element_Texture[i], P, textures);
		}
	}
	else
	{
		/* A tree exists --> step through the tree. */

		size = 0;

		Queue[size++] = Data->Tree;

		while (size > 0)
		{
			Tree = Queue[--size];

			/* Test if current node is a leaf. */

			if (Tree->Entries <= 0)
			{
				determine_element_texture(reinterpret_cast<Blob_Element *>(Tree->Node), Element_Texture[((Blob_Element *)Tree->Node)->index], P, textures);
			}
			else
			{
				/* Test all sub-nodes. */

				for (i = 0; i < (int)Tree->Entries; i++)
				{
					/* Insert sub-node if we are inside. */

					VSub(V1, P, Tree->Node[i]->C);

					VDot(rad2, V1, V1);

					if (rad2 <= Tree->Node[i]->r2)
						if (insert_node(Tree->Node[i], &size, Thread))
							Queue = reinterpret_cast<BSPHERE_TREE **>(Thread->Blob_Queue);
				}
			}
		}
	}

	/* Normalize weights so that their sum is 1. */

	if((textures.size() - firstinserted) > 0)
	{
		COLC sum = 0.0;

		for(size_t i = firstinserted; i < textures.size(); i++)
			sum += textures[i].weight;

		sum = 1.0 / sum;

		for(size_t i = firstinserted; i < textures.size(); i++)
			textures[i].weight *= sum;
	}
}


/*****************************************************************************
*
* FUNCTION
*
*   determine_element_texture
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Dieter Bayer
*
* DESCRIPTION
*
*   If the intersection point is inside the component calculate
*   the field density and store the element's texture and the field
*   value in the texture/weight list.
*
* CHANGES
*
*   Sep 1994 : Creation.
*
******************************************************************************/

void Blob::determine_element_texture(const Blob_Element *Element, TEXTURE *ElementTex, const VECTOR P, WeightedTextureVector& textures)
{
		DBL density = fabs(calculate_element_field(Element, P));

		if(density > 0.0)
				textures.push_back(WeightedTexture(density, ElementTex != NULL ? ElementTex : Texture));
}



/*****************************************************************************
*
* FUNCTION
*
*   Translate_Blob_Element
*
* INPUT
*
*   Element - Pointer to blob element
*   Vector  - Translation vector
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
*   Translate a blob element.
*
* CHANGES
*
*   Sep 1994 : Creation.
*
******************************************************************************/

void Blob::Translate_Blob_Element(Blob_Element *Element, const VECTOR Vector)
{
	TRANSFORM Trans;

	Compute_Translation_Transform(&Trans, Vector);

	if (Element->Trans == NULL)
	{
		/* This is a sphere component. */

		VAddEq(Element->O, Vector);
		Transform_Textures(Element->Texture, &Trans);
	}
	else
	{
		/* This is one of the other components. */

		Transform_Blob_Element(Element, &Trans);
	}
}



/*****************************************************************************
*
* FUNCTION
*
*   Rotate_Blob_Element
*
* INPUT
*
*   Element - Pointer to blob element
*   Vector  - Translation vector
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
*   Rotate a blob element.
*
* CHANGES
*
*   Sep 1994 : Creation.
*
******************************************************************************/

void Blob::Rotate_Blob_Element(Blob_Element *Element, const VECTOR Vector)
{
	TRANSFORM Trans;

	Compute_Rotation_Transform(&Trans, Vector);

	if (Element->Trans == NULL)
	{
		/* This is a sphere component. */

		MTransPoint(Element->O, Element->O, &Trans);
		Transform_Textures(Element->Texture, &Trans);
	}
	else
	{
		/* This is one of the other components. */

		Transform_Blob_Element(Element, &Trans);
	}
}



/*****************************************************************************
*
* FUNCTION
*
*   Scale_Blob_Element
*
* INPUT
*
*   Element - Pointer to blob element
*   Vector  - Translation vector
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
*   Scale a blob element.
*
* CHANGES
*
*   Sep 1994 : Creation.
*
******************************************************************************/

void Blob::Scale_Blob_Element(Blob_Element *Element, const VECTOR Vector)
{
	TRANSFORM Trans;

	if ((Vector[X] != Vector[Y]) || (Vector[X] != Vector[Z]))
	{
		if (Element->Trans == NULL)
		{
			/* This is a sphere component --> change to ellipsoid component. */

			Element->Type = BLOB_ELLIPSOID;

			Element->Trans = Create_Transform();
		}
	}

	Compute_Scaling_Transform(&Trans, Vector);

	if (Element->Trans == NULL)
	{
		/* This is a sphere component. */

		VScaleEq(Element->O, Vector[X]);

		Element->rad2 *= Sqr(Vector[X]);

		Transform_Textures(Element->Texture, &Trans);
	}
	else
	{
		/* This is one of the other components. */

		Transform_Blob_Element(Element, &Trans);
	}
}



/*****************************************************************************
*
* FUNCTION
*
*   Transform_Blob_Element
*
* INPUT
*
*   Element - Pointer to blob element
*   Trans   - Transformation
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
*   Transform a blob element.
*
* CHANGES
*
*   Sep 1994 : Creation.
*
******************************************************************************/

void Blob::Transform_Blob_Element(Blob_Element *Element, const TRANSFORM *Trans)
{
	if (Element->Trans == NULL)
	{
		/* This is a sphere component --> change to ellipsoid component. */

		Element->Type = BLOB_ELLIPSOID;

		Element->Trans = Create_Transform();
	}

	Compose_Transforms(Element->Trans, Trans);

	Transform_Textures(Element->Texture, Trans);
}



/*****************************************************************************
*
* FUNCTION
*
*   Invert_Blob_Element
*
* INPUT
*
*   Element - Pointer to blob element
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
*   Invert blob element by negating its strength.
*
* CHANGES
*
*   Sep 1994 : Creation.
*
******************************************************************************/

void Blob::Invert_Blob_Element(Blob_Element *Element)
{
	Element->c[2] *= -1.0;
}


/*****************************************************************************
*
* FUNCTION
*
*   insert_node
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Dieter Bayer
*
* DESCRIPTION
*
*   Insert a node into the node queue.
*
* CHANGES
*
*   Feb 1995 : Creation.
*
******************************************************************************/

bool Blob::insert_node(BSPHERE_TREE *Node, unsigned int *size, TraceThreadData *Thread)
{
	/* Resize queue if necessary. */
	bool rval = false ;
	BSPHERE_TREE **Queue = reinterpret_cast<BSPHERE_TREE **>(Thread->Blob_Queue);

	if (*size >= Thread->Max_Blob_Queue_Size)
	{
		Thread->Max_Blob_Queue_Size = (*size + 1) * 3 / 2;
		Queue = reinterpret_cast<BSPHERE_TREE **>(POV_REALLOC(Queue, Thread->Max_Blob_Queue_Size*sizeof(BSPHERE_TREE *), "blob queue"));
		Thread->Blob_Queue = reinterpret_cast<void **>(Queue);
		rval = true ;
	}

	Queue[(*size)++] = Node;
	return (rval) ;
}



/*****************************************************************************
*
* FUNCTION
*
*   Create_Blob_Element_Texture_List
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Dieter Bayer
*
* DESCRIPTION
*
*   Create a list of all textures in the blob.
*
*   The list actually contains copies of the textures not
*   just references to them.
*
* CHANGES
*
*   Mar 1996 : Created.
*
******************************************************************************/

void Blob::Create_Blob_Element_Texture_List(Blob_List_Struct *BlobList, int npoints)
{
	int i, count;
	Blob_List_Struct *bl;
	TEXTURE **et ;

	if (npoints < 1)
		throw POV_EXCEPTION_STRING("Need at least one component in a blob.");

	/* Figure out how many components there will be. */
	for (i = 0, count = npoints, bl = BlobList ; i < npoints; i++, bl = bl->next)
		if (bl->elem.Type & BLOB_CYLINDER)
			count += 2;

	/* Test for too many components. [DB 12/94] */
	if (count >= MAX_BLOB_COMPONENTS)
		throw POV_EXCEPTION_STRING("There are more than the maximum supported components in a blob.");

	Data = new Blob_Data (count) ;

	/* Allocate memory for list. */
	et = Element_Texture = reinterpret_cast<TEXTURE **>(POV_CALLOC(count,sizeof(TEXTURE *), "blob texture list"));
	for (i = 0, bl = BlobList; i < npoints; i++, bl = bl->next)
	{
		/*
		 * Copy texture into element texture list. This is neccessary
		 * because individual textures have to be transformed too if
		 * copies of the blob are transformed.
		 */
		*et++ = Copy_Textures(bl->elem.Texture);
		if (bl->elem.Type == BLOB_CYLINDER)
		{
			*et++ = Copy_Textures(bl->elem.Texture);
			*et++ = Copy_Textures(bl->elem.Texture);
		}
	}
}

/*****************************************************************************/

void Blob::getLocalIPoint(VECTOR lip, Intersection *isect) const
{
	if(isect->haveLocalIPoint == false)
	{
		if(Trans != NULL)
			MInvTransPoint(isect->LocalIPoint, isect->IPoint, Trans);
		else
			Assign_Vector(isect->LocalIPoint, isect->IPoint);

		isect->haveLocalIPoint = true;
	}

	Assign_Vector(lip, isect->LocalIPoint);
}

}
