/*******************************************************************************
 * sor.h
 *
 * This module contains all defines, typedefs, and prototypes for SOR.CPP.
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
 * POV-Ray is based on the popular DKB raytracer version 2.12.
 * DKBTrace was originally written by David K. Buck.
 * DKBTrace Ver 2.0-2.12 were written by David K. Buck & Aaron A. Collins.
 * ---------------------------------------------------------------------------
 * $File: //depot/public/povray/3.x/source/backend/shape/sor.h $
 * $Revision: #1 $
 * $Change: 6069 $
 * $DateTime: 2013/11/06 11:59:40 $
 * $Author: chrisc $
 *******************************************************************************/

#ifndef SOR_H
#define SOR_H

#include "backend/bounding/bcyl.h"

namespace pov
{

/*****************************************************************************
* Global preprocessor definitions
******************************************************************************/

#define SOR_OBJECT (STURM_OK_OBJECT)

/* Generate additional surface of revolution statistics. */

#define SOR_EXTRA_STATS 1



/*****************************************************************************
* Global typedefs
******************************************************************************/

typedef struct Sor_Spline_Entry_Struct SOR_SPLINE_ENTRY;
typedef struct Sor_Spline_Struct SOR_SPLINE;

struct Sor_Spline_Entry_Struct
{
	DBL A, B, C, D;
};

struct Sor_Spline_Struct
{
	int References;
	SOR_SPLINE_ENTRY *Entry;
	BCYL *BCyl;                 /* bounding cylinder.                  */
};

class Sor : public ObjectBase
{
	public:
		int Number;
		SOR_SPLINE *Spline;      /* List of spline segments     */
		DBL Height1, Height2;    /* Min./Max. height            */
		DBL Radius1, Radius2;    /* Min./Max. radius            */
		DBL Base_Radius_Squared; /* Radius**2 of the base plane */
		DBL Cap_Radius_Squared;  /* Radius**2 of the cap plane  */

		Sor();
		virtual ~Sor();

		virtual ObjectPtr Copy();

		virtual bool All_Intersections(const Ray&, IStack&, TraceThreadData *);
		virtual bool Inside(const VECTOR, TraceThreadData *) const;
		virtual void Normal(VECTOR, Intersection *, TraceThreadData *) const;
		virtual void UVCoord(UV_VECT, const Intersection *, TraceThreadData *) const;
		virtual void Translate(const VECTOR, const TRANSFORM *);
		virtual void Rotate(const VECTOR, const TRANSFORM *);
		virtual void Scale(const VECTOR, const TRANSFORM *);
		virtual void Transform(const TRANSFORM *);
		virtual void Invert();
		virtual void Compute_BBox();

		void Compute_Sor(UV_VECT *P, TraceThreadData *Thread);
	protected:
		bool Intersect(const Ray& ray, IStack& Depth_Stack, TraceThreadData *Thread);
		bool test_hit(const Ray&, IStack&, DBL, DBL, int, int, TraceThreadData *Thread);
};

}

#endif
