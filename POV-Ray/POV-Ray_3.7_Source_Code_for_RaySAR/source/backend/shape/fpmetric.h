/*******************************************************************************
 * fpmetric.h
 *
 * This module contains all defines, typedefs, and prototypes for fpmetric.cpp.
 *
 * This module was written by D.Skarda&T.Bily and modified by R.Suzuki.
 * Ported to POV-Ray 3.5 by Thorsten Froehlich.
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
 * $File: //depot/public/povray/3.x/source/backend/shape/fpmetric.h $
 * $Revision: #1 $
 * $Change: 6069 $
 * $DateTime: 2013/11/06 11:59:40 $
 * $Author: chrisc $
 *******************************************************************************/

#ifndef FPMETRIC_H
#define FPMETRIC_H

#include "backend/parser/parse.h"

namespace pov
{

/*****************************************************************************
* Global preprocessor defines
******************************************************************************/

#define PARAMETRIC_OBJECT        (PATCH_OBJECT)



/*****************************************************************************
* Global typedefs
******************************************************************************/

typedef struct PrecompParValues_Struct PRECOMP_PAR_DATA;

struct PrecompParValues_Struct
{
	int use, depth;
	char flags;
	DBL *Low[3], *Hi[3];     /*  X,Y,Z  */
};

class Parametric : public ObjectBase
{
	public:
		FunctionVM *vm;
		FUNCTION_PTR Function[3];
		DBL umin, umax, vmin, vmax;
		DBL accuracy;
		DBL max_gradient;
		int Inverted;

		int container_shape;
		union
		{
			struct
			{
				VECTOR center;
				DBL radius;
			} sphere;
			struct
			{
				VECTOR corner1;
				VECTOR corner2;
			} box;
		} container;

		Parametric();
		virtual ~Parametric();

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

		void Precompute_Parametric_Values(char flags, int depth, FPUContext *ctx);
	protected:
		void Precomp_Par_Int(int depth, DBL umin, DBL vmin, DBL umax, DBL vmax, FPUContext *ctx);
		PRECOMP_PAR_DATA *Copy_PrecompParVal();
		void Destroy_PrecompParVal();

		static inline DBL Evaluate_Function_UV(FPUContext *ctx, FUNCTION funct, const UV_VECT fnvec);
		static inline void Evaluate_Function_Interval_UV(FPUContext *ctx, FUNCTION funct, DBL threshold, const UV_VECT fnvec_low, const UV_VECT fnvec_hi, DBL max_gradient, DBL& low, DBL& hi);
		static void Interval(DBL dx, DBL a, DBL b, DBL max_gradient, DBL *Min, DBL *Max);
	private:
		PRECOMP_PAR_DATA *PData;
		int PrecompLastDepth;
};

}

#endif
