/*******************************************************************************
 * isosurf.h
 *
 * This module contains all defines, typedefs, and prototypes for isosurf.cpp.
 *
 * This module was written by D.Skarda & T.Bily and modified by R.Suzuki.
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
 * $File: //depot/public/povray/3.x/source/backend/shape/isosurf.h $
 * $Revision: #1 $
 * $Change: 6069 $
 * $DateTime: 2013/11/06 11:59:40 $
 * $Author: chrisc $
 *******************************************************************************/

#ifndef ISOSURF_H
#define ISOSURF_H

#include "backend/parser/parse.h"

namespace pov
{

/*****************************************************************************
* Global preprocessor defines
******************************************************************************/

#define ISOSURFACE_OBJECT      (BASIC_OBJECT)
#define ISOSURFACE_MAXTRACE    10


/*****************************************************************************
* Global typedefs
******************************************************************************/

class IsoSurface;
struct FPUContext;

struct ISO_Pair { DBL t,f; };

struct ISO_Max_Gradient
{
	unsigned int refcnt;
	DBL max_gradient, gradient;
	DBL eval_max, eval_cnt, eval_gradient_sum, eval_var;
};

struct ISO_ThreadData
{
	const IsoSurface *current;
	FPUContext *ctx;
	VECTOR Pglobal;
	VECTOR Dglobal;
	DBL Vlength;
	DBL tl;
	DBL fmax;
	bool cache;
	int Inv3;
};

class IsoSurface : public ObjectBase
{
	public:
		FunctionVM *vm;
		FUNCTION_PTR Function;
		volatile DBL max_gradient; // global in eval
		DBL gradient;
		DBL threshold;
		DBL accuracy;
		DBL eval_param[3];
		int max_trace;
		bool closed;
		bool eval;
		bool isCopy;

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

		IsoSurface();
		virtual ~IsoSurface();

		virtual ObjectPtr Copy();

		virtual bool All_Intersections(const Ray&, IStack&, TraceThreadData *);
		virtual bool Inside(const VECTOR, TraceThreadData *) const;
		virtual void Normal(VECTOR, Intersection *, TraceThreadData *) const;
		virtual void Translate(const VECTOR, const TRANSFORM *);
		virtual void Rotate(const VECTOR, const TRANSFORM *);
		virtual void Scale(const VECTOR, const TRANSFORM *);
		virtual void Transform(const TRANSFORM *);
		virtual void Invert();
		virtual void Compute_BBox();

		virtual void DispatchShutdownMessages(MessageFactory& messageFactory);

	protected:
		bool Function_Find_Root(ISO_ThreadData& itd, const VECTOR, const VECTOR, DBL*, DBL*, DBL& max_gradient, bool in_shadow_test);
		bool Function_Find_Root_R(ISO_ThreadData& itd, const ISO_Pair*, const ISO_Pair*, DBL, DBL, DBL, DBL& max_gradient);

		inline DBL Vector_Function(FPUContext *ctx, const VECTOR VPos) const;
		inline DBL Float_Function(ISO_ThreadData& itd, DBL t) const;
		static inline DBL Evaluate_Function(FPUContext *ctx, FUNCTION funct, const VECTOR fnvec);
	private:
		ISO_Max_Gradient *mginfo; // global, but just a statistic (read: not thread safe but we don't care) [trf]
};

}

#endif

