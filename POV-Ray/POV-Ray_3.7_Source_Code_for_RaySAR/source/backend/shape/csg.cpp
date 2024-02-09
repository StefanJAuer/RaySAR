/*******************************************************************************
 * csg.cpp
 *
 * This module implements routines for constructive solid geometry.
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
 * $File: //depot/public/povray/3.x/source/backend/shape/csg.cpp $
 * $Revision: #1 $
 * $Change: 6069 $
 * $DateTime: 2013/11/06 11:59:40 $
 * $Author: chrisc $
 *******************************************************************************/

// frame.h must always be the first POV file included (pulls in platform config)
#include "backend/frame.h"
#include "backend/math/vector.h"
#include "backend/bounding/bbox.h"
#include "backend/shape/csg.h"
#include "backend/math/matrices.h"
#include "backend/scene/objects.h"
#include "backend/shape/quadrics.h"
#include "backend/shape/hfield.h"
#include "backend/scene/threaddata.h"

#include "lightgrp.h" // TODO

#include <algorithm>

// this must be the last file included
#include "base/povdebug.h"

namespace pov
{

/*****************************************************************************
* Global preprocessor defines
******************************************************************************/

#define UNION_OBJECT        (IS_COMPOUND_OBJECT | IS_CSG_OBJECT)
#define MERGE_OBJECT        (IS_COMPOUND_OBJECT | IS_CSG_OBJECT)
#define INTERSECTION_OBJECT (IS_COMPOUND_OBJECT | IS_CSG_OBJECT)



inline bool Test_Ray_Flags(const Ray& ray, const ObjectBase* obj)
{
	// CJC 2005 if ray is primary ray ignore NO_IMAGE_FLAG to support the trace() SDL function
 	// TODO FIXME - I uess it would be better to have the trace() function use a different ray type [CLi]
	return ( ( !ray.IsPhotonRay() &&
	           (!Test_Flag(obj, NO_IMAGE_FLAG) || ray.IsImageRay() == false || ray.IsPrimaryRay() == true) &&
	           (!Test_Flag(obj, NO_REFLECTION_FLAG) || ray.IsReflectionRay() == false) &&
	           (!Test_Flag(obj, NO_RADIOSITY_FLAG) || ray.IsRadiosityRay() == false) ) ||
	         ( ray.IsPhotonRay() && !Test_Flag(obj, NO_SHADOW_FLAG) ) );
}

inline bool Test_Ray_Flags_Shadow(const Ray& ray, const ObjectBase* obj)
{
	// TODO CLARIFY - why does this function not ignore NO_IMAGE_FLAG for primary rays, as Test_Ray_Flags() does? [CLi]
	return ( ( !ray.IsPhotonRay() &&
	           (!Test_Flag(obj, NO_IMAGE_FLAG) || ray.IsImageRay() == false) &&
	           (!Test_Flag(obj, NO_REFLECTION_FLAG) || ray.IsReflectionRay() == false) &&
	           (!Test_Flag(obj, NO_RADIOSITY_FLAG) || ray.IsRadiosityRay() == false) ) ||
	         ( ray.IsPhotonRay() && !Test_Flag(obj, NO_SHADOW_FLAG) ) ||
	         ( ray.IsShadowTestRay() && !Test_Flag(obj, NO_SHADOW_FLAG) ) );
}

/*****************************************************************************
*
* FUNCTION
*
*   All_CSG_Union_Intersections
*
* INPUT
*   
* OUTPUT
*   
* RETURNS
*   
* AUTHOR
*
*   POV-Ray Team
*   
* DESCRIPTION
*
*   -
*
* CHANGES
*
*   Sep 1994 : Added code to count intersection tests. [DB]
*
******************************************************************************/

bool CSGUnion::All_Intersections(const Ray& ray, IStack& Depth_Stack, TraceThreadData *Thread)
{
	int Found;

	Thread->Stats()[Ray_CSG_Union_Tests]++;

	Found = false;

	// Use shortcut if no clip.

	if(Clip.empty())
	{
		for(vector<ObjectPtr>::const_iterator Current_Sib = children.begin(); Current_Sib != children.end(); Current_Sib++)
		{
			if(Test_Ray_Flags(ray, (*Current_Sib))) // TODO CLARIFY - why does CSGUnion use Test_Ray_Flags(), while CSGMerge uses Test_Ray_Flags_Shadow(), and CSGIntersection uses neither?
			{
				if((*Current_Sib)->Bound.empty() == true || Ray_In_Bound(ray, (*Current_Sib)->Bound, Thread))
				{
					if((*Current_Sib)->All_Intersections(ray, Depth_Stack, Thread))
						Found = true;
				}
			}
		}
	}
	else
	{
		IStack Local_Stack(Thread->stackPool);
		assert(Local_Stack->empty()); // verify that the IStack pulled from the pool is in a cleaned-up condition

		for(vector<ObjectPtr>::const_iterator Current_Sib = children.begin(); Current_Sib != children.end(); Current_Sib++)
		{
			if(Test_Ray_Flags(ray, (*Current_Sib))) // TODO CLARIFY - why does CSGUnion use Test_Ray_Flags(), while CSGMerge uses Test_Ray_Flags_Shadow(), and CSGIntersection uses neither?
			{
				if((*Current_Sib)->Bound.empty() == true || Ray_In_Bound(ray, (*Current_Sib)->Bound, Thread))
				{
					if((*Current_Sib)->All_Intersections (ray, Local_Stack, Thread))
					{
						while(Local_Stack->size() > 0)
						{
							if(Clip.empty() || Point_In_Clip(Local_Stack->top().IPoint, Clip, Thread))
							{
								Local_Stack->top().Csg = this;

								Depth_Stack->push(Local_Stack->top());

								Found = true;
							}

							Local_Stack->pop();
						}
					}
				}
			}
		}
		assert(Local_Stack->empty()); // verify that the IStack is in a cleaned-up condition (again)
	}

	if(Found)
		Thread->Stats()[Ray_CSG_Union_Tests_Succeeded]++;

	return (Found);
}



/*****************************************************************************
*
* FUNCTION
*
*   All_CSG_Intersection_Intersections
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   POV-Ray Team
*
* DESCRIPTION
*
*   -
*
* CHANGES
*
*   Sep 1994 : Added code to count intersection tests. [DB]
*
******************************************************************************/

bool CSGIntersection::All_Intersections(const Ray& ray, IStack& Depth_Stack, TraceThreadData *Thread)
{
	int Maybe_Found, Found;
	IStack Local_Stack(Thread->stackPool);
	assert(Local_Stack->empty()); // verify that the IStack pulled from the pool is in a cleaned-up condition

	Thread->Stats()[Ray_CSG_Intersection_Tests]++;

	Found = false;

	for(vector<ObjectPtr>::const_iterator Current_Sib = children.begin(); Current_Sib != children.end(); Current_Sib++)
	{
		if ((*Current_Sib)->Bound.empty() == true || Ray_In_Bound(ray, (*Current_Sib)->Bound, Thread))
		{
			if((*Current_Sib)->All_Intersections(ray, Local_Stack, Thread))
			{
				while(Local_Stack->size() > 0)
				{
					Maybe_Found = true;

					for(vector<ObjectPtr>::const_iterator Inside_Sib = children.begin(); Inside_Sib != children.end(); Inside_Sib++)
					{
						if(*Inside_Sib != *Current_Sib)
						{
							if(!((*Inside_Sib)->Type & LIGHT_SOURCE_OBJECT) || (!((LightSource *)(*Inside_Sib))->children.empty()))
							{
								if(!Inside_Object(Local_Stack->top().IPoint, *Inside_Sib, Thread))
								{
									Maybe_Found = false;
									break;
								}
							}
						}
					}

					if(Maybe_Found)
					{
						if(Clip.empty() || Point_In_Clip(Local_Stack->top().IPoint, Clip, Thread))
						{
							Local_Stack->top().Csg = this;

							Depth_Stack->push(Local_Stack->top());

							Found = true;
						}
					}

					Local_Stack->pop();
				}
			}
		}
	}

	if(Found)
		Thread->Stats()[Ray_CSG_Intersection_Tests_Succeeded]++;

	assert(Local_Stack->empty()); // verify that the IStack is in a cleaned-up condition (again)
	return (Found);
}



/*****************************************************************************
*
* FUNCTION
*
*   All_CSG_Merge_Intersections
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   POV-Ray Team
*
* DESCRIPTION
*
*   -
*
* CHANGES
*
*   Sep 1994 : Added code to count intersection tests. [DB]
*
******************************************************************************/

bool CSGMerge::All_Intersections(const Ray& ray, IStack& Depth_Stack, TraceThreadData *Thread)
{
	int Found;
	bool inside_flag;
	IStack Local_Stack(Thread->stackPool);
	assert(Local_Stack->empty()); // verify that the IStack pulled from the pool is in a cleaned-up condition

	Thread->Stats()[Ray_CSG_Merge_Tests]++;

	Found = false;

	// FIXME - though the name is misleading, the OPTIMISE_SHADOW_TEST flag can be used to
	//  determine if we're in a shadow ray, but it SHOULD be renamed.
	// We should probably change Optimization_Flags to a "ray-type" variable, that will tell
	// us if it is primary, reflection, refraction, shadow, primary photon, photon refleciton, or photon refraction ray.
	int shadow_flag = ray.IsShadowTestRay(); // TODO FIXME - why is this flag not used?!

	for(vector<ObjectPtr>::const_iterator Sib1 = children.begin(); Sib1 != children.end(); Sib1++)
	{
		if ( Test_Ray_Flags_Shadow(ray, (*Sib1)) )// TODO CLARIFY - why does CSGUnion use Test_Ray_Flags(), while CSGMerge uses Test_Ray_Flags_Shadow(), and CSGIntersection uses neither?
		{
			if ((*Sib1)->Bound.empty() == true || Ray_In_Bound (ray, (*Sib1)->Bound, Thread))
			{
				if ((*Sib1)->All_Intersections (ray, Local_Stack, Thread))
				{
					while (Local_Stack->size() > 0)
					{
						if (Clip.empty() || Point_In_Clip (Local_Stack->top().IPoint, Clip, Thread))
						{
							inside_flag = true;

							for(vector<ObjectPtr>::const_iterator Sib2 = children.begin(); (Sib2 != children.end()) && (inside_flag == true); Sib2++)
							{
								if (*Sib1 != *Sib2)
								{
									if (!((*Sib2)->Type & LIGHT_SOURCE_OBJECT) || (!((LightSource *)(*Sib2))->children.empty()))
									{
										if ( Test_Ray_Flags_Shadow(ray, (*Sib2)) )// TODO CLARIFY - why does CSGUnion use Test_Ray_Flags(), while CSGMerge uses Test_Ray_Flags_Shadow(), and CSGIntersection uses neither?
										{
											if (Inside_Object(Local_Stack->top().IPoint, *Sib2, Thread))
												inside_flag = false;
										}
									}
								}
							}

							if (inside_flag == true)
							{
								Local_Stack->top().Csg = this;

								Found = true;

								Depth_Stack->push(Local_Stack->top());
							}
						}

						Local_Stack->pop();
					}
				}
			}
		}
	}

	if (Found)
		Thread->Stats()[Ray_CSG_Merge_Tests_Succeeded]++;

	assert(Local_Stack->empty()); // verify that the IStack is in a cleaned-up condition (again)
	return (Found);
}



/*****************************************************************************
*
* FUNCTION
*
*   Inside_CSG_Union
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   POV-Ray Team
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

bool CSGUnion::Inside(const VECTOR IPoint, TraceThreadData *Thread) const
{
	for(vector<ObjectPtr>::const_iterator Current_Sib = children.begin(); Current_Sib != children.end(); Current_Sib++)
	{
		if(!((*Current_Sib)->Type & LIGHT_SOURCE_OBJECT) || (!((LightSource *)(*Current_Sib))->children.empty()))
		{
			if(Inside_Object(IPoint, *Current_Sib, Thread))
				return (true);
		}
	}

	return (false);
}



/*****************************************************************************
*
* FUNCTION
*
*   Inside_CSG_Intersection
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   POV-Ray Team
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

bool CSGIntersection::Inside(const VECTOR IPoint, TraceThreadData *Thread) const
{
	for(vector<ObjectPtr>::const_iterator Current_Sib = children.begin(); Current_Sib != children.end(); Current_Sib++)
		if(!((*Current_Sib)->Type & LIGHT_SOURCE_OBJECT) || (!((LightSource *)(*Current_Sib))->children.empty()))
			if(!Inside_Object(IPoint, (*Current_Sib), Thread))
				return (false);
	return (true);
}




/*****************************************************************************
*
* FUNCTION
*
*   Translate_CSG
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   POV-Ray Team
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


void CSG::Translate(const VECTOR Vector, const TRANSFORM *tr)
{
	for(vector<ObjectPtr>::iterator Current_Sib = children.begin(); Current_Sib != children.end(); Current_Sib++)
		Translate_Object (*Current_Sib, Vector, tr) ;

	Recompute_BBox(&BBox, tr);
}



/*****************************************************************************
*
* FUNCTION
*
*   Rotate_CSG
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   POV-Ray Team
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


void CSG::Rotate(const VECTOR Vector, const TRANSFORM *tr)
{
	for(vector<ObjectPtr>::iterator Current_Sib = children.begin(); Current_Sib != children.end(); Current_Sib++)
		Rotate_Object (*Current_Sib, Vector, tr) ;

	Recompute_BBox(&BBox, tr);
}



/*****************************************************************************
*
* FUNCTION
*
*   Scale_CSG
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   POV-Ray Team
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


void CSG::Scale(const VECTOR Vector, const TRANSFORM *tr)
{
	for(vector<ObjectPtr>::iterator Current_Sib = children.begin(); Current_Sib != children.end(); Current_Sib++)
		Scale_Object (*Current_Sib, Vector, tr) ;

	Recompute_BBox(&BBox, tr);
}



/*****************************************************************************
*
* FUNCTION
*
*   Transform_CSG
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   POV-Ray Team
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


void CSG::Transform(const TRANSFORM *tr)
{
	for(vector<ObjectPtr>::iterator Current_Sib = children.begin(); Current_Sib != children.end(); Current_Sib++)
		Transform_Object(*Current_Sib, tr);

	Recompute_BBox(&BBox, tr);
}

/*****************************************************************************
*
* FUNCTION
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   POV-Ray Team
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

void CSG::Invert()
{
	// REMINDER: Invert_Object will de-allocate the original object pointer and set it to NULL for any CSG children.
	for(vector<ObjectPtr>::iterator Current_Sib = children.begin(); Current_Sib != children.end(); Current_Sib++)
	{
		if (((*Current_Sib)->Type & IS_CSG_OBJECT) != 0)
			*Current_Sib = Invert_CSG_Object(*Current_Sib);
		else
			Invert_Object(*Current_Sib);
	}
	Invert_Flag(this, INVERTED_FLAG);
}

/*****************************************************************************
*
* FUNCTION
*
*   Create_CSG_Union
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   POV-Ray Team
*
* DESCRIPTION
*
*   -
*
* CHANGES
*
*   2000 : NK phmap
*
******************************************************************************/

CSGUnion::CSGUnion() : CSG(UNION_OBJECT)
{
	do_split = true;
}

CSGUnion::CSGUnion(int t) : CSG(t)
{
	do_split = true;
}



/*****************************************************************************
*
* FUNCTION
*
*   Create_CSG_Merge
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   POV-Ray Team
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

CSGMerge::CSGMerge() : CSGUnion(MERGE_OBJECT)
{
}

CSGMerge::CSGMerge(CompoundObject& o, bool transplant) : CSGUnion(MERGE_OBJECT, o, transplant)
{
}

/*****************************************************************************
*
* FUNCTION
*
*   Create_CSG_Intersection
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   POV-Ray Team
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

CSGIntersection::CSGIntersection(bool diff) : CSG(INTERSECTION_OBJECT), isDifference(diff)
{
	do_split = false; // TODO - not necessary but makes debugging clearer
}

CSGIntersection::CSGIntersection(bool diff, CompoundObject& o, bool transplant) : CSG(INTERSECTION_OBJECT, o, transplant), isDifference(diff)
{
	do_split = false; // TODO - not necessary but makes debugging clearer
}



/*****************************************************************************
*
* FUNCTION
*
*   Copy_CSG
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   POV-Ray Team
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


ObjectPtr CSGUnion::Copy()
{
	CSGUnion *New = new CSGUnion();
	Destroy_Transform(New->Trans);
	*New = *this;

	New->children.clear();
	New->children.reserve(children.size());
	for(vector<ObjectPtr>::iterator i(children.begin()); i != children.end(); i++)
		New->children.push_back(Copy_Object(*i));

	if(Type & LIGHT_GROUP_OBJECT)
	{
		New->LLights.clear();
		Promote_Local_Lights(New);
	}

	return (New);
}

ObjectPtr CSGMerge::Copy()
{
	CSGMerge *New = new CSGMerge();
	Destroy_Transform(New->Trans);
	*New = *this;

	New->children.clear();
	New->children.reserve(children.size());
	for(vector<ObjectPtr>::iterator i(children.begin()); i != children.end(); i++)
		New->children.push_back(Copy_Object(*i));

	if(Type & LIGHT_GROUP_OBJECT)
	{
		New->LLights.clear();
		Promote_Local_Lights(New);
	}

	return (New);
}

ObjectPtr CSGIntersection::Copy()
{
	CSGIntersection *New = new CSGIntersection(false);
	Destroy_Transform(New->Trans);
	*New = *this;

	New->children.clear();
	New->children.reserve(children.size());
	for(vector<ObjectPtr>::iterator i(children.begin()); i != children.end(); i++)
		New->children.push_back(Copy_Object(*i));

	if(Type & LIGHT_GROUP_OBJECT)
	{
		New->LLights.clear();
		Promote_Local_Lights(New);
	}

	return (New);
}

CSG *CSGMerge::Morph(void)
{
	CSGIntersection *New = new CSGIntersection(false, *this, true);
	delete this ;
	return (New);
}

CSG *CSGUnion::Morph(void)
{
	CSGIntersection *New = new CSGIntersection(false, *this, true);
	delete this ;
	return (New);
}

CSG *CSGIntersection::Morph(void)
{
	CSGMerge *New = new CSGMerge(*this, true);
	delete this ;
	return (New);
}

/*****************************************************************************
*
* FUNCTION
*
*   Compute_CSG_BBox
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   POV-Ray Team
*
* DESCRIPTION
*
*   -
*
* CHANGES
*
*   Sep 1994 : Improved bounding of quadrics used in CSG intersections. [DB]
*
******************************************************************************/

void CSG::Compute_BBox()
{
	DBL Old_Volume, New_Volume;
	VECTOR NewMin, NewMax, TmpMin, TmpMax, Min, Max;

	if(dynamic_cast<CSGIntersection *>(this) != NULL) // FIXME
	{
		/*
		 * Calculate the bounding box of a CSG intersection
		 * by intersecting the bounding boxes of all children.
		 */

		Make_Vector(NewMin, -BOUND_HUGE, -BOUND_HUGE, -BOUND_HUGE);
		Make_Vector(NewMax,  BOUND_HUGE,  BOUND_HUGE,  BOUND_HUGE);

		vector<Quadric *> Quadrics;

		/* Process all children. */

		for(vector<ObjectPtr>::iterator Current_Sib = children.begin(); Current_Sib != children.end(); Current_Sib++)
		{
			/* Inverted objects and height fields mustn't be considered */

			if(!Test_Flag((*Current_Sib), INVERTED_FLAG) && (dynamic_cast<HField *>(*Current_Sib) == NULL)) // FIXME
			{
				/* We store quadrics since they'll be processed last, to benefit from confining them to a certain range */
				if(dynamic_cast<Quadric *>(*Current_Sib) == NULL) // FIXME
				{
					if(dynamic_cast<Plane *>(*Current_Sib) != NULL) // FIXME
						Quadric::Compute_Plane_Min_Max((Plane *)(*Current_Sib), TmpMin, TmpMax);
					else
						Make_min_max_from_BBox(TmpMin, TmpMax, (*Current_Sib)->BBox);

					NewMin[X] = max(NewMin[X], TmpMin[X]);
					NewMin[Y] = max(NewMin[Y], TmpMin[Y]);
					NewMin[Z] = max(NewMin[Z], TmpMin[Z]);
					NewMax[X] = min(NewMax[X], TmpMax[X]);
					NewMax[Y] = min(NewMax[Y], TmpMax[Y]);
					NewMax[Z] = min(NewMax[Z], TmpMax[Z]);
				}
				else
					Quadrics.push_back(dynamic_cast<Quadric *>(*Current_Sib));
			}
		}

		/* Process any quadrics. */

		for(vector<Quadric *>::iterator i = Quadrics.begin(); i != Quadrics.end(); i++)
		{
			Quadric *q = *i;

			Assign_Vector(Min, NewMin);
			Assign_Vector(Max, NewMax);

			q->Compute_BBox(Min, Max);

			Make_min_max_from_BBox(TmpMin, TmpMax, q->BBox);

			NewMin[X] = max(NewMin[X], TmpMin[X]);
			NewMin[Y] = max(NewMin[Y], TmpMin[Y]);
			NewMin[Z] = max(NewMin[Z], TmpMin[Z]);
			NewMax[X] = min(NewMax[X], TmpMax[X]);
			NewMax[Y] = min(NewMax[Y], TmpMax[Y]);
			NewMax[Z] = min(NewMax[Z], TmpMax[Z]);
		}
	}
	else
	{
		/* Calculate the bounding box of a CSG merge/union object. */

		Make_Vector(NewMin,  BOUND_HUGE,  BOUND_HUGE,  BOUND_HUGE);
		Make_Vector(NewMax, -BOUND_HUGE, -BOUND_HUGE, -BOUND_HUGE);

		for(vector<ObjectPtr>::iterator Current_Sib = children.begin(); Current_Sib != children.end(); Current_Sib++)
		{
			Make_min_max_from_BBox(TmpMin, TmpMax, (*Current_Sib)->BBox);

			NewMin[X] = min(NewMin[X], TmpMin[X]);
			NewMin[Y] = min(NewMin[Y], TmpMin[Y]);
			NewMin[Z] = min(NewMin[Z], TmpMin[Z]);
			NewMax[X] = max(NewMax[X], TmpMax[X]);
			NewMax[Y] = max(NewMax[Y], TmpMax[Y]);
			NewMax[Z] = max(NewMax[Z], TmpMax[Z]);
		}
	}

	if((NewMin[X] > NewMax[X]) || (NewMin[Y] > NewMax[Y]) || (NewMin[Z] > NewMax[Z]))
		;// TODO MESSAGE    Warning(0, "Degenerate CSG bounding box (not used!).");
	else
	{
		New_Volume = (NewMax[X] - NewMin[X]) * (NewMax[Y] - NewMin[Y]) * (NewMax[Z] - NewMin[Z]);

		BOUNDS_VOLUME(Old_Volume, BBox);

		if(New_Volume < Old_Volume)
		{
			Make_BBox_from_min_max(BBox, NewMin, NewMax);

			/* Beware of bounding boxes too large. */

			if((BBox.Lengths[X] > CRITICAL_LENGTH) ||
			   (BBox.Lengths[Y] > CRITICAL_LENGTH) ||
			   (BBox.Lengths[Z] > CRITICAL_LENGTH))
				Make_BBox(BBox, -BOUND_HUGE/2, -BOUND_HUGE/2, -BOUND_HUGE/2, BOUND_HUGE, BOUND_HUGE, BOUND_HUGE);
		}
	}
}


/*****************************************************************************
*
* FUNCTION
*
*   Determine_CSG_Textures
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   POV-Ray Team
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

void CSG::Determine_Textures(Intersection *isect, bool hitinside, WeightedTextureVector& textures, TraceThreadData *threaddata)
{
	if(!children.empty())
	{
		if(Type & CSG_DIFFERENCE_OBJECT)
		{
			// For CSG Differences, use only the first object in the chain
			// (which is the first object in the POV file.  All other objects
			// are the ones that were "removed" from the first one, so their
			// textures should NOT be used.
			if(children[0]->Inside(isect->IPoint, threaddata))
			{
				if(children[0]->Type & IS_COMPOUND_OBJECT)
					children[0]->Determine_Textures(isect, hitinside, textures, threaddata);
				else if(children[0]->Texture != NULL)
					textures.push_back(WeightedTexture(1.0, children[0]->Texture));
			}
		}
		else
		{
			size_t firstinserted = textures.size();

			for(vector<ObjectPtr>::const_iterator Current_Sib = children.begin(); Current_Sib != children.end(); Current_Sib++)
			{
				if((*Current_Sib)->Inside(isect->IPoint, threaddata))
				{
					if((*Current_Sib)->Type & IS_COMPOUND_OBJECT)
						(*Current_Sib)->Determine_Textures(isect, hitinside, textures, threaddata);
					else if((*Current_Sib)->Texture != NULL)
						textures.push_back(WeightedTexture(1.0, (*Current_Sib)->Texture));
				}
			}

			COLC weight = 1.0f / max(COLC(textures.size() - firstinserted), 1.0f);

			for(size_t i = firstinserted; i < textures.size(); i++)
				textures[i].weight = weight;
		}
	}
}

}
