/*******************************************************************************
 * bezier.cpp
 *
 * This module implements the code for Bezier bicubic patch shapes
 *
 * This file was written by Alexander Enzmann.  He wrote the code for
 * bezier bicubic patches and generously provided us these enhancements.
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
 * $File: //depot/public/povray/3.x/source/backend/shape/bezier.cpp $
 * $Revision: #1 $
 * $Change: 6069 $
 * $DateTime: 2013/11/06 11:59:40 $
 * $Author: chrisc $
 *******************************************************************************/

// frame.h must always be the first POV file included (pulls in platform config)
#include "backend/frame.h"
#include "backend/math/vector.h"
#include "backend/shape/bezier.h"
#include "backend/math/matrices.h"
#include "backend/scene/objects.h"
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

const DBL BEZIER_EPSILON = 1.0e-10;
const DBL BEZIER_TOLERANCE = 1.0e-5;



/*****************************************************************************
*
* FUNCTION
*
*   create_new_bezier_node
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

BEZIER_NODE *BicubicPatch::create_new_bezier_node()
{
	BEZIER_NODE *Node = reinterpret_cast<BEZIER_NODE *>(POV_MALLOC(sizeof(BEZIER_NODE), "bezier node"));

	Node->Data_Ptr = NULL;

	return (Node);
}



/*****************************************************************************
*
* FUNCTION
*
*   create_bezier_vertex_block
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

BEZIER_VERTICES *BicubicPatch::create_bezier_vertex_block()
{
	BEZIER_VERTICES *Vertices;

	Vertices = reinterpret_cast<BEZIER_VERTICES *>(POV_MALLOC(sizeof(BEZIER_VERTICES), "bezier vertices"));

	return (Vertices);
}



/*****************************************************************************
*
* FUNCTION
*
*   create_bezier_child_block
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

BEZIER_CHILDREN *BicubicPatch::create_bezier_child_block()
{
	BEZIER_CHILDREN *Children;

	Children = reinterpret_cast<BEZIER_CHILDREN *>(POV_MALLOC(sizeof(BEZIER_CHILDREN), "bezier children"));

	return (Children);
}



/*****************************************************************************
*
* FUNCTION
*
*   bezier_tree_builder
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

BEZIER_NODE *BicubicPatch::bezier_tree_builder(const VECTOR (*Patch)[4][4], DBL u0, DBL u1, DBL v0, DBL v1, int depth, int& max_depth_reached)
{
	VECTOR Lower_Left[4][4], Lower_Right[4][4];
	VECTOR Upper_Left[4][4], Upper_Right[4][4];
	BEZIER_CHILDREN *Children;
	BEZIER_VERTICES *Vertices;
	BEZIER_NODE *Node = create_new_bezier_node();

	if (depth > max_depth_reached)
	{
		max_depth_reached = depth;
	}

	/* Build the bounding sphere for this subpatch. */

	bezier_bounding_sphere(Patch, Node->Center, &(Node->Radius_Squared));

	/*
	 * If the patch is close to being flat, then just perform
	 * a ray-plane intersection test.
	 */

	if (flat_enough(Patch))
	{
		/* The patch is now flat enough to simply store the corners. */

		Node->Node_Type = BEZIER_LEAF_NODE;

		Vertices = create_bezier_vertex_block();

		Assign_Vector(Vertices->Vertices[0], (*Patch)[0][0]);
		Assign_Vector(Vertices->Vertices[1], (*Patch)[0][3]);
		Assign_Vector(Vertices->Vertices[2], (*Patch)[3][3]);
		Assign_Vector(Vertices->Vertices[3], (*Patch)[3][0]);

		Vertices->uvbnds[0] = u0;
		Vertices->uvbnds[1] = u1;
		Vertices->uvbnds[2] = v0;
		Vertices->uvbnds[3] = v1;

		Node->Data_Ptr = reinterpret_cast<void *>(Vertices);
	}
	else
	{
		if (depth >= U_Steps)
		{
			if (depth >= V_Steps)
			{
				/* We are at the max recursion depth. Just store corners. */

				Node->Node_Type = BEZIER_LEAF_NODE;

				Vertices = create_bezier_vertex_block();

				Assign_Vector(Vertices->Vertices[0], (*Patch)[0][0]);
				Assign_Vector(Vertices->Vertices[1], (*Patch)[0][3]);
				Assign_Vector(Vertices->Vertices[2], (*Patch)[3][3]);
				Assign_Vector(Vertices->Vertices[3], (*Patch)[3][0]);

				Vertices->uvbnds[0] = u0;
				Vertices->uvbnds[1] = u1;
				Vertices->uvbnds[2] = v0;
				Vertices->uvbnds[3] = v1;

				Node->Data_Ptr = reinterpret_cast<void *>(Vertices);
			}
			else
			{
				bezier_split_up_down(Patch, &Lower_Left, &Upper_Left);

				Node->Node_Type = BEZIER_INTERIOR_NODE;

				Children = create_bezier_child_block();

				Children->Children[0] = bezier_tree_builder(&Lower_Left, u0, u1, v0, (v0 + v1) / 2.0, depth + 1, max_depth_reached);
				Children->Children[1] = bezier_tree_builder(&Upper_Left, u0, u1, (v0 + v1) / 2.0, v1, depth + 1, max_depth_reached);

				Node->Count = 2;

				Node->Data_Ptr = reinterpret_cast<void *>(Children);
			}
		}
		else
		{
			if (depth >= V_Steps)
			{
				bezier_split_left_right(Patch, &Lower_Left, &Lower_Right);

				Node->Node_Type = BEZIER_INTERIOR_NODE;

				Children = create_bezier_child_block();

				Children->Children[0] = bezier_tree_builder(&Lower_Left, u0, (u0 + u1) / 2.0, v0, v1, depth + 1, max_depth_reached);
				Children->Children[1] = bezier_tree_builder(&Lower_Right, (u0 + u1) / 2.0, u1, v0, v1, depth + 1, max_depth_reached);

				Node->Count = 2;

				Node->Data_Ptr = reinterpret_cast<void *>(Children);
			}
			else
			{
				bezier_split_left_right(Patch, &Lower_Left, &Lower_Right);

				bezier_split_up_down(&Lower_Left, &Lower_Left, &Upper_Left);

				bezier_split_up_down(&Lower_Right, &Lower_Right, &Upper_Right);

				Node->Node_Type = BEZIER_INTERIOR_NODE;

				Children = create_bezier_child_block();

				Children->Children[0] = bezier_tree_builder(&Lower_Left, u0, (u0 + u1) / 2.0, v0, (v0 + v1) / 2.0, depth + 1, max_depth_reached);
				Children->Children[1] = bezier_tree_builder(&Upper_Left, u0, (u0 + u1) / 2.0, (v0 + v1) / 2.0, v1, depth + 1, max_depth_reached);
				Children->Children[2] = bezier_tree_builder(&Lower_Right, (u0 + u1) / 2.0, u1, v0, (v0 + v1) / 2.0, depth + 1, max_depth_reached);
				Children->Children[3] = bezier_tree_builder(&Upper_Right, (u0 + u1) / 2.0, u1, (v0 + v1) / 2.0, v1, depth + 1, max_depth_reached);

				Node->Count = 4;

				Node->Data_Ptr = reinterpret_cast<void *>(Children);
			}
		}
	}

	return (Node);
}



/*****************************************************************************
*
* FUNCTION
*
*   bezier_value
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
*   Determine the position and normal at a single coordinate
*   point (u, v) on a Bezier patch.
*
* CHANGES
*
*   -
*
******************************************************************************/

void BicubicPatch::bezier_value(const VECTOR (*cp)[4][4], DBL u0, DBL  v0, VECTOR P, VECTOR  N)
{
	const DBL C[] = { 1.0, 3.0, 3.0, 1.0 };
	int i, j;
	DBL c, t, ut, vt;
	DBL u[4], uu[4], v[4], vv[4];
	DBL du[4], duu[4], dv[4], dvv[4];
	DBL squared_u1, squared_v1;
	VECTOR U1, V1;

	/* Calculate binomial coefficients times coordinate positions. */

	u[0] = 1.0; uu[0] = 1.0; du[0] = 0.0; duu[0] = 0.0;
	v[0] = 1.0; vv[0] = 1.0; dv[0] = 0.0; dvv[0] = 0.0;

	for (i = 1; i < 4; i++)
	{
		u[i] = u[i - 1] * u0;  uu[i] = uu[i - 1] * (1.0 - u0);
		v[i] = v[i - 1] * v0;  vv[i] = vv[i - 1] * (1.0 - v0);

		du[i] = i * u[i - 1];  duu[i] = -i * uu[i - 1];
		dv[i] = i * v[i - 1];  dvv[i] = -i * vv[i - 1];
	}

	/* Now evaluate position and tangents based on control points. */

	Make_Vector(P, 0, 0, 0);
	Make_Vector(U1, 0, 0, 0);
	Make_Vector(V1, 0, 0, 0);

	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++)
		{
			c = C[i] * C[j];

			ut = u[i] * uu[3 - i];
			vt = v[j] * vv[3 - j];

			t = c * ut * vt;

			VAddScaledEq(P, t, (*cp)[i][j]);

			t = c * vt * (du[i] * uu[3 - i] + u[i] * duu[3 - i]);

			VAddScaledEq(U1, t, (*cp)[i][j]);

			t = c * ut * (dv[j] * vv[3 - j] + v[j] * dvv[3 - j]);

			VAddScaledEq(V1, t, (*cp)[i][j]);
		}
	}

	/* Make the normal from the cross product of the tangents. */

	VCross(N, U1, V1);

	VDot(t, N, N);

	squared_u1 = VSumSqr(U1);
	squared_v1 = VSumSqr(V1);
	if (t > (BEZIER_EPSILON * squared_u1 * squared_v1))
	{
		t = 1.0 / sqrt(t);

		VScaleEq(N, t);
	}
	else
	{
		Make_Vector(N, 1, 0, 0);
	}
}



/*****************************************************************************
*
* FUNCTION
*
*   subpatch_normal
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
*   Calculate the normal to a subpatch (triangle) return the vector
*   <1.0 0.0 0.0> if the triangle is degenerate.
*
* CHANGES
*
*   -
*
******************************************************************************/

bool BicubicPatch::subpatch_normal(const VECTOR v1, const VECTOR v2, const VECTOR v3, VECTOR Result, DBL *d)
{
	VECTOR V1, V2;
	DBL squared_v1, squared_v2;
	DBL Length;

	VSub(V1, v1, v2);
	VSub(V2, v3, v2);

	VCross(Result, V1, V2);

	Length = VSumSqr(Result);
	squared_v1 = VSumSqr(V1);
	squared_v2 = VSumSqr(V2);

	if (Length <= (BEZIER_EPSILON * squared_v1 * squared_v2))
	{
		Make_Vector(Result, 1.0, 0.0, 0.0);

		*d = -1.0 * v1[X];

		return false;
	}
	else
	{
		Length = sqrt(Length);

		VInverseScale(Result, Result, Length);

		VDot(*d, Result, v1);

		*d = 0.0 - *d;

		return true;
	}
}

/*****************************************************************************
*
* FUNCTION
*
*   intersect_subpatch
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

bool BicubicPatch::intersect_subpatch(const Ray &ray, const VECTOR V1[3], const DBL uu[3], const DBL vv[3], DBL *Depth, VECTOR P, VECTOR N, DBL *u, DBL *v) const
{
	DBL squared_b0, squared_b1;
	DBL d, n, a, b, r;
	VECTOR Q, T1;
	VECTOR B[3], IB[3], NN[3];

	VSub(B[0], V1[1], V1[0]);
	VSub(B[1], V1[2], V1[0]);

	VCross(B[2], B[0], B[1]);

	VDot(d, B[2], B[2]);

	squared_b0 = VSumSqr(B[0]);
	squared_b1 = VSumSqr(B[1]);
	if (d <= (BEZIER_EPSILON * squared_b1 * squared_b0))
	{
		return false;
	}

	d = 1.0 / sqrt(d);

	VScaleEq(B[2], d);

	/* Degenerate triangle. */

	if (!MInvers3(B, IB))
	{
		return false;
	}

	VDot(d, ray.Direction, IB[2]);

	if (fabs(d) < BEZIER_EPSILON)
	{
		return false;
	}

	VSub(Q, V1[0], ray.Origin);

	VDot(n, Q, IB[2]);

	*Depth = n / d;

	if (*Depth < BEZIER_TOLERANCE)
	{
		return false;
	}

	VScale(T1, ray.Direction, *Depth);

	VAdd(P, ray.Origin, T1);

	VSub(Q, P, V1[0]);

	VDot(a, Q, IB[0]);
	VDot(b, Q, IB[1]);

	if ((a < 0.0) || (b < 0.0) || (a + b > 1.0))
	{
		return false;
	}

	r = 1.0 - a - b;

	Make_Vector(N, 0.0, 0.0, 0.0);

	bezier_value(&Control_Points, uu[0], vv[0], T1, NN[0]);
	bezier_value(&Control_Points, uu[1], vv[1], T1, NN[1]);
	bezier_value(&Control_Points, uu[2], vv[2], T1, NN[2]);

	VScale(T1, NN[0], r); VAddEq(N, T1);
	VScale(T1, NN[1], a); VAddEq(N, T1);
	VScale(T1, NN[2], b); VAddEq(N, T1);

	*u = r * uu[0] + a * uu[1] + b * uu[2];
	*v = r * vv[0] + a * vv[1] + b * vv[2];

	VDot(d, N, N);

	if (d > BEZIER_EPSILON)
	{
		d = 1.0 / sqrt(d);

		VScaleEq(N, d);
	}
	else
	{
		Make_Vector(N, 1, 0, 0);
	}

	return true;
}



/*****************************************************************************
*
* FUNCTION
*
*   find_average
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
*   Find a sphere that contains all of the points in the list "vectors".
*
* CHANGES
*
*   -
*
******************************************************************************/

void BicubicPatch::find_average(int vector_count, const VECTOR *vectors, VECTOR center, DBL *radius)
{
	int i;
	DBL r0, r1, xc = 0, yc = 0, zc = 0;
	DBL x0, y0, z0;

	for (i = 0; i < vector_count; i++)
	{
		xc += vectors[i][X];
		yc += vectors[i][Y];
		zc += vectors[i][Z];
	}

	xc /= (DBL)vector_count;
	yc /= (DBL)vector_count;
	zc /= (DBL)vector_count;

	r0 = 0.0;

	for (i = 0; i < vector_count; i++)
	{
		x0 = vectors[i][X] - xc;
		y0 = vectors[i][Y] - yc;
		z0 = vectors[i][Z] - zc;

		r1 = x0 * x0 + y0 * y0 + z0 * z0;

		if (r1 > r0)
		{
			r0 = r1;
		}
	}

	Make_Vector(center, xc, yc, zc);

	*radius = r0;
}



/*****************************************************************************
*
* FUNCTION
*
*   spherical_bounds_check
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

bool BicubicPatch::spherical_bounds_check(const Ray &ray, const VECTOR center, DBL radius)
{
	DBL x, y, z, dist1, dist2;

	x = center[X] - ray.Origin[X];
	y = center[Y] - ray.Origin[Y];
	z = center[Z] - ray.Origin[Z];

	dist1 = x * x + y * y + z * z;

	if (dist1 < radius)
	{
		/* ray starts inside sphere - assume it intersects. */

		return true;
	}
	else
	{
		dist2 = x*ray.Direction[X] + y*ray.Direction[Y] + z*ray.Direction[Z];

		dist2 *= dist2;

		if ((dist2 > 0) && ((dist1 - dist2) <= (radius + BEZIER_EPSILON) ))
		{
			return true;
		}
	}

	return false;
}



/*****************************************************************************
*
* FUNCTION
*
*   bezier_bounding_sphere
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
*   Find a sphere that bounds all of the control points of a Bezier patch.
*   The values returned are: the center of the bounding sphere, and the
*   square of the radius of the bounding sphere.
*
* CHANGES
*
*   -
*
******************************************************************************/

void BicubicPatch::bezier_bounding_sphere(const VECTOR (*Patch)[4][4], VECTOR center, DBL *radius)
{
	int i, j;
	DBL r0, r1, xc = 0, yc = 0, zc = 0;
	DBL x0, y0, z0;

	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++)
		{
			xc += (*Patch)[i][j][X];
			yc += (*Patch)[i][j][Y];
			zc += (*Patch)[i][j][Z];
		}
	}

	xc /= 16.0;
	yc /= 16.0;
	zc /= 16.0;

	r0 = 0.0;

	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++)
		{
			x0 = (*Patch)[i][j][X] - xc;
			y0 = (*Patch)[i][j][Y] - yc;
			z0 = (*Patch)[i][j][Z] - zc;

			r1 = x0 * x0 + y0 * y0 + z0 * z0;

			if (r1 > r0)
			{
				r0 = r1;
			}
		}
	}

	Make_Vector(center, xc, yc, zc);

	*radius = r0;
}



/*****************************************************************************
*
* FUNCTION
*
*   Precompute_Patch_Values
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
*   Precompute grid points and normals for a bezier patch.
*
* CHANGES
*
*   -
*
******************************************************************************/

void BicubicPatch::Precompute_Patch_Values()
{
	int i, j;
	VECTOR cp[16];
	VECTOR(*Patch_Ptr)[4][4] = &Control_Points;
	int max_depth_reached = 0;

	/* Calculate the bounding sphere for the entire patch. */

	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++)
		{
			Assign_Vector(cp[4*i + j], Control_Points[i][j]);
		}
	}

	find_average(16, cp, Bounding_Sphere_Center, &Bounding_Sphere_Radius);

	if (Patch_Type == 1)
	{
		if (Node_Tree != NULL)
		{
			bezier_tree_deleter(Node_Tree);
		}

		Node_Tree = bezier_tree_builder(Patch_Ptr, 0.0, 1.0, 0.0, 1.0, 0, max_depth_reached);
	}
}



/*****************************************************************************
*
* FUNCTION
*
*   point_plane_distance
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
*   Determine the distance from a point to a plane.
*
* CHANGES
*
*   -
*
******************************************************************************/

DBL BicubicPatch::point_plane_distance(const VECTOR p, const VECTOR n, DBL d)
{
	DBL temp1, temp2;

	VDot(temp1, p, n);

	temp1 += d;

	VLength(temp2, n);

	if (fabs(temp2) < EPSILON)
	{
		return (0.0);
	}

	temp1 /= temp2;

	return (temp1);
}



/*****************************************************************************
*
* FUNCTION
*
*   bezier_subpatch_intersect
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

int BicubicPatch::bezier_subpatch_intersect(const Ray &ray, const VECTOR (*Patch)[4][4], DBL u0, DBL  u1, DBL  v0, DBL  v1, IStack& Depth_Stack)
{
	int cnt = 0;
	VECTOR V1[3];
	DBL u, v, Depth;
	DBL uu[3], vv[3];
	VECTOR P, N;
	UV_VECT UV;
	DBL uv_point[2], tpoint[2];

	Assign_Vector(V1[0], (*Patch)[0][0]);
	Assign_Vector(V1[1], (*Patch)[0][3]);
	Assign_Vector(V1[2], (*Patch)[3][3]);

	uu[0] = u0; uu[1] = u0; uu[2] = u1;
	vv[0] = v0; vv[1] = v1; vv[2] = v1;

	if (intersect_subpatch(ray, V1, uu, vv, &Depth, P, N, &u, &v))
	{
		/* transform current point from uv space to texture space */
		uv_point[0] = v;
		uv_point[1] = u;
		Compute_Texture_UV(uv_point, ST, tpoint);

		UV[U] = tpoint[0];
		UV[V] = tpoint[1];
		Depth_Stack->push(Intersection(Depth, P, N, UV, this));

		cnt++;
	}

	Assign_Vector(V1[1], V1[2]);
	Assign_Vector(V1[2], (*Patch)[3][0]);

	uu[1] = uu[2]; uu[2] = u1;
	vv[1] = vv[2]; vv[2] = v0;

	if (intersect_subpatch(ray, V1, uu, vv, &Depth, P, N, &u, &v))
	{
		/* transform current point from uv space to texture space */
		uv_point[0] = v;
		uv_point[1] = u;
		Compute_Texture_UV(uv_point, ST, tpoint);

		UV[U] = tpoint[0];
		UV[V] = tpoint[1];
		Depth_Stack->push(Intersection(Depth, P, N, UV, this));

		cnt++;
	}

	return (cnt);
}




/*****************************************************************************
*
* FUNCTION
*
*   bezier_split_left_right
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

void BicubicPatch::bezier_split_left_right(const VECTOR (*Patch)[4][4], VECTOR (*Left_Patch)[4][4], VECTOR (*Right_Patch)[4][4])
{
	int i, j;
	VECTOR Half;
	VECTOR Temp1[4], Temp2[4];

	for (i = 0; i < 4; i++)
	{
		Assign_Vector(Temp1[0], (*Patch)[0][i]);

		VHalf(Temp1[1], (*Patch)[0][i], (*Patch)[1][i]);
		VHalf(Half, (*Patch)[1][i], (*Patch)[2][i]);
		VHalf(Temp1[2], Temp1[1], Half);
		VHalf(Temp2[2], (*Patch)[2][i], (*Patch)[3][i]);
		VHalf(Temp2[1], Half, Temp2[2]);
		VHalf(Temp1[3], Temp1[2], Temp2[1]);

		Assign_Vector(Temp2[0], Temp1[3]);
		Assign_Vector(Temp2[3], (*Patch)[3][i]);

		for (j = 0; j < 4; j++)
		{
			Assign_Vector((*Left_Patch)[j][i], Temp1[j]);
			Assign_Vector((*Right_Patch)[j][i], Temp2[j]);
		}
	}
}



/*****************************************************************************
*
* FUNCTION
*
*   bezier_split_up_down
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

void BicubicPatch::bezier_split_up_down(const VECTOR (*Patch)[4][4], VECTOR (*Bottom_Patch)[4][4], VECTOR (*Top_Patch)[4][4])
{
	int i, j;
	VECTOR Temp1[4], Temp2[4];
	VECTOR Half;

	for (i = 0; i < 4; i++)
	{
		Assign_Vector(Temp1[0], (*Patch)[i][0]);

		VHalf(Temp1[1], (*Patch)[i][0], (*Patch)[i][1]);
		VHalf(Half, (*Patch)[i][1], (*Patch)[i][2]);
		VHalf(Temp1[2], Temp1[1], Half);
		VHalf(Temp2[2], (*Patch)[i][2], (*Patch)[i][3]);
		VHalf(Temp2[1], Half, Temp2[2]);
		VHalf(Temp1[3], Temp1[2], Temp2[1]);

		Assign_Vector(Temp2[0], Temp1[3]);
		Assign_Vector(Temp2[3], (*Patch)[i][3]);

		for (j = 0; j < 4; j++)
		{
			Assign_Vector((*Bottom_Patch)[i][j], Temp1[j]);
			Assign_Vector((*Top_Patch)[i][j]   , Temp2[j]);
		}
	}
}



/*****************************************************************************
*
* FUNCTION
*
*   determine_subpatch_flatness
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
*   See how close to a plane a subpatch is, the patch must have at least
*   three distinct vertices. A negative result from this function indicates
*   that a degenerate value of some sort was encountered.
*
* CHANGES
*
*   -
*
******************************************************************************/

DBL BicubicPatch::determine_subpatch_flatness(const VECTOR (*Patch)[4][4])
{
	int i, j;
	DBL d, dist, temp1;
	VECTOR n, TempV;
	VECTOR vertices[4];

	Assign_Vector(vertices[0], (*Patch)[0][0]);
	Assign_Vector(vertices[1], (*Patch)[0][3]);

	VSub(TempV, vertices[0], vertices[1]);

	VLength(temp1, TempV);

	if (fabs(temp1) < EPSILON)
	{
		/*
		 * Degenerate in the V direction for U = 0. This is ok if the other
		 * two corners are distinct from the lower left corner - I'm sure there
		 * are cases where the corners coincide and the middle has good values,
		 * but that is somewhat pathalogical and won't be considered.
		 */

		Assign_Vector(vertices[1], (*Patch)[3][3]);

		VSub(TempV, vertices[0], vertices[1]);

		VLength(temp1, TempV);

		if (fabs(temp1) < EPSILON)
		{
			return (-1.0);
		}

		Assign_Vector(vertices[2], (*Patch)[3][0]);

		VSub(TempV, vertices[0], vertices[1]);

		VLength(temp1, TempV);

		if (fabs(temp1) < EPSILON)
		{
			return (-1.0);
		}

		VSub(TempV, vertices[1], vertices[2]);

		VLength(temp1, TempV);

		if (fabs(temp1) < EPSILON)
		{
			return (-1.0);
		}
	}
	else
	{
		Assign_Vector(vertices[2], (*Patch)[3][0]);

		VSub(TempV, vertices[0], vertices[1]);

		VLength(temp1, TempV);

		if (fabs(temp1) < EPSILON)
		{
			Assign_Vector(vertices[2], (*Patch)[3][3]);

			VSub(TempV, vertices[0], vertices[2]);

			VLength(temp1, TempV);

			if (fabs(temp1) < EPSILON)
			{
				return (-1.0);
			}

			VSub(TempV, vertices[1], vertices[2]);

			VLength(temp1, TempV);

			if (fabs(temp1) < EPSILON)
			{
				return (-1.0);
			}
		}
		else
		{
			VSub(TempV, vertices[1], vertices[2]);

			VLength(temp1, TempV);

			if (fabs(temp1) < EPSILON)
			{
				return (-1.0);
			}
		}
	}

	/*
	 * Now that a good set of candidate points has been found,
	 * find the plane equations for the patch.
	 */

	if (subpatch_normal(vertices[0], vertices[1], vertices[2], n, &d))
	{
		/*
		 * Step through all vertices and see what the maximum
		 * distance from the plane happens to be.
		 */

		dist = 0.0;

		for (i = 0; i < 4; i++)
		{
			for (j = 0; j < 4; j++)
			{
				temp1 = fabs(point_plane_distance(((*Patch)[i][j]), n, d));

				if (temp1 > dist)
				{
					dist = temp1;
				}
			}
		}

		return (dist);
	}
	else
	{
/*
		Debug_Info("Subpatch normal failed in determine_subpatch_flatness\n");
*/

		return (-1.0);
	}
}



/*****************************************************************************
*
* FUNCTION
*
*   flat_enough
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

bool BicubicPatch::flat_enough(const VECTOR (*Patch)[4][4]) const
{
	DBL Dist;

	Dist = determine_subpatch_flatness(Patch);

	if (Dist < 0.0)
	{
		return false;
	}
	else
	{
		if (Dist < Flatness_Value)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
}



/*****************************************************************************
*
* FUNCTION
*
*   bezier_subdivider
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

int BicubicPatch::bezier_subdivider(const Ray &ray, const VECTOR (*Patch)[4][4], DBL u0, DBL  u1, DBL  v0, DBL  v1, int recursion_depth, IStack& Depth_Stack)
{
	int cnt = 0;
	DBL ut, vt, radius;
	VECTOR Lower_Left[4][4], Lower_Right[4][4];
	VECTOR Upper_Left[4][4], Upper_Right[4][4];
	VECTOR center;

	/*
	 * Make sure the ray passes through a sphere bounding
	 * the control points of the patch.
	 */

	bezier_bounding_sphere(Patch, center, &radius);

	if (!spherical_bounds_check(ray, center, radius))
	{
		return (0);
	}

	/*
	 * If the patch is close to being flat, then just
	 * perform a ray-plane intersection test.
	 */

	if (flat_enough(Patch))
		return bezier_subpatch_intersect(ray, Patch, u0, u1, v0, v1, Depth_Stack);

	if (recursion_depth >= U_Steps)
	{
		if (recursion_depth >= V_Steps)
		{
			return bezier_subpatch_intersect(ray, Patch, u0, u1, v0, v1, Depth_Stack);
		}
		else
		{
			bezier_split_up_down(Patch, &Lower_Left, &Upper_Left);

			vt = (v1 + v0) / 2.0;

			cnt += bezier_subdivider(ray, &Lower_Left, u0, u1, v0, vt, recursion_depth + 1, Depth_Stack);
			cnt += bezier_subdivider(ray, &Upper_Left, u0, u1, vt, v1, recursion_depth + 1, Depth_Stack);
		}
	}
	else
	{
		if (recursion_depth >= V_Steps)
		{
			bezier_split_left_right(Patch, &Lower_Left, &Lower_Right);

			ut = (u1 + u0) / 2.0;

			cnt += bezier_subdivider(ray, &Lower_Left, u0, ut, v0, v1, recursion_depth + 1, Depth_Stack);
			cnt += bezier_subdivider(ray, &Lower_Right, ut, u1, v0, v1, recursion_depth + 1, Depth_Stack);
		}
		else
		{
			ut = (u1 + u0) / 2.0;
			vt = (v1 + v0) / 2.0;

			bezier_split_left_right(Patch, &Lower_Left, &Lower_Right);
			bezier_split_up_down(&Lower_Left, &Lower_Left, &Upper_Left) ;
			bezier_split_up_down(&Lower_Right, &Lower_Right, &Upper_Right);

			cnt += bezier_subdivider(ray, &Lower_Left, u0, ut, v0, vt, recursion_depth + 1, Depth_Stack);
			cnt += bezier_subdivider(ray, &Upper_Left, u0, ut, vt, v1, recursion_depth + 1, Depth_Stack);
			cnt += bezier_subdivider(ray, &Lower_Right, ut, u1, v0, vt, recursion_depth + 1, Depth_Stack);
			cnt += bezier_subdivider(ray, &Upper_Right, ut, u1, vt, v1, recursion_depth + 1, Depth_Stack);
		}
	}

	return (cnt);
}



/*****************************************************************************
*
* FUNCTION
*
*   bezier_tree_deleter
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

void BicubicPatch::bezier_tree_deleter(BEZIER_NODE *Node)
{
	int i;
	BEZIER_CHILDREN *Children;

	/* If this is an interior node then continue the descent. */

	if (Node->Node_Type == BEZIER_INTERIOR_NODE)
	{
		Children = reinterpret_cast<BEZIER_CHILDREN *>(Node->Data_Ptr);

		for (i = 0; i < Node->Count; i++)
		{
			bezier_tree_deleter(Children->Children[i]);
		}

		POV_FREE(Children);
	}
	else
	{
		if (Node->Node_Type == BEZIER_LEAF_NODE)
		{
			/* Free the memory used for the vertices. */

			POV_FREE(Node->Data_Ptr);
		}
	}

	/* Free the memory used for the node. */

	POV_FREE(Node);
}



/*****************************************************************************
*
* FUNCTION
*
*   bezier_tree_walker
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

int BicubicPatch::bezier_tree_walker(const Ray &ray, const BEZIER_NODE *Node, IStack& Depth_Stack)
{
	int i, cnt = 0;
	DBL Depth, u, v;
	DBL uu[3], vv[3];
	VECTOR N, P;
	VECTOR V1[3];
	UV_VECT UV;
	DBL uv_point[2], tpoint[2];
	const BEZIER_CHILDREN *Children;
	const BEZIER_VERTICES *Vertices;

	/*
	 * Make sure the ray passes through a sphere bounding
	 * the control points of the patch.
	 */

	if (!spherical_bounds_check(ray, Node->Center, Node->Radius_Squared))
	{
		return (0);
	}

	/*
	 * If this is an interior node then continue the descent,
	 * else do a check against the vertices.
	 */

	if (Node->Node_Type == BEZIER_INTERIOR_NODE)
	{
		Children = reinterpret_cast<const BEZIER_CHILDREN *>(Node->Data_Ptr);

		for (i = 0; i < Node->Count; i++)
		{
			cnt += bezier_tree_walker(ray, Children->Children[i], Depth_Stack);
		}
	}
	else if (Node->Node_Type == BEZIER_LEAF_NODE)
	{
		Vertices = reinterpret_cast<const BEZIER_VERTICES *>(Node->Data_Ptr);

		Assign_Vector(V1[0], Vertices->Vertices[0]);
		Assign_Vector(V1[1], Vertices->Vertices[1]);
		Assign_Vector(V1[2], Vertices->Vertices[2]);

		uu[0] = Vertices->uvbnds[0];
		uu[1] = Vertices->uvbnds[0];
		uu[2] = Vertices->uvbnds[1];
		vv[0] = Vertices->uvbnds[2];
		vv[1] = Vertices->uvbnds[3];
		vv[2] = Vertices->uvbnds[3];

		/*
		 * Triangulate this subpatch, then check for
		 * intersections in the triangles.
		 */

		if (intersect_subpatch( ray, V1, uu, vv, &Depth, P, N, &u, &v))
		{
			/* transform current point from uv space to texture space */
			uv_point[0] = v;
			uv_point[1] = u;
			Compute_Texture_UV(uv_point, ST, tpoint);

			UV[U] = tpoint[0];
			UV[V] = tpoint[1];
			Depth_Stack->push(Intersection(Depth, P, N, UV, this));

			cnt++;
		}

		Assign_Vector(V1[1], V1[2]);
		Assign_Vector(V1[2], Vertices->Vertices[3]);

		uu[1] = uu[2]; uu[2] = Vertices->uvbnds[1];
		vv[1] = vv[2]; vv[2] = Vertices->uvbnds[2];

		if (intersect_subpatch(ray, V1, uu, vv, &Depth, P, N, &u, &v))
		{
			/* transform current point from object space to texture space */
			uv_point[0] = v;
			uv_point[1] = u;
			Compute_Texture_UV(uv_point, ST, tpoint);

			UV[U] = tpoint[0];
			UV[V] = tpoint[1];
			Depth_Stack->push(Intersection(Depth, P, N, UV, this));

			cnt++;
		}
	}
	else
	{
		throw POV_EXCEPTION_STRING("Bad Node type in bezier_tree_walker().");
	}

	return (cnt);
}



/*****************************************************************************
*
* FUNCTION
*
*   intersect_bicubic_patch0
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

int BicubicPatch::intersect_bicubic_patch0(const Ray &ray, IStack& Depth_Stack)
{
	const VECTOR(*Patch)[4][4] = &Control_Points;

	return (bezier_subdivider(ray, Patch, 0.0, 1.0, 0.0, 1.0, 0, Depth_Stack));
}



/*****************************************************************************
*
* FUNCTION
*
*   All_Bicubic_Patch_Intersections
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

bool BicubicPatch::All_Intersections(const Ray& ray, IStack& Depth_Stack, TraceThreadData *Thread)
{
	int Found, cnt = 0;

	Found = false;

	Thread->Stats()[Ray_Bicubic_Tests]++;

	switch (Patch_Type)
	{
		case 0:

			cnt = intersect_bicubic_patch0(ray, Depth_Stack);

			break;

		case 1:

			cnt = bezier_tree_walker(ray, Node_Tree, Depth_Stack);

			break;

		default:

			throw POV_EXCEPTION_STRING("Bad patch type in All_Bicubic_Patch_Intersections.");
	}

	if (cnt > 0)
	{
		Thread->Stats()[Ray_Bicubic_Tests_Succeeded]++;

		Found = true;
	}

	return (Found);
}



/*****************************************************************************
*
* FUNCTION
*
*   Inside_Bicubic_Patch
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
*   A patch is not a solid, so an inside test doesn't make sense.
*
* CHANGES
*
*   -
*
******************************************************************************/

bool BicubicPatch::Inside(const VECTOR, TraceThreadData *) const
{
	return false;
}



/*****************************************************************************
*
* FUNCTION
*
*   Bicubic_Patch_Normal
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

void BicubicPatch::Normal(VECTOR Result, Intersection *Inter, TraceThreadData *Thread) const
{
	/* Use preocmputed normal. */

	Assign_Vector(Result, Inter->INormal);
}



/*****************************************************************************
*
* FUNCTION
*
*   Translate_Bicubic_Patch
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

void BicubicPatch::Translate(const VECTOR Vector, const TRANSFORM *)
{
	int i, j;

	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++)
		{
			VAdd(Control_Points[i][j], Control_Points[i][j], Vector);
		}
	}

	Precompute_Patch_Values();

	Compute_BBox();
}



/*****************************************************************************
*
* FUNCTION
*
*   Rotate_Bicubic_Patch
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

void BicubicPatch::Rotate(const VECTOR, const TRANSFORM *tr)
{
	Transform(tr);
}



/*****************************************************************************
*
* FUNCTION
*
*   Scale_Bicubic_Patch
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

void BicubicPatch::Scale(const VECTOR Vector, const TRANSFORM *)
{
	int i, j;

	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++)
		{
			VEvaluate(Control_Points[i][j], Control_Points[i][j], Vector);
		}
	}

	Precompute_Patch_Values();

	Compute_BBox();
}




/*****************************************************************************
*
* FUNCTION
*
*   Transform_Bicubic_Patch
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

void BicubicPatch::Transform(const TRANSFORM *tr)
{
	int i, j;

	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++)
		{
			MTransPoint(Control_Points[i][j], Control_Points[i][j], tr);
		}
	}

	Precompute_Patch_Values();

	Compute_BBox();
}



/*****************************************************************************
*
* FUNCTION
*
*   Invert_Bicubic_Patch
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
*   Inversion of a patch really doesn't make sense.
*
* CHANGES
*
*   -
*
******************************************************************************/

void BicubicPatch::Invert()
{
}



/*****************************************************************************
*
* FUNCTION
*
*   Create_Bicubic_Patch
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

BicubicPatch::BicubicPatch() : ObjectBase(BICUBIC_PATCH_OBJECT)
{
	Patch_Type = - 1;

	U_Steps = 0;
	V_Steps = 0;

	Flatness_Value = 0.0;
	accuracy = 0.01;

	Node_Tree = NULL;
	Weights = NULL;

	/*
	 * NOTE: Control_Points[4][4] is initialized in Parse_Bicubic_Patch.
	 * Bounding_Sphere_Center,Bounding_Sphere_Radius, Normal_Vector[], and
	 * IPoint[] are initialized in Precompute_Patch_Values.
	 */

	/* set the default uv-mapping coordinates */
	ST[0][U] = 0;
	ST[0][V] = 0;
	ST[1][U] = 1;
	ST[1][V] = 0;
	ST[2][U] = 1;
	ST[2][V] = 1;
	ST[3][U] = 0;
	ST[3][V] = 1;
}



/*****************************************************************************
*
* FUNCTION
*
*   Copy_Bicubic_Patch
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

ObjectPtr BicubicPatch::Copy()
{
	int i, j;
	BicubicPatch *New = new BicubicPatch();
	int m, h;

	/* Do not do *New = *Old so that Precompute works right */

	New->Patch_Type = Patch_Type;

	New->U_Steps = U_Steps;
	New->V_Steps = V_Steps;

	if ( Weights != NULL )
	{
		New->Weights = reinterpret_cast<WEIGHTS *>(POV_MALLOC( sizeof(WEIGHTS),"bicubic patch" ));
		POV_MEMCPY( New->Weights, Weights, sizeof(WEIGHTS) );
	}

	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++)
		{
			Assign_Vector(New->Control_Points[i][j], Control_Points[i][j]);
		}
	}

	New->Flatness_Value = Flatness_Value;

	New->Precompute_Patch_Values();

	/* copy the mapping */
	for (m = 0; m < 4; m++)
	{
		for (h = 0; h < 3; h++)
		{
			New->ST[m][h] = ST[m][h];
		}
	}

	return (New);
}



/*****************************************************************************
*
* FUNCTION
*
*   Destroy_Bicubic_Patch
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

BicubicPatch::~BicubicPatch()
{
	if (Patch_Type == 1)
	{
		if (Node_Tree != NULL)
		{
			bezier_tree_deleter(Node_Tree);
		}
	}

	if ( Weights != NULL ) POV_FREE(Weights);
}



/*****************************************************************************
*
* FUNCTION
*
*   Compute_Bicubic_Patch_BBox
*
* INPUT
*
*   Bicubic_Patch - Bicubic patch
*   
* OUTPUT
*
*   Bicubic_Patch
*   
* RETURNS
*   
* AUTHOR
*
*   Dieter Bayer
*   
* DESCRIPTION
*
*   Calculate the bounding box of a bicubic patch.
*
* CHANGES
*
*   Aug 1994 : Creation.
*
******************************************************************************/

void BicubicPatch::Compute_BBox()
{
	int i, j;
	VECTOR Min, Max;

	Make_Vector(Min, BOUND_HUGE, BOUND_HUGE, BOUND_HUGE);
	Make_Vector(Max, -BOUND_HUGE, -BOUND_HUGE, -BOUND_HUGE);

	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++)
		{
			Min[X] = min(Min[X], Control_Points[i][j][X]);
			Min[Y] = min(Min[Y], Control_Points[i][j][Y]);
			Min[Z] = min(Min[Z], Control_Points[i][j][Z]);
			Max[X] = max(Max[X], Control_Points[i][j][X]);
			Max[Y] = max(Max[Y], Control_Points[i][j][Y]);
			Max[Z] = max(Max[Z], Control_Points[i][j][Z]);
		}
	}

	Make_BBox_from_min_max(BBox, Min, Max);
}

/*****************************************************************************
*
* FUNCTION
*
*   Bicubic_Patch_UVCoord
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Nathan Kopp
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

void BicubicPatch::UVCoord(UV_VECT Result, const Intersection *Inter, TraceThreadData *Thread) const
{
	/* Use preocmputed uv coordinates. */

	Assign_UV_Vect(Result, Inter->Iuv);
}


/*****************************************************************************
*
* FUNCTION
*
*   Compute_UV_Point
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Mike Hough
*
* DESCRIPTION
*
*   Transform p from uv space to texture space (point t) using the
*   shape's ST mapping
*
* CHANGES
*
*   -
*
******************************************************************************/

void BicubicPatch::Compute_Texture_UV(const UV_VECT p, const UV_VECT st[4], UV_VECT t)
{
	UV_VECT u1, u2;

	u1[0] = st[0][0] + p[0] * (st[1][0] - st[0][0]);
	u1[1] = st[0][1] + p[0] * (st[1][1] - st[0][1]);

	u2[0] = st[3][0] + p[0] * (st[2][0] - st[3][0]);
	u2[1] = st[3][1] + p[0] * (st[2][1] - st[3][1]);

	t[0] = u1[0] + p[1] * (u2[0] - u1[0]);
	t[1] = u1[1] + p[1] * (u2[1] - u1[1]);
}

}
