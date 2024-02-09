/*******************************************************************************
 * camera.cpp
 *
 * This module implements methods for managing the viewpoint.
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
 * $File: //depot/public/povray/3.x/source/backend/scene/camera.cpp $
 * $Revision: #1 $
 * $Change: 6069 $
 * $DateTime: 2013/11/06 11:59:40 $
 * $Author: chrisc $
 *******************************************************************************/

// frame.h must always be the first POV file included (pulls in platform config)
#include "backend/frame.h"
#include "backend/scene/camera.h"
#include "backend/scene/objects.h"
#include "backend/texture/normal.h"
#include "backend/texture/pigment.h"
#include "backend/math/vector.h"
#include "backend/math/matrices.h"

// this must be the last file included
#include "base/povdebug.h"

namespace pov
{

/*****************************************************************************
*
* FUNCTION
*
*   Translate_Camera
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

void Camera::Translate(const VECTOR Vector)
{
	VAddEq(Location, Vector);
}



/*****************************************************************************
*
* FUNCTION
*
*   Rotate_Camera
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

void Camera::Rotate(const VECTOR Vector)
{
	TRANSFORM Trans;

	Compute_Rotation_Transform(&Trans, Vector);
	Transform(&Trans);
}



/*****************************************************************************
*
* FUNCTION
*
*   Scale_Camera
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

void Camera::Scale(const VECTOR Vector)
{
	TRANSFORM Trans;

	Compute_Scaling_Transform(&Trans, Vector);
	Transform(&Trans);
}



/*****************************************************************************
*
* FUNCTION
*
*   Transform_Camera
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

void Camera::Transform(const TRANSFORM *Trans)
{
	MTransPoint(Location, Location, Trans);
	MTransDirection(Direction, Direction, Trans);
	MTransDirection(Up, Up, Trans);
	MTransDirection(Right, Right, Trans);
}



/*****************************************************************************
*
* METHOD
*
*   Camera::Init
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

void Camera::Init()
{
	Make_Vector(Location,    0.0,  0.0, 0.0);
	Make_Vector(Direction,   0.0,  0.0, 1.0);
	Make_Vector(Up,          0.0,  1.0, 0.0);
	Make_Vector(Right,       1.33, 0.0, 0.0);
	Make_Vector(Sky,         0.0,  1.0, 0.0);
	Make_Vector(Look_At,     0.0,  0.0, 1.0);
	Make_Vector(Focal_Point, 0.0,  0.0, 1.0);

	/* Init focal blur stuff (not used by default). */
	Blur_Samples        = 0;
	Blur_Samples_Min    = 0;
	Confidence          = 0.9;
	Variance            = 1.0 / 10000.0;
	Aperture            = 0.0;
	Focal_Distance      = -1.0;

	/* Set default camera type and viewing angle. [DB 7/94] */
	Type = PERSPECTIVE_CAMERA;
	Angle = 90.0;

	/* Default view angle for spherical camera. [MH 6/99] */
	H_Angle = 360;
	V_Angle = 180;

	/* Do not perturb primary rays by default. [DB 7/94] */
	Tnormal = NULL;

	Bokeh = NULL; // no user-defined bokeh by default

	Trans = Create_Transform();

	Rays_Per_Pixel = 1;
	Face_Distribution_Method = 0;
	Smooth = false;
	Max_Ray_Distance = 0.0;
}

/*****************************************************************************
*
* FUNCTION
*
*   Create_Camera
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

Camera::Camera()
{
	Init();
}



/*****************************************************************************
*
* FUNCTION
*
*   Copy_Camera
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

Camera& Camera::operator=(const Camera& src)
{
	Assign_Vector(Location, src.Location);
	Assign_Vector(Direction, src.Direction);
	Assign_Vector(Up, src.Up);
	Assign_Vector(Right, src.Right);
	Assign_Vector(Sky, src.Sky);
	Assign_Vector(Look_At, src.Look_At);
	Assign_Vector(Focal_Point, src.Focal_Point);

	Focal_Distance = src.Focal_Distance;
	Aperture = src.Aperture;
	Blur_Samples = src.Blur_Samples;
	Blur_Samples_Min = src.Blur_Samples_Min;
	Confidence = src.Confidence;
	Variance = src.Variance;
	Type = src.Type;
	Angle = src.Angle;
	H_Angle = src.H_Angle;
	V_Angle = src.V_Angle;

	if (Tnormal != NULL)
		Destroy_Tnormal(Tnormal);
	Tnormal = src.Tnormal ? Copy_Tnormal(src.Tnormal) : NULL;
	if (Trans != NULL)
		Destroy_Transform(Trans);
	Trans = src.Trans ? Copy_Transform(src.Trans) : NULL;

	if (Bokeh != NULL)
		Destroy_Pigment(Bokeh);
	Bokeh = src.Bokeh ? Copy_Pigment(src.Bokeh) : NULL;

	for (std::vector<ObjectPtr>::iterator it = Meshes.begin(); it != Meshes.end(); it++)
		Destroy_Object(*it);
	Meshes.clear();
	for (std::vector<ObjectPtr>::const_iterator it = src.Meshes.begin(); it != src.Meshes.end(); it++)
		Meshes.push_back(Copy_Object(*it));
	Face_Distribution_Method = src.Face_Distribution_Method;
	Rays_Per_Pixel = src.Rays_Per_Pixel;
	Max_Ray_Distance = src.Max_Ray_Distance;
	Mesh_Index = src.Mesh_Index;
	for (int i = 0; i < 10; i++)
	{
		U_Xref[i] = src.U_Xref[i];
		V_Xref[i] = src.V_Xref[i];
	}
	Smooth = src.Smooth;

	return *this;
}

Camera::Camera(const Camera& src)
{
	Tnormal = NULL;
	Trans = NULL;
	Bokeh = NULL;
	operator=(src);
}

/*****************************************************************************
*
* FUNCTION
*
*   Destroy_Camera
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
*   -
*
* CHANGES
*
*   -
*
******************************************************************************/

Camera::~Camera()
{
	Destroy_Tnormal(Tnormal);
	Destroy_Transform(Trans);
	Destroy_Pigment(Bokeh);
	for (std::vector<ObjectPtr>::iterator it = Meshes.begin(); it != Meshes.end(); it++)
		Destroy_Object(*it);
	Meshes.clear();
}

}
