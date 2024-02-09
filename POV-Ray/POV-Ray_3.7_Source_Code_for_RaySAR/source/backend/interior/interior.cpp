/*******************************************************************************
 * interior.cpp
 *
 * This module contains all functions for interior stuff.
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
 * $File: //depot/public/povray/3.x/source/backend/interior/interior.cpp $
 * $Revision: #1 $
 * $Change: 6069 $
 * $DateTime: 2013/11/06 11:59:40 $
 * $Author: chrisc $
 *******************************************************************************/

// frame.h must always be the first POV file included (pulls in platform config)
#include "backend/frame.h"
#include "backend/interior/interior.h"
#include "backend/texture/texture.h"
#include "backend/lighting/subsurface.h"

// this must be the last file included
#include "base/povdebug.h"

namespace pov
{

/* How many subrays to trace for dispersive media */
#define DEFAULT_DISP_NELEMS  7

Interior::Interior()
{
	References = 1;

	IOR = 0.0;
	Old_Refract = 1.0;

	Dispersion  = 1.0;
	Disp_NElems = DEFAULT_DISP_NELEMS;

	Caustics = 0.0;

	Fade_Distance = 0.0;
	Fade_Power    = 0.0;

	hollow = false;

	subsurface = boost::shared_ptr<SubsurfaceInterior>();
}

Interior::Interior(const Interior& source)
{
	References = 1;
	Disp_NElems = source.Disp_NElems;
	Dispersion = source.Dispersion;
	Old_Refract = source.Old_Refract;
	Fade_Distance = source.Fade_Distance;
	Fade_Power = source.Fade_Power;
	Fade_Colour = source.Fade_Colour;
	media = source.media;
	hollow = source.hollow;
	IOR = source.IOR;
	subsurface = boost::shared_ptr<SubsurfaceInterior>(source.subsurface);
	Caustics = source.Caustics;
}

Interior::~Interior()
{
}

void Interior::Transform(const TRANSFORM *trans)
{
	for(vector<Media>::iterator i(media.begin());i != media.end(); i++)
		i->Transform(trans);
}

void Interior::PostProcess()
{
	for(vector<Media>::iterator i(media.begin());i != media.end(); i++)
		i->PostProcess();
}

/*****************************************************************************
*
* FUNCTION
*
*   Destroy_Interior
*
* INPUT
*
*   Interior - interior to destroy
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
*   Destroy an interior.
*
* CHANGES
*
*   Dec 1994 : Creation.
*
******************************************************************************/

void Destroy_Interior(Interior *interior)
{
	if((interior != NULL) && (--(interior->References) == 0))
		delete interior;
}

/*****************************************************************************
*
* FUNCTION
*
*   Copy_Interior_Pointer
*
* INPUT
*
*   Old - interior to copy
*
* OUTPUT
*
* RETURNS
*
*   INTERIOR * - new interior
*
* AUTHOR
*
*   Dieter Bayer
*
* DESCRIPTION
*
*   Copy an interior by increasing number of references.
*
* CHANGES
*
*   Dec 1994 : Creation.
*
******************************************************************************/

Interior *Copy_Interior_Pointer(Interior *Old)
{
	if (Old != NULL)
	{
		Old->References++;
	}

	return(Old);
}

MATERIAL *Create_Material()
{
	MATERIAL *New;

	New = (MATERIAL *)POV_MALLOC(sizeof(MATERIAL), "material");

	New->Texture  = NULL;
	New->Interior_Texture  = NULL;
	New->interior = NULL;

	return(New);
}

MATERIAL *Copy_Material(const MATERIAL *Old)
{
	MATERIAL *New;

	if (Old != NULL)
	{
		New = Create_Material();

		*New = *Old;

		New->Texture  = Copy_Textures(Old->Texture);
		New->Interior_Texture  = Copy_Textures(Old->Interior_Texture);
		if (Old->interior != NULL)
			New->interior = new Interior(*(Old->interior));

		return(New);
	}
	else
	{
		return(NULL);
	}
}

void Destroy_Material(MATERIAL *Material)
{
	if (Material != NULL)
	{
		Destroy_Textures(Material->Texture);
		Destroy_Textures(Material->Interior_Texture);
		Destroy_Interior(Material->interior);

		POV_FREE(Material);
	}
}

}
