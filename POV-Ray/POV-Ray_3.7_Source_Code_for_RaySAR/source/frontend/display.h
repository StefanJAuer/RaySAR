/*******************************************************************************
 * display.h
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
 * $File: //depot/public/povray/3.x/source/frontend/display.h $
 * $Revision: #1 $
 * $Change: 6069 $
 * $DateTime: 2013/11/06 11:59:40 $
 * $Author: chrisc $
 *******************************************************************************/

#ifndef POVRAY_FRONTEND_DISPLAY_H
#define POVRAY_FRONTEND_DISPLAY_H

#include "base/types.h"
#include "base/image/colourspace.h"

#include "frontend/configfrontend.h"

namespace pov_frontend
{

class Display
{
	public:
		struct RGBA8 { unsigned char red, green, blue, alpha; };

		Display(unsigned int w, unsigned int h, pov_base::GammaCurvePtr g);
		virtual ~Display();

		virtual void Initialise() = 0;

		unsigned int GetWidth();
		unsigned int GetHeight();

		pov_base::GammaCurvePtr GetGamma();

		virtual void DrawPixel(unsigned int x, unsigned int y, const RGBA8& colour) = 0;

		virtual void DrawRectangleFrame(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2, const RGBA8& colour);
		virtual void DrawFilledRectangle(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2, const RGBA8& colour);

		virtual void DrawPixelBlock(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2, const RGBA8 *colour);

		virtual void Clear();
	private:
		/// display width
		unsigned int width;
		/// display height
		unsigned int height;
		/// display gamma correction factor
		pov_base::GammaCurvePtr gamma;

		/// not available
		Display();
};

}

#endif // POVRAY_FRONTEND_DISPLAY_H
