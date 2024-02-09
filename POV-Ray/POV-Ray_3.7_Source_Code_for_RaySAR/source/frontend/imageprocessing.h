/*******************************************************************************
 * imageprocessing.h
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
 * $File: //depot/public/povray/3.x/source/frontend/imageprocessing.h $
 * $Revision: #1 $
 * $Change: 6069 $
 * $DateTime: 2013/11/06 11:59:40 $
 * $Author: chrisc $
 *******************************************************************************/

#ifndef POVRAY_FRONTEND_IMAGEPROCESSING_H
#define POVRAY_FRONTEND_IMAGEPROCESSING_H

#include "base/povmscpp.h"
#include "base/povmsgid.h"
#include "base/image/image.h"
#include "base/fileinputoutput.h"

#include "frontend/configfrontend.h"

#include <string>

#include <boost/scoped_ptr.hpp>

namespace pov_frontend
{

using namespace pov_base;

class ImageProcessing
{
	public:
		ImageProcessing(unsigned int width, unsigned int height);
		ImageProcessing(POVMS_Object& ropts);
		ImageProcessing(boost::shared_ptr<Image>& img);
		virtual ~ImageProcessing();

		UCS2String WriteImage(POVMS_Object& ropts, POVMSInt frame = 0, int digits = 0);

		boost::shared_ptr<Image>& GetImage();

		UCS2String GetOutputFilename(POVMS_Object& ropts, POVMSInt frame, int digits);
		bool OutputIsStdout(void) { return toStdout; }
		bool OutputIsStderr(void) { return toStderr; }
		virtual bool OutputIsStdout(POVMS_Object& ropts);
		virtual bool OutputIsStderr(POVMS_Object& ropts);
	protected:
		boost::shared_ptr<Image> image;
		bool toStdout;
		bool toStderr;

		void RGB2XYZ(const COLC *rgb, COLC *xyz);
		void XYZ2RGB(const COLC *xyz, COLC *rgb);
	private:
		ImageProcessing();
		ImageProcessing(const ImageProcessing&);
		ImageProcessing& operator=(const ImageProcessing&);
};

}

#endif // POVRAY_FRONTEND_IMAGEPROCESSING_H
