/*******************************************************************************
 * imageprocessing.cpp
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
 * $File: //depot/public/povray/3.x/source/frontend/imageprocessing.cpp $
 * $Revision: #1 $
 * $Change: 6069 $
 * $DateTime: 2013/11/06 11:59:40 $
 * $Author: chrisc $
 *******************************************************************************/

#include <string>
#include <cctype>

#include <boost/scoped_ptr.hpp>

// configbase.h must always be the first POV file included within base *.cpp files
#include "base/configbase.h"
#include "base/types.h"
#include "base/image/encoding.h"

#include "frontend/imageprocessing.h"

// this must be the last file included
#include "base/povdebug.h"

// TODO: update ImageProcessing with the means of accepting and caching
// blocks of pixels as opposed to individual ones, with a back-end that
// can serialize completed rows to the final image output file.

namespace pov_frontend
{

using namespace pov;

enum
{
	X = 0,
	Y = 1,
	Z = 2
};

ImageProcessing::ImageProcessing(unsigned int width, unsigned int height)
{
	image = boost::shared_ptr<Image>(Image::Create(width, height, Image::RGBFT_Float));
	toStderr = toStdout = false;

	// TODO FIXME - find a better place for this
	image->SetPremultiplied(true); // POV-Ray uses premultiplied opacity for its math, so that's what will end up in the image container
}

ImageProcessing::ImageProcessing(POVMS_Object& ropts)
{
	unsigned int width(ropts.TryGetInt(kPOVAttrib_Width, 160));
	unsigned int height(ropts.TryGetInt(kPOVAttrib_Height, 120));
	unsigned int blockSize(ropts.TryGetInt(kPOVAttrib_RenderBlockSize, 32));
	unsigned int maxBufferMem(ropts.TryGetInt(kPOVAttrib_MaxImageBufferMem, 128)); // number is megabytes

	image = boost::shared_ptr<Image>(Image::Create(width, height, Image::RGBFT_Float, maxBufferMem, blockSize * blockSize));
	toStdout = OutputIsStdout(ropts);
	toStderr = OutputIsStderr(ropts);

	// TODO FIXME - find a better place for this
	image->SetPremultiplied(true); // POV-Ray uses premultiplied opacity for its math, so that's what will end up in the image container
}

ImageProcessing::ImageProcessing(boost::shared_ptr<Image>& img)
{
	image = img;
	toStderr = toStdout = false;

	// TODO FIXME - find a better place for this
	image->SetPremultiplied(true); // POV-Ray uses premultiplied opacity for its math, so that's what will end up in the image container
}

ImageProcessing::~ImageProcessing()
{
}

UCS2String ImageProcessing::WriteImage(POVMS_Object& ropts, POVMSInt frame, int digits)
{
	if(ropts.TryGetBool(kPOVAttrib_OutputToFile, true) == true)
	{
		Image::WriteOptions wopts;
		Image::ImageFileType imagetype = Image::SYS;
		unsigned int filetype = POV_File_Image_System;

		wopts.bpcc = clip(ropts.TryGetInt(kPOVAttrib_BitsPerColor, 8), 5, 16);
		wopts.alphachannel = ropts.TryGetBool(kPOVAttrib_OutputAlpha, false);
		wopts.compress = clip(ropts.TryGetInt(kPOVAttrib_Compression, 0), 0, 255);
		wopts.grayscale = ropts.TryGetBool(kPOVAttrib_GrayscaleOutput, false);

		switch(ropts.TryGetInt(kPOVAttrib_OutputFileType, DEFAULT_OUTPUT_FORMAT))
		{
			case kPOVList_FileType_Targa:
				imagetype = Image::TGA;
				filetype = POV_File_Image_Targa;
				wopts.compress = 0;
				break;
			case kPOVList_FileType_CompressedTarga:
				imagetype = Image::TGA;
				filetype = POV_File_Image_Targa;
				wopts.compress = 1;
				break;
			case kPOVList_FileType_PNG:
				imagetype = Image::PNG;
				filetype = POV_File_Image_PNG;
				break;
			case kPOVList_FileType_JPEG:
				imagetype = Image::JPEG;
				filetype = POV_File_Image_JPEG;
				wopts.compress = clip(int(wopts.compress), 0, 100);
				break;
			case kPOVList_FileType_PPM:
				imagetype = Image::PPM;
				filetype = POV_File_Image_PPM;
				break;
			case kPOVList_FileType_BMP:
				imagetype = Image::BMP;
				filetype = POV_File_Image_BMP;
				break;
			case kPOVList_FileType_OpenEXR:
				imagetype = Image::EXR;
				filetype = POV_File_Image_EXR;
				break;
			case kPOVList_FileType_RadianceHDR:
				imagetype = Image::HDR;
				filetype = POV_File_Image_HDR;
				break;
			case kPOVList_FileType_System:
				imagetype = Image::SYS;
				filetype = POV_File_Image_System;
				break;
			default:
				throw POV_EXCEPTION_STRING("Invalid file type for output");
		}

		int gammaType = ropts.TryGetInt(kPOVAttrib_FileGammaType, DEFAULT_FILE_GAMMA_TYPE);
		float gamma = ropts.TryGetFloat(kPOVAttrib_FileGamma, DEFAULT_FILE_GAMMA);
		wopts.encodingGamma = GetGammaCurve(gammaType, gamma);
		// NB: RenderFrontend<...>::CreateView should have dealt with kPOVAttrib_LegacyGammaMode already and updated kPOVAttrib_WorkingGammaType and kPOVAttrib_WorkingGamma to fit.
		gammaType = ropts.TryGetInt(kPOVAttrib_WorkingGammaType, DEFAULT_WORKING_GAMMA_TYPE);
		gamma = ropts.TryGetFloat(kPOVAttrib_WorkingGamma, DEFAULT_WORKING_GAMMA);
		wopts.workingGamma = GetGammaCurve(gammaType, gamma);

		bool dither = ropts.TryGetBool(kPOVAttrib_Dither, false);
		int ditherMethod = kPOVList_DitherMethod_None;
		if (dither)
			ditherMethod = ropts.TryGetInt(kPOVAttrib_DitherMethod, kPOVList_DitherMethod_FloydSteinberg);
		wopts.dither = GetDitherHandler(ditherMethod, image->GetWidth());

		// in theory this should always return a filename since the frontend code
		// sets it via a call to GetOutputFilename() before the render starts.
		UCS2String filename = ropts.TryGetUCS2String(kPOVAttrib_OutputFile, "");
		if(filename.empty() == true)
			filename = GetOutputFilename(ropts, frame, digits);

		boost::scoped_ptr<OStream> imagefile(NewOStream(filename.c_str(), filetype, false)); // TODO - check file permissions somehow without macro [ttrf]
		if(imagefile == NULL)
			throw POV_EXCEPTION_CODE(kCannotOpenFileErr);

		Image::Write(imagetype, imagefile.get(), image.get(), wopts);

		return filename;
	}
	else
		return UCS2String();
}

boost::shared_ptr<Image>& ImageProcessing::GetImage()
{
	return image;
}

void ImageProcessing::RGB2XYZ(const COLC *rgb, COLC *xyz)
{
	// assumes D65 white point (slightly rounded sRGB)
	xyz[X] = (0.412453 * rgb[Colour::RED]) + (0.357580 * rgb[Colour::GREEN]) + (0.180423 * rgb[Colour::BLUE]);
	xyz[Y] = (0.212671 * rgb[Colour::RED]) + (0.715160 * rgb[Colour::GREEN]) + (0.072169 * rgb[Colour::BLUE]);
	xyz[Z] = (0.019334 * rgb[Colour::RED]) + (0.119193 * rgb[Colour::GREEN]) + (0.950227 * rgb[Colour::BLUE]);
}

void ImageProcessing::XYZ2RGB(const COLC *xyz, COLC *rgb)
{
	// assumes D65 white point (slightly rounded sRGB)
	rgb[Colour::RED] =    (3.240479 * xyz[X]) + (-1.537150 * xyz[X]) + (-0.498535 * xyz[X]);
	rgb[Colour::GREEN] = (-0.969256 * xyz[Y]) +  (1.875992 * xyz[Y]) +  (0.041556 * xyz[Y]);
	rgb[Colour::BLUE] =   (0.055648 * xyz[Z]) + (-0.204043 * xyz[Z]) +  (1.057311 * xyz[Z]);
}

bool ImageProcessing::OutputIsStdout(POVMS_Object& ropts)
{
	UCS2String path(ropts.TryGetUCS2String(kPOVAttrib_OutputFile, ""));

	toStdout = path == POVMS_ASCIItoUCS2String("-") || path == POVMS_ASCIItoUCS2String("stdout");
	toStderr = path == POVMS_ASCIItoUCS2String("stderr");
	return toStdout;
}

bool ImageProcessing::OutputIsStderr(POVMS_Object& ropts)
{
	OutputIsStdout(ropts);
	return toStderr;
}

UCS2String ImageProcessing::GetOutputFilename(POVMS_Object& ropts, POVMSInt frame, int digits)
{
	Path path(ropts.TryGetUCS2String(kPOVAttrib_OutputFile, ""));
	UCS2String filename = path.GetFile();
	UCS2String ext;
	Image::ImageFileType imagetype;

	switch(ropts.TryGetInt(kPOVAttrib_OutputFileType, DEFAULT_OUTPUT_FORMAT))
	{
		case kPOVList_FileType_Targa:
		case kPOVList_FileType_CompressedTarga:
			ext = ASCIItoUCS2String(".tga");
			imagetype = Image::TGA;
			break;

		case kPOVList_FileType_PNG:
			ext = ASCIItoUCS2String(".png");
			imagetype = Image::PNG;
			break;

		case kPOVList_FileType_JPEG:
			ext = ASCIItoUCS2String(".jpg");
			imagetype = Image::JPEG;
			break;

		case kPOVList_FileType_PPM:
			ext = ASCIItoUCS2String(".ppm"); // TODO FIXME - in case of greyscale output, extension should default to ".pgm"
			imagetype = Image::PPM;
			break;

		case kPOVList_FileType_BMP:
			ext = ASCIItoUCS2String(".bmp");
			imagetype = Image::BMP;
			break;

		case kPOVList_FileType_OpenEXR:
			ext = ASCIItoUCS2String(".exr");
			imagetype = Image::EXR;
			break;

		case kPOVList_FileType_RadianceHDR:
			ext = ASCIItoUCS2String(".hdr");
			imagetype = Image::HDR;
			break;

#ifdef SYS_TO_STANDARD
		case kPOVList_FileType_System:
			ext = ASCIItoUCS2String(POV_SYS_FILE_EXTENSION);
			imagetype = Image::SYS_TO_STANDARD;
			break;
#endif

		default:
			throw POV_EXCEPTION_STRING("Invalid file type for output");
	}

	if (OutputIsStdout(ropts) || OutputIsStderr())
	{
		switch (imagetype)
		{
			case Image::HDR:
			case Image::PNG:
			case Image::TGA:
			case Image::PPM:
			case Image::BMP:
				break;

			default:
				throw POV_EXCEPTION_STRING("Output to STDOUT/STDERR not supported for selected file format");
		}
		return POVMS_ASCIItoUCS2String(OutputIsStdout() ? "stdout" : "stderr");
	}

	// we disallow an output filename that consists purely of the default extension
	// (e.g. Output_File_Name=".png").
	if((filename == ext) || (filename.empty() == true))
	{
		// get the input file name and merge the existing path if need be.
		if (path.Empty() == true)
		{
			path = ropts.TryGetUCS2String(kPOVAttrib_InputFile, "object.pov");
			filename = path.GetFile();
		}
		else
			filename = Path(ropts.TryGetUCS2String(kPOVAttrib_InputFile, "object.pov")).GetFile();

		// if the input file name ends with '.' or '.anything', we remove it
		UCS2String::size_type pos = filename.find_last_of('.');
		if(pos != string::npos)
			filename.erase(pos);
	}
	else if ((path.HasVolume() == false) && (path.Empty() == false))
	{
		// to get here, path must be a relative path with filename
		// if the filename ends with a '.' or with the default extension (case-sensitive),
		// we remove it.
		UCS2String::size_type pos = filename.find_last_of('.');
		if((pos != UCS2String::npos) && ((pos == filename.size() - 1) || (filename.substr(pos) == ext)))
			filename.erase(pos);
	}
	else
	{
		// if there is no path component already, get it from the input file.
		if (path.Empty() == true)
			path = ropts.TryGetUCS2String(kPOVAttrib_InputFile, "object.pov");

		// if the filename ends with a '.' or with the default extension (case-sensitive),
		// we remove it.
		UCS2String::size_type pos = filename.find_last_of('.');
		if((pos != UCS2String::npos) && ((pos == filename.size() - 1) || (filename.substr(pos) == ext)))
			filename.erase(pos);
	}

	if (digits > 0)
	{
		for(int i = 0; i < digits; i++)
			filename += '0';
		for(UCS2String::size_type i = filename.length() - 1; frame > 0; i--, frame /= 10)
			filename[i] = '0' + (frame % 10);
	}

	path.SetFile(filename + ext);

	return path();
}

}
