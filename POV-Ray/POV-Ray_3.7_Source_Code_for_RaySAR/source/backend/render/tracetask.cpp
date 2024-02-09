/*******************************************************************************
 * tracetask.cpp
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
 * $File: //depot/public/povray/3.x/source/backend/render/tracetask.cpp $
 * $Revision: #1 $
 * $Change: 6069 $
 * $DateTime: 2013/11/06 11:59:40 $
 * $Author: chrisc $
 *******************************************************************************/

#include <vector>

#include <boost/thread.hpp>

// frame.h must always be the first POV file included (pulls in platform config)
#include "backend/frame.h"
#include "backend/colour/colour.h"
#include "backend/math/vector.h"
#include "backend/math/matrices.h"
#include "backend/render/trace.h"
#include "backend/render/tracetask.h"
#include "backend/support/jitter.h"
#include "backend/texture/normal.h"
#include "backend/math/chi2.h"

// new for RaySAR simulator
#include "backend/parser/parse.h"
#include "backend/parser/reswords.h"

#ifdef PROFILE_INTERSECTIONS
#include "base/image/image.h"
#endif

// this must be the last file included
#include "base/povdebug.h"

namespace pov
{

// Global variables
// New for RaySAR simulator
	/*
DBL Intensity_Radar = 0.0;
DBL Distance_Radar = 0.0;
DBL Azimuth_Initial = 0.0;
DBL Elevation_Initial = 0.0;
DBL Reflectivity_Bounce = 0.0;
DBL Reflectivity_Bounce_Temp;
DBL Total_Depth = 0.0;
int Trace_Level = 0;
int mark_specular;
VECTOR Pixel_Initial_Radar;
VECTOR Direction_Primary_Ray;
VECTOR Direction_Current_Ray;
*/
raysar::RefContribWriter contribWriter("Contributions.txt");
// New for RaySAR simulator

#ifdef PROFILE_INTERSECTIONS
	bool gDoneBSP;
	bool gDoneBVH;
	POV_ULONG gMinVal = ULLONG_MAX ;
	POV_ULONG gMaxVal = 0;
	POV_ULONG gIntersectionTime;
	vector<vector<POV_ULONG> > gBSPIntersectionTimes;
	vector<vector<POV_ULONG> > gBVHIntersectionTimes;
	vector <vector<POV_ULONG> > *gIntersectionTimes;
#endif

class SmartBlock
{
	public:
		SmartBlock(int ox, int oy, int bw, int bh);

		bool GetFlag(int x, int y) const;
		void SetFlag(int x, int y, bool f);

		Colour& operator()(int x, int y);
		const Colour& operator()(int x, int y) const;

		vector<Colour>& GetPixels();
	private:
		vector<Colour> framepixels;
		vector<Colour> pixels;
		vector<bool> frameflags;
		vector<bool> flags;
		int offsetx;
		int offsety;
		int blockwidth;
		int blockheight;

		int GetOffset(int x, int y) const;
};

SmartBlock::SmartBlock(int ox, int oy, int bw, int bh) :
	offsetx(ox),
	offsety(oy),
	blockwidth(bw),
	blockheight(bh)
{
	framepixels.resize((blockwidth * 2) + (blockheight * 2) + 4);
	pixels.resize(blockwidth * blockheight);
	frameflags.resize((blockwidth * 2) + (blockheight * 2) + 4);
	flags.resize(blockwidth * blockheight);
}

bool SmartBlock::GetFlag(int x, int y) const
{
	int offset = GetOffset(x, y);

	if(offset < 0)
		return frameflags[-1 - offset];
	else
		return flags[offset];
}

void SmartBlock::SetFlag(int x, int y, bool f)
{
	int offset = GetOffset(x, y);

	if(offset < 0)
		frameflags[-1 - offset] = f;
	else
		flags[offset] = f;
}

Colour& SmartBlock::operator()(int x, int y)
{
	int offset = GetOffset(x, y);

	if(offset < 0)
		return framepixels[-1 - offset];
	else
		return pixels[offset];
}

const Colour& SmartBlock::operator()(int x, int y) const
{
	int offset = GetOffset(x, y);

	if(offset < 0)
		return framepixels[-1 - offset];
	else
		return pixels[offset];
}

vector<Colour>& SmartBlock::GetPixels()
{
	return pixels;
}

int SmartBlock::GetOffset(int x, int y) const
{
	x -= offsetx;
	y -= offsety;

	if(x < 0)
		x = -1;
	else if(x >= blockwidth)
		x = blockwidth;

	if(y < 0)
		y = -1;
	else if(y >= blockheight)
		y = blockheight;

	if((x < 0) && (y < 0))
		return -1;
	else if((x >= blockwidth) && (y < 0))
		return -2;
	else if((x < 0) && (y >= blockheight))
		return -3;
	else if((x >= blockwidth) && (y >= blockheight))
		return -4;
	else if(x < 0)
		return -(4 + y);
	else if(y < 0)
		return -(4 + x + blockheight);
	else if(x >= blockwidth)
		return -(4 + y + blockheight + blockwidth);
	else if(y >= blockheight)
		return -(4 + x + blockheight + blockwidth + blockheight);
	else
		return (x + (y * blockwidth));
}

TraceTask::CooperateFunction::CooperateFunction(Task& t) :
	task(t)
{

}

void TraceTask::CooperateFunction::operator()()
{
	task.Cooperate();
}

TraceTask::SubdivisionBuffer::SubdivisionBuffer(size_t s) :
	colors(s * s),
	sampled(s * s),
	size(s)
{
	Clear();
}

void TraceTask::SubdivisionBuffer::SetSample(size_t x, size_t y, const Colour& col)
{
	colors[x + (y * size)] = col;
	sampled[x + (y * size)] = true;
}

bool TraceTask::SubdivisionBuffer::Sampled(size_t x, size_t y)
{
	return sampled[x + (y * size)];
}

Colour& TraceTask::SubdivisionBuffer::operator()(size_t x, size_t y)
{
	return  colors[x + (y * size)];
}

void TraceTask::SubdivisionBuffer::Clear()
{
	for(vector<bool>::iterator i(sampled.begin()); i != sampled.end(); i++)
		*i = false;
}

TraceTask::TraceTask(ViewData *vd, unsigned int tm, DBL js, DBL aat, unsigned int aad, GammaCurvePtr& aag, unsigned int ps, bool psc, bool final, bool hr) :
	RenderTask(vd),
	trace(vd, GetViewDataPtr(), vd->GetSceneData()->parsedMaxTraceLevel, vd->GetSceneData()->parsedAdcBailout,
	      vd->GetQualityFeatureFlags(), cooperate, media, radiosity),
	cooperate(*this),
	tracingMethod(tm),
	jitterScale(js),
	aaThreshold(aat),
	aaDepth(aad),
	aaGamma(aag),
	previewSize(ps),
	previewSkipCorner(psc),
	finalTrace(final),
	highReproducibility(hr),
	media(GetViewDataPtr(), &trace, &photonGatherer),
	radiosity(vd->GetSceneData(), GetViewDataPtr(),
	          vd->GetSceneData()->radiositySettings, vd->GetRadiosityCache(), cooperate, final, Vector3d(vd->GetCamera().Location)),
	photonGatherer(&vd->GetSceneData()->mediaPhotonMap, vd->GetSceneData()->photonSettings)	
{
#ifdef PROFILE_INTERSECTIONS
	Rectangle ra = vd->GetRenderArea();
	if (vd->GetSceneData()->boundingMethod == 2)
	{
		gBSPIntersectionTimes.clear();
		gBSPIntersectionTimes.resize(ra.bottom + 1);
		for (int i = 0; i < ra.bottom + 1; i++)
			gBSPIntersectionTimes[i].resize(ra.right + 1);
		gIntersectionTimes = &gBSPIntersectionTimes;
		gDoneBSP = true;
	}
	else
	{
		gBVHIntersectionTimes.clear();
		gBVHIntersectionTimes.resize(ra.bottom + 1);
		for (int i = 0; i < ra.bottom + 1; i++)
			gBVHIntersectionTimes[i].resize(ra.right + 1);
		gIntersectionTimes = &gBVHIntersectionTimes;
		gDoneBVH = true;
	}
#endif
	// TODO: this could be initialised someplace more suitable
	GetViewDataPtr()->qualityFlags = vd->GetQualityFeatureFlags();
}

TraceTask::~TraceTask()
{
}

void TraceTask::Run()
{

#ifdef RTR_HACK
	bool forever = GetViewData()->GetRealTimeRaytracing();
	do
	{
#endif
		switch(tracingMethod)
		{
			case 0:
				if(previewSize > 0)
					SimpleSamplingM0P();
				else
					SimpleSamplingM0();	

				break;
			case 1:
				NonAdaptiveSupersamplingM1();
				break;
			case 2:
				AdaptiveSupersamplingM2();
				break;
		}

#ifdef RTR_HACK
		if(forever)
		{
			const Camera *camera = GetViewData()->GetRTRData()->CompletedFrame();
			Cooperate();
			if(camera != NULL)
				trace.SetupCamera(*camera);
		}
	} while(forever);
#endif

	GetViewData()->SetHighestTraceLevel(trace.GetHighestTraceLevel());

}

void TraceTask::Stopped()
{
	// nothing to do for now [trf]
}

void TraceTask::Finish()
{
	GetViewDataPtr()->timeType = SceneThreadData::kRenderTime;
	GetViewDataPtr()->realTime = ConsumedRealTime();
	GetViewDataPtr()->cpuTime = ConsumedCPUTime();

	
#ifdef PROFILE_INTERSECTIONS
	if (gDoneBSP && gDoneBVH)
	{
		int width = gBSPIntersectionTimes[0].size();
		int height = gBSPIntersectionTimes.size();
		if (width == gBVHIntersectionTimes[0].size() && height == gBVHIntersectionTimes.size())
		{
			SNGL scale = 1.0 / (gMaxVal - gMinVal);
			Image::WriteOptions opts;
			opts.bpcc = 16;
			Image *img = Image::Create(width, height, Image::Gray_Int16, false);
			for (int y = 0 ; y < height ; y++)
				for (int x = 0 ; x < width ; x++)
					img->SetGrayValue(x, y, (gBSPIntersectionTimes[y][x] - gMinVal) * scale);
			OStream *imagefile(NewOStream("bspprofile.png", 0, false));
			Image::Write(Image::PNG, imagefile, img, opts);
			delete imagefile;
			delete img;

			img = Image::Create(width, height, Image::Gray_Int16, false);
			imagefile = NewOStream("bvhprofile.png", 0, false);
			for (int y = 0 ; y < height ; y++)
				for (int x = 0 ; x < width ; x++)
					img->SetGrayValue(x, y, (gBVHIntersectionTimes[y][x] - gMinVal) * scale);
			Image::Write(Image::PNG, imagefile, img, opts);
			delete imagefile;
			delete img;

			img = Image::Create(width, height, Image::Gray_Int16, false);
			imagefile = NewOStream("summedprofile.png", 0, false);
			for (int y = 0 ; y < height ; y++)
				for (int x = 0 ; x < width ; x++)
					img->SetGrayValue(x, y, 0.5f + ((((gBSPIntersectionTimes[y][x] - gMinVal) - (gBVHIntersectionTimes[y][x] - gMinVal)) * scale) / 2));
			Image::Write(Image::PNG, imagefile, img, opts);
			delete imagefile;
			delete img;

			img = Image::Create(width, height, Image::RGBFT_Float, false);
			imagefile = NewOStream("rgbprofile.png", 0, false);
			for (int y = 0 ; y < height ; y++)
			{
				for (int x = 0 ; x < width ; x++)
				{
					Colour col;
					float bspval = (gBSPIntersectionTimes[y][x] - gMinVal) * scale ;
					float bvhval = (gBVHIntersectionTimes[y][x] - gMinVal) * scale ;
					float diff = bspval - bvhval ;
					if (diff > 0.0)
						col.blue() += diff ;
					else
						col.red() -= diff ;
					img->SetRGBFTValue(x, y, col);
				}
			}
			Image::Write(Image::PNG, imagefile, img, opts);
			delete imagefile;
			delete img;
		}
		gDoneBSP = gDoneBVH = false;
		gMinVal = ULLONG_MAX;
		gMaxVal = 0;
	}
#endif
}

void TraceTask::SimpleSamplingM0()
{
	POVRect rect;
	vector<Colour> pixels;
	unsigned int serial;

	while(GetViewData()->GetNextRectangle(rect, serial) == true)
	{
		radiosity.BeforeTile(highReproducibility? serial : 0);

		pixels.clear();
		pixels.reserve(rect.GetArea());

		for(DBL y = DBL(rect.top); y <= DBL(rect.bottom); y++)
		{
			for(DBL x = DBL(rect.left); x <= DBL(rect.right); x++)
			{
#ifdef PROFILE_INTERSECTIONS
				POV_LONG it = ULLONG_MAX;
				for (int i = 0 ; i < 3 ; i++)
				{
					Colour c;
					gIntersectionTime = 0;
					trace(x, y, GetViewData()->GetWidth(), GetViewData()->GetHeight(), c);
					if (gIntersectionTime < it)
						it = gIntersectionTime;
				}
				(*gIntersectionTimes)[(int) y] [(int) x] = it;
				if (it < gMinVal)
					gMinVal = it;
				if (it > gMaxVal)
					gMaxVal = it;
#endif
				Colour col;

				// New for RaySAR Simulator
				/*
				Reflectivity_Bounce = 1;  // Reset
				Trace_Level = 1; // Reset
				Total_Depth = 0.0; // Reset
				*/
				// New for RaySAR Simulator
				trace.ResetDepth();

				trace(x, y, GetViewData()->GetWidth(), GetViewData()->GetHeight(), col);
				GetViewDataPtr()->Stats()[Number_Of_Pixels]++;

				pixels.push_back(col);

				Cooperate();
			}
		}

		radiosity.AfterTile();

		GetViewDataPtr()->AfterTile();
		GetViewData()->CompletedRectangle(rect, serial, pixels, 1, finalTrace);
		
		Cooperate();

	}
}

void TraceTask::SimpleSamplingM0P()
{
	DBL stepsize(previewSize);
	POVRect rect;
	vector<Vector2d> pixelpositions;
	vector<Colour> pixelcolors;
	unsigned int serial;

	while(GetViewData()->GetNextRectangle(rect, serial) == true)
	{
		radiosity.BeforeTile(highReproducibility? serial : 0);

		unsigned int px = (rect.GetWidth() + previewSize - 1) / previewSize;
		unsigned int py = (rect.GetHeight() + previewSize - 1) / previewSize;

		pixelpositions.clear();
		pixelpositions.reserve(px * py);
		pixelcolors.clear();
		pixelcolors.reserve(px * py);

		for(DBL y = DBL(rect.top); y <= DBL(rect.bottom); y += stepsize)
		{
			for(DBL x = DBL(rect.left); x <= DBL(rect.right); x += stepsize)
			{
				if((previewSkipCorner == true) && (fmod(x, stepsize * 2.0) < EPSILON) && (fmod(y, stepsize * 2.0) < EPSILON))
					continue;

#ifdef PROFILE_INTERSECTIONS
				POV_LONG it = ULLONG_MAX;
				for (int i = 0 ; i < 3 ; i++)
				{
					Colour c;
					gIntersectionTime = 0;
					trace(x, y, GetViewData()->GetWidth(), GetViewData()->GetHeight(), c);
					if (gIntersectionTime < it)
						it = gIntersectionTime;
				}
				(*gIntersectionTimes)[(int) y] [(int) x] = it;
				if (it < gMinVal)
					gMinVal = it;
				if (it > gMaxVal)
					gMaxVal = it;
#endif
				Colour col;

				// New for RaySAR Simulator
				// Reflectivity_Bounce = 1;  // Reset

				trace(x, y, GetViewData()->GetWidth(), GetViewData()->GetHeight(), col);
				GetViewDataPtr()->Stats()[Number_Of_Pixels]++;

				pixelpositions.push_back(Vector2d(x, y));
				pixelcolors.push_back(col);

				Cooperate();
			}
		}

		radiosity.AfterTile();

		GetViewDataPtr()->AfterTile();
		if(pixelpositions.size() > 0)
			GetViewData()->CompletedRectangle(rect, serial, pixelpositions, pixelcolors, previewSize, finalTrace);

		Cooperate();
	}
}

void TraceTask::NonAdaptiveSupersamplingM1()
{
	POVRect rect;
	unsigned int serial;

	jitterScale = jitterScale / DBL(aaDepth);

	while(GetViewData()->GetNextRectangle(rect, serial) == true)
	{
		radiosity.BeforeTile(highReproducibility? serial : 0);

		SmartBlock pixels(rect.left, rect.top, rect.GetWidth(), rect.GetHeight());

		// sample line above current block
		for(int x = rect.left; x <= rect.right; x++)
		{
			trace(DBL(x), DBL(rect.top) - 1.0, GetViewData()->GetWidth(), GetViewData()->GetHeight(), pixels(x, rect.top - 1));
			GetViewDataPtr()->Stats()[Number_Of_Pixels]++;

			// Cannot supersample this pixel, so just claim it was already supersampled! [trf]
			// [CJC] see comment for leftmost pixels below; similar situation applies here
			pixels.SetFlag(x, rect.top - 1, true);

			Cooperate();
		}

		for(int y = rect.top; y <= rect.bottom; y++)
		{
			trace(DBL(rect.left) - 1.0, y, GetViewData()->GetWidth(), GetViewData()->GetHeight(), pixels(rect.left - 1, y)); // sample pixel left of current line in block
			GetViewDataPtr()->Stats()[Number_Of_Pixels]++;

			// Cannot supersample this pixel, so just claim it was already supersampled! [trf]

			// [CJC] NB this could in some circumstances cause an artifact at a block boundary,
			// if the leftmost pixel on this blockline ends up not being supersampled because the
			// difference between it and the rightmost pixel on the same line in the last block
			// was insufficient to trigger the supersample right now, *BUT* if the rightmost pixel
			// *had* been supersampled, AND the difference was then enough to trigger the call
			// to supersample the current pixel, AND when the block on the left was/is rendered,
			// the abovementioned rightmost pixel *does* get supersampled due to the logic applied
			// when the code rendered *that* block ... [a long set of preconditions but possible].

			// there's no easy solution to this because if we *do* supersample right now, the
			// reverse situation could apply if the rightmost pixel in the last block ends up
			// not being supersampled ...
			pixels.SetFlag(rect.left - 1, y, true);

			Cooperate();

			for(int x = rect.left; x <= rect.right; x++)
			{

				// trace current pixel
				trace(DBL(x), DBL(y), GetViewData()->GetWidth(), GetViewData()->GetHeight(), pixels(x, y));
				GetViewDataPtr()->Stats()[Number_Of_Pixels]++;

				Cooperate();

				bool sampleleft = (pixels.GetFlag(x - 1, y) == false);
				bool sampletop = (pixels.GetFlag(x, y - 1) == false);
				bool samplecurrent = true;

				// perform antialiasing
				NonAdaptiveSupersamplingForOnePixel(DBL(x), DBL(y), pixels(x - 1, y), pixels(x, y - 1), pixels(x, y), sampleleft, sampletop, samplecurrent);

				// if these pixels have been supersampled, set their supersampling flag
				if(sampleleft == true)
					pixels.SetFlag(x - 1, y, true);
				if(sampletop == true)
					pixels.SetFlag(x, y - 1, true);
				if(samplecurrent == true)
					pixels.SetFlag(x, y, true);
			}
		}

		radiosity.AfterTile();

		GetViewDataPtr()->AfterTile();
		GetViewData()->CompletedRectangle(rect, serial, pixels.GetPixels(), 1, finalTrace);

		Cooperate();
	}
}

void TraceTask::AdaptiveSupersamplingM2()
{
	POVRect rect;
	unsigned int serial;
	size_t subsize = (1 << aaDepth);
	SubdivisionBuffer buffer(subsize + 1);

	jitterScale = jitterScale / DBL((1 << aaDepth) + 1);

	while(GetViewData()->GetNextRectangle(rect, serial) == true)
	{
		radiosity.BeforeTile(highReproducibility? serial : 0);

		SmartBlock pixels(rect.left, rect.top, rect.GetWidth(), rect.GetHeight());

		for(int y = rect.top; y <= rect.bottom + 1; y++)
		{
			for(int x = rect.left; x <= rect.right + 1; x++)
			{
				// trace upper-left corners of all pixels
				trace(DBL(x) - 0.5, DBL(y) - 0.5, GetViewData()->GetWidth(), GetViewData()->GetHeight(), pixels(x, y));
				GetViewDataPtr()->Stats()[Number_Of_Pixels]++;

				Cooperate();
			}
		}

		// note that the bottom and/or right corner are the
		// upper-left corner of the bottom and/or right pixels
		for(int y = rect.top; y <= rect.bottom; y++)
		{
			for(int x = rect.left; x <= rect.right; x++)
			{
				buffer.Clear();

				buffer.SetSample(0, 0, pixels(x, y));
				buffer.SetSample(0, subsize, pixels(x, y + 1));
				buffer.SetSample(subsize, 0, pixels(x + 1, y));
				buffer.SetSample(subsize, subsize, pixels(x + 1, y + 1));

				SubdivideOnePixel(DBL(x), DBL(y), 0.5, 0, 0, subsize, buffer, pixels(x, y), aaDepth - 1);

				Cooperate();
			}
		}

		radiosity.AfterTile();

		GetViewDataPtr()->AfterTile();
		GetViewData()->CompletedRectangle(rect, serial, pixels.GetPixels(), 1, finalTrace);

		Cooperate();
	}
}

void TraceTask::NonAdaptiveSupersamplingForOnePixel(DBL x, DBL y, Colour& leftcol, Colour& topcol, Colour& curcol, bool& sampleleft, bool& sampletop, bool& samplecurrent)
{
	Colour gcLeft = GammaCurve::Encode(aaGamma, leftcol);
	Colour gcTop  = GammaCurve::Encode(aaGamma, topcol);
	Colour gcCur  = GammaCurve::Encode(aaGamma, curcol);

	bool leftdiff = (Colour_Distance_RGBT(gcLeft, gcCur) >= aaThreshold);
	bool topdiff  = (Colour_Distance_RGBT(gcTop,  gcCur) >= aaThreshold);

	sampleleft = sampleleft && leftdiff;
	sampletop = sampletop && topdiff;
	samplecurrent = ((leftdiff == true) || (topdiff == true));

	if(sampleleft == true)
		SupersampleOnePixel(x - 1.0, y, leftcol);

	if(sampletop == true)
		SupersampleOnePixel(x, y - 1.0, topcol);

	if(samplecurrent == true)
		SupersampleOnePixel(x, y, curcol);
}

void TraceTask::SupersampleOnePixel(DBL x, DBL y, Colour& col)
{
	DBL step(1.0 / DBL(aaDepth));
	DBL range(0.5 - (step * 0.5));
	DBL rx, ry;
	Colour tempcol;

	GetViewDataPtr()->Stats()[Number_Of_Pixels_Supersampled]++;

	for(DBL yy = -range; yy <= (range + EPSILON); yy += step)
	{
		for(DBL xx = -range; xx <= (range + EPSILON); xx += step)
		{
			if (jitterScale > 0.0)
			{
				Jitter2d(x + xx, y + yy, rx, ry);
				trace(x + xx + (rx * jitterScale), y + yy + (ry * jitterScale), GetViewData()->GetWidth(), GetViewData()->GetHeight(), tempcol);
			}
			else
				trace(x + xx, y + yy, GetViewData()->GetWidth(), GetViewData()->GetHeight(), tempcol);

			col += tempcol;
			GetViewDataPtr()->Stats()[Number_Of_Samples]++;

			Cooperate();
		}
	}

	col /= (aaDepth * aaDepth + 1);
}

void TraceTask::SubdivideOnePixel(DBL x, DBL y, DBL d, size_t bx, size_t by, size_t bstep, SubdivisionBuffer& buffer, Colour& result, int level)
{
	Colour& cx0y0 = buffer(bx, by);
	Colour& cx0y2 = buffer(bx, by + bstep);
	Colour& cx2y0 = buffer(bx + bstep, by);
	Colour& cx2y2 = buffer(bx + bstep, by + bstep);
	size_t bstephalf = bstep / 2;

	// o = no operation, + = input, * = output

	// Input:
	// +o+
	// ooo
	// +o+

	Colour cx0y0g = GammaCurve::Encode(aaGamma, cx0y0);
	Colour cx0y2g = GammaCurve::Encode(aaGamma, cx0y2);
	Colour cx2y0g = GammaCurve::Encode(aaGamma, cx2y0);
	Colour cx2y2g = GammaCurve::Encode(aaGamma, cx2y2);

	if((level > 0) &&
	   ((Colour_Distance_RGBT(cx0y0g, cx0y2g) >= aaThreshold) ||
	    (Colour_Distance_RGBT(cx0y0g, cx2y0g) >= aaThreshold) ||
	    (Colour_Distance_RGBT(cx0y0g, cx2y2g) >= aaThreshold) ||
	    (Colour_Distance_RGBT(cx0y2g, cx2y0g) >= aaThreshold) ||
	    (Colour_Distance_RGBT(cx0y2g, cx2y2g) >= aaThreshold) ||
	    (Colour_Distance_RGBT(cx2y0g, cx2y2g) >= aaThreshold)))
	{
		Colour rcx0y0;
		Colour rcx0y1;
		Colour rcx1y0;
		Colour rcx1y1;
		Colour col;
		DBL rxcx0y1, rycx0y1;
		DBL rxcx1y0, rycx1y0;
		DBL rxcx2y1, rycx2y1;
		DBL rxcx1y2, rycx1y2;
		DBL rxcx1y1, rycx1y1;
		DBL d2 = d * 0.5;

		// Trace:
		// ooo
		// *oo
		// ooo
		if(buffer.Sampled(bx, by + bstephalf) == false)
		{
			if (jitterScale > 0.0)
			{
				Jitter2d(x - d, y, rxcx0y1, rycx0y1);
				trace(x - d + (rxcx0y1 * jitterScale), y + (rycx0y1 * jitterScale), GetViewData()->GetWidth(), GetViewData()->GetHeight(), col);
			}
			else
				trace(x - d, y, GetViewData()->GetWidth(), GetViewData()->GetHeight(), col);

			buffer.SetSample(bx, by + bstephalf, col);

			GetViewDataPtr()->Stats()[Number_Of_Samples]++;
			Cooperate();
		}

		// Trace:
		// o*o
		// ooo
		// ooo
		if(buffer.Sampled(bx + bstephalf, by) == false)
		{
			if (jitterScale > 0.0)
			{
				Jitter2d(x, y - d, rxcx1y0, rycx1y0);
				trace(x + (rxcx1y0 * jitterScale), y - d + (rycx1y0 * jitterScale), GetViewData()->GetWidth(), GetViewData()->GetHeight(), col);
			}
			else
				trace(x, y - d, GetViewData()->GetWidth(), GetViewData()->GetHeight(), col);

			buffer.SetSample(bx + bstephalf, by, col);

			GetViewDataPtr()->Stats()[Number_Of_Samples]++;
			Cooperate();
		}

		// Trace:
		// ooo
		// oo*
		// ooo
		if(buffer.Sampled(bx + bstep, by + bstephalf) == false)
		{
			if (jitterScale > 0.0)
			{
				Jitter2d(x + d, y, rxcx2y1, rycx2y1);
				trace(x + d + (rxcx2y1 * jitterScale), y + (rycx2y1 * jitterScale), GetViewData()->GetWidth(), GetViewData()->GetHeight(), col);
			}
			else
				trace(x + d, y, GetViewData()->GetWidth(), GetViewData()->GetHeight(), col);

			buffer.SetSample(bx + bstep, by + bstephalf, col);

			GetViewDataPtr()->Stats()[Number_Of_Samples]++;
			Cooperate();
		}

		// Trace:
		// ooo
		// ooo
		// o*o
		if(buffer.Sampled(bx + bstephalf, by + bstep) == false)
		{
			if (jitterScale > 0.0)
			{
				Jitter2d(x, y + d, rxcx1y2, rycx1y2);
				trace(x + (rxcx1y2 * jitterScale), y + d + (rycx1y2 * jitterScale), GetViewData()->GetWidth(), GetViewData()->GetHeight(), col);
			}
			else
				trace(x, y + d, GetViewData()->GetWidth(), GetViewData()->GetHeight(), col);

			buffer.SetSample(bx + bstephalf, by + bstep, col);

			GetViewDataPtr()->Stats()[Number_Of_Samples]++;
			Cooperate();
		}

		// Trace:
		// ooo
		// o*o
		// ooo
		if(buffer.Sampled(bx + bstephalf, by + bstephalf) == false)
		{
			if (jitterScale > 0.0)
			{
				Jitter2d(x, y, rxcx1y1, rycx1y1);
				trace(x + (rxcx1y1 * jitterScale), y + (rycx1y1 * jitterScale), GetViewData()->GetWidth(), GetViewData()->GetHeight(), col);
			}
			else
				trace(x, y, GetViewData()->GetWidth(), GetViewData()->GetHeight(), col);

			buffer.SetSample(bx + bstephalf, by + bstephalf, col);

			GetViewDataPtr()->Stats()[Number_Of_Samples]++;
			Cooperate();
		}

		// Subdivide Input:
		// ++o
		// ++o
		// ooo
		// Subdivide Output:
		// *o
		// oo
		SubdivideOnePixel(x - d2, y - d2, d2, bx, by, bstephalf, buffer, rcx0y0, level - 1);

		// Subdivide Input:
		// ooo
		// ++o
		// ++o
		// Subdivide Output:
		// oo
		// *o
		SubdivideOnePixel(x - d2, y + d2, d2, bx, by + bstephalf, bstephalf, buffer, rcx0y1, level - 1);

		// Subdivide Input:
		// o++
		// o++
		// ooo
		// Subdivide Output:
		// o*
		// oo
		SubdivideOnePixel(x + d2, y - d2, d2, bx + bstephalf, by, bstephalf, buffer, rcx1y0, level - 1);

		// Subdivide Input:
		// ooo
		// o++
		// o++
		// Subdivide Output:
		// oo
		// o*
		SubdivideOnePixel(x + d2, y + d2, d2, bx + bstephalf, by + bstephalf, bstephalf, buffer, rcx1y1, level - 1);

		result = (rcx0y0 + rcx0y1 + rcx1y0 + rcx1y1) / 4.0;
	}
	else
	{
		result = (cx0y0 + cx0y2 + cx2y0 + cx2y2) / 4.0;
	}
}

}
