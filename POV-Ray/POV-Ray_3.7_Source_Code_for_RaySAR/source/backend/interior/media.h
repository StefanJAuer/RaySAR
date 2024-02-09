/*******************************************************************************
 * media.h
 *
 * This module contains all defines, typedefs, and prototypes for MEDIA.CPP.
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
 * $File: //depot/public/povray/3.x/source/backend/interior/media.h $
 * $Revision: #1 $
 * $Change: 6069 $
 * $DateTime: 2013/11/06 11:59:40 $
 * $Author: chrisc $
 *******************************************************************************/

#ifndef MEDIA_H
#define MEDIA_H

#include "backend/render/trace.h"

namespace pov
{

// Scattering types.
enum
{
	ISOTROPIC_SCATTERING            = 1,
	MIE_HAZY_SCATTERING             = 2,
	MIE_MURKY_SCATTERING            = 3,
	RAYLEIGH_SCATTERING             = 4,
	HENYEY_GREENSTEIN_SCATTERING    = 5,
	SCATTERING_TYPES                = 5
};

void Transform_Density(PIGMENT *Density, const TRANSFORM *Trans);

class MediaFunction : public Trace::MediaFunctor
{
	public:
		MediaFunction(TraceThreadData *td, Trace *t, PhotonGatherer *pg);

		virtual void ComputeMedia(vector<Media>& mediasource, const Ray& ray, Intersection& isect, Colour& colour, Trace::TraceTicket& ticket);
		virtual void ComputeMedia(const RayInteriorVector& mediasource, const Ray& ray, Intersection& isect, Colour& colour, Trace::TraceTicket& ticket);
		virtual void ComputeMedia(MediaVector& medias, const Ray& ray, Intersection& isect, Colour& colour, Trace::TraceTicket& ticket);
	protected:
		/// pseudo-random number sequence
		RandomDoubleSequence randomNumbers;
		/// pseudo-random number generator based on random number sequence
		RandomDoubleSequence::Generator randomNumberGenerator;
		/// thread data
		TraceThreadData *threadData;
		/// tracing functions
		Trace *trace;
		/// photon gather functions
		PhotonGatherer *photonGatherer;

		void ComputeMediaRegularSampling(MediaVector& medias, LightSourceEntryVector& lights, MediaIntervalVector& mediaintervals,
		                                 const Ray& ray, const Media *IMedia, int minsamples, bool ignore_photons, bool use_scattering,
		                                 bool all_constant_and_light_ray, Trace::TraceTicket& ticket);
		void ComputeMediaAdaptiveSampling(MediaVector& medias, LightSourceEntryVector& lights, MediaIntervalVector& mediaintervals,
		                                  const Ray& ray, const Media *IMedia, DBL aa_threshold, int minsamples, bool ignore_photons, bool use_scattering, Trace::TraceTicket& ticket);
		void ComputeMediaColour(MediaIntervalVector& mediaintervals, Colour& colour);
		void ComputeMediaSampleInterval(LitIntervalVector& litintervals, MediaIntervalVector& mediaintervals, const Media *media);
		void ComputeMediaLightInterval(LightSourceEntryVector& lights, LitIntervalVector& litintervals, const Ray& ray, const Intersection& isect);
		void ComputeOneMediaLightInterval(LightSource *light, LightSourceEntryVector&lights, const Ray& ray, const Intersection& isect);
		bool ComputeSpotLightInterval(const Ray &ray, const LightSource *Light, DBL *d1, DBL *d2);
		bool ComputeCylinderLightInterval(const Ray &ray, const LightSource *Light, DBL *d1, DBL *d2);
		void ComputeOneMediaSample(MediaVector& medias, LightSourceEntryVector& lights, MediaInterval& mediainterval, const Ray &ray, DBL d0, RGBColour& SampCol,
		                           RGBColour& SampOptDepth, int sample_method, bool ignore_photons, bool use_scattering, bool photonPass, Trace::TraceTicket& ticket);
		void ComputeOneMediaSampleRecursive(MediaVector& medias, LightSourceEntryVector& lights, MediaInterval& mediainterval, const Ray& ray,
		                                    DBL d1, DBL d3, RGBColour& Result, const RGBColour& C1, const RGBColour& C3, RGBColour& ODResult, const RGBColour& od1, const RGBColour& od3,
		                                    int depth, DBL Jitter, DBL aa_threshold, bool ignore_photons, bool use_scattering, bool photonPass, Trace::TraceTicket& ticket);
		void ComputeMediaPhotons(MediaVector& medias, RGBColour& Te, const RGBColour& Sc, const Ray& ray, const VECTOR H);
		void ComputeMediaScatteringAttenuation(MediaVector& medias, RGBColour& OutputColor, const RGBColour& Sc, const RGBColour& Light_Colour, const Ray &ray, const Ray &Light_Ray);
};

}

#endif
