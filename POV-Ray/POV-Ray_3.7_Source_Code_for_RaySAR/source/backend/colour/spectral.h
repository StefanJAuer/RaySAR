/*******************************************************************************
 * spectral.h
 *
 * This module contains all defines, typedefs, and prototypes for spectral.cpp.
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
 * $File: //depot/public/povray/3.x/source/backend/colour/spectral.h $
 * $Revision: #1 $
 * $Change: 6069 $
 * $DateTime: 2013/11/06 11:59:40 $
 * $Author: chrisc $
 *******************************************************************************/

#ifndef SPECTRAL_H
#define SPECTRAL_H

namespace pov
{

using namespace pov_base;

#define SPECTRAL_VIOLET         380.0   // extreme violet
#define SPECTRAL_RED            730.0   // extreme red
#define SPECTRAL_BANDWIDTH      (SPECTRAL_RED-SPECTRAL_VIOLET)
#define SPECTRAL_CENTER         ((SPECTRAL_VIOLET + SPECTRAL_RED)/2)    // TODO - maybe should define this as yellow

/// Class representing a spectral band.
class SpectralBand
{
	public:
		/// Default Constructor.
		SpectralBand():
			wavelength  ( (SPECTRAL_VIOLET + SPECTRAL_RED)/2 ), // deliberately NOT using SPECTRAL_CENTER
			bandwidth   ( SPECTRAL_BANDWIDTH )
		{}

		/// Construct by abstract band index (ranging from 0 to N-1).
		SpectralBand(unsigned int band, unsigned int bands):
			wavelength  ( SPECTRAL_VIOLET + SPECTRAL_BANDWIDTH*( ((float)band+0.5)/(float)bands ) ),
			bandwidth   ( SPECTRAL_BANDWIDTH/(float)bands )
		{}

		/// Construct by spectral band and sub-band index (ranging from 1 to N).
		SpectralBand(const SpectralBand& super, unsigned int subBand, unsigned int subBands):
			wavelength  ( super.wavelength + super.bandwidth*( ((float)subBand+0.5)/(float)subBands - 0.5 ) ),
			bandwidth   ( super.bandwidth/(float)subBands )
		{}

		/// Construct by physical parameters.
		SpectralBand(float wl, float bw):
			wavelength(wl),
			bandwidth(bw)
		{}

		/// Get peak wavelength.
		float GetWavelength() const
		{
			return wavelength;
		}

		/// Get corrected IOR
		double GetDispersionIOR(double nominalIOR, double nominalDispersion) const
		{
			return nominalIOR * pow(nominalDispersion, -(wavelength-SPECTRAL_CENTER)/SPECTRAL_BANDWIDTH );
		}

		/// Get hue
		RGBColour GetHue() const
		{
			return (GetHueIntegral(wavelength+bandwidth/2) - GetHueIntegral(wavelength-bandwidth/2)) * (SPECTRAL_BANDWIDTH/bandwidth);
		}

	private:
		/// Peak wavelength.
		float wavelength;
		/// Nominal bandwidth.
		float bandwidth;

		static RGBColour GetHueIntegral(double wavelength);
};

}

#endif // SPECTRAL_H
