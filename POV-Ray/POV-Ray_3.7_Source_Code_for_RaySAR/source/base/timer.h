/*******************************************************************************
 * timer.h
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
 * $File: //depot/public/povray/3.x/source/base/timer.h $
 * $Revision: #1 $
 * $Change: 6069 $
 * $DateTime: 2013/11/06 11:59:40 $
 * $Author: chrisc $
 *******************************************************************************/

#ifndef POVRAY_BASE_TIMER_H
#define POVRAY_BASE_TIMER_H

#include <boost/thread.hpp>

#include "base/configbase.h"

namespace pov_base
{

/**
 *  Wait for the specified time.
 *  @param  msec  Milliseconds to wait.
 */
void Delay(unsigned int msec);


#ifndef POV_TIMER

/// Use default time if no platform specific implementation is provided.
#define POV_TIMER TimerDefault

/**
 *  Default class for millisecond-precision timers.
 */
class TimerDefault
{
	public:
		/**
		 *  Create a new timer and start it.
		 *  @param  tctime  CPU time shall only be that consumed
		 *                  by the current thread (if supported).
		 */
		TimerDefault(bool tctime = false);

		/**
		 *  Destructor.
		 */
		~TimerDefault();

		/**
		 *  Determine elapsed real time since creation or last reset.
		 *  @return         Elapsed real time in milliseconds.
		 */
		POV_LONG ElapsedRealTime() const;

		/**
		 *  Determine elapsed CPU time since creation or last reset.
		 *  If not supported, return elapsed real time instead!
		 *  @return         Elapsed CPU time in milliseconds.
		 */
		POV_LONG ElapsedCPUTime() const;

		/**
		 *  Reset the timer.
		 */
		void Reset();

		/**
		 *  Determine if CPU time is supported for the current settings.
		 *  in particular this has to return false if keeping the CPU time
		 *  for the current thread is not supported but was requsted!
		 *  @return         True if the CPU time is validfor the current
		 *                  requested settings, false otherwise.
		 */
		bool HasValidCPUTime() const;
	private:
		/// thread CPU time flag
		bool threadCPUTimeOnly;
		/// real time at last reset
		boost::xtime realTimeStart;
		/// CPU time at last reset
		boost::xtime cpuTimeStart;
};

#endif

}

#ifdef USE_SYSPROTO
#include "syspovprotobase.h"
#endif

namespace pov_base
{

/// Millisecond-precision timer
typedef POV_TIMER Timer;

}

#endif // POVRAY_BASE_TIMER_H
