/*******************************************************************************
 * console.h
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
 * $File: //depot/public/povray/3.x/source/frontend/console.h $
 * $Revision: #1 $
 * $Change: 6069 $
 * $DateTime: 2013/11/06 11:59:40 $
 * $Author: chrisc $
 *******************************************************************************/

#ifndef POVRAY_FRONTEND_CONSOLE_H
#define POVRAY_FRONTEND_CONSOLE_H

#include "base/textstreambuffer.h"

#include "frontend/configfrontend.h"

#include <string>

namespace pov_frontend
{

using namespace pov_base;

class Console : public TextStreamBuffer
{
	public:
		Console(unsigned int wrapwidth = 80);
		virtual ~Console();

		virtual void Initialise() = 0;
		virtual void Output(const string&) = 0;
	private:
		void lineoutput(const char *str, unsigned int chars);
};

}

#endif // POVRAY_FRONTEND_CONSOLE_H
