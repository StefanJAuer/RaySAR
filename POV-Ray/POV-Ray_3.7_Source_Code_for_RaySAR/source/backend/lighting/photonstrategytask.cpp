/*******************************************************************************
 * photonstrategytask.cpp
 *
 * This module implements Photon Mapping.
 *
 * Author: Nathan Kopp
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
 * $File: //depot/public/povray/3.x/source/backend/lighting/photonstrategytask.cpp $
 * $Revision: #1 $
 * $Change: 6069 $
 * $DateTime: 2013/11/06 11:59:40 $
 * $Author: chrisc $
 *******************************************************************************/

// frame.h must always be the first POV file included (pulls in platform config)
#include "backend/frame.h"
#include "base/povms.h"
#include "base/povmsgid.h"
#include "backend/math/vector.h"
#include "backend/math/matrices.h"
#include "backend/scene/objects.h"
#include "backend/shape/csg.h"
#include "backend/support/octree.h"
#include "backend/bounding/bbox.h"
#include "backend/scene/threaddata.h"
#include "backend/scene/scene.h"
#include "backend/scene/view.h"
#include "backend/support/msgutil.h"
#include "backend/lighting/point.h"
#include "backend/lighting/photonshootingstrategy.h"
#include "backend/lighting/photonstrategytask.h"
#include "lightgrp.h"

#include <algorithm>

// this must be the last file included
#include "base/povdebug.h"

namespace pov
{

PhotonStrategyTask::PhotonStrategyTask(ViewData *vd, PhotonShootingStrategy* strategy) :
	RenderTask(vd),
	strategy(strategy),
	messageFactory(10, 370, "Photon", vd->GetSceneData()->backendAddress, vd->GetSceneData()->frontendAddress, vd->GetSceneData()->sceneId, 0), // TODO FIXME - Values need to come from the correct place!
	cooperate(*this)
{
	// do nothing
}

PhotonStrategyTask::~PhotonStrategyTask()
{
	//delete strategy;
}

void PhotonStrategyTask::SendProgress(void)
{
	if (timer.ElapsedRealTime() > 1000)
	{
		timer.Reset();
		POVMS_Object obj(kPOVObjectClass_PhotonProgress);
		obj.SetInt(kPOVAttrib_CurrentPhotonCount, GetSceneData()->surfacePhotonMap.numPhotons + GetSceneData()->mediaPhotonMap.numPhotons);
		RenderBackend::SendViewOutput(GetViewData()->GetViewId(), GetSceneData()->frontendAddress, kPOVMsgIdent_Progress, obj);
	}
}

void PhotonStrategyTask::Run()
{
	// quit right away if photons not enabled
	if (!GetSceneData()->photonSettings.photonsEnabled) return;

	Cooperate();

	/*  loop through global light sources  */
	GetViewDataPtr()->Light_Is_Global = true;
	for(vector<LightSource *>::iterator Light = GetSceneData()->lightSources.begin(); Light != GetSceneData()->lightSources.end(); Light++)
	{
		if ((*Light)->Light_Type != FILL_LIGHT_SOURCE)
		{
			if ((*Light)->Light_Type == CYLINDER_SOURCE && !(*Light)->Parallel)
				messageFactory.Warning(0,"Cylinder lights should be parallel when used with photons.");

			/* do object-specific lighting */
			SearchThroughObjectsCreateUnits(GetSceneData()->objects, (*Light));
		}

		Cooperate();
	}

	// loop through light_group light sources
/*
	TODO
	GetSceneData()->photonSettings.Light_Is_Global = false;
	for(vector<LightSource *>::iterator Light_Group_Light = GetSceneData()->lightGroupLights.begin(); Light_Group_Light != GetSceneData()->lightGroupLights.end(); Light_Group_Light++)
	{
		Light = Light_Group_Light->Light;

		if (Light->Light_Type == CYLINDER_SOURCE && !Light->Parallel)
			messageFactory.Warning(0,"Cylinder lights should be parallel when used with photons.");

		// do object-specific lighting
		SearchThroughObjectsCreateUnits(GetSceneData()->objects, Light);

		Cooperate();
		SendProgress();
	}
*/
	// good idea to make sure all warnings and errors arrive frontend now [trf]
	Cooperate();

	strategy->start();
}

void PhotonStrategyTask::Stopped()
{
	// nothing to do for now [trf]
}

void PhotonStrategyTask::Finish()
{
	GetViewDataPtr()->timeType = SceneThreadData::kPhotonTime;
	GetViewDataPtr()->realTime = ConsumedRealTime();
	GetViewDataPtr()->cpuTime = ConsumedCPUTime();
}


/*****************************************************************************

 FUNCTION

  SearchThroughObjectsCreateUnits()

  Searches through 'object' and all siblings  and children of 'object' to
  locate objects with PH_TARGET_FLAG set.  This flag means that the object
  receives photons.

  Preconditions:
    Photon mapping initialized (InitBacktraceEverything() called)
    'Object' is a object (with or without siblings)
    'Light' is a light source in the scene

  Postconditions:
    Work units (combination of light + target) have been added to the
    strategy.

******************************************************************************/

void PhotonStrategyTask::SearchThroughObjectsCreateUnits(vector<ObjectPtr>& Objects, LightSource *Light)
{
	boost::shared_ptr<SceneData> sceneData = GetSceneData();

	/* check this object and all siblings */
	for(vector<ObjectPtr>::iterator Sib = Objects.begin(); Sib != Objects.end(); Sib++)
	{
		if(Test_Flag((*Sib), PH_TARGET_FLAG) &&
		   !((*Sib)->Type & LIGHT_SOURCE_OBJECT))
		{
			/* do not shoot photons if global lights are turned off for ObjectPtr */
			if(!Test_Flag((*Sib), NO_GLOBAL_LIGHTS_FLAG))
			{
				strategy->createUnitsForCombo((*Sib), Light, sceneData);
			}

			Cooperate();
			SendProgress();
		}
		/* if it has children, check them too */
		else if(((*Sib)->Type & IS_COMPOUND_OBJECT))
		{
			SearchThroughObjectsCreateUnits(((CSG *)(*Sib))->children, Light);
		}
	}
}


}
