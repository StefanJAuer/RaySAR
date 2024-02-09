/*******************************************************************************
 * boundingtask.cpp
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
 * $File: //depot/public/povray/3.x/source/backend/bounding/boundingtask.cpp $
 * $Revision: #1 $
 * $Change: 6069 $
 * $DateTime: 2013/11/06 11:59:40 $
 * $Author: chrisc $
 *******************************************************************************/

#include <set>

#include <boost/thread.hpp>
#include <boost/bind.hpp>

// frame.h must always be the first POV file included (pulls in platform config)
#include "backend/frame.h"
#include "backend/math/vector.h"
#include "backend/math/matrices.h"
#include "backend/scene/objects.h"
#include "backend/support/bsptree.h"
#include "backend/bounding/boundingtask.h"
#include "backend/shape/cones.h"
#include "backend/texture/texture.h"
#include "backend/texture/pigment.h"

// this must be the last file included
#include "base/povdebug.h"

namespace pov
{

class SceneObjects : public BSPTree::Objects
{
	public:
		vector<ObjectPtr> infinite;
		vector<ObjectPtr> finite;
		unsigned int numLights;

		SceneObjects(vector<ObjectPtr>& objects)
		{
			numLights = 0;
			for(vector<ObjectPtr>::iterator i(objects.begin()); i != objects.end(); i++)
			{
				if(Test_Flag((*i), INFINITE_FLAG))
				{
					infinite.push_back(*i);
					if (((*i)->Type & LIGHT_SOURCE_OBJECT) != 0)
						numLights++;
				}
				else
					finite.push_back(*i);
			}
		}

		virtual ~SceneObjects()
		{
			// nothing to do
		}

		virtual unsigned int size() const
		{
			return finite.size();
		}

		virtual float GetMin(unsigned int axis, unsigned int i) const
		{
			return finite[i]->BBox.lowerleft[axis];
		}

		virtual float GetMax(unsigned int axis, unsigned int i) const
		{
			return (finite[i]->BBox.lowerleft[axis] + finite[i]->BBox.length[axis]);
		}
};

class BSPProgress : public BSPTree::Progress
{
	public:
		BSPProgress(RenderBackend::SceneId sid, POVMSAddress addr) :
			sceneId(sid),
			frontendAddress(addr),
			lastProgressTime(0)
		{
		}

		virtual void operator()(unsigned int nodes) const
		{
			if((timer.ElapsedRealTime() - lastProgressTime) > 1000) // update progress at most every second
			{
				POVMS_Object obj(kPOVObjectClass_BoundingProgress);
				obj.SetLong(kPOVAttrib_RealTime, timer.ElapsedRealTime());
				obj.SetLong(kPOVAttrib_CurrentNodeCount, nodes);
				RenderBackend::SendSceneOutput(sceneId, frontendAddress, kPOVMsgIdent_Progress, obj);

				Task::CurrentTaskCooperate();

				lastProgressTime = timer.ElapsedRealTime();
			}
		}
	private:
		RenderBackend::SceneId sceneId;
		POVMSAddress frontendAddress;
		Timer timer;
		mutable POV_LONG lastProgressTime;

		BSPProgress();
};

BoundingTask::BoundingTask(boost::shared_ptr<SceneData> sd, unsigned int bt) :
	Task(new SceneThreadData(sd), boost::bind(&BoundingTask::SendFatalError, this, _1)),
	sceneData(sd),
	boundingThreshold(bt)
{
}

BoundingTask::~BoundingTask()
{
}

void BoundingTask::AppendObject(ObjectPtr p)
{
	sceneData->objects.push_back(p);
}

void BoundingTask::Run()
{
	if((sceneData->objects.size() < boundingThreshold) || (sceneData->boundingMethod == 0))
	{
		SceneObjects objects(sceneData->objects);
		sceneData->boundingMethod = 0;
		sceneData->numberOfFiniteObjects = objects.finite.size();
		sceneData->numberOfInfiniteObjects = objects.infinite.size() - objects.numLights;
		return;
	}

	switch(sceneData->boundingMethod)
	{
		case 2:
		{
			// new BSP tree code
			SceneObjects objects(sceneData->objects);
			BSPProgress progress(sceneData->sceneId, sceneData->frontendAddress);

			sceneData->objects.clear();
			sceneData->objects.insert(sceneData->objects.end(), objects.finite.begin(), objects.finite.end());
			sceneData->objects.insert(sceneData->objects.end(), objects.infinite.begin(), objects.infinite.end());
			sceneData->numberOfFiniteObjects = objects.finite.size();
			sceneData->numberOfInfiniteObjects = objects.infinite.size() - objects.numLights;
			sceneData->tree = new BSPTree(sceneData->bspMaxDepth, sceneData->bspObjectIsectCost, sceneData->bspBaseAccessCost, sceneData->bspChildAccessCost, sceneData->bspMissChance);
			sceneData->tree->build(progress, objects,
			                       sceneData->nodes, sceneData->splitNodes, sceneData->objectNodes, sceneData->emptyNodes,
			                       sceneData->maxObjects, sceneData->averageObjects, sceneData->maxDepth, sceneData->averageDepth,
			                       sceneData->aborts, sceneData->averageAborts, sceneData->averageAbortObjects, sceneData->inputFile);
			break;
		}
		case 1:
		{
			// old bounding box code
			unsigned int numberOfLightSources;

			Build_Bounding_Slabs(&(sceneData->boundingSlabs), sceneData->objects, sceneData->numberOfFiniteObjects,
			                     sceneData->numberOfInfiniteObjects, numberOfLightSources);
			break;
		}
	}
}

void BoundingTask::Stopped()
{
}

void BoundingTask::Finish()
{
	GetSceneDataPtr()->timeType = SceneThreadData::kBoundingTime;
	GetSceneDataPtr()->realTime = ConsumedRealTime();
	GetSceneDataPtr()->cpuTime = ConsumedCPUTime();
}

void BoundingTask::SendFatalError(Exception& e)
{
	// if the front-end has been told about this exception already, we don't tell it again
	if (e.frontendnotified(true))
		return;

	POVMS_Message msg(kPOVObjectClass_ControlData, kPOVMsgClass_SceneOutput, kPOVMsgIdent_Error);

	msg.SetString(kPOVAttrib_EnglishText, e.what());
	msg.SetInt(kPOVAttrib_Error, 0);
	msg.SetInt(kPOVAttrib_SceneId, sceneData->sceneId);
	msg.SetSourceAddress(sceneData->backendAddress);
	msg.SetDestinationAddress(sceneData->frontendAddress);

	POVMS_SendMessage(msg);
}

}
