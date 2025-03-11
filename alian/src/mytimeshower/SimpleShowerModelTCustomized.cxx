// ShowerModel.cc is a part of the PYTHIA event generator.
// Copyright (C) 2023 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the
// ShowerModel class.

#include "SimpleShowerModelTCustomized.hh"
#include "SimpleTimeShowerCustomized.hh"
#include "Pythia8/SimpleSpaceShower.h"
#include "Pythia8/Merging.h"
#include "Pythia8/MergingHooks.h"

namespace Pythia8 {

//==========================================================================

// Initialize the SimpleShowerModel.
void SimpleShowerModelTCustomized::addParameters(Settings &settings)
{
		SimpleTimeShowerCustomized::addParameters(settings);
}

bool SimpleShowerModelTCustomized::init(MergingPtr mergPtrIn,
																					MergingHooksPtr mergHooksPtrIn, PartonVertexPtr,
																					WeightContainer *)
	{
		std::cout << "Initializing SimpleShowerModelTCustomized" << std::endl;
		subObjects.clear();
		mergingPtr = mergPtrIn;
		if (mergingPtr)
			registerSubObject(*mergingPtr);
		mergingHooksPtr = mergHooksPtrIn;
		if (mergingHooksPtr)
			registerSubObject(*mergingHooksPtr);
		timesPtr = timesDecPtr = make_shared<SimpleTimeShowerCustomized>();
		registerSubObject(*timesPtr);
		spacePtr = make_shared<SimpleSpaceShower>();
		registerSubObject(*spacePtr);
		return true;
	}

	//==========================================================================

} // end namespace Pythia8
