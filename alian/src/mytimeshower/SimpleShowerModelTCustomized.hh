#ifndef Pythia8_SimpleShowerModelTCustomized_H
#define Pythia8_SimpleShowerModelTCustomized_H

#include "Pythia8/ShowerModel.h"

namespace Pythia8
{
	//==========================================================================

	// The shower model class handling the default Pythia shower model
	// with SimpleTimeShower and SimpleSpaceShower classes.

	class SimpleShowerModelTCustomized : public ShowerModel
	{

	public:
		static void addParameters(Settings &settings);
		
		// Empty constructor.
		SimpleShowerModelTCustomized() = default;

		// Empty virtual destructor
		virtual ~SimpleShowerModelTCustomized() override {}

		// Function called from Pythia after the basic pointers has been set.
		virtual bool init(MergingPtr mergPtrIn, MergingHooksPtr mergHooksPtrIn,
											PartonVertexPtr partonVertexPtrIn,
											WeightContainer *weightContainerPtrIn) override;

		// Function called from Pythia after the beam particles have been set up,
		// so that showers may be initialized after the beams are initialized.
		// Currently only dummy dunction.
		virtual bool initAfterBeams() override { return true; }
	};

	//==========================================================================
}

#endif // Pythia8_SimpleShowerModelTCustomized_H
