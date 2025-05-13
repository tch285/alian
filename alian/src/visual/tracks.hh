#ifndef ALIAN_DEMO_TRACKS_HH
#define ALIAN_DEMO_TRACKS_HH

#include <ROOT/REveManager.hxx>
#include <ROOT/REveElement.hxx>

namespace alian
{
	void makeTracks(int N_Tracks, ROOT::Experimental::REveElement *trackHolder);
	ROOT::Experimental::REveManager *tracks();
}

#endif