#ifndef HEPPYY_DEMO_HH
#define HEPPYY_DEMO_HH

#include <vector>
#include <list>
#include <fstream>
#include <string>

#include<fastjet/PseudoJet.hh>
#include <TDatabasePDG.h>

namespace alian
{
	class DemoClass
	{
	public:
		DemoClass();
		std::vector<fastjet::PseudoJet> divideby(const std::vector<fastjet::PseudoJet> &input, double scalar); // just to show fastjet functionality ok
		virtual ~DemoClass() {;}
	private:
		Double_t _scalar; // just to show root functionality ok
	};
};
#endif
