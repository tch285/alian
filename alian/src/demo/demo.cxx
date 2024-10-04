#include "demo.hh"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <limits>
#include <cmath>

#include <fastjet/PseudoJet.hh>

namespace alian
{
	DemoClass::DemoClass()
		: _scalar(0.0)
	{
	}

	std::vector<fastjet::PseudoJet> DemoClass::divideby(const std::vector<fastjet::PseudoJet> &input, double scalar)
	{
		_scalar = scalar;
		std::vector<fastjet::PseudoJet> output;
		for (std::vector<fastjet::PseudoJet>::const_iterator it = input.begin(); it != input.end(); ++it)
		{
			fastjet::PseudoJet psj = fastjet::PseudoJet(it->px(), it->py(), it->pz(), it->E());
			psj.set_user_index(it->user_index());
			psj /= _scalar;
			output.push_back(psj);
		}
		return output;
	}
};