// ***************************************************************
// This file was created using the ./CreateProject.sh script
// for project CrossSectionFinder
// ./CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#ifndef __COUNTINGEXP__H
#define __COUNTINGEXP__H

#include <BAT/BCModel.h>

// This is a CrossSectionFinder header file.
// Model source code is located in file CrossSectionFinder/CrossSectionFinder.cxx

// ---------------------------------------------------------
class CrossSectionFinder : public BCModel
{
	public:

		// Constructors and destructor
		CrossSectionFinder();
		CrossSectionFinder(const char * name);
		~CrossSectionFinder();

		// Methods to overload, see file CrossSectionFinder.cxx
		void DefineParameters();
		double LogAPrioriProbability(std::vector <double> parameters);
		double LogLikelihood(std::vector <double> parameters);
};
// ---------------------------------------------------------

#endif

