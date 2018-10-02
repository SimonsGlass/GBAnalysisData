#ifndef POTENTIALS_H
#define POTENTIALS_H

#include "../Resources/std_include.h"
#include "BasePotential.h"
#include "HarmonicPotential.h"
#include "HertzianPotential.h"
#include "SoftPotential.h"
#include "LJPotential.h"

namespace LiuJamming
{



/*
 *	For each potential type, create an instance using the default constructor.
 *	Each instance then gets put in a map so that they can be referenced by name.
 *	All future potential objects can then be created by calling the Clone() method
 *	on one of these instances. For example
 *
 *		string name = CHarmonicPotential::GetName();
 *		CPotential *pot = m[name]->Clone();
 *
 */

static std::map<string,CPotential*> CreatePotentialMap()
{
	std::map<string,CPotential*> m;
	CPotential *p;
	
	//For each potential...
	p = new CHarmonicPotential(); m[CHarmonicPotential::GetName()] = p;
	p = new CHertzianPotential(); m[CHertzianPotential::GetName()] = p;
	p = new CSoftPotential();     m[CSoftPotential    ::GetName()] = p;
	p = new CLJPotential(); m[CLJPotential::GetName()] = p;
	return m;
}

std::map<string,CPotential*> CPotential::PotentialTypes = CreatePotentialMap();


};

#endif //POTENTIALS_H




