#ifndef BOX_H
#define BOX_H


#include "../Resources/std_include.h"
#include "BaseBox.h"
#include "PeriodicBox.h"

namespace LiuJamming
{

template <int Dim>
static std::map< string, CBox<Dim>* > CreateBoxMap()
{
	std::map<string,CBox<Dim>*> m;
	CBox<Dim> *p;

	//For each box...
	p = new CPeriodicBox<Dim,0>();   m[CPeriodicBox<Dim,0>::GetName()] = p;
	p = new CPeriodicBox<Dim,1>();   m[CPeriodicBox<Dim,1>::GetName()] = p;
	p = new CPeriodicBox<Dim,Dim>(); m[CPeriodicBox<Dim,Dim>::GetName()] = p;

	return m;

}

template<> std::map<string,CBox<2>*> CBox<2>::BoxTypes = CreateBoxMap<2>();
template<> std::map<string,CBox<3>*> CBox<3>::BoxTypes = CreateBoxMap<3>();











}

#endif //BOX_H

