/*
 * RNG_taus.cc
 */

#include "RNG_taus.h"
#include <cmath>                       // this is only here for the
                                       // random_normal functions

/* Originally rng/taus.c, then RNG_taus.h, now split into RNG_taus.h,
 * RNG_taus.cc
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 James Theiler, Brian Gough
 * Modified A. Alan Middleton with minor mods to fit it into a small
 * portable C++ class. See original code for comparison.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
 * USA.
 */



#define MASK 0xffffffffUL
#define TAUSWORTHE(s,a,b,c,d) (((s &c) <<d) &MASK) ^ ((((s <<a) &MASK)^s) >>b)
#define LCG(n) ((69069 * n) & 0xffffffffUL)


RNG_taus::RNG_taus()
{
  set(0);
}
RNG_taus::RNG_taus (unsigned long int s)
{
  set(s);
}
unsigned long RNG_taus::get_UL ()
{
  s1 = TAUSWORTHE (s1, 13, 19, 4294967294UL, 12);
  s2 = TAUSWORTHE (s2, 2, 25, 4294967288UL, 4);
  s3 = TAUSWORTHE (s3, 3, 11, 4294967280UL, 17);
  return (s1 ^ s2 ^ s3);
}

double RNG_taus::get_double ()
{
  return get_UL () / ((double)MAX_UL + 1);
}

void RNG_taus::set (unsigned long int s)
{
  if (s == 0)
    s = 1;      /* default seed is 1 */

  s1 = LCG (s);
  if (s1 < 2) s1 += 2UL;
    s2 = LCG (s1);
  if (s2 < 8) s2 += 8UL;
    s3 = LCG (s2);
  if (s3 < 16) s3 += 16UL;
    
  /* "warm it up" */
  get_UL ();
  get_UL ();
  get_UL ();
  get_UL ();
  get_UL ();
  get_UL ();
  return;
}



double random_normal(RNG_taus& rng)
{
  static double x;
  static double y;
  static int use_or_make = 0;

  double return_value;

  if (use_or_make == 0)
  {    
    x = ((double)rng.get_UL()+1.) / ((double)MAX_UL + 2.);
                                       // increase in this case to exclude 0
    y = (double)rng.get_UL() / ((double)MAX_UL + 1.);
                                       // want to include 0 in this variable

    return_value = sqrt(-2*log(x)) * sin(2*M_PI*y);
  }
  else
    return_value = sqrt(-2*log(x)) * cos(2*M_PI*y);
  use_or_make = (use_or_make==1)?0:1;
  return return_value;
}


/*
 * I copied the algorithm from pmyrand. I don't understand why this is
 * necessary, and I didn't bother keeping a separate generator so the random
 * number generation is unchanged by switching from rantom_normal to
 * perturbed_random_normal.  I can do this in the future if necessary.
 *
 * Creighton
 */
double perturbed_random_normal(RNG_taus& rng, double delta)
{
  double pfac = sqrt(1.0 + delta * delta);
  double x;
  double y;

  double return_value;

  x = ((double)rng.get_UL()+1.) / ((double)MAX_UL + 2.);
                                       // increase in this case to exclude 0
  y = (double)rng.get_UL() / ((double)MAX_UL + 1.);
                                       // want to include 0 in this variable

  return_value = sqrt(-2*log(x)) * ( sin(2*M_PI*y) + delta*cos(2*M_PI*y) ) / pfac;
  return return_value;
}

