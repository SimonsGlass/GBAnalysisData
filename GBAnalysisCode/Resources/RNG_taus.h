/*
 * RNG_taus.h
 */

#ifndef RNG_TAUS
#define RNG_TAUS
#include <cstdlib>

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

/* Modified A. Alan Middleton with minor mods to fit it into a small
 * portable C++ class. See original code for comparison.
 * Changes: 1) made the RNG a class and the 3 integer state
 *             state part of the class, so don't need state in argument.
 *          2) removed coupling to standard GNU RNG.
 *          3) renamed methods
 * At least replicates x_{10007} cited below (see taus_test.cpp)
 */
/* This is a maximally equidistributed combined Tausworthe
   generator. The sequence is,

   x_n = (s1_n ^ s2_n ^ s3_n) 

   s1_{n+1} = (((s1_n & 4294967294) <<12) ^ (((s1_n <<13) ^ s1_n) >>19))
   s2_{n+1} = (((s2_n & 4294967288) << 4) ^ (((s2_n << 2) ^ s2_n) >>25))
   s3_{n+1} = (((s3_n & 4294967280) <<17) ^ (((s3_n << 3) ^ s3_n) >>11))

   computed modulo 2^32. In the three formulas above '^' means
   exclusive-or (C-notation), not exponentiation. Note that the
   algorithm relies on the properties of 32-bit unsigned integers (it
   is formally defined on bit-vectors of length 32). I have added a
   bitmask to make it work on 64 bit machines.

   We initialize the generator with s1_1 .. s3_1 = s_n MOD m, where
   s_n = (69069 * s_{n-1}) mod 2^32, and s_0 = s is the user-supplied
   seed.

   The theoretical value of x_{10007} is 2733957125. The subscript
   10007 means (1) seed the generator with s=1 (2) do six warm-up
   iterations, (3) then do 10000 actual iterations.

   The period of this generator is about 2^88.

   From: P. L'Ecuyer, "Maximally Equidistributed Combined Tausworthe
   Generators", Mathematics of Computation, 65, 213 (1996), 203--213.

   This is available on the net from L'Ecuyer's home page,

   http://www.iro.umontreal.ca/~lecuyer/myftp/papers/tausme.ps
   ftp://ftp.iro.umontreal.ca/pub/simulation/lecuyer/papers/tausme.ps 

   Update: April 2002

   There is an erratum in the paper "Tables of Maximally
   Equidistributed Combined LFSR Generators", Mathematics of
   Computation, 68, 225 (1999), 261--269:
   http://www.iro.umontreal.ca/~lecuyer/myftp/papers/tausme2.ps

        ... the k_j most significant bits of z_j must be non-
        zero, for each j. (Note: this restriction also applies to the 
        computer code given in [4], but was mistakenly not mentioned in
        that paper.)
   
   This affects the seeding procedure by imposing the requirement
   s1 > 1, s2 > 7, s3 > 15.

   The generator taus2 has been added to satisfy this requirement.
   The original taus generator is unchanged.

   Update: November 2002

   There was a bug in the correction to the seeding procedure for s2.
   It affected the following seeds 254679140 1264751179 1519430319
   2274823218 2529502358 3284895257 3539574397 (s2 < 8).

*/

class RNG_taus {
public:
  RNG_taus ();                         // default constructor: seed of 1
  RNG_taus (unsigned long int s);      // seeding constructor: seed of s
  unsigned long get_UL (); 
  double get_double ();                // return a random double in [0,1)
  void set (unsigned long int s);
private:
  unsigned long int s1, s2, s3;
};

const unsigned long int MAX_UL = 4294967295U;
double random_normal(RNG_taus& rng);
double perturbed_random_normal(RNG_taus& rng, double delta);
#endif
