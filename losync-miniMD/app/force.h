/* ----------------------------------------------------------------------
   miniMD is a simple, parallel molecular dynamics (MD) code.   miniMD is
   an MD microapplication in the Mantevo project at Sandia National
   Laboratories ( http://www.mantevo.org ). The primary
   authors of miniMD are Steve Plimpton (sjplimp@sandia.gov) , Paul Crozier
   (pscrozi@sandia.gov) and Christian Trott (crtrott@sandia.gov).

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This library is free software; you
   can redistribute it and/or modify it under the terms of the GNU Lesser
   General Public License as published by the Free Software Foundation;
   either version 3 of the License, or (at your option) any later
   version.

   This library is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this software; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA.  See also: http://www.gnu.org/licenses/lgpl.txt .

   For questions, contact Paul S. Crozier (pscrozi@sandia.gov) or
   Christian Trott (crtrott@sandia.gov).

   Please read the accompanying README and LICENSE files.
---------------------------------------------------------------------- */

#ifndef FORCE_H_
#define FORCE_H_

#include "ljs.h"
#include "atom.h"
#include "neighbor.h"
#include "comm.h"

class Force
{
  public:
    MMD_float cutforce;     // DSM Force cutoff. Constant, input parameter.
    MMD_float* cutforcesq;  // DSM Set to constant cutforce * cutforce in setup() for all atom types
    MMD_float mass;         // DSM Seems unused: only ever written to in EAM branch, no instances of reads found.
    int ntypes;             // DSM Number of atom types. Constant set on object creation

    // DSM Multibox change: these used to be single variables but extended to arrays (one per box)
    MMD_float* eng_vdwl;    // DSM One of the outputs of compute(). Used in energy()
    MMD_float* virial;      // DSM One of the outputs of compute(). Used in pressure()
    MMD_int* evflag;        // DSM Controls whether eng_vdwl and virial are set in a compute() call or not
    int boxes_per_process;  // DSM: Multibox change

    Force() {};
    virtual ~Force() {};
    virtual void setup() {};
    virtual void finalise() {};
    virtual void compute(Atom &, Neighbor &, Comm &, int) {};

    int use_sse;            // DSM Flag/constant in main() - signals use of vectorised compute function
    int use_oldcompute;     // DSM Flag/constant in main() - signals use of unoptimised compute function
    ThreadData* threads;    // DSM Object containing thread and rank IDs (all constants)
    MMD_int reneigh;        // DSM Constant set in ForceLJ constructor
    Timer* timer;           // DSM For timestamping and profiling

    // DSM Lennard-Jones potential input parameters. Constants calculated in main().
    // DSM Arrays as one entry per atom type but all are set to the same value in main() (in.epsilon and in.sigma)
    MMD_float *epsilon, *sigma6, *sigma; //Parameters for LJ only

    ForceStyle style;
  protected:

    MMD_int me;             // DSM Duplicate information, already in threads object (only used in EAM in any case).
};

#endif
