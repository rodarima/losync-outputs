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

#include "atom.h"
#include "force.h"
#include "neighbor.h"
#include "comm.h"
#include "thermo.h"
#include "timer.h"
#include "threadData.h"

class Integrate
{
  public:
    MMD_float dt;
    MMD_float dtforce;
    MMD_int ntimes;
    // DSM Multibox: These converted to local variables in the functions where they are used.
    /*MMD_int nlocal, nmax;
    MMD_float* x, *v, *f, *xold;*/
    MMD_float mass;

    MMD_int sort_every;

    Integrate();
    ~Integrate();
    void setup();
    // DSM Multibox: These functions signatures changes as x, v, f and xold are no longer class attributes
    void initialIntegrate(MMD_float* x, MMD_float* v, MMD_float* f, MMD_int nlocal);
    void finalIntegrate(MMD_float* v, MMD_float* f, MMD_int nlocal);
    // DSM Multibox change: Atom now an array of atoms[] and Neighbor now an attribute of atoms
    void run(Atom* atoms[], Force*, Comm &, Thermo &, Timer &);

    ThreadData* threads;
};
