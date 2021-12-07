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
//#define PRINTDEBUG(a) a
#define PRINTDEBUG(a)
#include "stdio.h"
#include "integrate.h"
#include "math.h"

Integrate::Integrate() {sort_every=20;}
Integrate::~Integrate() {}

void Integrate::setup()
{
  dtforce = 0.5 * dt;
}

void Integrate::initialIntegrate(MMD_float* x, MMD_float* v, MMD_float* f, MMD_int nlocal)
{
  OMPFORSCHEDULE
  for(MMD_int i = 0; i < nlocal; i++) {
    v[i * PAD + 0] += dtforce * f[i * PAD + 0];
    v[i * PAD + 1] += dtforce * f[i * PAD + 1];
    v[i * PAD + 2] += dtforce * f[i * PAD + 2];
    x[i * PAD + 0] += dt * v[i * PAD + 0];
    x[i * PAD + 1] += dt * v[i * PAD + 1];
    x[i * PAD + 2] += dt * v[i * PAD + 2];
  }
}

// DSM Velocity update
void Integrate::finalIntegrate(MMD_float* v, MMD_float* f, MMD_int nlocal)
{
  OMPFORSCHEDULE
  for(MMD_int i = 0; i < nlocal; i++) {
    v[i * PAD + 0] += dtforce * f[i * PAD + 0];
    v[i * PAD + 1] += dtforce * f[i * PAD + 1];
    v[i * PAD + 2] += dtforce * f[i * PAD + 2];
  }

}

void Integrate::run(Atom* atoms[], Force* force,
                    Comm &comm, Thermo &thermo, Timer &timer)
{
  int i, n;

  comm.timer = &timer;
  timer.array[TIME_TEST] = 0.0;

  int check_safeexchange = comm.check_safeexchange;
  const int every = (*atoms[0]->neighbor).every; // DSM Multibox: "every" is constant across all neighbors and we can assume at least 1 box
  const int boxes_per_process = atoms[0]->boxes_per_process;

  char *initialIntegrateSentinels = comm.initialIntegrateSentinels;
  char *sortSentinels = comm.sortSentinels;
  char *exchangePBCSentinels = comm.exchangePBCSentinels;
  char *communicateSentinels = comm.communicateSentinels;
  char **communicatePackSentinels = comm.communicatePackSentinels;
  char *communicateInternalPackSentinels = comm.communicateInternalPackSentinels;
  char *communicateInternalUnpackSentinels = comm.communicateInternalUnpackSentinels;
  char *bordersInternalSentinels = comm.bordersInternalSentinels;
  char *forceComputeSentinels = comm.forceComputeSentinels;
  char *neighbourBuildSentinels = comm.neighbourBuildSentinels;
  char *bordersUnpackSentinels = comm.bordersUnpackSentinels;
  char *exchangePackSentinels = comm.exchangePackSentinels;
  char *bordersPackSentinels = comm.bordersPackSentinels;

  mass = atoms[0]->mass; // DSM Multibox change: assuming all atoms have the same mass value here and we have at least 1
  dtforce = dtforce / mass;
  //#pragma omp parallel private(i,n)
  {
    // Replaced with an array (one per box)
    int next_sort[atoms[0]->boxes_per_process];
    for (int box_index = 0; box_index < atoms[0]->boxes_per_process; ++box_index) {
      next_sort[box_index] = sort_every > 0 ? sort_every : ntimes + 1;
    }

    // DSM: The main loop over all timesteps. Split into:
    //   * initialIntegrate()
    //   * Almost always:
    //       - comm.communicate()
    //     but every neighbor.every timesteps instead do:
    //       - comm.exchange()
    //       - comm.borders()
    //       - neighbor.build()
    //   * force->compute()
    //   * finalIntegrate()
    //   * thermo.compute() (once, gives final output)
    //
    // comm.communicate is a cheaper function that is called most of the time. Updates atom positions.
    // comm.exchange+comm.borders are more expensive and called occasionally to update neighbours and change atom boxes
    for(n = 0; n < ntimes; n++) {
      
      //if (threads->mpi_me == 0) {
      //  printf("Starting iteration %d\n", n);
      //}

      // --- initialIntegrate ---
      for(int box_index = 0; box_index < atoms[0]->boxes_per_process; ++box_index) {
        // in-dependencies on all packs for this box. 26 way communication but
        // 3 (box layers) * 10 (neighbours) elements in sentinels array to simplify allocation.
        // Packs do not need to complete until this point; force calculation and finalIntegrate can be done before or
        // after sends are posted.
        // Recvs from neighbours must have completed, however.
        #pragma oss task label(initialIntegrate) \
                         in(forceComputeSentinels[box_index]) \
                         in(communicatePackSentinels[box_index][0:30]) \
                         in(communicateInternalPackSentinels[box_index]) \
                         in(exchangePackSentinels[box_index]) \
                         in(bordersPackSentinels[box_index]) \
                         out(initialIntegrateSentinels[box_index]) firstprivate(box_index)
        {
          initialIntegrate(atoms[box_index]->x, atoms[box_index]->v, atoms[box_index]->f, atoms[box_index]->nlocal);
        }
      }

      // Timers broken in tasked version - only used when built with tasking disabled
#ifndef USE_TASKS
      timer.stamp();
#endif

      // --- comm.communicate() ---
      // DSM: communicate() almost every timestep. neighbor.every number of timesteps do more expensive exchange().
      if((n + 1) % every) {
        comm.communicate(atoms); // DSM: Multibox
#ifndef USE_TASKS
        timer.stamp(TIME_COMM);
#endif
      }

      // --- comm.exchange() ---
      if( !((n + 1) % every) ) {
#ifndef USE_TASKS
        timer.stamp_extra_start();
#endif
        comm.exchange(atoms);
      }

      // --- atom.sort() ---
      for(int box_index = 0; box_index < atoms[0]->boxes_per_process; ++box_index) {
        if( !((n + 1) % every) ) {
          if(n+1>=next_sort[box_index]) {
            #pragma oss task label(atom->sort) in(exchangePBCSentinels[box_index]) out(sortSentinels[box_index]) firstprivate(box_index)
            atoms[box_index]->sort(*atoms[box_index]->neighbor);

            next_sort[box_index] +=  sort_every;
          }
        }
      }

      // --- comm.borders() ---
      if( !((n + 1) % every) ) {
        comm.borders(atoms);
#ifndef USE_TASKS
        timer.stamp_extra_stop(TIME_TEST);
        timer.stamp(TIME_COMM);
#endif
      }

      // Only perform neighbour rebuild "every" iterations
      if( !((n + 1) % every) ) {
        // --- neighbor.build() ---
        // --- force->compute() ---
        // --- finalIntegrate --- combined tasks
        for(int box_index = 0; box_index < boxes_per_process; ++box_index) {
          // Depends on all borders unpack tasks completing, i.e. must have full
          // knowledge of ghost atoms before rebuilding neighbour list.
          // Does not modify x, no pack depedencies.
          // Dependencies on communicate tasks to prevent this task running before all non-rebuild iterations are complete
          #pragma oss task label("neighbor->build & force->compute & finalIntegrate") \
                           in(bordersUnpackSentinels[box_index]) \
                           in(bordersInternalSentinels[box_index]) \
                           in(communicateSentinels[box_index]) \
                           in(communicateInternalUnpackSentinels[box_index]) \
                           out(forceComputeSentinels[box_index]) \
                           firstprivate(box_index, n)
          {
            atoms[box_index]->neighbor->build(*atoms[box_index]);
#ifndef USE_TASKS
            timer.stamp(TIME_NEIGH);
#endif
            // DSM: thermo.nstat is a constant, an input file parameter fixed at initial setup
            force->evflag[box_index] = (n + 1) % thermo.nstat == 0; // Controls whether eng_vdwl and virial are set this compute call or not
            force->compute(*atoms[box_index], *atoms[box_index]->neighbor, *(Comm*)0, NULL); // last 2 arguments (comm & comm.me) are not used in force_lj implementation. Replace with nulls
#ifndef USE_TASKS
            timer.stamp(TIME_FORCE);
#endif
            finalIntegrate(atoms[box_index]->v, atoms[box_index]->f, atoms[box_index]->nlocal);
          }
        }
      } else {
        // --- force->compute() ---
        // --- finalIntegrate --- combined tasks
        for(int box_index = 0; box_index < boxes_per_process; ++box_index) {
          // No need for borders dependencies, borders tasks run only in reneighbouring branch
          #pragma oss task label("force->compute & finalIntegrate") \
                           in(communicateSentinels[box_index]) \
                           in(communicateInternalUnpackSentinels[box_index]) \
                           out(forceComputeSentinels[box_index]) \
                           firstprivate(box_index, n)
          {
            // DSM: thermo.nstat is a constant, an input file parameter fixed at initial setup
            force->evflag[box_index] = (n + 1) % thermo.nstat == 0; // Controls whether eng_vdwl and virial are set this compute call or not
            force->compute(*atoms[box_index], *atoms[box_index]->neighbor, *(Comm*)0, NULL); // last 2 arguments (comm & comm.me) are not used in force_lj implementation. Replace with nulls
#ifndef USE_TASKS
            timer.stamp(TIME_FORCE);
#endif
            finalIntegrate(atoms[box_index]->v, atoms[box_index]->f, atoms[box_index]->nlocal);
          }
        }
      }

      // --- thermo.compute() ---
      // DSM Multibox: This function requires output of force->compute() on the timesteps it is called.
      // Consists of energy(), temperature() and pressure() calculations = 3x MPI_Allreduces
      // Loops over boxes_per_process occur within these functions.

      if (!((n+1) % thermo.nstat)) {
        // Wait for all tasks to complete before doing final thermo calculation
        #pragma oss taskwait
        thermo.compute(n + 1, atoms, force, timer);
      }

    } // End main loop over timesteps
  } //end OpenMP parallel
} // End run() function
