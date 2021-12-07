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

#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "mpi.h"
#include "atom.h"
#include "neighbor.h"

#ifdef USE_TAMPI
#include <TAMPI.h>
#endif

#define DELTA 20000

Atom::Atom(int ntypes_, int boxes_per_process_)
{
  natoms = 0;
  nlocal = 0;
  nghost = 0;
  nmax = 0;
  copy_size = 0;
  boxes_per_process = boxes_per_process_;

  x = v = f = xold = x_copy = v_copy = NULL;
  type = type_copy = NULL;
  comm_size = 3;
  reverse_size = 3;
  border_size = 4;

  mass = 1;

  ntypes = ntypes_;

  // DSM: Multibox change - need one neighborlist per box rather than per process
  Neighbor _neighbor(ntypes_);
  neighbor = &_neighbor;
}

Atom::~Atom()
{
  if(nmax) {
    destroy_2d_MMD_float_array(x);
    destroy_2d_MMD_float_array(v);
    destroy_2d_MMD_float_array(f);
    destroy_2d_MMD_float_array(xold);
    destroy_1d_int_array(type);
  }
}

void Atom::growarray()
{
  int nold = nmax;
  nmax += DELTA;
  x = (MMD_float*) realloc_2d_MMD_float_array(x, nmax, PAD, PAD * nold);
  v = (MMD_float*) realloc_2d_MMD_float_array(v, nmax, PAD, PAD * nold);
  f = (MMD_float*) realloc_2d_MMD_float_array(f, nmax, PAD, PAD * nold);
  type = realloc_1d_int_array(type, nmax, nold);
  xold = (MMD_float*) realloc_2d_MMD_float_array(xold, nmax, PAD, PAD * nold);

  if(x == NULL || v == NULL || f == NULL || xold == NULL) {
    printf("ERROR: No memory for atoms\n");
  }
}

// DSM: Called in loop at initial setup from setup.cpp's create_atoms()
// DSM: x_in, y_in, z_in are atom coordinates. vx_in, vy_in, vz_in are atom velocities
void Atom::addatom(MMD_float x_in, MMD_float y_in, MMD_float z_in,
                   MMD_float vx_in, MMD_float vy_in, MMD_float vz_in)
{
  // DSM: If we're out of memory for new atoms, allocate arrays that are "DELTA" elements larger, copy existing contents
  // over and free the originals. On first call, nlocal and nmax are 0 - initial allocation is done without any copies.
  if(nlocal == nmax) growarray();

  // DSM: x-array is atom position, v-array is atom velocity
  x[nlocal*PAD + 0] = x_in;
  x[nlocal*PAD + 1] = y_in;
  x[nlocal*PAD + 2] = z_in;
  v[nlocal*PAD + 0] = vx_in;
  v[nlocal*PAD + 1] = vy_in;
  v[nlocal*PAD + 2] = vz_in;
  type[nlocal] = rand()%ntypes; // DSM: Counts for each atom type randomly decided here, not in input file.

  nlocal++;
}

/* enforce PBC
   order of 2 tests is important to insure lo-bound <= coord < hi-bound
   even with round-off errors where (coord +/- epsilon) +/- period = bound */

void Atom::pbc()
{
  //#pragma omp for
  for(int i = 0; i < nlocal; i++) {
    if(x[i*PAD + 0] < 0.0) x[i * PAD + 0] += box.xprd;

    if(x[i * PAD + 0] >= box.xprd) x[i * PAD + 0] -= box.xprd;

    if(x[i * PAD + 1] < 0.0) x[i * PAD + 1] += box.yprd;

    if(x[i * PAD + 1] >= box.yprd) x[i * PAD + 1] -= box.yprd;

    if(x[i * PAD + 2] < 0.0) x[i * PAD + 2] += box.zprd;

    if(x[i * PAD + 2] >= box.zprd) x[i * PAD + 2] -= box.zprd;
  }
}

void Atom::copy(int i, int j)
{
  x[j * PAD + 0] = x[i * PAD + 0];
  x[j * PAD + 1] = x[i * PAD + 1];
  x[j * PAD + 2] = x[i * PAD + 2];
  v[j * PAD + 0] = v[i * PAD + 0];
  v[j * PAD + 1] = v[i * PAD + 1];
  v[j * PAD + 2] = v[i * PAD + 2];
  type[j] = type[i];
}

// DSM: Packs box's position array x into buf for use in MPI_Send. PBC corrections applied here.
// DSM: Parameters:
//   - n          number of atoms/size of list
//   - list       array of indices/atom IDs into the x/atom position array. These are the atoms to be sent.
//   - buf        buffer atoms are packed into
//   - pbc_flags  periodic boundary conditions array
void Atom::pack_comm(int n, int* list, MMD_float* buf, int pbc_any, MMD_float pbc_x, MMD_float pbc_y, MMD_float pbc_z)
{
  int i, j;

  // DSM: Why does this branch exist? pbc_flags[1-3] can be -1, 1 or 0 depending on whether the box length is to be
  // subtracted, added or not included in the calculation (i.e. atom is within the box). Top and bottom branches give
  // identical answers when pbc_flags[1-3] == 0, i.e. when pbc_flags[0] == 0.
  if(pbc_any == 0) {

	//#pragma omp for schedule(static)
    for(i = 0; i < n; i++) {
      j = list[i];
      buf[3 * i] = x[j * PAD + 0];
      buf[3 * i + 1] = x[j * PAD + 1];
      buf[3 * i + 2] = x[j * PAD + 2];
    }
  } else {

    //#pragma omp for schedule(static)
    for(i = 0; i < n; i++) {
      j = list[i];
      buf[3 * i] = x[j * PAD + 0] + pbc_x;
      buf[3 * i + 1] = x[j * PAD + 1] + pbc_y;
      buf[3 * i + 2] = x[j * PAD + 2] + pbc_z;
    }
  }
}

void Atom::unpack_comm(int n, int first, MMD_float* buf)
{
  int i;

  //#pragma omp for schedule(static)
  for(i = 0; i < n; i++) {
    x[(first + i) * PAD + 0] = buf[3 * i];
    x[(first + i) * PAD + 1] = buf[3 * i + 1];
    x[(first + i) * PAD + 2] = buf[3 * i + 2];
  }
}

void Atom::pack_reverse(int n, int first, MMD_float* buf)
{
  int i;

  //#pragma omp for schedule(static)
  for(i = 0; i < n; i++) {
    buf[3 * i] = f[(first + i) * PAD + 0];
    buf[3 * i + 1] = f[(first + i) * PAD + 1];
    buf[3 * i + 2] = f[(first + i) * PAD + 2];
  }
}

void Atom::unpack_reverse(int n, int* list, MMD_float* buf)
{
  int i, j;

  //#pragma omp for schedule(static)
  for(i = 0; i < n; i++) {
    j = list[i];
    f[j * PAD + 0] += buf[3 * i];
    f[j * PAD + 1] += buf[3 * i + 1];
    f[j * PAD + 2] += buf[3 * i + 2];
  }
}

int Atom::pack_border(int i, MMD_float* buf, int* pbc_flags)
{
  int m = 0;

  // DSM: Is this branch necessary? pbc_flags are set to 0 if not periodic in a dimension - are the branches not equivalent?
  if(pbc_flags[0] == 0) {
    buf[m++] = x[i * PAD + 0];
    buf[m++] = x[i * PAD + 1];
    buf[m++] = x[i * PAD + 2];
    buf[m++] = type[i];
  } else {
    buf[m++] = x[i * PAD + 0] + pbc_flags[1] * box.xprd;
    buf[m++] = x[i * PAD + 1] + pbc_flags[2] * box.yprd;
    buf[m++] = x[i * PAD + 2] + pbc_flags[3] * box.zprd;
    buf[m++] = type[i];
  }

  return m;
}

int Atom::unpack_border(int i, MMD_float* buf)
{
  if(i == nmax) growarray(); // DSM: No change needed for multibox, this grows this box's x,v,f arrays

  int m = 0;
  x[i * PAD + 0] = buf[m++];
  x[i * PAD + 1] = buf[m++];
  x[i * PAD + 2] = buf[m++];
  type[i] = buf[m++];
  return m;
}

int Atom::pack_exchange(int i, MMD_float* buf)
{
  int m = 0;
  buf[m++] = x[i * PAD + 0];
  buf[m++] = x[i * PAD + 1];
  buf[m++] = x[i * PAD + 2];
  buf[m++] = v[i * PAD + 0];
  buf[m++] = v[i * PAD + 1];
  buf[m++] = v[i * PAD + 2];
  buf[m++] = type[i];
  return m;
}

int Atom::unpack_exchange(int i, MMD_float* buf)
{
  if(i == nmax) growarray();

  int m = 0;
  x[i * PAD + 0] = buf[m++];
  x[i * PAD + 1] = buf[m++];
  x[i * PAD + 2] = buf[m++];
  v[i * PAD + 0] = buf[m++];
  v[i * PAD + 1] = buf[m++];
  v[i * PAD + 2] = buf[m++];
  type[i] = buf[m++];
  return m;
}

int Atom::skip_exchange(MMD_float* buf)
{
  return 7;
}

/* realloc a 2-d MMD_float array */

MMD_float* Atom::realloc_2d_MMD_float_array(MMD_float* array,
    int n1, int n2, int nold)

{
  MMD_float* newarray;

  newarray = create_2d_MMD_float_array(n1, n2);

  if(nold) memcpy(newarray, array, nold * sizeof(MMD_float));

  destroy_2d_MMD_float_array(array);

  return newarray;
}

/* create a 2-d MMD_float array */

MMD_float* Atom::create_2d_MMD_float_array(int n1, int n2)
{
  MMD_float* array;

  if(n1 * n2 == 0) return NULL;

  #ifdef ALIGNMALLOC
    array = (MMD_float*) _mm_malloc((n1 * n2 + 1024 + 1) * sizeof(MMD_float), ALIGNMALLOC);
  #else
    array = (MMD_float*) malloc((n1 * n2 + 1024 + 1) * sizeof(MMD_float));
  #endif

  return array;
}

/* free memory of a 2-d MMD_float array */

void Atom::destroy_2d_MMD_float_array(MMD_float* array)
{
  if(array != NULL) {
  #ifdef ALIGNMALLOC
	_mm_free(array);
  #else
      free(array);
  #endif
  }
}

int * Atom::realloc_1d_int_array(int* array,
    int n1, int nold)

{
  int* newarray;

  newarray = create_1d_int_array(n1);

  if(nold) memcpy(newarray, array, nold * sizeof(int));

  destroy_1d_int_array(array);

  return newarray;
}

/* create a 2-d MMD_float array */

int * Atom::create_1d_int_array(int n1)
{
  int ALIGN = 16;
  int* data;
  int i, n;

  if(n1 == 0) return NULL;

  #ifdef ALIGNMALLOC
    data = (int*) _mm_malloc((n1 + 1024 + 1) * sizeof(int), ALIGNMALLOC);
  #else
    data = (int*) malloc((n1) * sizeof(int));
  #endif

  return data;
}

/* free memory of a 2-d MMD_float array */

void Atom::destroy_1d_int_array(int* array)
{
  if(array != NULL) {
  #ifdef ALIGNMALLOC
    _mm_free(array);
  #else
    free(array);
  #endif
  }
}

void Atom::sort(Neighbor &neighbor)
{
  neighbor.binatoms(*this,nlocal);
  //#pragma omp barrier

  binpos = neighbor.bincount;
  bins = neighbor.bins;

  const int mbins = neighbor.mbins;
  const int atoms_per_bin = neighbor.atoms_per_bin;

  //#pragma omp master
  {
    for(int i=1; i<mbins; i++)
	  binpos[i] += binpos[i-1];
    if(copy_size<nmax) {
	    destroy_2d_MMD_float_array(x_copy);
	    destroy_2d_MMD_float_array(v_copy);
	    destroy_1d_int_array(type_copy);
      x_copy = (MMD_float*) create_2d_MMD_float_array(nmax, PAD);
      v_copy = (MMD_float*) create_2d_MMD_float_array(nmax, PAD);
      type_copy = create_1d_int_array(nmax);
      copy_size = nmax;
    }
  }

  //#pragma omp barrier
  MMD_float* new_x = x_copy;
  MMD_float* new_v = v_copy;
  int* new_type = type_copy;
  MMD_float* old_x = x;
  MMD_float* old_v = v;
  int* old_type = type;

  //#pragma omp for
  for(int mybin = 0; mybin<mbins; mybin++) {
    const int start = mybin>0?binpos[mybin-1]:0;
    const int count = binpos[mybin] - start;
    for(int k=0; k<count; k++) {
	  const int new_i = start+k;
	  const int old_i = bins[mybin*atoms_per_bin+k];

	  new_x[new_i*PAD+0] = old_x[old_i*PAD+0];
	  new_x[new_i*PAD+1] = old_x[old_i*PAD+1];
	  new_x[new_i*PAD+2] = old_x[old_i*PAD+2];
	  new_v[new_i*PAD+0] = old_v[old_i*PAD+0];
	  new_v[new_i*PAD+1] = old_v[old_i*PAD+1];
	  new_v[new_i*PAD+2] = old_v[old_i*PAD+2];
	  new_type[new_i] = old_type[old_i];
    }
  }

  //#pragma omp master
  {
    MMD_float* x_tmp = x;
    MMD_float* v_tmp = v;
    int* type_tmp = type;

    x = x_copy;
    v = v_copy;
    type = type_copy;
    x_copy = x_tmp;
    v_copy = v_tmp;
    type_copy = type_tmp;
  }
  //#pragma omp barrier
}
