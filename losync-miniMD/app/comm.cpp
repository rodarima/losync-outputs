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
#include "stdlib.h"
#include "mpi.h"
#include "comm.h"

#ifdef USE_TAMPI
#include <TAMPI.h>
#else
// Dummy definitions for tasked versions of communicate, borders, exchange
int TAMPI_Iwait(MPI_Request*, MPI_Status*) { return 0; }
int TAMPI_Iwaitall(int, MPI_Request*, MPI_Status*) { return 0; }
#endif

// DSM Added for memcpy
#include <string.h>

#define BUFFACTOR 1.5
#define BUFMIN 10000000
#define BUFEXTRA 100
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

// DSM Helpers to keep neighbour directions consistent
// Box Layers:
#define SAME_LAYER 0
#define LOWER_LAYER 1
#define UPPER_LAYER 2
// Special case layer: used to terminate loop over atoms in borders()
#define NO_LAYER -1
// Neighbours
#define UP 0
#define TOP_RIGHT 1
#define RIGHT 2
#define BOTTOM_RIGHT 3
#define DOWN 4
#define BOTTOM_LEFT 5
#define LEFT 6
#define TOP_LEFT 7
// Special case neighbours: for internal send in the y-direction
#define ME -1

// DSM Helper for selecting send/recv buffers and communicators for each function
#define COMMUNICATE_FUNCTION 0
#define EXCHANGE_FUNCTION 1
#define BORDERS_FUNCTION 2

Comm::Comm(int boxes_per_process_)
{
  // DSM: Constant per box
  check_safeexchange = 0;
  do_safeexchange = 0;

  // DSM Multibox
  boxes_per_process = boxes_per_process_;
  boxBufs = (BoxBufs*)malloc(boxes_per_process * sizeof(BoxBufs));

  // DSM Allocate initial send and recv buffers per box
  // TODO: Remove all this boilerplate. Replace with initialiseBuffer() method?
  for (int i = 0; i < boxes_per_process; ++i) {
    // 3 sets of buffers per function that performs communication
    for (int function_index = 0; function_index < 3; ++function_index) {
      // Internal buffers for communication between boxes on this process
      boxBufs[i].internal_buf_send_up[function_index].buf = (MMD_float *) malloc((BUFMIN + BUFMIN) * sizeof(MMD_float));
      boxBufs[i].internal_buf_send_up[function_index].maxsize = BUFMIN;
      boxBufs[i].internal_buf_send_up[function_index].natoms = 0;
      boxBufs[i].internal_buf_send_up[function_index].pbc_any = 0;
      boxBufs[i].internal_buf_send_up[function_index].pbc_x = 0.0;
      boxBufs[i].internal_buf_send_up[function_index].pbc_y = 0.0;
      boxBufs[i].internal_buf_send_up[function_index].pbc_z = 0.0;

      boxBufs[i].internal_buf_send_down[function_index].buf = (MMD_float *) malloc((BUFMIN + BUFMIN) * sizeof(MMD_float));
      boxBufs[i].internal_buf_send_down[function_index].maxsize = BUFMIN;
      boxBufs[i].internal_buf_send_down[function_index].natoms = 0;
      boxBufs[i].internal_buf_send_down[function_index].pbc_any = 0;
      boxBufs[i].internal_buf_send_down[function_index].pbc_x = 0.0;
      boxBufs[i].internal_buf_send_down[function_index].pbc_y = 0.0;
      boxBufs[i].internal_buf_send_down[function_index].pbc_z = 0.0;

      boxBufs[i].internal_buf_recv_up[function_index].buf = (MMD_float *) malloc((BUFMIN + BUFMIN) * sizeof(MMD_float));
      boxBufs[i].internal_buf_recv_up[function_index].maxsize = BUFMIN;
      boxBufs[i].internal_buf_recv_up[function_index].natoms = 0;
      boxBufs[i].internal_buf_recv_up[function_index].pbc_any = 0;
      boxBufs[i].internal_buf_recv_up[function_index].pbc_x = 0.0;
      boxBufs[i].internal_buf_recv_up[function_index].pbc_y = 0.0;
      boxBufs[i].internal_buf_recv_up[function_index].pbc_z = 0.0;

      boxBufs[i].internal_buf_recv_down[function_index].buf = (MMD_float *) malloc((BUFMIN + BUFMIN) * sizeof(MMD_float));
      boxBufs[i].internal_buf_recv_down[function_index].maxsize = BUFMIN;
      boxBufs[i].internal_buf_recv_down[function_index].natoms = 0;
      boxBufs[i].internal_buf_recv_down[function_index].pbc_any = 0;
      boxBufs[i].internal_buf_recv_down[function_index].pbc_x = 0.0;
      boxBufs[i].internal_buf_recv_down[function_index].pbc_y = 0.0;
      boxBufs[i].internal_buf_recv_down[function_index].pbc_z = 0.0;

      // 3 layers of boxes * 8 neighbours per layer
      for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {
        for (int box_neigh_index = 0; box_neigh_index < 8; ++box_neigh_index) {
          boxBufs[i].bufs_send[function_index][box_layer_index][box_neigh_index].buf = (MMD_float *) malloc(
                  (BUFMIN + BUFMIN) * sizeof(MMD_float));
          boxBufs[i].bufs_send[function_index][box_layer_index][box_neigh_index].maxsize = BUFMIN;
          boxBufs[i].bufs_send[function_index][box_layer_index][box_neigh_index].natoms = 0;
          boxBufs[i].bufs_send[function_index][box_layer_index][box_neigh_index].pbc_any = 0;
          boxBufs[i].bufs_send[function_index][box_layer_index][box_neigh_index].pbc_x = 0.0;
          boxBufs[i].bufs_send[function_index][box_layer_index][box_neigh_index].pbc_y = 0.0;
          boxBufs[i].bufs_send[function_index][box_layer_index][box_neigh_index].pbc_z = 0.0;

          boxBufs[i].bufs_recv[function_index][box_layer_index][box_neigh_index].buf = (MMD_float *) malloc(BUFMIN * sizeof(MMD_float));
          boxBufs[i].bufs_recv[function_index][box_layer_index][box_neigh_index].maxsize = BUFMIN;
          boxBufs[i].bufs_recv[function_index][box_layer_index][box_neigh_index].natoms = 0;
          boxBufs[i].bufs_recv[function_index][box_layer_index][box_neigh_index].pbc_any = 0;
          boxBufs[i].bufs_recv[function_index][box_layer_index][box_neigh_index].pbc_x = 0.0;
          boxBufs[i].bufs_recv[function_index][box_layer_index][box_neigh_index].pbc_y = 0.0;
          boxBufs[i].bufs_recv[function_index][box_layer_index][box_neigh_index].pbc_z = 0.0;
        }
      }
    }

    // Buffers used exclusively in communicate meighbour tasks implmentation using a single message per neighbour
    // Send
    boxBufs[i].buf_send_single.buf = (MMD_float *) malloc((BUFMIN + BUFMIN) * sizeof(MMD_float));
    boxBufs[i].buf_send_single.maxsize = BUFMIN;
    boxBufs[i].buf_send_single.natoms = 0;
    boxBufs[i].buf_send_single.pbc_any = 0;
    boxBufs[i].buf_send_single.pbc_x = 0.0;
    boxBufs[i].buf_send_single.pbc_y = 0.0;
    boxBufs[i].buf_send_single.pbc_z = 0.0;
    // Recv
    boxBufs[i].buf_recv_single.buf = (MMD_float *) malloc((BUFMIN + BUFMIN) * sizeof(MMD_float));
    boxBufs[i].buf_recv_single.maxsize = BUFMIN;
    boxBufs[i].buf_recv_single.natoms = 0;
    boxBufs[i].buf_recv_single.pbc_any = 0;
    boxBufs[i].buf_recv_single.pbc_x = 0.0;
    boxBufs[i].buf_recv_single.pbc_y = 0.0;
    boxBufs[i].buf_recv_single.pbc_z = 0.0;

    // DSM set initial values for attributes
    boxBufs[i].maxthreads = 0;
    boxBufs[i].maxnlocal = 0;
  }
}

Comm::~Comm() {
  // DSM Multibox: Free any memory allocated for communicator attributes
  for (int i = 0; i < boxes_per_process; ++i) {
    for (int function_index = 0; function_index < 3; ++function_index) {
      free(boxBufs[i].internal_buf_send_up[function_index].buf);
      free(boxBufs[i].internal_buf_send_down[function_index].buf);
      free(boxBufs[i].internal_buf_recv_up[function_index].buf);
      free(boxBufs[i].internal_buf_recv_down[function_index].buf);

      for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {
        for (int box_neigh_index = 0; box_neigh_index < 8; ++box_neigh_index) {
          free(boxBufs[i].bufs_send[function_index][box_layer_index][box_neigh_index].buf);
          free(boxBufs[i].bufs_recv[function_index][box_layer_index][box_neigh_index].buf);
        }
      }
    }
    free(boxBufs[i].buf_send_single.buf);
    free(boxBufs[i].buf_recv_single.buf);
  }
  free(boxBufs);

  // Task sentinels
  free(initialIntegrateSentinels);
  free(sortSentinels);
  free(forceComputeSentinels);
  free(communicateSentinels);
  free(communicateInternalPackSentinels);
  free(communicateInternalUnpackSentinels);
  free(exchangeSentinels);
  free(exchangePackSentinels);
  free(exchangeInternalSendSentinels);
  free(exchangeInternalRecvSentinels);
  free(exchangePBCSentinels);
  free(bordersSentinels);
  free(bordersPackSentinels);
  free(bordersUnpackSentinels);
  free(bordersSendSentinels);
  free(bordersRecvSentinels);
  free(bordersInternalSentinels);
  free(bordersInternalSendSentinels);
  free(neighbourBuildSentinels);
  for (int i = 0; i < boxes_per_process; ++i) {
    free(communicatePackSentinels[i]);
    free(communicateSendSentinel[i]);
    free(communicateRecvSentinels[i]);
  }
  free(communicatePackSentinels);
  free(communicateSendSentinel);
  free(communicateRecvSentinels);
  for (int i = 0; i < boxes_per_process; ++i) {
    free(exchangeSend1Sentinels[i]);
    free(exchangeRecv1Sentinels[i]);
    free(exchangeSend2Sentinels[i]);
    free(exchangeRecv2Sentinels[i]);
  }
  free(exchangeSend1Sentinels);
  free(exchangeRecv1Sentinels);
  free(exchangeSend2Sentinels);
  free(exchangeRecv2Sentinels);
}

// DSM new function for freeing box communicators. Cannot add this to destructor as would result in MPI_Comm_free being
// called after MPI_Finalize.
/*void Comm::free_box_comms() {
  for (int function_index=0; function_index < 3; ++function_index) {
    if (box_comms[function_index]) {
      for (int i = 0; i < boxes_per_process; ++i) {
        MPI_Comm_free(&box_comms[function_index][i]);
      }
      free(box_comms[function_index]);
    }
  }
}
*/

/* setup spatial-decomposition communication patterns */

// Helper function for mapping dimension & swap number to neighbour direction
// i.e. Maps from idmin and swap 0/1 to UP, TOP_RIGHT, etc. for using in exchange communication
int dim_to_neigh(int idim, int swapnum) {
  // idim can be 0:4 indicating: x, y, z, diagonal-right, diagonal-left
  // swapnum can be 0:1 indicating -ve direction neighbour or +ve direction neighbour
  switch(idim) {
    // x direction
    case 0:
      return swapnum? RIGHT : LEFT;
    // y direction. Special case: internal copy between boxes
    case 1:
      return swapnum? ME : ME;
    // z direction
    case 2:
      return swapnum? UP : DOWN;
    // diagonal-right
    case 3:
      return swapnum? TOP_RIGHT : BOTTOM_LEFT;
    // diagonal-left
    case 4:
      return swapnum? TOP_LEFT : BOTTOM_RIGHT;
  }

  // Something went wrong, we should only have 5 dimensions
  printf("ERROR: Attempted to map dimension %d to a neighbour direction and failed!", idim);
  MPI_Abort(MPI_COMM_WORLD, 1);
  return -1;
}

// Helper function for determining box IDs to use at each layer of communication
void layer_to_targets(Atom* atom, int box_layer_index, int* send_target_box_id, int* recv_target_box_id) {
  // Always the same box IDs for each swap on a layer
  if (box_layer_index == 0) {
    // Own layer - send/recv from boxes with the same ID
    *send_target_box_id = atom->box_id;
    *recv_target_box_id = atom->box_id;
  } else if (box_layer_index == 1) {
    // Lower layer - send to boxes below and receive from above
    *send_target_box_id = atom->boxneigh_negative;
    *recv_target_box_id = atom->boxneigh_positive;
  } else {
    // Upper layer - send to boxes above and receive from below
    *send_target_box_id = atom->boxneigh_positive;
    *recv_target_box_id = atom->boxneigh_negative;
  }
}

int Comm::setup(MMD_float cutneigh, int nprocsx, int nprocsz, Atom* atoms[], int nonblocking_enabled)
{
  int i;
  int nprocs;
  int periods[3];
  MMD_float prd[3];
  int myloc[3];
  MPI_Comm cartesian;
  int ineed, idim, nbox;

  // Sentinels for tasking
  // Integrate
  initialIntegrateSentinels = (char*)malloc(atoms[0]->boxes_per_process);
  sortSentinels = (char*)malloc(atoms[0]->boxes_per_process);
  neighbourBuildSentinels = (char*)malloc(atoms[0]->boxes_per_process);
  forceComputeSentinels = (char*)malloc(atoms[0]->boxes_per_process);
  // Communicate
  communicateSentinels = (char*)malloc(atoms[0]->boxes_per_process);
  communicateInternalPackSentinels = (char*)malloc(atoms[0]->boxes_per_process);
  communicateInternalUnpackSentinels = (char*)malloc(atoms[0]->boxes_per_process);
  exchangeSentinels = (char*)malloc(atoms[0]->boxes_per_process);
  exchangePackSentinels = (char*)malloc(atoms[0]->boxes_per_process);
  exchangeInternalSendSentinels = (char*)malloc(atoms[0]->boxes_per_process);
  exchangeInternalRecvSentinels = (char*)malloc(atoms[0]->boxes_per_process);
  exchangePBCSentinels = (char*)malloc(atoms[0]->boxes_per_process);
  bordersSentinels = (char*)malloc(atoms[0]->boxes_per_process);
  bordersPackSentinels = (char*)malloc(atoms[0]->boxes_per_process);
  bordersSendSentinels = (char*)malloc(atoms[0]->boxes_per_process);
  bordersUnpackSentinels = (char*)malloc(atoms[0]->boxes_per_process);
  bordersInternalSentinels = (char*)malloc(atoms[0]->boxes_per_process);
  bordersInternalSendSentinels = (char*)malloc(atoms[0]->boxes_per_process);
  // For fully-tasked version of communicate (separate pack, send, recv, unpack tasks per loop iteration, per box)
  communicatePackSentinels = (char**)malloc(atoms[0]->boxes_per_process * sizeof(char*));
  communicateSendSentinel = (char**)malloc(atoms[0]->boxes_per_process * sizeof(char*));
  communicateRecvSentinels = (char**)malloc(atoms[0]->boxes_per_process * sizeof(char*));
  for (int i=0; i<atoms[0]->boxes_per_process; ++i) {
    // 3 layers of boxes * 10 neighbours per layer, including internal communication between boxes in y-dimension (these are skipped)
    communicatePackSentinels[i] = (char*)malloc(3 * 10);
    communicateSendSentinel[i] = (char*)malloc(3 * 10);
    communicateRecvSentinels[i] = (char*)malloc(3 * 10);
  }
  // For fully-tasked version of exchange
  exchangeSend1Sentinels = (char**)malloc(atoms[0]->boxes_per_process * sizeof(char*));
  exchangeRecv1Sentinels = (char**)malloc(atoms[0]->boxes_per_process * sizeof(char*));
  exchangeSend2Sentinels = (char**)malloc(atoms[0]->boxes_per_process * sizeof(char*));
  exchangeRecv2Sentinels = (char**)malloc(atoms[0]->boxes_per_process * sizeof(char*));
  for (int i=0; i<atoms[0]->boxes_per_process; ++i) {
    // Loop over 5 dimensions in exchange
    exchangeSend1Sentinels[i] = (char*)malloc(5);
    exchangeRecv1Sentinels[i] = (char*)malloc(5);
    exchangeSend2Sentinels[i] = (char*)malloc(5);
    exchangeRecv2Sentinels[i] = (char*)malloc(5);
  }
  // Fully-tasked borders
  // 3d array: 5 dimensions in borders, 2 swaps per dimension, all per box
  bordersRecvSentinels = (char(*)[5][2])(malloc(atoms[0]->boxes_per_process * sizeof(*bordersRecvSentinels)));

  // Enable/disable non-blocking mode
  this->nonblocking_enabled = nonblocking_enabled;

  // Orderings of 24 buffers per box (one send and one recv buffer per exchange)
  int box_layers[] = {SAME_LAYER, LOWER_LAYER, UPPER_LAYER};
  int neighbours[] = {UP, TOP_RIGHT, RIGHT, BOTTOM_RIGHT, DOWN, BOTTOM_LEFT, LEFT, TOP_LEFT};

  // DSM Box dimensions. Product of units and lattice in create_box(
  // DSM: Multibox change. All boxes have the same prd dimensions and we assume we have at least one box.
  prd[0] = atoms[0]->box.xprd;
  prd[1] = atoms[0]->box.yprd;
  prd[2] = atoms[0]->box.zprd;

  /* setup 3-d grid of procs */
  // DSM: Could pull this from threads.mpi_me and threads.mpi_num_threads to save an MPI call.
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  MMD_float area[3];

  // DSM: entire simulation box dimensions: x*y, x*z, y*z
  area[0] = prd[0] * prd[1];
  area[1] = prd[0] * prd[2];
  area[2] = prd[1] * prd[2];

  MMD_float bestsurf = 2.0 * (area[0] + area[1] + area[2]);

  // loop thru all possible factorizations of nprocs
  // surf = surface area of a proc sub-domain
  // for 2d, insure ipz = 1

  int ipx, ipy, ipz, nremain;
  MMD_float surf;

  ipx = 1;

  // DSM: TODO: This generates the process grid. Need to change for multibox as we're ignoring the outcome
  // DSM: Looks like it tries every possible process grid "loop thru all possible factorizations of nprocs" until it
  // finds the smallest surface (bestsurf), where "surf = surface area of a proc sub-domain"
  while(ipx <= nprocs) {
    if(nprocs % ipx == 0) {
      nremain = nprocs / ipx;
      ipy = 1;

      while(ipy <= nremain) {
        if(nremain % ipy == 0) {
          ipz = nremain / ipy;
          surf = area[0] / ipx / ipy + area[1] / ipx / ipz + area[2] / ipy / ipz;

          if(surf < bestsurf) {
            bestsurf = surf;
            procgrid[0] = ipx;
            procgrid[1] = ipy;
            procgrid[2] = ipz;
          }
        }

        ipy++;
      }
    }

    ipx++;
  }

  // DSM: Multiple boxes change. Ignoring results from above factorisation and hardcoding a grid with one column of box
  // space per process.
  // TODO: Fix factorisation to avoid need for user to manually specify process grid in input file
  procgrid[0] = nprocsx; procgrid[1] = 1; procgrid[2] = nprocsz;

  // DSM: Could this happen? Is this a redundant assert?
  if(procgrid[0]*procgrid[1]*procgrid[2] != nprocs) {
    printf("ERROR: Bad grid of processors\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // DSM Set equivalent boxgrid (# boxes per dimension) based on procgrid.
  // DSM This will always be equal to procgrid in x and z dimensions as we are splitting the grid over boxes in y.
  boxgrid[0] = procgrid[0];
  boxgrid[1] = boxes_per_process;
  boxgrid[2] = procgrid[2];

  /* determine where I am and my neighboring procs in 3d grid of procs */

  // DSM: COMM_WORLD rank order preserved in new communicator. All dimensions are periodic, i.e. coordinate 0 is a
  // neighbour of coordinate nmax ("it circles around").
  int reorder = 0;
  periods[0] = periods[1] = periods[2] = 1;

  MPI_Cart_create(MPI_COMM_WORLD, 3, procgrid, periods, reorder, &cartesian);
  // DSM: Retrieving myloc: array of calling process's coordinates in the new communicator
  MPI_Cart_get(cartesian, 3, procgrid, periods, myloc);
  // DSM: Getting ranks of neighbours in all dimensions
  // We have 8 process neighbours * 3 layers of boxes + 2 internal communications between boxes on this rank = 26-way
  // Process neighbours are constant for all boxes on this rank. Neighbouring box IDs tracked in atom->boxneigh
  // Sides:
  // Due to procgrid[1] being fixed to 1, shift in direction 1 will always == me.
  MPI_Cart_shift(cartesian, 0, 1, &procneigh[0][0], &procneigh[0][1]); // DSM: Left/Right
  MPI_Cart_shift(cartesian, 1, 1, &procneigh[1][0], &procneigh[1][1]); // DSM: Down/Up
  MPI_Cart_shift(cartesian, 2, 1, &procneigh[2][0], &procneigh[2][1]); // DSM: Front/Back
  // Corners:
  int cornerLoc[3];
  cornerLoc[0] = myloc[0]-1; cornerLoc[1] = myloc[1]; cornerLoc[2] = myloc[2]-1;
  MPI_Cart_rank(cartesian, cornerLoc, &procneigh[3][0]); // DSM: "bottom-left"
  cornerLoc[0] = myloc[0]+1; cornerLoc[1] = myloc[1]; cornerLoc[2] = myloc[2]+1;
  MPI_Cart_rank(cartesian, cornerLoc, &procneigh[3][1]); // DSM: "top-right"
  cornerLoc[0] = myloc[0]+1; cornerLoc[1] = myloc[1]; cornerLoc[2] = myloc[2]-1;
  MPI_Cart_rank(cartesian, cornerLoc, &procneigh[4][0]); // DSM: "bottom-right"
  cornerLoc[0] = myloc[0]-1; cornerLoc[1] = myloc[1]; cornerLoc[2] = myloc[2]+1;
  MPI_Cart_rank(cartesian, cornerLoc, &procneigh[4][1]); // DSM: "top-left"

  /* lo/hi = my local box bounds */
  // DSM: Decomposing the experiment box across processes
  // DSM: Multibox change. Different lo/hi entries for each box on this process
  MMD_float perprocess_lo, perprocess_hi;
  for (i=0; i < boxes_per_process; ++i) {
    // DSM: We're only splitting the y-dimension into boxes. In x and z, we keep the dimensions of the current process
    atoms[i]->box.xlo = myloc[0] * prd[0] / procgrid[0];
    atoms[i]->box.xhi = (myloc[0] + 1) * prd[0] / procgrid[0];

    // DSM: Divide each "side" of the box between processes...
    perprocess_lo = myloc[1] * prd[1] / procgrid[1];
    perprocess_hi = (myloc[1] + 1) * prd[1] / procgrid[1];
    // ...then further divide into multiple boxes on each process.
    // DSM TODO: Simplify this.
    atoms[i]->box.ylo = perprocess_lo + (i * (perprocess_hi-perprocess_lo) / boxes_per_process);
    atoms[i]->box.yhi = perprocess_lo + ((i + 1) * (perprocess_hi-perprocess_lo) / boxes_per_process);

    atoms[i]->box.zlo = myloc[2] * prd[2] / procgrid[2];
    atoms[i]->box.zhi = (myloc[2] + 1) * prd[2] / procgrid[2];

    // DSM: Set own local box ID and IDs of neighbouring boxes (processes owning these boxes are tracked in procneigh)
    atoms[i]->box_id = i;
    // DSM Neighbouring box IDs are in range 0->boxes_per_process-1 and are periodic
    int positivedir, negativedir;
    // Box below and box above
    negativedir = atoms[i]->box_id - 1;
    positivedir = atoms[i]->box_id + 1;
    // Fix boundaries if outside range of boxes per process
    if (negativedir < 0) { negativedir = boxes_per_process-1; }
    if (positivedir >= boxes_per_process) { positivedir = 0; }
    atoms[i]->boxneigh_negative = negativedir;
    atoms[i]->boxneigh_positive = positivedir;
  }

  /* need = # of boxes I need atoms from in each dimension */

  need[0] = static_cast<int>(cutneigh * procgrid[0] / prd[0] + 1);
  //need[1] = static_cast<int>(cutneigh * procgrid[1] / prd[1] + 1);
  need[1] = static_cast<int>(cutneigh * boxes_per_process / prd[1] + 1);
  need[2] = static_cast<int>(cutneigh * procgrid[2] / prd[2] + 1);
  // Added for multibox (corner dimensions
  need[3] = 1;
  need[4] = 1;

  // DSM The "need" calculations above allow for >3^3 communications (i.e. >=2 communications in each dimension).
  // Choosing to support nearest neighbour communications only, i.e. need[0] == need[1] == need[2] == 1.
  if (need[0] != 1 || need[1] != 1 || need[2] != 1) {
      printf("ERROR: # boxes Rank %d needs atoms from in each dimension calculated as %d, %d, %d.\n"
             "Only communication with nearest neighbours supported in taskified version.\n"
             "Atoms cannot move >=2 neighbours away in a single round of communication.\n", me, need[0], need[1], need[2]);
      MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // DSM Multibox: When swapping atoms, we need to know both which process rank we are communicating with, which box we
  // are sending into and which box we are receiving into:
  //   Process rank tracked by procneigh.
  //   Box we are sending to: atom->boxneigh.
  //   Box we are receiving into: atom->box_id (we receive into our own box)
  // Could use tag to encode both sending-to and receiving-into box IDs but would need an encoding scheme and would
  // potentially run into maximum tag value limitations.
  // Solution:
  //   Create a separate MPI communicator per box (multiple clones of COMM_WORLD in practice)
  //   Every MPI_Send communication uses own box_id as tag and send target box_id as index into communicator duplicates array
  //   Every MPI_Recv uses receiver box_id as tag and own box_id as communicator index
  // Tag value now only in range 0->boxes_per_process-1, i.e. only an issue when >32768 boxes per process (upper bound
  // guaranteed by the MPI standard)
  // Need 3 independent sets of communicators (one per function) to avoid matching borders() messages with exchange(),
  // for example
  for (int box_index=0; box_index < boxes_per_process; ++box_index) {
    for (int function_index=0; function_index < 3; ++function_index) {
      MPI_Comm_dup(MPI_COMM_WORLD, &boxBufs[box_index].comm[function_index]);
    }
  }

  /* alloc comm memory */
  // DSM: maxswap 2* as we need to communicate in +ve and -ve directions per dimension
  //int maxswap = 2 * (need[0] + need[1] + need[2]); // DSM: Will almost always be 2*(3) = 6

  // DSM: Hardcoding this as each box does same number of communications:
  // Need to communicate with 8 neighbours + 2 internal communications in y-axis (box above and box below) = 10 swaps
  // 8 neighbours * 3 box layers (box_id-1, box_id, box_id+1) + 2 internal swaps = each process doing 26-way communication
  int maxswap = 10;

  // DSM Multibox: Duplicate all these per box (could possibly share some between boxes in a future optimisation)
  for (i=0; i < boxes_per_process; ++i) {
    // DSM TODO: These don't appear to be free()-ed at any point. Add to destructor?

    // Slab arrays replaced by xneg_slab_lo, etc. constants in atom.box
    //boxBufs[i].slablo = (MMD_float*) malloc(maxswap * sizeof(MMD_float));  // bounds of slabs to send to other procs
    //boxBufs[i].slabhi = (MMD_float*) malloc(maxswap * sizeof(MMD_float));  // as above

    // PBC arrays replaced by per buffer flags
    //boxBufs[i].pbc_any = (int*) malloc(maxswap * sizeof(int));             // whether any PBC on this swap (Periodic Boundary Condition)
    //boxBufs[i].pbc_flagx = (int*) malloc(maxswap * sizeof(int));           // PBC correction in x for this swap
    //boxBufs[i].pbc_flagy = (int*) malloc(maxswap * sizeof(int));           // same in y
    //boxBufs[i].pbc_flagz = (int*) malloc(maxswap * sizeof(int));           // same in z

    boxBufs[i].sendproc = (int*) malloc(maxswap * sizeof(int));        // proc to send to at each swap
    boxBufs[i].recvproc = (int*) malloc(maxswap * sizeof(int));        // proc to recv from at each swap
    boxBufs[i].sendproc_exc = (int*) malloc(maxswap * sizeof(int));    // proc to send to at each swap for safe exchange DSM: TODO: Ignoring this for multibox at the moment
    boxBufs[i].recvproc_exc = (int*) malloc(maxswap * sizeof(int));    // proc to recv from at each swap for safe exchange DSM: TODO: Ignoring this for multibox at the moment
    boxBufs[i].sendneigh = (int *) malloc(maxswap * sizeof(int));      // neighbour/direction to send to at each swap
    boxBufs[i].recvneigh = (int *) malloc(maxswap * sizeof(int));      // neighbour/direction to recv from at each swap

    // Duplicate these per box layer
    for (int box_layer_index = 0; box_layer_index<3; ++box_layer_index) {
      boxBufs[i].sendnum[box_layer_index] = (int *) malloc(maxswap * sizeof(int));             // # of atoms to send in each swap
      boxBufs[i].recvnum[box_layer_index] = (int *) malloc(maxswap * sizeof(int));             // # of atoms to recv in each swap
      boxBufs[i].firstrecv[box_layer_index] = (int *) malloc(maxswap * sizeof(int));           // where to put 1st recv atom in each swap
      boxBufs[i].comm_send_size[box_layer_index] = (int *) malloc(maxswap * sizeof(int));      // # of values to send in each comm
      boxBufs[i].comm_recv_size[box_layer_index] = (int *) malloc(maxswap * sizeof(int));      // # of values to recv in each comm
      boxBufs[i].reverse_send_size[box_layer_index] = (int *) malloc(maxswap * sizeof(int));   // # of values to send in each reverse
      boxBufs[i].reverse_recv_size[box_layer_index] = (int *) malloc(maxswap * sizeof(int));   // # of values to recv in each reverse
    }
  }
  int iswap = 0;

  // DSM: TODO: isn't this just exactly the same as procneigh but flattened into 1d?
  // DSM: TODO: Broken, but we can ignore this for multibox for now. sendproc_exc/recvproc_exc are only used in exchange_all, which is only called if safeexchange is enabled.
  for(int idim = 0; idim < 3; idim++) {
    // DSM: Change to get multibox to compile. This branch still broken.
    for (int box_index=0; box_index < boxes_per_process; ++box_index) {
      iswap = 0;
      for (int i = 1; i <= need[idim]; i++, iswap += 2) {
        MPI_Cart_shift(cartesian, idim, i, &boxBufs[box_index].sendproc_exc[iswap], &boxBufs[box_index].sendproc_exc[iswap + 1]);
        MPI_Cart_shift(cartesian, idim, i, &boxBufs[box_index].recvproc_exc[iswap + 1], &boxBufs[box_index].recvproc_exc[iswap]);
      }
    }
  }

  MPI_Comm_free(&cartesian);

  // DSM Multibox: Duplicate these per box, per buffer
  // Definitions:
  //   firstrecv         index of the first atom in the x array in the block we're receiving into
  //   sendlist          list of indices into x array representing the atoms to send
  for (int box_index=0; box_index < boxes_per_process; ++box_index) {
    for (int layer_index = 0; layer_index < 3; ++layer_index) {
      for (int neighbour_index = 0; neighbour_index < 8; ++neighbour_index) {
        // Which buffer are we modifying?
        int layer = box_layers[layer_index];
        int neighbour = neighbours[neighbour_index];
        boxBufs[box_index].sendlists[layer][neighbour] = (int *) malloc(BUFMIN * sizeof(int));
        boxBufs[box_index].maxsendlists[layer][neighbour] = BUFMIN;
      } // End of loop over neighbours
    } // End of loop over layers

    // Allocate lists for special case buffers for internal communication between boxes on this process
    // Send
    boxBufs[box_index].internal_sendlist_up = (int *) malloc(BUFMIN * sizeof(int));
    boxBufs[box_index].internal_maxsendlist_up = BUFMIN;
    boxBufs[box_index].internal_sendlist_down = (int *) malloc(BUFMIN * sizeof(int));
    boxBufs[box_index].internal_maxsendlist_down = BUFMIN;
  } // End of loop over boxes

  /* setup 4 parameters for each exchange: (spart,rpart,slablo,slabhi)
     sendproc(nswap) = proc to send to at each swap
     recvproc(nswap) = proc to recv from at each swap
     slablo/slabhi(nswap) = slab boundaries (in correct dimension) of atoms
                            to send at each swap
     1st part of if statement is sending to the west/south/down
     2nd part of if statement is sending to the east/north/up
     nbox = atoms I send originated in this box */

  /* set commflag if atoms are being exchanged across a box boundary
     commflag(idim,nswap) =  0 -> not across a boundary
                          =  1 -> add box-length to position when sending
                          = -1 -> subtract box-length from pos when sending */

  // DSM: Variables named above seem to have been renamed in the actual code:
  // spart/rpart -> sendproc/recvproc
  // commflag -> pbc_any and pbc_flag[x/y/z]

  // DSM: Remove slablo/slabhi and replaced with atom.box.xneg_slab_lo, etc. to facilitate 26-way communication

  for (int box_index = 0; box_index < boxes_per_process; ++box_index) {

    // DSM Set up this box's slab boundary constants for use in comm.borders()
    // "Slab boundaries" are regions of the box: lo to lo+cutneigh, and hi-cutneigh to hi.
    // Any atoms within these regions have their positions updated with neighbours each timestep.
    // Relevant functions: comm.borders() and comm.communicate()
    // "Left": X -ve direction
    atoms[box_index]->box.xneg_slab_lo = atoms[box_index]->box.xlo;
    atoms[box_index]->box.xneg_slab_hi = atoms[box_index]->box.xlo + cutneigh;
    // "Right": X +ve direction:
    atoms[box_index]->box.xpos_slab_lo = atoms[box_index]->box.xhi - cutneigh;
    atoms[box_index]->box.xpos_slab_hi = atoms[box_index]->box.xhi;
    // "Down": Y -ve direction
    atoms[box_index]->box.yneg_slab_lo = atoms[box_index]->box.ylo;
    atoms[box_index]->box.yneg_slab_hi = atoms[box_index]->box.ylo + cutneigh;
    // "Up": Y +ve direction
    atoms[box_index]->box.ypos_slab_lo = atoms[box_index]->box.yhi - cutneigh;
    atoms[box_index]->box.ypos_slab_hi = atoms[box_index]->box.yhi;
    // "Back": Z -ve direction
    atoms[box_index]->box.zneg_slab_lo = atoms[box_index]->box.zlo;
    atoms[box_index]->box.zneg_slab_hi = atoms[box_index]->box.zlo + cutneigh;
    // "Front": Z +ve direction
    atoms[box_index]->box.zpos_slab_lo = atoms[box_index]->box.zhi - cutneigh;
    atoms[box_index]->box.zpos_slab_hi = atoms[box_index]->box.zhi;

    // Create a set of buffers for each of the 3 functions that perform communication
    for (int function_index = 0; function_index < 3; ++function_index) {
      // Set Periodic Boundary Conditions (PBCs) on buffers
      // PBC corrections are done on the sender's side => we only need to set for send buffers
      // Testing: am I a periphery process on the boundaries of the process grid in these dimensions?
      // pbc_flags replaced by storing the value for the correction. No longer have to compute the product at every
      // packing
      for (int layer_index = 0; layer_index < 3; ++layer_index) {
        int layer = box_layers[layer_index];
        // -ve X boundary
        if (myloc[0] == 0) {
          // Left
          boxBufs[box_index].bufs_send[function_index][layer][LEFT].pbc_any = 1;
          boxBufs[box_index].bufs_send[function_index][layer][LEFT].pbc_x = atoms[box_index]->box.xprd;
          // Diagonal Down
          boxBufs[box_index].bufs_send[function_index][layer][BOTTOM_LEFT].pbc_any = 1;
          boxBufs[box_index].bufs_send[function_index][layer][BOTTOM_LEFT].pbc_x = atoms[box_index]->box.xprd;
          // Diagonal Up
          boxBufs[box_index].bufs_send[function_index][layer][TOP_LEFT].pbc_any = 1;
          boxBufs[box_index].bufs_send[function_index][layer][TOP_LEFT].pbc_x = atoms[box_index]->box.xprd;
        }

        // +ve X boundary
        if (myloc[0] == procgrid[0] - 1) {
          // Right
          boxBufs[box_index].bufs_send[function_index][layer][RIGHT].pbc_any = 1;
          boxBufs[box_index].bufs_send[function_index][layer][RIGHT].pbc_x = -atoms[box_index]->box.xprd;
          // Diagonal Down
          boxBufs[box_index].bufs_send[function_index][layer][BOTTOM_RIGHT].pbc_any = 1;
          boxBufs[box_index].bufs_send[function_index][layer][BOTTOM_RIGHT].pbc_x = -atoms[box_index]->box.xprd;
          // Diagonal Up
          boxBufs[box_index].bufs_send[function_index][layer][TOP_RIGHT].pbc_any = 1;
          boxBufs[box_index].bufs_send[function_index][layer][TOP_RIGHT].pbc_x = -atoms[box_index]->box.xprd;
        }

        // Y dimension special case: we are periodic if this is a periphery *box*
        // -ve Y boundary
        if (atoms[box_index]->box_id == 0 && layer == LOWER_LAYER) {
          // Set all buffers periodic
          for (int neighbour_index=0; neighbour_index < 8; ++neighbour_index) {
            int neighbour = neighbours[neighbour_index];
            boxBufs[box_index].bufs_send[function_index][layer][neighbour].pbc_any = 1;
            boxBufs[box_index].bufs_send[function_index][layer][neighbour].pbc_y = atoms[box_index]->box.yprd;
          }
        }

        // +ve Y boundary
        if (atoms[box_index]->box_id == boxes_per_process - 1 && layer == UPPER_LAYER) {
          for (int neighbour_index=0; neighbour_index < 8; ++neighbour_index) {
            int neighbour = neighbours[neighbour_index];
            boxBufs[box_index].bufs_send[function_index][layer][neighbour].pbc_any = 1;
            boxBufs[box_index].bufs_send[function_index][layer][neighbour].pbc_y = -atoms[box_index]->box.yprd;
          }
        }

        // -ve Z boundary
        if (myloc[2] == 0) {
          // Down
          boxBufs[box_index].bufs_send[function_index][layer][DOWN].pbc_any = 1;
          boxBufs[box_index].bufs_send[function_index][layer][DOWN].pbc_z = atoms[box_index]->box.zprd;
          // Diagonal Left
          boxBufs[box_index].bufs_send[function_index][layer][BOTTOM_LEFT].pbc_any = 1;
          boxBufs[box_index].bufs_send[function_index][layer][BOTTOM_LEFT].pbc_z = atoms[box_index]->box.zprd;
          // Diagonal Right
          boxBufs[box_index].bufs_send[function_index][layer][BOTTOM_RIGHT].pbc_any = 1;
          boxBufs[box_index].bufs_send[function_index][layer][BOTTOM_RIGHT].pbc_z = atoms[box_index]->box.zprd;
        }

        // +ve Z Boundary
        if (myloc[2] == procgrid[2] - 1) {
          // Up
          boxBufs[box_index].bufs_send[function_index][layer][UP].pbc_any = 1;
          boxBufs[box_index].bufs_send[function_index][layer][UP].pbc_z = -atoms[box_index]->box.zprd;
          // Diagonal Left
          boxBufs[box_index].bufs_send[function_index][layer][TOP_LEFT].pbc_any = 1;
          boxBufs[box_index].bufs_send[function_index][layer][TOP_LEFT].pbc_z = -atoms[box_index]->box.zprd;
          // Diagonal Right
          boxBufs[box_index].bufs_send[function_index][layer][TOP_RIGHT].pbc_any = 1;
          boxBufs[box_index].bufs_send[function_index][layer][TOP_RIGHT].pbc_z = -atoms[box_index]->box.zprd;
        }

      } // End of loop over buffer layers

      // Set PBCs for internal communication between boxes on this process - Y boundaries
      // No need to check myloc[1]==0 or myloc[1]==procgrid[1]-1 to confirm if this is a periphery process as we have
      // hardcoded procgrid[1]=1.
      // -ve Y boundary
      if (atoms[box_index]->box_id == 0) {
        boxBufs[box_index].internal_buf_send_down[function_index].pbc_any = 1;
        boxBufs[box_index].internal_buf_send_down[function_index].pbc_y = atoms[box_index]->box.yprd;
      }

      // +ve Y boundary
      if (atoms[box_index]->box_id == boxes_per_process - 1) {
        boxBufs[box_index].internal_buf_send_up[function_index].pbc_any = 1;
        boxBufs[box_index].internal_buf_send_up[function_index].pbc_y = -atoms[box_index]->box.yprd;
      }

    } // End of loop over functions

    // DSM Reference all variables/arrays for this box instance
    int& nswap = boxBufs[box_index].nswap;
    int*& sendproc = boxBufs[box_index].sendproc;
    int*& recvproc = boxBufs[box_index].recvproc;
    int*& sendneigh = boxBufs[box_index].sendneigh;
    int*& recvneigh = boxBufs[box_index].recvneigh;
    nswap = 0;

    // DSM: Rewritten loop to be over 5 dimensions rather than 3:
    // x, y, z, diagonal-right, diagonal-left
    // Loop will be executed 10 times, once for each swap: 5 dimensions * 2 swaps per dimension = 10 swaps overall.
    // Dimension y is a special case - an internal data copy between boxes on the same process. Functions comm.borders(),
    // comm.exchange() and comm.communicate() skip the y dimension and perform internal copies outside their main loops
    for(idim = 0; idim < 5; idim++) {
      //int swaps_this_dim = 2 * need[idim];
      int swaps_this_dim = 2; // DSM Multibox: Removed "need" parameter, we have the same number of swaps (always 2) in each dimension (per box)
      for(ineed = 0; ineed < swaps_this_dim; ineed++) { // DSM Multibox: swaps_this_dim guaranteed to be 2 => nswap 0->6 (idim*2)

        if (ineed % 2 == 0) { // DSM: This is the "west/south/down" branch. Now also "diagonal-top"
          // DSM: Send to -ve in this dimension and receive from +ve
          // DSM: No change for multibox - we send to the same processes, just differing box IDs
          sendproc[nswap] = procneigh[idim][0];
          recvproc[nswap] = procneigh[idim][1];
          // DSM Multibox: Changed PBC calculations to be outside swaps loop, flags are now per buffer
        } else { // DSM: "east/north/up". Also "diagonal-bottom"
          // DSM: Send to +ve in this dimension and receive from -ve
          sendproc[nswap] = procneigh[idim][1];
          recvproc[nswap] = procneigh[idim][0];
        }

        // Determine which buffers we are using in each communication
        // Map from idim/ineed to UP, TOP_RIGHT, etc. neighbour ID
        // Value of ineed swaps between 0 and 1:
        //  ineed=0: we want -ve send-neighbour and +ve recv-neighbour
        //  ineed=1: we want +ve send-neighbour -ve recv-neighbour
        // Second parameter of mapping function is 0 for -ve neighbour and 1 for +ve neighbour, hence 0+ineed/1-ineed
        sendneigh[nswap] = dim_to_neigh(idim, 0 + ineed);
        recvneigh[nswap] = dim_to_neigh(idim, 1 - ineed);

        nswap++;
      }  // DSM: End of loop over swaps per dimension
    } // DSM: End of loop over dimensions
  } // DSM: End of loop over boxes per process

  return 0;
}

/* communication of atom info every timestep */

void Comm::communicate(Atom* atoms[]) {

  // Select between blocking and non-blocking communication modes
  if (nonblocking_enabled) {
    communicate_nonblocking(atoms);
    // Wait for all other communications to finish before doing internal swaps between boxes on same process
    for(int box_index = 0; box_index < atoms[0]->boxes_per_process; ++box_index) {
      communicate_internal(*atoms[box_index], box_index);
    }
  } else {
    // Call once per box instance
    for(int box_index = 0; box_index < atoms[0]->boxes_per_process; ++box_index) {
      //#pragma oss task label(communicate_blocking) in(initialIntegrateSentinels[box_index]) out(communicateSentinels[box_index]) firstprivate(box_index)
      //communicate_blocking_isend(*atoms[box_index], box_index);

      // Tasks created within this functions
      //communicate_blocking_alltasks(atoms[box_index], box_index);
      //communicate_blocking_alltasks_recvfirst(atoms[box_index], box_index);
      //communicate_nonblocking_neighbourtasks(atoms[box_index], box_index);
      //communicate_blocking_neighbourtasks(atoms[box_index], box_index);
      //communicate_blocking_single_message_neighbourtasks(atoms[box_index], box_index); // TODO: BROKEN
      //communicate_nonblocking_alltasks_tampi_iwait(atoms[box_index], box_index);
      communicate_nonblocking_neighbourtasks_tampi_iwaitall(atoms[box_index], box_index);

      // Generate tasks for internal swaps between boxes only if we have some atoms to send
      if (boxBufs[box_index].sendnum[0][2] > 0 && boxBufs[box_index].sendnum[0][3] > 0) {
        // 2 sets of tasks here: the pushing of data to both of this box's neighbours (i.e. "sends") and the unpacking of
        // data done in loop following (i.e. "recvs")
        #pragma oss task label(communicate_internal_pack) \
                       in(initialIntegrateSentinels[box_index]) \
                       out(communicateInternalPackSentinels[box_index]) firstprivate(box_index)
        {
          // Lower layer (usually box_id-1)
          int neighbour = atoms[box_index]->boxneigh_negative;
          communicate_internal_send(atoms[box_index],
                                    &boxBufs[box_index].internal_buf_send_down[COMMUNICATE_FUNCTION],
                                    &boxBufs[neighbour].internal_buf_recv_down[COMMUNICATE_FUNCTION],
                                    boxBufs[box_index].internal_sendlist_down,
                                    boxBufs[box_index].comm_send_size[0][2], boxBufs[box_index].sendnum[0][2]);
          // Upper layer (usually box_id+1)
          neighbour = atoms[box_index]->boxneigh_positive;
          communicate_internal_send(atoms[box_index],
                                    &boxBufs[box_index].internal_buf_send_up[COMMUNICATE_FUNCTION],
                                    &boxBufs[neighbour].internal_buf_recv_up[COMMUNICATE_FUNCTION],
                                    boxBufs[box_index].internal_sendlist_up,
                                    boxBufs[box_index].comm_send_size[0][3], boxBufs[box_index].sendnum[0][3]);
        }
      }
    }

    // Receive/unpack half of internal swaps between boxes
    for(int box_index = 0; box_index < atoms[0]->boxes_per_process; ++box_index) {
      // Wait for both of my neighbours' pack tasks to finish before unpacking (if there is data to unpack)
      if (boxBufs[box_index].recvnum[0][2] > 0 && boxBufs[box_index].recvnum[0][3] > 0) {
        #pragma oss task label(communicate_internal_unpack) \
                     in(communicateInternalPackSentinels[atoms[box_index]->boxneigh_positive]) \
                     in(communicateInternalPackSentinels[atoms[box_index]->boxneigh_negative]) \
                     out(communicateInternalUnpackSentinels[box_index]) firstprivate(box_index)
        {
          // Values of firstrecv are set in borders function: [2] is always y-ve swap and [3] always y+ve.
          // Data has been pushed into my internal buffers. Unpack into end of atom list
          // Unlike other unpack routines, unpack_comm has an internal loop over atoms
          atoms[box_index]->unpack_comm(boxBufs[box_index].recvnum[0][2], boxBufs[box_index].firstrecv[0][2],
                                        boxBufs[box_index].internal_buf_recv_down[COMMUNICATE_FUNCTION].buf);
          // As above but y+ve direction ([3] always y+ve)
          atoms[box_index]->unpack_comm(boxBufs[box_index].recvnum[0][3], boxBufs[box_index].firstrecv[0][3],
                                        boxBufs[box_index].internal_buf_recv_up[COMMUNICATE_FUNCTION].buf);
        }
      } // End of recvnum branch
    } // End of loop over boxes
  } // End of blocking branch
} // End of function

// DSM:
// Calculated for each timestep in borders()+exchange():
// sendnum           number of atoms to send
// sendlist          list of indices into x array representing the atoms to send
// comm_recv_size    number of doubles to receive in a swap (will be 3x sending proc's sendnum)
// comm_send_size    number of doubles to send in a swap (will be 3x this proc's sendnum)
// firstrecv         index of the first atom in the x array in the region we're receiving into (used as offset to x in unpack_comm())
//
// Calculated in setup(), constant for all timesteps:
// nswap             number of swaps to perform
// sendproc          list of process ranks to send to
// recvproc          list of process ranks to receive from
// pbc_any           whether this process+box combination is periodic in any dimension (flag == 0 or 1)
// pbc_flagx         whether this process+box combination is periodic in the x dimension (flag == 0 or 1 or -1)
// pbc_flagy         whether this process+box combination is periodic in the y dimension (flag == 0 or 1 or -1)
// pbc_flagz         whether this process+box combination is periodic in the z dimension (flag == 0 or 1 or -1)
void Comm::communicate_nonblocking(Atom* atoms[])
{

  // Arrays for isend/irecv requests
  // nswaps communications per box layer, 3 box layers per box on this process.
  // Internal swaps in y-dimension are skipped => subtract 6 per box (2 y-swaps per box layer)
  // nswaps and boxes_per_process constant per box
  int nrequests = boxBufs[0].nswap * 3 * atoms[0]->boxes_per_process - 6 * atoms[0]->boxes_per_process;
  MPI_Request recv_requests[nrequests];
  MPI_Request send_requests[nrequests];
  int recv_req_i = 0;
  int send_req_i = 0;

  for(int box_id = 0; box_id < atoms[0]->boxes_per_process; ++box_id) {

    int nsend;
    AtomBuffer *buf_send = NULL, *buf_recv = NULL;
    // ID of boxes we're sending to and receiving from
    int send_target_box_id, recv_target_box_id;
    // Directions of buffers we're using to send to/recieve into
    int sendNeighbour, recvNeighbour;
    // Message tags
    int sendtag, recvtag;
    // Communicators to use
    MPI_Comm *send_comm, *recv_comm;
    // List of atoms to be sent
    int *sendlist;
    // Swap counter for this box
    int swapnum = 0;

    // DSM Multibox: Reference buffers/variables specific to this box.
    Atom& atom = *atoms[box_id];
    int &nswap = boxBufs[box_id].nswap;
    int *(&sendnum)[3] = boxBufs[box_id].sendnum;
    int *(&recvnum)[3] = boxBufs[box_id].recvnum;
    int *(&firstrecv)[3] = boxBufs[box_id].firstrecv;
    int *&sendproc = boxBufs[box_id].sendproc;
    int *&recvproc = boxBufs[box_id].recvproc;
    int *&sendneigh = boxBufs[box_id].sendneigh;
    int *&recvneigh = boxBufs[box_id].recvneigh;
    int *(&comm_send_size)[3] = boxBufs[box_id].comm_send_size;
    int *(&comm_recv_size)[3] = boxBufs[box_id].comm_recv_size;
    AtomBuffer (&bufs_send)[3][8] = boxBufs[box_id].bufs_send[COMMUNICATE_FUNCTION];
    AtomBuffer (&bufs_recv)[3][8] = boxBufs[box_id].bufs_recv[COMMUNICATE_FUNCTION];
    int *(&sendlists)[3][8] = boxBufs[box_id].sendlists;
    AtomBuffer &internal_buf_send_up = boxBufs[box_id].internal_buf_send_up[COMMUNICATE_FUNCTION];
    AtomBuffer &internal_buf_send_down = boxBufs[box_id].internal_buf_send_down[COMMUNICATE_FUNCTION];
    AtomBuffer &internal_buf_recv_up = boxBufs[box_id].internal_buf_recv_up[COMMUNICATE_FUNCTION];
    AtomBuffer &internal_buf_recv_down = boxBufs[box_id].internal_buf_recv_down[COMMUNICATE_FUNCTION];
    int *&internal_sendlist_up = boxBufs[box_id].internal_sendlist_up;
    int *&internal_sendlist_down = boxBufs[box_id].internal_sendlist_down;

    for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {

      // Determine box IDs to send to/receive from on this layer
      layer_to_targets(&atom, box_layer_index, &send_target_box_id, &recv_target_box_id);

      // DSM nswap (number of swaps) is a constant calculated in setup()
      for (int iswap = 0; iswap < nswap; iswap++) {

        // Check if this is an internal memory copy between boxes on same proc. Skip MPI communication if so
        if (iswap == 2) {         // iswap=2 guaranteed to be swap in y -ve direction
          // Only perform the internal swap once, do not repeat for every box layer
          if (box_layer_index > 0) { continue; }
          // First exchange to lower layer (usually box_id-1 but periodic)
          int boxNeighbour = atom.boxneigh_negative;
          buf_send = &internal_buf_send_down;
          buf_recv = &(boxBufs[boxNeighbour].internal_buf_recv_down[COMMUNICATE_FUNCTION]);
          sendlist = internal_sendlist_down;
          // Increase size of receive buffers if packing/receiving more than can be held
          nsend = comm_send_size[box_layer_index][iswap];
          // TODO: Could optimise by packing directly into neighbour's buffer (PBC would complicate things)
          if(nsend > buf_send->maxsize) {
            buf_send->growrecv(nsend);
          }
          if(nsend > buf_recv->maxsize) {
            buf_recv->growrecv(nsend);
          }
          // Pack data into send buffer
          atom.pack_comm(sendnum[box_layer_index][iswap], sendlist, buf_send->buf, buf_send->pbc_any, buf_send->pbc_x,
                         buf_send->pbc_y, buf_send->pbc_z);
          // Push to neighbour.
          memcpy(buf_recv->buf, buf_send->buf, nsend * sizeof(MMD_float));
          continue;
        } else if (iswap == 3) {  // iswap=3 guaranteed to be swap in y +ve direction
          if (box_layer_index > 0) { continue; }
          // As above but to upper layer (usually box_id+1)
          int boxNeighbour = atom.boxneigh_positive;
          buf_send = &internal_buf_send_up;
          buf_recv = &(boxBufs[boxNeighbour].internal_buf_recv_up[COMMUNICATE_FUNCTION]);
          sendlist = internal_sendlist_up;
          nsend = comm_send_size[box_layer_index][iswap];
          if(nsend > buf_send->maxsize) {
            buf_send->growrecv(nsend);
          }
          if(nsend > buf_recv->maxsize) {
            buf_recv->growrecv(nsend);
          }
          atom.pack_comm(sendnum[box_layer_index][iswap], sendlist, buf_send->buf, buf_send->pbc_any, buf_send->pbc_x,
                         buf_send->pbc_y, buf_send->pbc_z);
          memcpy(buf_recv->buf, buf_send->buf, nsend * sizeof(MMD_float));
          continue;
        }

        // Determine which buffers we are using in each communication
        sendNeighbour = sendneigh[iswap];
        recvNeighbour = recvneigh[iswap];
        buf_send = &bufs_send[box_layer_index][sendNeighbour];
        buf_recv = &bufs_recv[box_layer_index][recvNeighbour];
        send_comm = &boxBufs[send_target_box_id].comm[COMMUNICATE_FUNCTION];
        recv_comm = &boxBufs[atom.box_id].comm[COMMUNICATE_FUNCTION];
        sendtag = swapnum;
        recvtag = sendtag;
        ++swapnum;
        sendlist = sendlists[box_layer_index][sendNeighbour];

        /* pack buffer */

        // DSM Added buffer size checks. Previously weren't necessary as buffers were reused from borders() (which ensured
        // ensured sufficient size) but now independent buffers used for each function
        if (comm_send_size[box_layer_index][iswap] > buf_send->maxsize) {
          buf_send->growsend(comm_send_size[box_layer_index][iswap]);
        }
        if (comm_recv_size[box_layer_index][iswap] > buf_recv->maxsize) {
          buf_recv->growrecv(comm_recv_size[box_layer_index][iswap]);
        }

        //#pragma omp barrier
        atom.pack_comm(sendnum[box_layer_index][iswap], sendlist, buf_send->buf, buf_send->pbc_any, buf_send->pbc_x,
                       buf_send->pbc_y, buf_send->pbc_z);

        /* exchange with another proc
           if self, set recv buffer to send buffer */

        // DSM Multibox: Reverse order of this if statement and include box id
        if (sendproc[iswap] == me && send_target_box_id == atom.box_id) {
          // Skip copies to myself, data will be unpacked directly from send buffer in next step
          continue;
        } else {
        //#pragma omp master
          {
            if (sizeof(MMD_float) == 4) {
              MPI_Irecv(buf_recv->buf, comm_recv_size[box_layer_index][iswap], MPI_FLOAT,
                        recvproc[iswap], recvtag, *recv_comm, &recv_requests[recv_req_i]);
              MPI_Isend(buf_send->buf, comm_send_size[box_layer_index][iswap], MPI_FLOAT,
                       sendproc[iswap], sendtag, *send_comm, &send_requests[send_req_i]);
            } else {
              MPI_Irecv(buf_recv->buf, comm_recv_size[box_layer_index][iswap], MPI_DOUBLE,
                        recvproc[iswap], recvtag, *recv_comm, &recv_requests[recv_req_i]);
              MPI_Isend(buf_send->buf, comm_send_size[box_layer_index][iswap], MPI_DOUBLE,
                       sendproc[iswap], sendtag, *send_comm, &send_requests[send_req_i]);
            }

            recv_req_i++;
            send_req_i++;
          }
        }

      } // End of loop over number of swap
    } // End of loop over box layers
  } // End of loop over boxes per processs

  // All swaps posted, wait for completion
  // req_is rather than nrequests as may have skipped requests if sending to self
  MPI_Waitall(recv_req_i, recv_requests, MPI_STATUSES_IGNORE);
  MPI_Waitall(send_req_i, send_requests, MPI_STATUSES_IGNORE);

  /* unpack buffer */

  // Do all unpacks unless this was an internal swap
  for(int box_id = 0; box_id < atoms[0]->boxes_per_process; ++box_id) {

    int neighbour, send_target_box_id, recv_target_box_id;
    // Variables for this box
    AtomBuffer *buf = NULL;
    AtomBuffer (&bufs_send)[3][8] = boxBufs[box_id].bufs_send[COMMUNICATE_FUNCTION];
    AtomBuffer (&bufs_recv)[3][8] = boxBufs[box_id].bufs_recv[COMMUNICATE_FUNCTION];
    Atom &atom = *atoms[box_id];
    int &nswap = boxBufs[box_id].nswap;
    int *&sendproc = boxBufs[box_id].sendproc;
    int *(&recvnum)[3] = boxBufs[box_id].recvnum;
    int *(&firstrecv)[3] = boxBufs[box_id].firstrecv;
    int *&sendneigh = boxBufs[box_id].sendneigh;
    int *&recvneigh = boxBufs[box_id].recvneigh;

    for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {
      // Determine box IDs to send to/receive from on this layer
      layer_to_targets(&atom, box_layer_index, &send_target_box_id, &recv_target_box_id);
      for (int iswap = 0; iswap < nswap; iswap++) {
        // Skip internal swaps - 2/3 are always y-ve and y+ve directions
        if (iswap == 2 || iswap == 3) {
          continue;
        }

        // Which buffer are we unpacking?
        if (sendproc[iswap] == me && send_target_box_id == atom.box_id) {
          // Sending to self. Just unpack the send buffer
          neighbour = sendneigh[iswap];
          buf = &bufs_send[box_layer_index][neighbour];
        } else {
          // Otherwise, unpack the receive buffer as normal
          neighbour = recvneigh[iswap];
          buf = &bufs_recv[box_layer_index][neighbour];
        }

        atom.unpack_comm(recvnum[box_layer_index][iswap], firstrecv[box_layer_index][iswap], buf->buf);
      } // End of loop over number of swap
    } // End of loop over box layers
  } // End of loop over boxes per processs
} // End of function

void Comm::communicate_blocking(Atom &atom, int box_id)
{

  MPI_Request request;
  MPI_Status status;

  int nsend;
  AtomBuffer *buf_send = NULL, *buf_recv = NULL, *buf = NULL;
  // ID of box we're sending to and receiving from
  int send_target_box_id, recv_target_box_id;
  // Directions of buffers we're using to send to/recieve into
  int sendNeighbour, recvNeighbour;
  // Message tags
  int sendtag, recvtag;
  // Communicators to use
  MPI_Comm *send_comm, *recv_comm;
  // List of atoms to be sent
  int *sendlist;
  // Swap counter for this box
  int swapnum = 0;

  // DSM Multibox: Reference buffers/variables specific to this box.
  int& nswap = boxBufs[box_id].nswap;
  int* (&sendnum)[3] = boxBufs[box_id].sendnum;
  int* (&recvnum)[3] = boxBufs[box_id].recvnum;
  int* (&firstrecv)[3] = boxBufs[box_id].firstrecv;
  int*& sendproc = boxBufs[box_id].sendproc;
  int*& recvproc = boxBufs[box_id].recvproc;
  int*& sendneigh = boxBufs[box_id].sendneigh;
  int*& recvneigh = boxBufs[box_id].recvneigh;
  int* (&comm_send_size)[3] = boxBufs[box_id].comm_send_size;
  int* (&comm_recv_size)[3] = boxBufs[box_id].comm_recv_size;
  AtomBuffer (&bufs_send)[3][8] = boxBufs[box_id].bufs_send[COMMUNICATE_FUNCTION];
  AtomBuffer (&bufs_recv)[3][8] = boxBufs[box_id].bufs_recv[COMMUNICATE_FUNCTION];
  int* (&sendlists)[3][8] = boxBufs[box_id].sendlists;
  AtomBuffer& internal_buf_send_up = boxBufs[box_id].internal_buf_send_up[COMMUNICATE_FUNCTION];
  AtomBuffer& internal_buf_send_down = boxBufs[box_id].internal_buf_send_down[COMMUNICATE_FUNCTION];
  AtomBuffer& internal_buf_recv_up = boxBufs[box_id].internal_buf_recv_up[COMMUNICATE_FUNCTION];
  AtomBuffer& internal_buf_recv_down = boxBufs[box_id].internal_buf_recv_down[COMMUNICATE_FUNCTION];
  int*& internal_sendlist_up = boxBufs[box_id].internal_sendlist_up;
  int*& internal_sendlist_down = boxBufs[box_id].internal_sendlist_down;

  for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {

    // Determine box IDs to send to/receive from on this layer
    layer_to_targets(&atom, box_layer_index, &send_target_box_id, &recv_target_box_id);

    // DSM nswap (number of swaps) is a constant calculated in setup()
    for(int iswap = 0; iswap < nswap; iswap++) {

      // Skip MPI communication in y dimension. Internal memory copies between boxes on same proc done in separate task
      if (iswap == 2 || iswap == 3) {         // iswap=2 swap in y -ve direction, iswap=3 swap in y +ve direction
        continue;
      }

      // Determine which buffers we are using in each communication
      sendNeighbour = sendneigh[iswap];
      recvNeighbour = recvneigh[iswap];
      buf_send = &bufs_send[box_layer_index][sendNeighbour];
      buf_recv = &bufs_recv[box_layer_index][recvNeighbour];
      send_comm = &boxBufs[send_target_box_id].comm[COMMUNICATE_FUNCTION];
      recv_comm = &boxBufs[atom.box_id].comm[COMMUNICATE_FUNCTION];
      sendtag = swapnum;
      recvtag = sendtag;
      ++swapnum;
      sendlist = sendlists[box_layer_index][sendNeighbour];

      /* pack buffer */

      // Added buffer size checks. Previously weren't necessary as buffers were reused from borders() (which ensured
      // ensured sufficient size) but now independent buffers used for each function
      if (comm_send_size[box_layer_index][iswap] > buf_send->maxsize) {
        buf_send->growsend(comm_send_size[box_layer_index][iswap]);
      }
      if (comm_recv_size[box_layer_index][iswap] > buf_recv->maxsize) {
        buf_recv->growrecv(comm_recv_size[box_layer_index][iswap]);
      }

      //#pragma omp barrier
      atom.pack_comm(sendnum[box_layer_index][iswap], sendlist, buf_send->buf, buf_send->pbc_any, buf_send->pbc_x,
                     buf_send->pbc_y, buf_send->pbc_z);

      //#pragma omp barrier

      /* exchange with another proc
         if self, set recv buffer to send buffer */

      // DSM Multibox: Reverse order of this if statement and include box id
      if (sendproc[iswap] == me && send_target_box_id == atom.box_id) {
        // Skip MPI communication to myself, data will be unpacked directly from send buffer in next step
        buf = buf_send;
      } else {
        //#pragma omp master
        {
          if(sizeof(MMD_float) == 4) {
            MPI_Irecv(buf_recv->buf, comm_recv_size[box_layer_index][iswap], MPI_FLOAT,
                      recvproc[iswap], recvtag, *recv_comm, &request);
            MPI_Send(buf_send->buf, comm_send_size[box_layer_index][iswap], MPI_FLOAT,
                     sendproc[iswap], sendtag, *send_comm);
          } else {
            MPI_Irecv(buf_recv->buf, comm_recv_size[box_layer_index][iswap], MPI_DOUBLE,
                      recvproc[iswap], recvtag, *recv_comm, &request);
            MPI_Send(buf_send->buf, comm_send_size[box_layer_index][iswap], MPI_DOUBLE,
                     sendproc[iswap], sendtag, *send_comm);
          }

          MPI_Wait(&request, &status);
          buf = buf_recv;
        }
      }

      /* unpack buffer */

      atom.unpack_comm(recvnum[box_layer_index][iswap], firstrecv[box_layer_index][iswap], buf->buf);

    } // End of loop over number of swap
  } // End of loop over box layers
}

void Comm::communicate_blocking_isend(Atom &atom, int box_id)
{
  int nsend;
  AtomBuffer *buf_send = NULL, *buf_recv = NULL, *buf = NULL;
  // ID of box we're sending to and receiving from
  int send_target_box_id, recv_target_box_id;
  // Directions of buffers we're using to send to/recieve into
  int sendNeighbour, recvNeighbour;
  // Message tags
  int sendtag, recvtag;
  // Communicators to use
  MPI_Comm *send_comm, *recv_comm;
  // List of atoms to be sent
  int *sendlist;
  // Swap counter for this box
  int swapnum = 0;

  // DSM Multibox: Reference buffers/variables specific to this box.
  int& nswap = boxBufs[box_id].nswap;
  int* (&sendnum)[3] = boxBufs[box_id].sendnum;
  int* (&recvnum)[3] = boxBufs[box_id].recvnum;
  int* (&firstrecv)[3] = boxBufs[box_id].firstrecv;
  int*& sendproc = boxBufs[box_id].sendproc;
  int*& recvproc = boxBufs[box_id].recvproc;
  int*& sendneigh = boxBufs[box_id].sendneigh;
  int*& recvneigh = boxBufs[box_id].recvneigh;
  int* (&comm_send_size)[3] = boxBufs[box_id].comm_send_size;
  int* (&comm_recv_size)[3] = boxBufs[box_id].comm_recv_size;
  AtomBuffer (&bufs_send)[3][8] = boxBufs[box_id].bufs_send[COMMUNICATE_FUNCTION];
  AtomBuffer (&bufs_recv)[3][8] = boxBufs[box_id].bufs_recv[COMMUNICATE_FUNCTION];
  int* (&sendlists)[3][8] = boxBufs[box_id].sendlists;
  AtomBuffer& internal_buf_send_up = boxBufs[box_id].internal_buf_send_up[COMMUNICATE_FUNCTION];
  AtomBuffer& internal_buf_send_down = boxBufs[box_id].internal_buf_send_down[COMMUNICATE_FUNCTION];
  AtomBuffer& internal_buf_recv_up = boxBufs[box_id].internal_buf_recv_up[COMMUNICATE_FUNCTION];
  AtomBuffer& internal_buf_recv_down = boxBufs[box_id].internal_buf_recv_down[COMMUNICATE_FUNCTION];
  int*& internal_sendlist_up = boxBufs[box_id].internal_sendlist_up;
  int*& internal_sendlist_down = boxBufs[box_id].internal_sendlist_down;

  // 3 layers * nswaps per layer * 2 (sends and recvs)
  MPI_Request requests[3*nswap*2];
  int reqi = 0;

  for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {

    // Determine box IDs to send to/receive from on this layer
    layer_to_targets(&atom, box_layer_index, &send_target_box_id, &recv_target_box_id);

    // DSM nswap (number of swaps) is a constant calculated in setup()
    for(int iswap = 0; iswap < nswap; iswap++) {

      // Skip MPI communication in y dimension. Internal memory copies between boxes on same proc done in separate task
      if (iswap == 2 || iswap == 3) {         // iswap=2 swap in y -ve direction, iswap=3 swap in y +ve direction
        continue;
      }

      // Determine which buffers we are using in each communication
      sendNeighbour = sendneigh[iswap];
      recvNeighbour = recvneigh[iswap];
      buf_send = &bufs_send[box_layer_index][sendNeighbour];
      buf_recv = &bufs_recv[box_layer_index][recvNeighbour];
      send_comm = &boxBufs[send_target_box_id].comm[COMMUNICATE_FUNCTION];
      recv_comm = &boxBufs[atom.box_id].comm[COMMUNICATE_FUNCTION];
      sendtag = swapnum;
      recvtag = sendtag;
      ++swapnum;
      sendlist = sendlists[box_layer_index][sendNeighbour];

      /* pack buffer */

      // Added buffer size checks. Previously weren't necessary as buffers were reused from borders() (which ensured
      // ensured sufficient size) but now independent buffers used for each function
      if (comm_send_size[box_layer_index][iswap] > buf_send->maxsize) {
        buf_send->growsend(comm_send_size[box_layer_index][iswap]);
      }
      if (comm_recv_size[box_layer_index][iswap] > buf_recv->maxsize) {
        buf_recv->growrecv(comm_recv_size[box_layer_index][iswap]);
      }

      //#pragma omp barrier
      atom.pack_comm(sendnum[box_layer_index][iswap], sendlist, buf_send->buf, buf_send->pbc_any, buf_send->pbc_x,
                     buf_send->pbc_y, buf_send->pbc_z);

      //#pragma omp barrier

      /* exchange with another proc
         if self, set recv buffer to send buffer */

      // DSM Multibox: Reverse order of this if statement and include box id
      if (sendproc[iswap] == me && send_target_box_id == atom.box_id) {
        // Skip copies to myself, data will be unpacked directly from send buffer in next step
        continue;
      } else {
        //#pragma omp master
        {
          if(sizeof(MMD_float) == 4) {
            MPI_Irecv(buf_recv->buf, comm_recv_size[box_layer_index][iswap], MPI_FLOAT,
                      recvproc[iswap], recvtag, *recv_comm, &requests[reqi++]);
            MPI_Isend(buf_send->buf, comm_send_size[box_layer_index][iswap], MPI_FLOAT,
                     sendproc[iswap], sendtag, *send_comm, &requests[reqi++]);
          } else {
            MPI_Irecv(buf_recv->buf, comm_recv_size[box_layer_index][iswap], MPI_DOUBLE,
                      recvproc[iswap], recvtag, *recv_comm, &requests[reqi++]);
            MPI_Isend(buf_send->buf, comm_send_size[box_layer_index][iswap], MPI_DOUBLE,
                     sendproc[iswap], sendtag, *send_comm, &requests[reqi++]);
          }
        }
      }
    } // End of loop over number of swap
  } // End of loop over box layers

  // Wait for all send and recv requests to complete
  MPI_Waitall(reqi, requests, MPI_STATUSES_IGNORE);

  // Now do unpacks
  for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {
    // Determine box IDs to send to/receive from on this layer
    layer_to_targets(&atom, box_layer_index, &send_target_box_id, &recv_target_box_id);
    for (int iswap = 0; iswap < nswap; iswap++) {
      // Skip internal swaps - 2/3 are always y-ve and y+ve directions
      if (iswap == 2 || iswap == 3) {
        continue;
      }

      // Which buffer are we unpacking?
      if (sendproc[iswap] == me && send_target_box_id == atom.box_id) {
        // Sending to self. Just unpack the send buffer
        sendNeighbour = sendneigh[iswap];
        buf = &bufs_send[box_layer_index][sendNeighbour];
      } else {
        // Otherwise, unpack the receive buffer as normal
        recvNeighbour = recvneigh[iswap];
        buf = &bufs_recv[box_layer_index][recvNeighbour];
      }

      atom.unpack_comm(recvnum[box_layer_index][iswap], firstrecv[box_layer_index][iswap], buf->buf);
    } // End of loop over number of swap
  } // End of loop over box layers
}

// All steps in an iteration are individual tasks per box
// Buffer pack, send to neighbour, recv from neighbour and buffer unpack are their own tasks
void Comm::communicate_blocking_alltasks(Atom *atom, int box_id)
{
  // ID of box we're sending to and receiving from
  int send_target_box_id, recv_target_box_id;
  // Message tags
  int sendtag, recvtag;
  // Swap counter for this box
  int swapnum = 0;

  for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {

    // Determine box IDs to send to/receive from on this layer
    layer_to_targets(atom, box_layer_index, &send_target_box_id, &recv_target_box_id);

    // DSM nswap (number of swaps) is a constant calculated in setup()
    for(int iswap = 0; iswap < boxBufs[box_id].nswap; iswap++) {

      // Skip MPI communication in y dimension. Internal memory copies between boxes on same proc done in separate task
      if (iswap == 2 || iswap == 3) {         // iswap=2 swap in y -ve direction, iswap=3 swap in y +ve direction
        continue;
      }

      /* pack buffer */

      // swapnum as index to communicate sentinels as depend on individual iterations of this loop
      // No longer have 1:1 relationship between boxes and tasks in this implementation
      #pragma oss task label (communicate_pack) \
                             in(initialIntegrateSentinels[box_id]) \
                             out(communicatePackSentinels[box_id][swapnum]) \
                             firstprivate(atom, box_id, box_layer_index, iswap)
      {
        // Determine which buffers we are using in each communication
        int sendNeighbour = boxBufs[box_id].sendneigh[iswap];
        AtomBuffer *buf_send = &(boxBufs[box_id].bufs_send[COMMUNICATE_FUNCTION][box_layer_index][sendNeighbour]);

        // Added buffer size checks. Previously weren't necessary as buffers were reused from borders() (which ensured
        // ensured sufficient size) but now independent buffers used for each function
        if (boxBufs[box_id].comm_send_size[box_layer_index][iswap] > buf_send->maxsize) {
          buf_send->growsend(boxBufs[box_id].comm_send_size[box_layer_index][iswap]);
        }

        atom->pack_comm(boxBufs[box_id].sendnum[box_layer_index][iswap], boxBufs[box_id].sendlists[box_layer_index][sendNeighbour],
                        buf_send->buf, buf_send->pbc_any, buf_send->pbc_x, buf_send->pbc_y, buf_send->pbc_z);
      }

      // Set tags
      sendtag = swapnum;
      recvtag = sendtag;

      /* exchange with another proc
         if self, set recv buffer to send buffer */

      // DSM Multibox: Reverse order of this if statement and include box id
      if (boxBufs[box_id].sendproc[iswap] == me && send_target_box_id == atom->box_id) {

        // Skip MPI communication to myself
        // This unpack task has dependencies on the packing task as we unpack our own send buffer.
        // Unpack task in other branch depends only on recv tasks
        #pragma oss task label (communicate_unpack_self) \
                         in(communicatePackSentinels[box_id][swapnum]) \
                         out(communicateSentinels[box_id]) \
                         firstprivate(atom, box_id, box_layer_index, iswap)
        {
          // Sending to self. Just unpack the send buffer
          int neighbour = boxBufs[box_id].sendneigh[iswap];
          atom->unpack_comm(boxBufs[box_id].recvnum[box_layer_index][iswap],
                            boxBufs[box_id].firstrecv[box_layer_index][iswap],
                            boxBufs[box_id].bufs_send[COMMUNICATE_FUNCTION][box_layer_index][neighbour].buf);
        }

      } else {
        #pragma oss task label (communicate_send) \
                         in(communicatePackSentinels[box_id][swapnum]) \
                         out(communicateSendSentinel[box_id][swapnum]) \
                         firstprivate(box_id, box_layer_index, iswap, sendtag, send_target_box_id)
        {
          int sendNeighbour = boxBufs[box_id].sendneigh[iswap];
          AtomBuffer *buf_send = &boxBufs[box_id].bufs_send[COMMUNICATE_FUNCTION][box_layer_index][sendNeighbour];
          MPI_Comm *send_comm = &boxBufs[send_target_box_id].comm[COMMUNICATE_FUNCTION];

          if (sizeof(MMD_float) == 4) {
            MPI_Send(buf_send->buf, boxBufs[box_id].comm_send_size[box_layer_index][iswap], MPI_FLOAT,
                     boxBufs[box_id].sendproc[iswap], sendtag, *send_comm);
          } else {
            MPI_Send(buf_send->buf, boxBufs[box_id].comm_send_size[box_layer_index][iswap], MPI_DOUBLE,
                     boxBufs[box_id].sendproc[iswap], sendtag, *send_comm);
          }
        }

        // No in-dependencies on pack or send tasks. Recvs can be posted as soon as communicate starts
        #pragma oss task label (communicate_recv) \
                         in(initialIntegrateSentinels[box_id]) \
                         out(communicateRecvSentinels[box_id][swapnum]) \
                         firstprivate(atom, box_id, box_layer_index, iswap, recvtag)
        {
          int recvNeighbour = boxBufs[box_id].recvneigh[iswap];
          AtomBuffer *buf_recv = &boxBufs[box_id].bufs_recv[COMMUNICATE_FUNCTION][box_layer_index][recvNeighbour];
          MPI_Comm *recv_comm = &boxBufs[atom->box_id].comm[COMMUNICATE_FUNCTION];

          // Ensure recv buffer is large enough. Equivalent check for send buffer done in packing task
          if (boxBufs[box_id].comm_recv_size[box_layer_index][iswap] > buf_recv->maxsize) {
            buf_recv->growrecv(boxBufs[box_id].comm_recv_size[box_layer_index][iswap]);
          }

          if (sizeof(MMD_float) == 4) {
            MPI_Recv(buf_recv->buf, boxBufs[box_id].comm_recv_size[box_layer_index][iswap], MPI_FLOAT,
                     boxBufs[box_id].recvproc[iswap], recvtag, *recv_comm, MPI_STATUS_IGNORE);
          } else {
            MPI_Recv(buf_recv->buf, boxBufs[box_id].comm_recv_size[box_layer_index][iswap], MPI_DOUBLE,
                     boxBufs[box_id].recvproc[iswap], recvtag, *recv_comm, MPI_STATUS_IGNORE);
          }
        }

        // Now do unpacks
        // We only depend on corresponding recvs here; no dependencies on send tasks.
        // Sends don't need to complete until after force/before initialIntegrate calculations
        #pragma oss task label (communicate_unpack) \
                         in(communicateRecvSentinels[box_id][swapnum]) \
                         out(communicateSentinels[box_id]) \
                         firstprivate(atom, box_id, box_layer_index, iswap)
        {
          int neighbour = boxBufs[box_id].recvneigh[iswap];
          atom->unpack_comm(boxBufs[box_id].recvnum[box_layer_index][iswap],
                            boxBufs[box_id].firstrecv[box_layer_index][iswap],
                            boxBufs[box_id].bufs_recv[COMMUNICATE_FUNCTION][box_layer_index][neighbour].buf);
        }

      } // End of skip-MPI branch

      // Counts total number of swaps over all box layers. i.e. value is 0 to 30 (3 box layers * 10 swaps per layer)
      ++swapnum;

    } // End of loop over number of swaps
  } // End of loop over box layers
} // End of function

// Alltasks but recv tasks posted first
void Comm::communicate_blocking_alltasks_recvfirst(Atom *atom, int box_id)
{
  // ID of box we're sending to and receiving from
  int send_target_box_id, recv_target_box_id;
  // Message tags
  int sendtag, recvtag;
  // Swap counter for this box
  int swapnum = 0;

  for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {

    // Determine box IDs to send to/receive from on this layer
    layer_to_targets(atom, box_layer_index, &send_target_box_id, &recv_target_box_id);

    // No in-dependencies on pack or send tasks. Recvs posted as soon as communicate starts to avoid late receiver
    // performance problems
    for (int iswap = 0; iswap < boxBufs[box_id].nswap; iswap++) {

      // Skip MPI communication in y dimension. Internal memory copies between boxes on same proc done in separate task
      if (iswap == 2 || iswap == 3) {         // iswap=2 swap in y -ve direction, iswap=3 swap in y +ve direction
        continue;
      }

      // Set tags
      recvtag = swapnum;

      // Only generate recv if not sending to self
      if (!(boxBufs[box_id].sendproc[iswap] == me && send_target_box_id == atom->box_id)) {
        #pragma oss task label (communicate_recv) \
                         in(initialIntegrateSentinels[box_id]) \
                         out(communicateRecvSentinels[box_id][swapnum]) \
                         firstprivate(atom, box_id, box_layer_index, iswap, recvtag)
        {
          int recvNeighbour = boxBufs[box_id].recvneigh[iswap];
          AtomBuffer *buf_recv = &boxBufs[box_id].bufs_recv[COMMUNICATE_FUNCTION][box_layer_index][recvNeighbour];
          MPI_Comm *recv_comm = &boxBufs[atom->box_id].comm[COMMUNICATE_FUNCTION];

          // Ensure recv buffer is large enough. Equivalent check for send buffer done in packing task
          if (boxBufs[box_id].comm_recv_size[box_layer_index][iswap] > buf_recv->maxsize) {
            buf_recv->growrecv(boxBufs[box_id].comm_recv_size[box_layer_index][iswap]);
          }

          if (sizeof(MMD_float) == 4) {
            MPI_Recv(buf_recv->buf, boxBufs[box_id].comm_recv_size[box_layer_index][iswap], MPI_FLOAT,
                     boxBufs[box_id].recvproc[iswap], recvtag, *recv_comm, MPI_STATUS_IGNORE);
          } else {
            MPI_Recv(buf_recv->buf, boxBufs[box_id].comm_recv_size[box_layer_index][iswap], MPI_DOUBLE,
                     boxBufs[box_id].recvproc[iswap], recvtag, *recv_comm, MPI_STATUS_IGNORE);
          }
        }
      }

      ++swapnum;
    }
  }

  // Reset swapnum for next loop to generate remaining pack, send and unpack tasks
  swapnum = 0;

  for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {

    // Determine box IDs to send to/receive from on this layer
    layer_to_targets(atom, box_layer_index, &send_target_box_id, &recv_target_box_id);

    // DSM nswap (number of swaps) is a constant calculated in setup()
    for(int iswap = 0; iswap < boxBufs[box_id].nswap; iswap++) {

      // Skip MPI communication in y dimension. Internal memory copies between boxes on same proc done in separate task
      if (iswap == 2 || iswap == 3) {         // iswap=2 swap in y -ve direction, iswap=3 swap in y +ve direction
        continue;
      }

      /* pack buffer */

      // swapnum as index to communicate sentinels as depend on individual iterations of this loop
      // No longer have 1:1 relationship between boxes and tasks in this implementation
      #pragma oss task label (communicate_pack) \
                             in(initialIntegrateSentinels[box_id]) \
                             out(communicatePackSentinels[box_id][swapnum]) \
                             firstprivate(atom, box_id, box_layer_index, iswap)
      {
        // Determine which buffers we are using in each communication
        int sendNeighbour = boxBufs[box_id].sendneigh[iswap];
        AtomBuffer *buf_send = &(boxBufs[box_id].bufs_send[COMMUNICATE_FUNCTION][box_layer_index][sendNeighbour]);

        // Added buffer size checks. Previously weren't necessary as buffers were reused from borders() (which ensured
        // ensured sufficient size) but now independent buffers used for each function
        if (boxBufs[box_id].comm_send_size[box_layer_index][iswap] > buf_send->maxsize) {
          buf_send->growsend(boxBufs[box_id].comm_send_size[box_layer_index][iswap]);
        }

        atom->pack_comm(boxBufs[box_id].sendnum[box_layer_index][iswap], boxBufs[box_id].sendlists[box_layer_index][sendNeighbour],
                        buf_send->buf, buf_send->pbc_any, buf_send->pbc_x, buf_send->pbc_y, buf_send->pbc_z);
      }

      // Set tags
      sendtag = swapnum;

      /* exchange with another proc
         if self, set recv buffer to send buffer */

      // DSM Multibox: Reverse order of this if statement and include box id
      if (boxBufs[box_id].sendproc[iswap] == me && send_target_box_id == atom->box_id) {

        // Skip MPI communication to myself
        // This unpack task has dependencies on the packing task as we unpack our own send buffer.
        // Unpack task in other branch depends only on recv tasks
        #pragma oss task label (communicate_unpack_self) \
                         in(communicatePackSentinels[box_id][swapnum]) \
                         out(communicateSentinels[box_id]) \
                         firstprivate(atom, box_id, box_layer_index, iswap)
        {
          // Sending to self. Just unpack the send buffer
          int neighbour = boxBufs[box_id].sendneigh[iswap];
          atom->unpack_comm(boxBufs[box_id].recvnum[box_layer_index][iswap],
                            boxBufs[box_id].firstrecv[box_layer_index][iswap],
                            boxBufs[box_id].bufs_send[COMMUNICATE_FUNCTION][box_layer_index][neighbour].buf);
        }

      } else {
        #pragma oss task label (communicate_send) \
                         in(communicatePackSentinels[box_id][swapnum]) \
                         out(communicateSendSentinel[box_id][swapnum]) \
                         firstprivate(box_id, box_layer_index, iswap, sendtag, send_target_box_id)
        {
          int sendNeighbour = boxBufs[box_id].sendneigh[iswap];
          AtomBuffer *buf_send = &boxBufs[box_id].bufs_send[COMMUNICATE_FUNCTION][box_layer_index][sendNeighbour];
          MPI_Comm *send_comm = &boxBufs[send_target_box_id].comm[COMMUNICATE_FUNCTION];

          if (sizeof(MMD_float) == 4) {
            MPI_Send(buf_send->buf, boxBufs[box_id].comm_send_size[box_layer_index][iswap], MPI_FLOAT,
                     boxBufs[box_id].sendproc[iswap], sendtag, *send_comm);
          } else {
            MPI_Send(buf_send->buf, boxBufs[box_id].comm_send_size[box_layer_index][iswap], MPI_DOUBLE,
                     boxBufs[box_id].sendproc[iswap], sendtag, *send_comm);
          }
        }

        // Recv tasks already generated at start of function

        // Now do unpacks
        // We only depend on corresponding recvs here; no dependencies on send tasks.
        // Sends don't need to complete until after force/before initialIntegrate calculations
        #pragma oss task label (communicate_unpack) \
                         in(communicateRecvSentinels[box_id][swapnum]) \
                         out(communicateSentinels[box_id]) \
                         firstprivate(atom, box_id, box_layer_index, iswap)
        {
          int neighbour = boxBufs[box_id].recvneigh[iswap];
          atom->unpack_comm(boxBufs[box_id].recvnum[box_layer_index][iswap],
                            boxBufs[box_id].firstrecv[box_layer_index][iswap],
                            boxBufs[box_id].bufs_recv[COMMUNICATE_FUNCTION][box_layer_index][neighbour].buf);
        }

      } // End of skip-MPI branch

      // Counts total number of swaps over all box layers. i.e. value is 0 to 30 (3 box layers * 10 swaps per layer)
      ++swapnum;

    } // End of loop over number of swaps
  } // End of loop over box layers
} // End of function

// Alltasks but using nonblocking communication to reduce task overhead
void Comm::communicate_nonblocking_alltasks_tampi_iwait(Atom *atom, int box_id)
{
  // ID of box we're sending to and receiving from
  int send_target_box_id, recv_target_box_id;
  // Message tags
  int sendtag, recvtag;
  // Swap counter for this box
  int swapnum = 0;

  for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {

    // Determine box IDs to send to/receive from on this layer
    layer_to_targets(atom, box_layer_index, &send_target_box_id, &recv_target_box_id);

    // DSM nswap (number of swaps) is a constant calculated in setup()
    for(int iswap = 0; iswap < boxBufs[box_id].nswap; iswap++) {

      // Skip MPI communication in y dimension. Internal memory copies between boxes on same proc done in separate task
      if (iswap == 2 || iswap == 3) {         // iswap=2 swap in y -ve direction, iswap=3 swap in y +ve direction
        continue;
      }

      /* pack buffer */

      // swapnum as index to communicate sentinels as depend on individual iterations of this loop
      // No longer have 1:1 relationship between boxes and tasks in this implementation
      #pragma oss task label (communicate_pack) \
                             in(initialIntegrateSentinels[box_id]) \
                             out(communicatePackSentinels[box_id][swapnum]) \
                             firstprivate(atom, box_id, box_layer_index, iswap)
      {
        // Determine which buffers we are using in each communication
        int sendNeighbour = boxBufs[box_id].sendneigh[iswap];
        AtomBuffer *buf_send = &(boxBufs[box_id].bufs_send[COMMUNICATE_FUNCTION][box_layer_index][sendNeighbour]);

        // Added buffer size checks. Previously weren't necessary as buffers were reused from borders() (which ensured
        // ensured sufficient size) but now independent buffers used for each function
        if (boxBufs[box_id].comm_send_size[box_layer_index][iswap] > buf_send->maxsize) {
          buf_send->growsend(boxBufs[box_id].comm_send_size[box_layer_index][iswap]);
        }

        atom->pack_comm(boxBufs[box_id].sendnum[box_layer_index][iswap], boxBufs[box_id].sendlists[box_layer_index][sendNeighbour],
                        buf_send->buf, buf_send->pbc_any, buf_send->pbc_x, buf_send->pbc_y, buf_send->pbc_z);
      }

      // Set tags
      sendtag = swapnum;
      recvtag = sendtag;

      /* exchange with another proc
         if self, set recv buffer to send buffer */

      // DSM Multibox: Reverse order of this if statement and include box id
      if (boxBufs[box_id].sendproc[iswap] == me && send_target_box_id == atom->box_id) {

        // Skip MPI communication to myself
        // This unpack task has dependencies on the packing task as we unpack our own send buffer.
        // Unpack task in other branch depends only on recv tasks
        #pragma oss task label (communicate_unpack_self) \
                         in(communicatePackSentinels[box_id][swapnum]) \
                         out(communicateSentinels[box_id]) \
                         firstprivate(atom, box_id, box_layer_index, iswap)
        {
          // Sending to self. Just unpack the send buffer
          int neighbour = boxBufs[box_id].sendneigh[iswap];
          atom->unpack_comm(boxBufs[box_id].recvnum[box_layer_index][iswap],
                            boxBufs[box_id].firstrecv[box_layer_index][iswap],
                            boxBufs[box_id].bufs_send[COMMUNICATE_FUNCTION][box_layer_index][neighbour].buf);
        }

      } else {
        #pragma oss task label (communicate_send) \
                         in(communicatePackSentinels[box_id][swapnum]) \
                         out(communicateSendSentinel[box_id][swapnum]) \
                         firstprivate(box_id, box_layer_index, iswap, sendtag, send_target_box_id)
        {
          MPI_Request request;
          int sendNeighbour = boxBufs[box_id].sendneigh[iswap];
          AtomBuffer *buf_send = &boxBufs[box_id].bufs_send[COMMUNICATE_FUNCTION][box_layer_index][sendNeighbour];
          MPI_Comm *send_comm = &boxBufs[send_target_box_id].comm[COMMUNICATE_FUNCTION];

          if (sizeof(MMD_float) == 4) {
            MPI_Isend(buf_send->buf, boxBufs[box_id].comm_send_size[box_layer_index][iswap], MPI_FLOAT,
                     boxBufs[box_id].sendproc[iswap], sendtag, *send_comm, &request);
          } else {
            MPI_Isend(buf_send->buf, boxBufs[box_id].comm_send_size[box_layer_index][iswap], MPI_DOUBLE,
                     boxBufs[box_id].sendproc[iswap], sendtag, *send_comm, &request);
          }

          // Allows task to return but dependencies are not released until MPI operations complete
          TAMPI_Iwait(&request, MPI_STATUS_IGNORE);
        }

        // No in-dependencies on pack or send tasks. Recvs can be posted as soon as communicate starts
        #pragma oss task label (communicate_recv) \
                         in(initialIntegrateSentinels[box_id]) \
                         out(communicateRecvSentinels[box_id][swapnum]) \
                         firstprivate(atom, box_id, box_layer_index, iswap, recvtag)
        {
          MPI_Request request;
          int recvNeighbour = boxBufs[box_id].recvneigh[iswap];
          AtomBuffer *buf_recv = &boxBufs[box_id].bufs_recv[COMMUNICATE_FUNCTION][box_layer_index][recvNeighbour];
          MPI_Comm *recv_comm = &boxBufs[atom->box_id].comm[COMMUNICATE_FUNCTION];

          // Ensure recv buffer is large enough. Equivalent check for send buffer done in packing task
          if (boxBufs[box_id].comm_recv_size[box_layer_index][iswap] > buf_recv->maxsize) {
            buf_recv->growrecv(boxBufs[box_id].comm_recv_size[box_layer_index][iswap]);
          }

          if (sizeof(MMD_float) == 4) {
            MPI_Irecv(buf_recv->buf, boxBufs[box_id].comm_recv_size[box_layer_index][iswap], MPI_FLOAT,
                     boxBufs[box_id].recvproc[iswap], recvtag, *recv_comm, &request);
          } else {
            MPI_Irecv(buf_recv->buf, boxBufs[box_id].comm_recv_size[box_layer_index][iswap], MPI_DOUBLE,
                     boxBufs[box_id].recvproc[iswap], recvtag, *recv_comm, &request);
          }

          TAMPI_Iwait(&request, MPI_STATUS_IGNORE);
        }

        // Now do unpacks
        // We only depend on corresponding recvs here; no dependencies on send tasks.
        // Sends don't need to complete until after force/before initialIntegrate calculations
        #pragma oss task label (communicate_unpack) \
                         in(communicateRecvSentinels[box_id][swapnum]) \
                         out(communicateSentinels[box_id]) \
                         firstprivate(atom, box_id, box_layer_index, iswap)
        {
          int neighbour = boxBufs[box_id].recvneigh[iswap];
          atom->unpack_comm(boxBufs[box_id].recvnum[box_layer_index][iswap],
                            boxBufs[box_id].firstrecv[box_layer_index][iswap],
                            boxBufs[box_id].bufs_recv[COMMUNICATE_FUNCTION][box_layer_index][neighbour].buf);
        }

      } // End of skip-MPI branch

      // Counts total number of swaps over all box layers. i.e. value is 0 to 30 (3 box layers * 10 swaps per layer)
      ++swapnum;

    } // End of loop over number of swaps
  } // End of loop over box layers
} // End of function

// Buffer pack, send to neighbour, recv from neighbour and buffer unpack are their own tasks per neighbour
void Comm::communicate_nonblocking_neighbourtasks(Atom *atom, int box_id)
{
  // ID of box we're sending to and receiving from
  int send_target_box_id[3], recv_target_box_id[3];

  // No outer loop over box layers this implementation. All 3 layers captured in each task

  // Determine box IDs to send to/receive from on all layers
  layer_to_targets(atom, 0, &send_target_box_id[0], &recv_target_box_id[0]);
  layer_to_targets(atom, 1, &send_target_box_id[1], &recv_target_box_id[1]);
  layer_to_targets(atom, 2, &send_target_box_id[2], &recv_target_box_id[2]);

  // DSM nswap (number of swaps) is a constant calculated in setup()
  for(int iswap = 0; iswap < boxBufs[box_id].nswap; iswap++) {

    // Skip MPI communication in y dimension. Internal memory copies between boxes on same proc done in separate task
    if (iswap == 2 || iswap == 3) {         // iswap=2 swap in y -ve direction, iswap=3 swap in y +ve direction
      continue;
    }

    /* pack buffer */

    // swapnum as index to communicate sentinels as depend on individual iterations of this loop
    // No longer have 1:1 relationship between boxes and tasks in this implementation
    #pragma oss task label (communicate_pack) \
                     in(initialIntegrateSentinels[box_id]) \
                     out(communicatePackSentinels[box_id][iswap]) \
                     firstprivate(atom, box_id, iswap)
    {
      int sendNeighbour = boxBufs[box_id].sendneigh[iswap];
      for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {
        // Determine which buffers we are using in each communication
        AtomBuffer *buf_send = &(boxBufs[box_id].bufs_send[COMMUNICATE_FUNCTION][box_layer_index][sendNeighbour]);

        // Added buffer size checks. Previously weren't necessary as buffers were reused from borders() (which ensured
        // ensured sufficient size) but now independent buffers used for each function
        if (boxBufs[box_id].comm_send_size[box_layer_index][iswap] > buf_send->maxsize) {
          buf_send->growsend(boxBufs[box_id].comm_send_size[box_layer_index][iswap]);
        }

        atom->pack_comm(boxBufs[box_id].sendnum[box_layer_index][iswap],
                        boxBufs[box_id].sendlists[box_layer_index][sendNeighbour], buf_send->buf, buf_send->pbc_any,
                        buf_send->pbc_x, buf_send->pbc_y, buf_send->pbc_z);
      }
    }

    /* exchange with another proc
       if self, set recv buffer to send buffer */

    #pragma oss task label (communicate_send) \
                     in(communicatePackSentinels[box_id][iswap]) \
                     out(communicateSendSentinel[box_id][iswap]) \
                     firstprivate(box_id, iswap)
    {
      MPI_Request requests[] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL};
      int sendNeighbour = boxBufs[box_id].sendneigh[iswap];
      for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {
        AtomBuffer *buf_send = &boxBufs[box_id].bufs_send[COMMUNICATE_FUNCTION][box_layer_index][sendNeighbour];
        MPI_Comm *send_comm = &boxBufs[send_target_box_id[box_layer_index]].comm[COMMUNICATE_FUNCTION];

        // Calculate send tag. In range 0 to 29, equivalent to sendnum from other implementations
        int sendtag = box_layer_index * 10 + iswap;

        if (sizeof(MMD_float) == 4) {
          MPI_Isend(buf_send->buf, boxBufs[box_id].comm_send_size[box_layer_index][iswap], MPI_FLOAT,
                   boxBufs[box_id].sendproc[iswap], sendtag, *send_comm, &requests[box_layer_index]);
        } else {
          MPI_Isend(buf_send->buf, boxBufs[box_id].comm_send_size[box_layer_index][iswap], MPI_DOUBLE,
                   boxBufs[box_id].sendproc[iswap], sendtag, *send_comm, &requests[box_layer_index]);
        }
      }

      MPI_Waitall(3, requests, MPI_STATUSES_IGNORE);
    }

    // No in-dependencies on pack or send tasks. Recvs can be posted as soon as communicate starts
    #pragma oss task label (communicate_recv) \
                     in(initialIntegrateSentinels[box_id]) \
                     out(communicateRecvSentinels[box_id][iswap]) \
                     firstprivate(atom, box_id, iswap)
    {
      MPI_Request requests[] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL};
      int recvNeighbour = boxBufs[box_id].recvneigh[iswap];
      for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {
        AtomBuffer *buf_recv = &boxBufs[box_id].bufs_recv[COMMUNICATE_FUNCTION][box_layer_index][recvNeighbour];
        MPI_Comm *recv_comm = &boxBufs[atom->box_id].comm[COMMUNICATE_FUNCTION];

        // Calculate recv tag. In range 0 to 29, equivalent to sendnum from other implementations
        int recvtag = box_layer_index * 10 + iswap;

        // Ensure recv buffer is large enough. Equivalent check for send buffer done in packing task
        if (boxBufs[box_id].comm_recv_size[box_layer_index][iswap] > buf_recv->maxsize) {
          buf_recv->growrecv(boxBufs[box_id].comm_recv_size[box_layer_index][iswap]);
        }

        if (sizeof(MMD_float) == 4) {
          MPI_Irecv(buf_recv->buf, boxBufs[box_id].comm_recv_size[box_layer_index][iswap], MPI_FLOAT,
                   boxBufs[box_id].recvproc[iswap], recvtag, *recv_comm,&requests[box_layer_index]);
        } else {
          MPI_Irecv(buf_recv->buf, boxBufs[box_id].comm_recv_size[box_layer_index][iswap], MPI_DOUBLE,
                   boxBufs[box_id].recvproc[iswap], recvtag, *recv_comm, &requests[box_layer_index]);
        }
      }

      MPI_Waitall(3, requests, MPI_STATUSES_IGNORE);
    }

    // Now do unpacks
    #pragma oss task label (communicate_unpack) \
                     in(communicateRecvSentinels[box_id][iswap]) \
                     out(communicateSentinels[box_id]) \
                     firstprivate(atom, box_id, iswap)
    {
      int recvNeighbour = boxBufs[box_id].recvneigh[iswap];
      for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {
        atom->unpack_comm(boxBufs[box_id].recvnum[box_layer_index][iswap],
                          boxBufs[box_id].firstrecv[box_layer_index][iswap],
                          boxBufs[box_id].bufs_recv[COMMUNICATE_FUNCTION][box_layer_index][recvNeighbour].buf);
      }
    }

  } // End of loop over number of swaps
} // End of function

// One message per neighbour containing all three box layers' data
// TODO: This approach is wrong. We're sending our single buffer to only one of the boxes on the neighbour. Would need a
// task that unpacks data into multiple boxes, rather than just its own box_id
void Comm::communicate_blocking_single_message_neighbourtasks(Atom *atom, int box_id)
{
  // ID of box we're sending to and receiving from
  int send_target_box_id[3], recv_target_box_id[3];

  // No outer loop over box layers this implementation. All 3 layers captured in each task

  // Determine box IDs to send to/receive from on all layers
  layer_to_targets(atom, 0, &send_target_box_id[0], &recv_target_box_id[0]);
  layer_to_targets(atom, 1, &send_target_box_id[1], &recv_target_box_id[1]);
  layer_to_targets(atom, 2, &send_target_box_id[2], &recv_target_box_id[2]);

  // nswap (number of swaps) is a constant calculated in setup()
  for(int iswap = 0; iswap < boxBufs[box_id].nswap; iswap++) {

    // Skip MPI communication in y dimension. Internal memory copies between boxes on same proc done in separate task
    if (iswap == 2 || iswap == 3) {         // iswap=2 swap in y -ve direction, iswap=3 swap in y +ve direction
      continue;
    }

    /* pack buffer */

    // swapnum as index to communicate sentinels as depend on individual iterations of this loop
    // No longer have 1:1 relationship between boxes and tasks in this implementation
    #pragma oss task label("communicate_pack") \
                     in(initialIntegrateSentinels[box_id]) \
                     out(communicatePackSentinels[box_id][iswap]) \
                     firstprivate(atom, box_id, iswap)
    {
      int sendNeighbour = boxBufs[box_id].sendneigh[iswap];

      // Calculate the total we're sending over all 3 boxes layers and grow array as appropriate
      int totalSend = boxBufs[box_id].comm_send_size[0][iswap] +
                      boxBufs[box_id].comm_send_size[1][iswap] +
                      boxBufs[box_id].comm_send_size[2][iswap];
      // Use layer 0's buffer for all 3 layers' data
      AtomBuffer *buf_send = &(boxBufs[box_id].bufs_send[COMMUNICATE_FUNCTION][0][sendNeighbour]);
      // Added buffer size checks. Previously weren't necessary as buffers were reused from borders() (which ensured
      // ensured sufficient size) but now independent buffers used for each function
      if (totalSend > buf_send->maxsize) {
        buf_send->growsend(totalSend);
      }

      // Pack all data into single buffer: box layer 0 to 2, end-to-end.
      // Update buffer offset every loop iteration by size of data packed
      int bufferOffset = 0;
      for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {
        // Get PBC information from buffer we *would* be using in 3 messages-per-neighbour implementation
        AtomBuffer *buf_send_pbc = &(boxBufs[box_id].bufs_send[COMMUNICATE_FUNCTION][box_layer_index][sendNeighbour]);
        atom->pack_comm(boxBufs[box_id].sendnum[box_layer_index][iswap],
                        boxBufs[box_id].sendlists[box_layer_index][sendNeighbour], &(buf_send->buf[bufferOffset]),
                        buf_send_pbc->pbc_any, buf_send_pbc->pbc_x, buf_send_pbc->pbc_y, buf_send_pbc->pbc_z);
        bufferOffset += boxBufs[box_id].comm_send_size[box_layer_index][iswap];
      }
    }

    #pragma oss task label("communicate_send") \
                     in(communicatePackSentinels[box_id][iswap]) \
                     out(communicateSendSentinel[box_id][iswap]) \
                     firstprivate(box_id, iswap)
    {
      int sendNeighbour = boxBufs[box_id].sendneigh[iswap];

      // Calculate the total we're sending over all 3 boxes layers and pass to MPI send call
      int totalSend = boxBufs[box_id].comm_send_size[0][iswap] +
                      boxBufs[box_id].comm_send_size[1][iswap] +
                      boxBufs[box_id].comm_send_size[2][iswap];

      // Calculate send tag. In range 0 to 8, equivalent to sendnum from other implementations
      int sendtag = iswap;

      // TODO: This is wrong. We're sending our single buffer to only one of the boxes on the neighbour. Would need a
      // task that unpacks data into multiple boxes, rather than just its own box_id

      // Use layer 0's buffer for all 3 layers' data
      AtomBuffer *buf_send = &boxBufs[box_id].bufs_send[COMMUNICATE_FUNCTION][0][sendNeighbour];
      MPI_Comm *send_comm = &boxBufs[send_target_box_id[0]].comm[COMMUNICATE_FUNCTION];

      // Perform the (single) send for this neighbour
      if (sizeof(MMD_float) == 4) {
        MPI_Send(buf_send->buf, totalSend, MPI_FLOAT, boxBufs[box_id].sendproc[iswap], sendtag, *send_comm);
      } else {
        MPI_Send(buf_send->buf, totalSend, MPI_DOUBLE, boxBufs[box_id].sendproc[iswap], sendtag, *send_comm);
      }
    }

    // No in-dependencies on pack or send tasks. Recvs can be posted as soon as any previous recvs are completed and the
    // buffers become free again
    #pragma oss task label("communicate_recv & unpack") \
                     inout(communicateSentinels[box_id]) \
                     firstprivate(atom, box_id, iswap)
    {
      int recvNeighbour = boxBufs[box_id].recvneigh[iswap];

      // Calculate the total we're receiving over all 3 boxes layers and pass to MPI recv call
      int totalRecv = boxBufs[box_id].comm_send_size[0][iswap] +
                      boxBufs[box_id].comm_send_size[1][iswap] +
                      boxBufs[box_id].comm_send_size[2][iswap];

      // Calculate recv tag. In range 0 to 8, equivalent to sendnum from other implementations
      int recvtag = iswap;

      // Use layer 0's buffer for all 3 layers' data
      AtomBuffer *buf_recv = &boxBufs[box_id].bufs_recv[COMMUNICATE_FUNCTION][0][recvNeighbour];
      MPI_Comm *recv_comm = &boxBufs[atom->box_id].comm[COMMUNICATE_FUNCTION];

      // Ensure recv buffer is large enough. Equivalent check for send buffer done in packing task
      if (totalRecv > buf_recv->maxsize) {
        buf_recv->growrecv(totalRecv);
      }

      // Perform the (single) receive
      if (sizeof(MMD_float) == 4) {
        MPI_Recv(buf_recv->buf, totalRecv, MPI_FLOAT, boxBufs[box_id].recvproc[iswap], recvtag, *recv_comm, MPI_STATUS_IGNORE);
      } else {
        MPI_Recv(buf_recv->buf, totalRecv, MPI_DOUBLE, boxBufs[box_id].recvproc[iswap], recvtag, *recv_comm, MPI_STATUS_IGNORE);
      }

      // Immediately unpack received buffer.
      // Data is in single buffer: box layer 0 to 2, packed end-to-end.
      // Update buffer offset every loop iteration by size of data unpacked
      int bufferOffset = 0;
      for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {
        atom->unpack_comm(boxBufs[box_id].recvnum[box_layer_index][iswap],
                          boxBufs[box_id].firstrecv[box_layer_index][iswap],
                          &(buf_recv->buf[bufferOffset]));
        bufferOffset += boxBufs[box_id].comm_recv_size[box_layer_index][iswap];
      }
    }

  } // End of loop over number of swaps
} // End of function

// Nonblocking_neighbourtasks but using TAMPI_Iwaitall
// Combined recv+unpack tasks. Only sends are non-blocking now
void Comm::communicate_nonblocking_neighbourtasks_tampi_iwaitall(Atom *atom, int box_id)
{
  // ID of box we're sending to and receiving from
  int send_target_box_id[3], recv_target_box_id[3];

  // No outer loop over box layers this implementation. All 3 layers captured in each task

  // Determine box IDs to send to/receive from on all layers
  layer_to_targets(atom, 0, &send_target_box_id[0], &recv_target_box_id[0]);
  layer_to_targets(atom, 1, &send_target_box_id[1], &recv_target_box_id[1]);
  layer_to_targets(atom, 2, &send_target_box_id[2], &recv_target_box_id[2]);

  // DSM nswap (number of swaps) is a constant calculated in setup()
  for(int iswap = 0; iswap < boxBufs[box_id].nswap; iswap++) {

    // Skip MPI communication in y dimension. Internal memory copies between boxes on same proc done in separate task
    if (iswap == 2 || iswap == 3) {         // iswap=2 swap in y -ve direction, iswap=3 swap in y +ve direction
      continue;
    }

    /* pack buffer */

    // swapnum as index to communicate sentinels as depend on individual iterations of this loop
    // No longer have 1:1 relationship between boxes and tasks in this implementation
    #pragma oss task label (communicate_pack) \
                     in(initialIntegrateSentinels[box_id]) \
                     out(communicatePackSentinels[box_id][iswap]) \
                     firstprivate(atom, box_id, iswap)
    {
      int sendNeighbour = boxBufs[box_id].sendneigh[iswap];
      for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {
        // Determine which buffers we are using in each communication
        AtomBuffer *buf_send = &(boxBufs[box_id].bufs_send[COMMUNICATE_FUNCTION][box_layer_index][sendNeighbour]);

        // Added buffer size checks. Previously weren't necessary as buffers were reused from borders() (which ensured
        // ensured sufficient size) but now independent buffers used for each function
        if (boxBufs[box_id].comm_send_size[box_layer_index][iswap] > buf_send->maxsize) {
          buf_send->growsend(boxBufs[box_id].comm_send_size[box_layer_index][iswap]);
        }

        atom->pack_comm(boxBufs[box_id].sendnum[box_layer_index][iswap],
                        boxBufs[box_id].sendlists[box_layer_index][sendNeighbour], buf_send->buf, buf_send->pbc_any,
                        buf_send->pbc_x, buf_send->pbc_y, buf_send->pbc_z);
      }
    }

    /* exchange with another proc
       if self, set recv buffer to send buffer */

    #pragma oss task label (communicate_send) \
                     in(communicatePackSentinels[box_id][iswap]) \
                     out(communicateSendSentinel[box_id][iswap]) \
                     firstprivate(box_id, iswap)
    {
      MPI_Request requests[] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL};
      int sendNeighbour = boxBufs[box_id].sendneigh[iswap];
      for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {
        AtomBuffer *buf_send = &boxBufs[box_id].bufs_send[COMMUNICATE_FUNCTION][box_layer_index][sendNeighbour];
        MPI_Comm *send_comm = &boxBufs[send_target_box_id[box_layer_index]].comm[COMMUNICATE_FUNCTION];

        // Calculate send tag. In range 0 to 29, equivalent to sendnum from other implementations
        int sendtag = box_layer_index * 10 + iswap;

        if (sizeof(MMD_float) == 4) {
          MPI_Isend(buf_send->buf, boxBufs[box_id].comm_send_size[box_layer_index][iswap], MPI_FLOAT,
                    boxBufs[box_id].sendproc[iswap], sendtag, *send_comm, &requests[box_layer_index]);
        } else {
          MPI_Isend(buf_send->buf, boxBufs[box_id].comm_send_size[box_layer_index][iswap], MPI_DOUBLE,
                    boxBufs[box_id].sendproc[iswap], sendtag, *send_comm, &requests[box_layer_index]);
        }
      }

      TAMPI_Iwaitall(3, requests, MPI_STATUSES_IGNORE);
    }

    // No in-dependencies on pack or send tasks. Recvs can be posted as soon as any previous recvs are completed and the
    // buffers become free again
    #pragma oss task label ("communicate_recv & unpack") \
                     inout(communicateSentinels[box_id]) \
                     firstprivate(atom, box_id, iswap)
    {
      MPI_Request requests[] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL};
      int recvNeighbour = boxBufs[box_id].recvneigh[iswap];
      for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {
        AtomBuffer *buf_recv = &boxBufs[box_id].bufs_recv[COMMUNICATE_FUNCTION][box_layer_index][recvNeighbour];
        MPI_Comm *recv_comm = &boxBufs[atom->box_id].comm[COMMUNICATE_FUNCTION];

        // Calculate recv tag. In range 0 to 29, equivalent to sendnum from other implementations
        int recvtag = box_layer_index * 10 + iswap;

        // Ensure recv buffer is large enough. Equivalent check for send buffer done in packing task
        if (boxBufs[box_id].comm_recv_size[box_layer_index][iswap] > buf_recv->maxsize) {
          buf_recv->growrecv(boxBufs[box_id].comm_recv_size[box_layer_index][iswap]);
        }

        if (sizeof(MMD_float) == 4) {
          MPI_Recv(buf_recv->buf, boxBufs[box_id].comm_recv_size[box_layer_index][iswap], MPI_FLOAT,
                    boxBufs[box_id].recvproc[iswap], recvtag, *recv_comm, MPI_STATUS_IGNORE);
        } else {
          MPI_Recv(buf_recv->buf, boxBufs[box_id].comm_recv_size[box_layer_index][iswap], MPI_DOUBLE,
                    boxBufs[box_id].recvproc[iswap], recvtag, *recv_comm, MPI_STATUS_IGNORE);
        }

        // Immediately unpack received buffer.
        atom->unpack_comm(boxBufs[box_id].recvnum[box_layer_index][iswap], boxBufs[box_id].firstrecv[box_layer_index][iswap], buf_recv->buf);
      }
    }

  } // End of loop over number of swaps
} // End of function

// Buffer pack, send to neighbour, recv from neighbour and buffer unpack are their own tasks per neighbour
// Communication is blocking
void Comm::communicate_blocking_neighbourtasks(Atom *atom, int box_id)
{
  // ID of box we're sending to and receiving from
  int send_target_box_id[3], recv_target_box_id[3];

  // No outer loop over box layers this implementation. All 3 layers captured in each task

  // Determine box IDs to send to/receive from on all layers
  layer_to_targets(atom, 0, &send_target_box_id[0], &recv_target_box_id[0]);
  layer_to_targets(atom, 1, &send_target_box_id[1], &recv_target_box_id[1]);
  layer_to_targets(atom, 2, &send_target_box_id[2], &recv_target_box_id[2]);

  // DSM nswap (number of swaps) is a constant calculated in setup()
  for(int iswap = 0; iswap < boxBufs[box_id].nswap; iswap++) {

    // Skip MPI communication in y dimension. Internal memory copies between boxes on same proc done in separate task
    if (iswap == 2 || iswap == 3) {         // iswap=2 swap in y -ve direction, iswap=3 swap in y +ve direction
      continue;
    }

    /* pack buffer */

    // swapnum as index to communicate sentinels as depend on individual iterations of this loop
    // No longer have 1:1 relationship between boxes and tasks in this implementation
    #pragma oss task label (communicate_pack) \
                     in(initialIntegrateSentinels[box_id]) \
                     out(communicatePackSentinels[box_id][iswap]) \
                     firstprivate(atom, box_id, iswap)
    {
      for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {
        // Determine which buffers we are using in each communication
        int sendNeighbour = boxBufs[box_id].sendneigh[iswap];
        AtomBuffer *buf_send = &(boxBufs[box_id].bufs_send[COMMUNICATE_FUNCTION][box_layer_index][sendNeighbour]);

        // Added buffer size checks. Previously weren't necessary as buffers were reused from borders() (which ensured
        // ensured sufficient size) but now independent buffers used for each function
        if (boxBufs[box_id].comm_send_size[box_layer_index][iswap] > buf_send->maxsize) {
          buf_send->growsend(boxBufs[box_id].comm_send_size[box_layer_index][iswap]);
        }

        atom->pack_comm(boxBufs[box_id].sendnum[box_layer_index][iswap],
                        boxBufs[box_id].sendlists[box_layer_index][sendNeighbour], buf_send->buf, buf_send->pbc_any,
                        buf_send->pbc_x, buf_send->pbc_y, buf_send->pbc_z);
      }
    }

    // Always use MPI to send in this implementation, even in case where sending to self
    // Simplifies dependencies, e.g. unpack task does not have to depend on pack task
    #pragma oss task label (communicate_send) \
                     in(communicatePackSentinels[box_id][iswap]) \
                     out(communicateSendSentinel[box_id][iswap]) \
                     firstprivate(box_id, iswap)
    {
      for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {
        int sendNeighbour = boxBufs[box_id].sendneigh[iswap];
        AtomBuffer *buf_send = &boxBufs[box_id].bufs_send[COMMUNICATE_FUNCTION][box_layer_index][sendNeighbour];
        MPI_Comm *send_comm = &boxBufs[send_target_box_id[box_layer_index]].comm[COMMUNICATE_FUNCTION];

        // Calculate send tag. In range 0 to 29, equivalent to sendnum from other implementations
        int sendtag = box_layer_index * 10 + iswap;

        if (sizeof(MMD_float) == 4) {
          MPI_Send(buf_send->buf, boxBufs[box_id].comm_send_size[box_layer_index][iswap], MPI_FLOAT,
                    boxBufs[box_id].sendproc[iswap], sendtag, *send_comm);
        } else {
          MPI_Send(buf_send->buf, boxBufs[box_id].comm_send_size[box_layer_index][iswap], MPI_DOUBLE,
                    boxBufs[box_id].sendproc[iswap], sendtag, *send_comm);
        }
      }
    }

    // No in-dependencies on pack or send tasks. Recvs can be posted as soon as communicate starts
    #pragma oss task label (communicate_recv) \
                     in(initialIntegrateSentinels[box_id]) \
                     out(communicateRecvSentinels[box_id][iswap]) \
                     firstprivate(atom, box_id, iswap)
    {
      for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {
        int recvNeighbour = boxBufs[box_id].recvneigh[iswap];
        AtomBuffer *buf_recv = &boxBufs[box_id].bufs_recv[COMMUNICATE_FUNCTION][box_layer_index][recvNeighbour];
        MPI_Comm *recv_comm = &boxBufs[atom->box_id].comm[COMMUNICATE_FUNCTION];

        // Calculate recv tag. In range 0 to 29, equivalent to sendnum from other implementations
        int recvtag = box_layer_index * 10 + iswap;

        // Ensure recv buffer is large enough. Equivalent check for send buffer done in packing task
        if (boxBufs[box_id].comm_recv_size[box_layer_index][iswap] > buf_recv->maxsize) {
          buf_recv->growrecv(boxBufs[box_id].comm_recv_size[box_layer_index][iswap]);
        }

        if (sizeof(MMD_float) == 4) {
          MPI_Recv(buf_recv->buf, boxBufs[box_id].comm_recv_size[box_layer_index][iswap], MPI_FLOAT,
                    boxBufs[box_id].recvproc[iswap], recvtag, *recv_comm, MPI_STATUS_IGNORE);
        } else {
          MPI_Recv(buf_recv->buf, boxBufs[box_id].comm_recv_size[box_layer_index][iswap], MPI_DOUBLE,
                    boxBufs[box_id].recvproc[iswap], recvtag, *recv_comm, MPI_STATUS_IGNORE);
        }
      }
    }

    // Now do unpacks. Unlike in borders(), unpack function cannot realloc coordinates array. Packs and unpacks
    // operate on entirely independent sections of the array:
    //   - packs on x[0:nlocal-1]
    //   - unpacks on x[nlocal:nghost-1]
    // Every unpack also uses an independent section (calculated in borders), hence unpack can be merged with recv
    // task, no need for separate commutative tasks
    // TODO: Merging greatly reduces performance. Why?
    #pragma oss task label (communicate_unpack) \
                     in(communicateRecvSentinels[box_id][iswap]) \
                     out(communicateSentinels[box_id]) \
                     firstprivate(atom, box_id, iswap)
    {
      int recvNeighbour = boxBufs[box_id].recvneigh[iswap];
      for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {
        atom->unpack_comm(boxBufs[box_id].recvnum[box_layer_index][iswap],
                          boxBufs[box_id].firstrecv[box_layer_index][iswap],
                          boxBufs[box_id].bufs_recv[COMMUNICATE_FUNCTION][box_layer_index][recvNeighbour].buf);
      }
    }

  } // End of loop over number of swaps
} // End of function


// As above but single MPI_send/MPI_recv operation per neighbour.
// Data from each box layer is copied to single intermediate buffer which is used in the MPI_send.
// The recv task uses another intermediate buffer which is copied from to ensure individual recv buffers contain the
// data expected by unpack tasks.
//void Comm::communicate_blocking_single_message_neighbourtasks(Atom *atom, int box_id)
//{
//  // ID of box we're sending to and receiving from
//  int send_target_box_id[3], recv_target_box_id[3];
//
//  // No outer loop over box layers this implementation. All 3 layers captured in each task
//
//  // Determine box IDs to send to/receive from on all layers
//  layer_to_targets(atom, 0, &send_target_box_id[0], &recv_target_box_id[0]);
//  layer_to_targets(atom, 1, &send_target_box_id[1], &recv_target_box_id[1]);
//  layer_to_targets(atom, 2, &send_target_box_id[2], &recv_target_box_id[2]);
//
//  // DSM nswap (number of swaps) is a constant calculated in setup()
//  for(int iswap = 0; iswap < boxBufs[box_id].nswap; iswap++) {
//
//    // Skip MPI communication in y dimension. Internal memory copies between boxes on same proc done in separate task
//    if (iswap == 2 || iswap == 3) {         // iswap=2 swap in y -ve direction, iswap=3 swap in y +ve direction
//      continue;
//    }
//
//    /* pack buffer */
//
//    // swapnum as index to communicate sentinels as depend on individual iterations of this loop
//    // No longer have 1:1 relationship between boxes and tasks in this implementation
//    #pragma oss task label (communicate_pack) \
//                     in(initialIntegrateSentinels[box_id]) \
//                     out(communicatePackSentinels[box_id][iswap]) \
//                     firstprivate(atom, box_id, iswap)
//    {
//      for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {
//        // Determine which buffers we are using in each communication
//        int sendNeighbour = boxBufs[box_id].sendneigh[iswap];
//        AtomBuffer *buf_send = &(boxBufs[box_id].bufs_send[COMMUNICATE_FUNCTION][box_layer_index][sendNeighbour]);
//
//        // Added buffer size checks. Previously weren't necessary as buffers were reused from borders() (which ensured
//        // ensured sufficient size) but now independent buffers used for each function
//        if (boxBufs[box_id].comm_send_size[box_layer_index][iswap] > buf_send->maxsize) {
//          buf_send->growsend(boxBufs[box_id].comm_send_size[box_layer_index][iswap]);
//        }
//
//        atom->pack_comm(boxBufs[box_id].sendnum[box_layer_index][iswap],
//                        boxBufs[box_id].sendlists[box_layer_index][sendNeighbour], buf_send->buf, buf_send->pbc_any,
//                        buf_send->pbc_x, buf_send->pbc_y, buf_send->pbc_z);
//      }
//    }
//
//    // Memory copy atoms from all box layers into a single large send buffer.
//    // Only 1 MPI_send call per task rather than 3
//    #pragma oss task label (communicate_send) \
//                     in(communicatePackSentinels[box_id][iswap]) \
//                     out(communicateSendSentinel[box_id][iswap]) \
//                     firstprivate(box_id, iswap)
//    {
//      int sendNeighbour = boxBufs[box_id].sendneigh[iswap];
//
//      // Send size for all 3 box layers
//      int total_send_size = boxBufs[box_id].comm_send_size[0][iswap] +
//                            boxBufs[box_id].comm_send_size[1][iswap] +
//                            boxBufs[box_id].comm_send_size[2][iswap];
//
//      // Resize send buffer if needed
//      if (total_send_size > boxBufs[box_id].buf_send_single.maxsize) {
//        boxBufs[box_id].buf_send_single.growsend(total_send_size);
//      }
//
//      // Copy data into single send buffer
//      int send_count = 0;
//      for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {
//        AtomBuffer *buf_send = &boxBufs[box_id].bufs_send[COMMUNICATE_FUNCTION][box_layer_index][sendNeighbour];
//        memcpy(&boxBufs[box_id].buf_send_single.buf[send_count], buf_send, boxBufs[box_id].comm_send_size[box_layer_index][iswap]);
//        send_count += boxBufs[box_id].comm_send_size[box_layer_index][iswap];
//      }
//
//      // Now perform single send
//      MPI_Comm *send_comm = &boxBufs[0].comm[COMMUNICATE_FUNCTION];
//      int sendtag = iswap; // Only single send/recv pair, do not need to include box_layer_index in calculation
//
//      if (sizeof(MMD_float) == 4) {
//        MPI_Send(boxBufs[box_id].buf_send_single.buf, total_send_size, MPI_FLOAT,
//                 boxBufs[box_id].sendproc[iswap], sendtag, *send_comm);
//      } else {
//        MPI_Send(boxBufs[box_id].buf_send_single.buf, total_send_size, MPI_DOUBLE,
//                 boxBufs[box_id].sendproc[iswap], sendtag, *send_comm);
//      }
//    }
//
//    // No in-dependencies on pack or send tasks. Recvs can be posted as soon as communicate starts
//    #pragma oss task label (communicate_recv) \
//                     in(initialIntegrateSentinels[box_id]) \
//                     out(communicateRecvSentinels[box_id][iswap]) \
//                     firstprivate(atom, box_id, iswap)
//    {
//      // Perform single recv
//      int recvNeighbour = boxBufs[box_id].recvneigh[iswap];
//      MPI_Comm *recv_comm = &boxBufs[0].comm[COMMUNICATE_FUNCTION];
//      int recvtag = iswap;
//
//      // Recv size for all 3 box layers
//      int total_recv_size = boxBufs[box_id].comm_recv_size[0][iswap] +
//                            boxBufs[box_id].comm_recv_size[1][iswap] +
//                            boxBufs[box_id].comm_recv_size[2][iswap];
//
//      // Resize recv buffer if needed
//      if (total_recv_size > boxBufs[box_id].buf_recv_single.maxsize) {
//        boxBufs[box_id].buf_recv_single.growrecv(total_recv_size);
//      }
//
//      if (sizeof(MMD_float) == 4) {
//        MPI_Recv(boxBufs[box_id].buf_recv_single.buf, total_recv_size, MPI_FLOAT,
//                 boxBufs[box_id].recvproc[iswap], recvtag, *recv_comm, MPI_STATUS_IGNORE);
//      } else {
//        MPI_Recv(boxBufs[box_id].buf_recv_single.buf, total_recv_size, MPI_DOUBLE,
//                 boxBufs[box_id].recvproc[iswap], recvtag, *recv_comm, MPI_STATUS_IGNORE);
//      }
//
//      // Copy recv buffer into individual box layer buffers expected by the unpack task
//      int recv_count = 0;
//      for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {
//        AtomBuffer *buf_recv = &boxBufs[box_id].bufs_recv[COMMUNICATE_FUNCTION][box_layer_index][recvNeighbour];
//
//        // Ensure recv buffer is large enough. Equivalent check for send buffer done in packing task
//        if (boxBufs[box_id].comm_recv_size[box_layer_index][iswap] > buf_recv->maxsize) {
//          buf_recv->growrecv(boxBufs[box_id].comm_recv_size[box_layer_index][iswap]);
//        }
//
//        memcpy(buf_recv->buf, &boxBufs[box_id].buf_recv_single.buf[recv_count], boxBufs[box_id].comm_recv_size[box_layer_index][iswap]);
//        recv_count += boxBufs[box_id].comm_recv_size[box_layer_index][iswap];
//      }
//    }
//
//    // Now do unpacks. Unlike in borders(), unpack function cannot realloc coordinates array. Packs and unpacks
//    // operate on entirely independent sections of the array:
//    //   - packs on x[0:nlocal-1]
//    //   - unpacks on x[nlocal:nghost-1]
//    // Every unpack also uses an independent section (calculated in borders), hence unpack can be merged with recv
//    // task, no need for separate commutative tasks
//    // TODO: Merging greatly reduces performance. Why?
//    #pragma oss task label (communicate_unpack) \
//                     in(communicateRecvSentinels[box_id][iswap]) \
//                     out(communicateSentinels[box_id]) \
//                     firstprivate(atom, box_id, iswap)
//    {
//      int recvNeighbour = boxBufs[box_id].recvneigh[iswap];
//      for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {
//        atom->unpack_comm(boxBufs[box_id].recvnum[box_layer_index][iswap],
//                          boxBufs[box_id].firstrecv[box_layer_index][iswap],
//                          boxBufs[box_id].bufs_recv[COMMUNICATE_FUNCTION][box_layer_index][recvNeighbour].buf);
//      }
//    }
//
//  } // End of loop over number of swaps
//} // End of function

/* reverse communication of atom info every timestep */

void Comm::reverse_communicate(Atom &atom)
{
  // DSM Disable for multibox
  printf("Rank %d called comm.reverse_communicate(). Aborting.", me);
  MPI_Abort(MPI_COMM_WORLD, 1);

//  int iswap;
//  MMD_float* buf;
//  MPI_Request request;
//  MPI_Status status;
//
//  for(iswap = nswap - 1; iswap >= 0; iswap--) {
//
//    /* pack buffer */
//
//    // #pragma omp barrier
//    atom.pack_reverse(recvnum[iswap], firstrecv[iswap], buf_send);
//
//    // #pragma omp barrier
//    /* exchange with another proc
//       if self, set recv buffer to send buffer */
//
//    if(sendproc[iswap] != me) {
//
//      #pragma omp master
//      {
//        if(sizeof(MMD_float) == 4) {
//          MPI_Irecv(buf_recv, reverse_recv_size[iswap], MPI_FLOAT,
//          sendproc[iswap], 0, MPI_COMM_WORLD, &request);
//          MPI_Send(buf_send, reverse_send_size[iswap], MPI_FLOAT,
//          recvproc[iswap], 0, MPI_COMM_WORLD);
//        } else {
//          MPI_Irecv(buf_recv, reverse_recv_size[iswap], MPI_DOUBLE,
//          sendproc[iswap], 0, MPI_COMM_WORLD, &request);
//          MPI_Send(buf_send, reverse_send_size[iswap], MPI_DOUBLE,
//          recvproc[iswap], 0, MPI_COMM_WORLD);
//        }
//        MPI_Wait(&request, &status);
//      }
//      buf = buf_recv;
//    } else buf = buf_send;
//
//    /* unpack buffer */
//
//    #pragma omp barrier
//    atom.unpack_reverse(sendnum[iswap], sendlist[iswap], buf);
//    // #pragma omp barrier
//  }
}

/* exchange:
   move atoms to correct proc boxes
   send out atoms that have left my box, receive ones entering my box
   this routine called before every reneighboring
   atoms exchanged with all 6 stencil neighbors
*/

void Comm::exchange_pack(Atom* atom, AtomBuffer bufs_send[3][8], AtomBuffer* internal_buf_send_up, AtomBuffer* internal_buf_send_down) {
  // Skip exchange entirely if running with only one process and one box
  // Not just an optimisation; algorithm breaks otherwise.
  if (threads->mpi_num_threads == 1 && boxes_per_process == 1) {
    return;
  }

  // DSM: This is the atom position array
  MMD_float* x = atom->x;

  // DSM: Do single loop over all atoms placing each in the appropriate send buffer. Replaces existing code's approach
  // of multiple loops over all atoms (one per dimension/direction)
  // DSM Checking whether each local atom is in our local box (not a slab boundary check as in borders()). Mark for
  // exchange with neighbour if outside box
  // Save current count of local atoms. Decrement in loop as atoms are marked for sending and deleted from local box.
  // A while loop as i is not incremented on every iteration, only when an atom is _not_ marked for sending (follows
  // the implementation in exchange_all)
  int nlocal = atom->nlocal;
  int i = 0, send_flag = 0;
  while (i < nlocal) {
    // Get atom's current coordinates
    MMD_float xcoord, ycoord, zcoord;
    xcoord = x[i * PAD];
    ycoord = x[i * PAD + 1];
    zcoord = x[i * PAD + 2];
    // DSM Initially mark atom as not to be sent (0 indicates atom _is_ sent. 1 indicates it is _not_).
    send_flag = 1;
    // Indices into array of send/recv buffers
    int atom_layer = SAME_LAYER;         // box_id (same layer), box_id-1 (below), box_id+1 (above)
    int neighbour = UP;                  // neighbours in clockwise direction from neighbour directly to the "back" or "front" (z -ve/+ve directions)

    // 26 exchanges = 8 process neighbours * 3 box layers + 2 internal exchanges between boxes on same process

    // First, check which layer this atom is on
    if (ycoord < atom->box.ylo) {
      // Outside in -ve y direction: box_id-1 / "lower" layer.
      atom_layer = LOWER_LAYER;
      // Outside box, necessarily mark as to-be-sent
      send_flag = 0;
      // At first, assume sending to self (internal box below)
      neighbour = ME;
    } else if (ycoord >= atom->box.yhi) {
      // Outside in +ve y direction: box_id+1 / "upper" layer
      atom_layer = UPPER_LAYER;
      // Outside box, necessarily mark as to-be-sent
      send_flag = 0;
      // At first, assume sending to self (internal box above)
      neighbour = ME;
    } else {
      // Not outside in either y direction: box_id / "same" layer
      atom_layer = SAME_LAYER;
    }

    if (xcoord < atom->box.xlo) {  // Is atom outside box in x -ve direction?
      // Yes, mark as to-be-sent and decide which left-hand side neighbour it is going to.
      send_flag = 0;

      // Outside in z -ve direction?
      if (zcoord < atom->box.zlo) {
        // x-ve and z-ve: send to "bottom-left"
        neighbour = BOTTOM_LEFT;
      } else if (zcoord >= atom->box.zhi) {
        // x-ve and z+ve: send to "top-left"
        neighbour = TOP_LEFT;
      } else {
        // x-ve but within z: send "left"
        neighbour = LEFT;
      }

    } else if (xcoord >= atom->box.xhi) { // Is atom outside box in x +ve direction?
      // Yes, mark as to-be-sent and decide which right-hand side neighbour it is going to.
      send_flag = 0;

      // Outside in z -ve direction?
      if (zcoord < atom->box.zlo) {
        // x+ve and z-ve: send to "bottom-right"
        neighbour = BOTTOM_RIGHT;
      } else if (zcoord >= atom->box.zhi) {
        // x+ve and z+ve: send to "top-right"
        neighbour = TOP_RIGHT;
      } else {
        // x-ve but within z: send "left"
        neighbour = RIGHT;
      }

    } else if (zcoord < atom->box.zlo) { // Is atom outside box in z -ve direction?
      // Within x but outside z-ve. Send down
      send_flag = 0;
      neighbour = DOWN;

    } else if (zcoord >= atom->box.zhi) { // Is atom outside box in z +ve direction?
      // Within x but outside z+ve. Send up
      send_flag = 0;
      neighbour = UP;

    }

    // All directions accounted for. If this atom marked for sending, copy into appropriate buffer
    AtomBuffer *buf_send;
    if (send_flag == 0) {

      // Special case: use separate internal buffers if this is a copy between boxes on the same process
      if (neighbour == ME) {
        if (atom_layer == LOWER_LAYER) {
          buf_send = internal_buf_send_down;
        } else if (atom_layer == UPPER_LAYER) {
          buf_send = internal_buf_send_up;
        } else {
          printf("ERROR: Rank %d attempting to exchange with its own box layer. Aborting.", me);
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
      }
        // Otherwise use the buffer allocated for this specific off-process neighbour
      else {
        buf_send = &bufs_send[atom_layer][neighbour];
      }

      // Increment number of atoms to be sent by this buffer
      buf_send->natoms++;
      // Increase size of buffer if cannot hold the current number of atoms
      // *7 as sending 7 elements per atom: x,y,z coordinates; x,y,z velocities; and flag for atom type
      if (buf_send->natoms * 7 > buf_send->maxsize) {
        buf_send->growsend(buf_send->natoms * 7);
      }

      // Now perform the packing
      int buf_send_index = (buf_send->natoms - 1) * 7;
      atom->pack_exchange(i, &(buf_send->buf[buf_send_index]));
      // Delete atom from local box by copying over it with last atom. Do not increment i in this branch - immediately
      // check the copied over atom in next iteration. Same implementation as in exchange_all()
      atom->copy(nlocal - 1, i);
      --nlocal;

    } else {
      // Atom was not marked for sending. Check the next one.
      ++i;
    }

  } // end of loop over atoms

  // Update count of atoms in this box with new total following removal of atoms marked for exchange
  atom->nlocal = nlocal;
}

// DSM: exchange() is always called before borders()
void Comm::exchange(Atom* atoms[]) {
  // Select between blocking and non-blocking communication modes
  if (nonblocking_enabled) {
    exchange_nonblocking(atoms);

    // After exchanges over all boxes completed with their off-process neighbours,
    // perform internal exchanges with all boxes on this process.
    for(int box_index = 0; box_index < atoms[0]->boxes_per_process; ++box_index) {
      exchange_internal(*atoms[box_index], box_index);
      atoms[box_index]->pbc();
    }
  } else {
    // Call once per box instance
    for(int box_index = 0; box_index < atoms[0]->boxes_per_process; ++box_index) {
      // Pack all atoms into relevant buffers to be sent to new owners and update local atom count
      #pragma oss task label(exchange_pack) in(initialIntegrateSentinels[box_index]) out(exchangePackSentinels[box_index]) firstprivate(box_index)
      exchange_pack(atoms[box_index],
                    boxBufs[box_index].bufs_send[EXCHANGE_FUNCTION],
                    &boxBufs[box_index].internal_buf_send_up[EXCHANGE_FUNCTION],
                    &boxBufs[box_index].internal_buf_send_down[EXCHANGE_FUNCTION]);
    }

    for(int box_index = 0; box_index < atoms[0]->boxes_per_process; ++box_index) {

        //#pragma oss task label(exchange_blocking) in(exchangePackSentinels[box_index]) out(exchangeSentinels[box_index]) firstprivate(box_index)
        //exchange_blocking(*atoms[box_index], box_index);

        // Tasks created within these functions
        //exchange_blocking_neighbourtasks(atoms[box_index], box_index);
        exchange_nonblocking_neighbourtasks_tampi_iwaitall(atoms[box_index], box_index);
    }

    for(int box_index = 0; box_index < atoms[0]->boxes_per_process; ++box_index) {
      // Generate tasks for internal swaps between boxes
      // Represent the pushing of data to both of this box's neighbours (i.e. "sends")
      // Unpacking of data done in separate loop following (i.e. "recvs")
      #pragma oss task label(exchange_internal_send) \
                     in(exchangePackSentinels[box_index]) \
                     out(exchangeInternalSendSentinels[box_index]) firstprivate(box_index)
      {
        // To lower layer (usually box_id-1)
        int box_neighbour = atoms[box_index]->boxneigh_negative;
        exchange_internal_send(&boxBufs[box_index].internal_buf_send_down[EXCHANGE_FUNCTION],
                               &boxBufs[box_neighbour].internal_buf_recv_down[EXCHANGE_FUNCTION]);

        // To upper layer (usually box_id+1)
        box_neighbour = atoms[box_index]->boxneigh_positive;
        exchange_internal_send(&boxBufs[box_index].internal_buf_send_up[EXCHANGE_FUNCTION],
                               &boxBufs[box_neighbour].internal_buf_recv_up[EXCHANGE_FUNCTION]);
      }
    }

    // Receive half of internal swaps between boxes
    for(int box_index = 0; box_index < atoms[0]->boxes_per_process; ++box_index) {
      // Wait for my neighbours' send tasks to finish before unpacking
      // Also wait for my box's exchange task to complete. Avoids race condition with two tasks attempting to unpack
      // into position array concurrently
      // TODO: Break depedencies on non-internal sends
      // Is dependency on pack needed? Yes, what if unpack tries to grow array while a pack task is running
      // commutative with other unpack tasks to avoid concurrent unpacks into the same buffer
      #pragma oss task label(exchange_internal_recv) commutative(atoms[box_index]->nlocal) \
                     in(exchangePackSentinels[box_index]) \
                     in(exchangeInternalSendSentinels[atoms[box_index]->boxneigh_positive]) \
                     in(exchangeInternalSendSentinels[atoms[box_index]->boxneigh_negative]) \
                     in(exchangeSentinels[box_index]) \
                     out(exchangeInternalRecvSentinels[box_index]) firstprivate(box_index)
      {
        exchange_internal_recv(atoms[box_index], &boxBufs[box_index].internal_buf_recv_down[EXCHANGE_FUNCTION]);
        exchange_internal_recv(atoms[box_index], &boxBufs[box_index].internal_buf_recv_up[EXCHANGE_FUNCTION]);
      }
    }

    // TODO: Merge with sort task?
    /* enforce PBC */
    // Moved this to from the beginning of original exchange to here as would break new checks for 26-way communication
    // e.g. An x coordinate of -1 would get corrected to +ve box length and erroneously be sent *RIGHT* instead of *LEFT*
    // Correcting after all atoms have been received and unpacked addresses this
    for(int box_index = 0; box_index < atoms[0]->boxes_per_process; ++box_index) {
      #pragma oss task label(atom->pbc) \
                       in(exchangeSentinels[box_index]) \
                       in(exchangeInternalRecvSentinels[box_index]) \
                       out(exchangePBCSentinels[box_index]) \
                       firstprivate(box_index)
      atoms[box_index]->pbc();
    }

  } // End of blocking branch
} // End of function

void Comm::exchange_nonblocking(Atom* atoms[])
{
  // Skip exchange entirely if running with only one process and one box
  // Not just an optimisation; algorithm breaks otherwise.
  if (threads->mpi_num_threads == 1 && atoms[0]->boxes_per_process == 1) {
    return;
  }

  // Arrays for isend/irecv requests
  // nswaps communications per box layer, 3 box layers per box on this process.
  // Internal swaps in y-dimension are skipped => subtract 6 per box (2 y-swaps per box layer)
  // /2 as swaps are split between -ve and +ve communications
  // nswaps and boxes_per_process constant per box
  int nrequests = (boxBufs[0].nswap * 3 * atoms[0]->boxes_per_process - 6 * atoms[0]->boxes_per_process)/2;
  int req_i = 0;
  MPI_Request recv1_requests[nrequests];
  MPI_Request recv2_requests[nrequests];
  MPI_Request send1_requests[nrequests];
  MPI_Request send2_requests[nrequests];
  // Buffer sizes
  int nrecvs1[nrequests];
  int nrecvs2[nrequests];
  int nsends1[nrequests];
  int nsends2[nrequests];

  MPI_Comm *send_comm, *recv_comm;

  for(int box_id = 0; box_id < atoms[0]->boxes_per_process; ++box_id) {

    Atom& atom = *atoms[box_id];

    if (do_safeexchange)
      return exchange_all(atom);

    int i, idim, nlocal, send_flag;
    MMD_float *x;
    AtomBuffer *buf_send = NULL, *buf_recv1 = NULL, *buf_recv2 = NULL;
    int swapnum = 200; // to ensure tag is unique accross communicate, borders and exchange

    // DSM Multibox: Reference buffers/variables specific to this box.
    int &maxnlocal = boxBufs[box_id].maxnlocal;
    int &maxthreads = boxBufs[box_id].maxthreads;
    int *&nsend_thread = boxBufs[box_id].nsend_thread;
    int *&nrecv_thread = boxBufs[box_id].nrecv_thread;
    int *&nholes_thread = boxBufs[box_id].nholes_thread;
    int *&maxsend_thread = boxBufs[box_id].maxsend_thread;
    int **&exc_sendlist_thread = boxBufs[box_id].exc_sendlist_thread;
    AtomBuffer (&bufs_send)[3][8] = boxBufs[box_id].bufs_send[EXCHANGE_FUNCTION];
    AtomBuffer (&bufs_recv)[3][8] = boxBufs[box_id].bufs_recv[EXCHANGE_FUNCTION];
    AtomBuffer &internal_buf_send_up = boxBufs[box_id].internal_buf_send_up[EXCHANGE_FUNCTION];
    AtomBuffer &internal_buf_send_down = boxBufs[box_id].internal_buf_send_down[EXCHANGE_FUNCTION];
    AtomBuffer &internal_buf_recv_up = boxBufs[box_id].internal_buf_recv_up[EXCHANGE_FUNCTION];
    AtomBuffer &internal_buf_recv_down = boxBufs[box_id].internal_buf_recv_down[EXCHANGE_FUNCTION];

    // DSM: This is the atom position array
    x = atom.x;

    // DSM: Do single loop over all atoms placing each in the appropriate send buffer. Replaces existing code's approach
    // of multiple loops over all atoms (one per dimension/direction)
    // DSM Checking whether each local atom is in our local box (not a slab boundary check as in borders()). Mark for
    // exchange with neighbour if outside box
    // Save current count of local atoms. Decrement in loop as atoms are marked for sending and deleted from local box.
    // A while loop as i is not incremented on every iteration, only when an atom is _not_ marked for sending (follows
    // the implementation in exchange_all)
    nlocal = atom.nlocal;
    i = 0;
    while (i < nlocal) {
      // Get atom's current coordinates
      MMD_float xcoord, ycoord, zcoord;
      xcoord = x[i * PAD];
      ycoord = x[i * PAD + 1];
      zcoord = x[i * PAD + 2];
      // DSM Initially mark atom as not to be sent (0 indicates atom _is_ sent. 1 indicates it is _not_).
      send_flag = 1;
      // Indices into array of send/recv buffers
      int atom_layer = SAME_LAYER;         // box_id (same layer), box_id-1 (below), box_id+1 (above)
      int neighbour = UP;                  // neighbours in clockwise direction from neighbour directly to the "back" or "front" (z -ve/+ve directions)

      // 26 exchanges = 8 process neighbours * 3 box layers + 2 internal exchanges between boxes on same process

      // First, check which layer this atom is on
      if (ycoord < atom.box.ylo) {
        // Outside in -ve y direction: box_id-1 / "lower" layer.
        atom_layer = LOWER_LAYER;
        // Outside box, necessarily mark as to-be-sent
        send_flag = 0;
        // At first, assume sending to self (internal box below)
        neighbour = ME;
      } else if (ycoord >= atom.box.yhi) {
        // Outside in +ve y direction: box_id+1 / "upper" layer
        atom_layer = UPPER_LAYER;
        // Outside box, necessarily mark as to-be-sent
        send_flag = 0;
        // At first, assume sending to self (internal box above)
        neighbour = ME;
      } else {
        // Not outside in either y direction: box_id / "same" layer
        atom_layer = SAME_LAYER;
      }

      if (xcoord < atom.box.xlo) {  // Is atom outside box in x -ve direction?
        // Yes, mark as to-be-sent and decide which left-hand side neighbour it is going to.
        send_flag = 0;

        // Outside in z -ve direction?
        if (zcoord < atom.box.zlo) {
          // x-ve and z-ve: send to "bottom-left"
          neighbour = BOTTOM_LEFT;
        } else if (zcoord >= atom.box.zhi) {
          // x-ve and z+ve: send to "top-left"
          neighbour = TOP_LEFT;
        } else {
          // x-ve but within z: send "left"
          neighbour = LEFT;
        }

      } else if (xcoord >= atom.box.xhi) { // Is atom outside box in x +ve direction?
        // Yes, mark as to-be-sent and decide which right-hand side neighbour it is going to.
        send_flag = 0;

        // Outside in z -ve direction?
        if (zcoord < atom.box.zlo) {
          // x+ve and z-ve: send to "bottom-right"
          neighbour = BOTTOM_RIGHT;
        } else if (zcoord >= atom.box.zhi) {
          // x+ve and z+ve: send to "top-right"
          neighbour = TOP_RIGHT;
        } else {
          // x-ve but within z: send "left"
          neighbour = RIGHT;
        }

      } else if (zcoord < atom.box.zlo) { // Is atom outside box in z -ve direction?
        // Within x but outside z-ve. Send down
        send_flag = 0;
        neighbour = DOWN;

      } else if (zcoord >= atom.box.zhi) { // Is atom outside box in z +ve direction?
        // Within x but outside z+ve. Send up
        send_flag = 0;
        neighbour = UP;

      }

      // All directions accounted for. If this atom marked for sending, copy into appropriate buffer
      AtomBuffer *buf_send;
      if (send_flag == 0) {

        // Special case: use separate internal buffers if this is a copy between boxes on the same process
        if (neighbour == ME) {
          if (atom_layer == LOWER_LAYER) {
            buf_send = &internal_buf_send_down;
          } else if (atom_layer == UPPER_LAYER) {
            buf_send = &internal_buf_send_up;
          } else {
            printf("ERROR: Rank %d attempting to exchange with its own box layer. Aborting.", me);
            MPI_Abort(MPI_COMM_WORLD, 1);
          }
        }
          // Otherwise use the buffer allocated for this specific off-process neighbour
        else {
          buf_send = &bufs_send[atom_layer][neighbour];
        }

        // Increment number of atoms to be sent by this buffer
        buf_send->natoms++;
        // Increase size of buffer if cannot hold the current number of atoms
        // *7 as sending 7 elements per atom: x,y,z coordinates; x,y,z velocities; and flag for atom type
        if (buf_send->natoms * 7 > buf_send->maxsize) {
          buf_send->growsend(buf_send->natoms * 7);
        }

        // Now perform the packing
        int buf_send_index = (buf_send->natoms - 1) * 7;
        atom.pack_exchange(i, &(buf_send->buf[buf_send_index]));
        // Delete atom from local box by copying over it with last atom. Do not increment i in this branch - immediately
        // check the copied over atom in next iteration. Same implementation as in exchange_all()
        atom.copy(nlocal - 1, i);
        --nlocal;

      } else {
        // Atom was not marked for sending. Check the next one.
        ++i;
      }

    } // end of loop over atoms

    // Update count of atoms in this box with new total following removal of atoms marked for exchange
    atom.nlocal = nlocal;

    // Loop over box layers: box_id, box_id-1, box_id+1
    for (int box_layer_index = 0; box_layer_index < 3; box_layer_index++) {
      // DSM: Additionally, now loop over 5 dimensions (sides + corners) rather than only 3 (sides)
      for (idim = 0; idim < 5; idim++) {

        // DSM: Do not perform MPI communication in the y / idim==1 dimension
        // Instead, perform memory copy into neighbouring box on this proc's recv buffer.
        // The exchange_internal function, called following completion of all exchange operations, will then pull this
        // data into each box's atoms list.
        if (idim == 1) {
          // Only perform this copy step once
          if (box_layer_index == 0) {

            // First exchange to lower layer (usually box_id-1 but periodic)
            int box_neighbour = atom.boxneigh_negative;
            buf_send = &internal_buf_send_down;
            buf_recv1 = &(boxBufs[box_neighbour].internal_buf_recv_down[EXCHANGE_FUNCTION]);
            int nsend = buf_send->natoms * 7; // *7 as 7 MMD_floats per atom
            // Increase size of receive buffers if receiving more than can be held
            if (nsend > buf_recv1->maxsize) {
              buf_recv1->growrecv(nsend);
            }
            // Push to neighbour.
            memcpy(buf_recv1->buf, buf_send->buf, nsend * sizeof(MMD_float));
            // Set receiving buffer sizes
            buf_recv1->natoms = buf_send->natoms;
            // Reset atom counter for next use of this buffer
            buf_send->natoms = 0;

            // Second exchange to upper layer (usually box_id+1)
            box_neighbour = atom.boxneigh_positive;
            buf_send = &internal_buf_send_up;
            buf_recv2 = &(boxBufs[box_neighbour].internal_buf_recv_up[EXCHANGE_FUNCTION]);
            nsend = buf_send->natoms * 7;
            if (nsend > buf_recv2->maxsize) {
              buf_recv2->growrecv(nsend);
            }
            memcpy(buf_recv2->buf, buf_send->buf, nsend * sizeof(MMD_float));
            buf_recv2->natoms = buf_send->natoms;
            buf_send->natoms = 0;

          }
          continue; // Skip to next dimension
        }

        /* only exchange if more than one proc in this dimension */
        /*if (idim < 3) { // Only applies to non-corner, non-y dimensions
          if (procgrid[idim] == 1) {
            continue;
          }
        }*/

        //#pragma omp barrier
        //#pragma omp master
        {
          int sendDirection, send_target_box_id, recv_target_box_id, sendtag, recvtag;
          layer_to_targets(&atom, box_layer_index, &send_target_box_id, &recv_target_box_id);

          // Determine which buffers we are using in first communication
          // Map from idim to UP, TOP_RIGHT, etc. neighbour ID
          sendDirection = dim_to_neigh(idim, 0);
          buf_send = &bufs_send[box_layer_index][sendDirection];
          nsends1[req_i] = buf_send->natoms * 7; // *7 to represent how many MMD_floats are due to be sent
          buf_send->natoms = 0; // Reset atom counter for next use of this buffer
          send_comm = &boxBufs[send_target_box_id].comm[EXCHANGE_FUNCTION];
          recv_comm = &boxBufs[atom.box_id].comm[EXCHANGE_FUNCTION];
          sendtag = swapnum;
          recvtag = sendtag;
          ++swapnum;

          /* send/recv atoms in both directions
             only if neighboring procs are different */

          nrecvs1[req_i] = 0;
          nrecvs2[req_i] = 0;

          // DSM First send/recv number of elements we're expecting to/from every neighbour
          // Differs from blocking implementation in that data is only sent after size exchanges with every neighbour
          MPI_Irecv(&nrecvs1[req_i], 1, MPI_INT, procneigh[idim][1], recvtag,
                    *recv_comm, &recv1_requests[req_i]);
          MPI_Isend(&nsends1[req_i], 1, MPI_INT, procneigh[idim][0], sendtag,
                    *send_comm, &send1_requests[req_i]);

          // Determine which buffers we are using in second communication
          sendDirection = dim_to_neigh(idim, 1);
          buf_send = &bufs_send[box_layer_index][sendDirection];
          nsends2[req_i] = buf_send->natoms * 7; // *7 to represent how many MMD_floats are due to be sent
          buf_send->natoms = 0; // Reset atom counter for next use of this buffer
          sendtag = swapnum;
          recvtag = sendtag;
          ++swapnum;
          // TODO: Do we need to swap send/recv target box IDs and communicators here?

          // Update: Unconditionally do second swap now. Required by 26-way communication
          MPI_Irecv(&nrecvs2[req_i], 1, MPI_INT, procneigh[idim][0], recvtag,
                    *recv_comm, &recv2_requests[req_i]);
          MPI_Isend(&nsends2[req_i], 1, MPI_INT, procneigh[idim][1], sendtag,
                    *send_comm, &send2_requests[req_i]);

          req_i++;
        } // end of omp master region
      } // end of loop over dims
    } // end of loop over box layers
  } // end of loop over boxes per process

  // All buffer size messages posted to all neighbours, wait for completion
  // Using req_i instead of nrequests here as we might not be using the maximum length of the array, e.g. if skipped
  // any exchanges due to procgrid==1. req_i has the actual number of exchanges
  MPI_Waitall(req_i, recv1_requests, MPI_STATUSES_IGNORE);
  MPI_Waitall(req_i, send1_requests, MPI_STATUSES_IGNORE);
  MPI_Waitall(req_i, recv2_requests, MPI_STATUSES_IGNORE);
  MPI_Waitall(req_i, send2_requests, MPI_STATUSES_IGNORE);

  // Now post the actual data
  req_i = 0;
  // Individual indices for buffers as sends/recv are conditional (do not post messages if buffer sizes are 0)
  int recv1_req_i = 0, send1_req_i = 0, recv2_req_i = 0, send2_req_i = 0;
  for(int box_id = 0; box_id < atoms[0]->boxes_per_process; ++box_id) {

    Atom& atom = *atoms[box_id];
    AtomBuffer (&bufs_send)[3][8] = boxBufs[box_id].bufs_send[EXCHANGE_FUNCTION];
    AtomBuffer (&bufs_recv)[3][8] = boxBufs[box_id].bufs_recv[EXCHANGE_FUNCTION];
    int swapnum = 200; // to ensure tag is unique accross communicate, borders and exchange

    for (int box_layer_index = 0; box_layer_index < 3; box_layer_index++) {
      for (int idim = 0; idim < 5; idim++) {

        AtomBuffer *buf_send = NULL, *buf_recv1 = NULL, *buf_recv2 = NULL;

        // Skip y-dimension and any other non-corner dimensions with only a single process
        if (idim == 1) { continue; }
        //if (idim < 3 && procgrid[idim] == 1) { continue; }

        // Redetermine buffers and directions
        int sendDirection, recvDirection, send_target_box_id, recv_target_box_id, sendtag, recvtag;
        layer_to_targets(&atom, box_layer_index, &send_target_box_id, &recv_target_box_id);

        sendDirection = dim_to_neigh(idim, 0);
        buf_send = &bufs_send[box_layer_index][sendDirection];
        recvDirection = dim_to_neigh(idim, 1);
        buf_recv1 = &bufs_recv[box_layer_index][recvDirection];
        send_comm = &boxBufs[send_target_box_id].comm[EXCHANGE_FUNCTION];
        recv_comm = &boxBufs[atom.box_id].comm[EXCHANGE_FUNCTION];
        sendtag = swapnum;
        recvtag = sendtag;
        ++swapnum;

        // Increase size of receive buffer if receiving more than it can hold
        if (nrecvs1[req_i] > buf_recv1->maxsize) {
          buf_recv1->growrecv(nrecvs1[req_i]);
        }

        // Perform exchange with first neighbour
        // DSM TODO: Could we not make this a static check? #define MMD_MPI_FLOAT as MPI_FLOAT or MPI_DOUBLE?
        if (sizeof(MMD_float) == 4) {
          if (nrecvs1[req_i] != 0) {
            MPI_Irecv(buf_recv1->buf, nrecvs1[req_i], MPI_FLOAT, procneigh[idim][1], recvtag,
                      *recv_comm, &recv1_requests[recv1_req_i]);
            recv1_req_i++;
          }

          if (nsends1[req_i] != 0) {
            MPI_Isend(buf_send->buf, nsends1[req_i], MPI_FLOAT, procneigh[idim][0], sendtag,
                      *send_comm, &send1_requests[send1_req_i]);
            send1_req_i++;
          }
        } else {
          if (nrecvs1[req_i] != 0) {
            MPI_Irecv(buf_recv1->buf, nrecvs1[req_i], MPI_DOUBLE, procneigh[idim][1], recvtag,
                      *recv_comm, &recv1_requests[recv1_req_i]);
            recv1_req_i++;
          }

          if (nsends1[req_i] != 0) {
            MPI_Isend(buf_send->buf, nsends1[req_i], MPI_DOUBLE, procneigh[idim][0], sendtag,
                      *send_comm, &send1_requests[send1_req_i]);
            send1_req_i++;
          }
        }

        // The second exchange
        sendDirection = dim_to_neigh(idim, 1);
        buf_send = &bufs_send[box_layer_index][sendDirection];
        recvDirection = dim_to_neigh(idim, 0);
        buf_recv2 = &bufs_recv[box_layer_index][recvDirection];
        sendtag = swapnum;
        recvtag = sendtag;
        ++swapnum;
        // TODO: Do we need to swap send/recv target box IDs and communicators here?

        if (nrecvs2[req_i] > buf_recv2->maxsize) {
          buf_recv2->growrecv(nrecvs2[req_i]);
        }

        if (sizeof(MMD_float) == 4) {
          if (nrecvs2[req_i] != 0) {
            MPI_Irecv(buf_recv2->buf, nrecvs2[req_i], MPI_FLOAT, procneigh[idim][0], recvtag,
                      *recv_comm, &recv2_requests[recv2_req_i]);
            recv2_req_i++;
          }

          if (nsends2[req_i] != 0) {
            MPI_Isend(buf_send->buf, nsends2[req_i], MPI_FLOAT, procneigh[idim][1], sendtag,
                      *send_comm, &send2_requests[send2_req_i]);
            send2_req_i++;
          }
        } else {
          if (nrecvs2[req_i] != 0) {
            MPI_Irecv(buf_recv2->buf, nrecvs2[req_i], MPI_DOUBLE, procneigh[idim][0], recvtag,
                      *recv_comm, &recv2_requests[recv2_req_i]);
            recv2_req_i++;
          }

          if (nsends2[req_i] != 0) {
            MPI_Isend(buf_send->buf, nsends2[req_i], MPI_DOUBLE, procneigh[idim][1], sendtag,
                      *send_comm, &send2_requests[send2_req_i]);
            send2_req_i++;
          }
        }

        req_i++;
      } // end of loop over dims
    } // end of loop over box layers
  } // end of loop over boxes per process

  // All buffers posted to all neighbours, wait for completion
  // Using individual req_is as we might not be using the maximum length of the array, e.g. if skipped any exchanges due
  // to buffer size == 0
  MPI_Waitall(recv1_req_i, recv1_requests, MPI_STATUSES_IGNORE);
  MPI_Waitall(send1_req_i, send1_requests, MPI_STATUSES_IGNORE);
  MPI_Waitall(recv2_req_i, recv2_requests, MPI_STATUSES_IGNORE);
  MPI_Waitall(send2_req_i, send2_requests, MPI_STATUSES_IGNORE);

  // Now unpack all buffers
  req_i = 0;
  for(int box_id = 0; box_id < atoms[0]->boxes_per_process; ++box_id) {

    Atom &atom = *atoms[box_id];
    AtomBuffer (&bufs_recv)[3][8] = boxBufs[box_id].bufs_recv[EXCHANGE_FUNCTION];

    for (int box_layer_index = 0; box_layer_index < 3; box_layer_index++) {
      for (int idim = 0; idim < 5; idim++) {

        // Skip y-dimension and any other non-corner dimensions with only a single process
        if (idim == 1) { continue; }
        //if (idim < 3 && procgrid[idim] == 1) { continue; }

        // Again, /7 as 7 MMD_floats per atom received (position, velocity and type)
        nrecvs1[req_i] /= 7;
        nrecvs2[req_i] /= 7;

        // Determine buffer pointers
        int recvDirection1 = dim_to_neigh(idim, 1);
        int recvDirection2 = dim_to_neigh(idim, 0);
        AtomBuffer *buf_recv1 = &bufs_recv[box_layer_index][recvDirection1];
        AtomBuffer *buf_recv2 = &bufs_recv[box_layer_index][recvDirection2];

        // Add new atoms from each buffer to end of list. Update local count of atoms after each unpack loop.
        for (int i = 0; i < nrecvs1[req_i]; ++i) {
          atom.unpack_exchange(atom.nlocal + i, &buf_recv1->buf[i * 7]);
        }
        atom.nlocal += nrecvs1[req_i];

        for (int i = 0; i < nrecvs2[req_i]; ++i) {
          atom.unpack_exchange(atom.nlocal + i, &buf_recv2->buf[i * 7]);
        }
        atom.nlocal += nrecvs2[req_i];

        req_i++;
      } // end of loop over dims
    } // end of loop over box layers
  } // end of loop over boxes per process
} // end of function

void Comm::exchange_blocking(Atom &atom, int box_id)
{
  // Skip exchange entirely if running with only one process and one box
  // Not just an optimisation; algorithm breaks otherwise.
  if (threads->mpi_num_threads == 1 && atom.boxes_per_process == 1) {
    return;
  }

  if(do_safeexchange)
    return exchange_all(atom);

  int i, idim, nlocal, nsend, nrecv1, nrecv2;
  MMD_float* x;
  AtomBuffer *buf_send = NULL, *buf_recv1 = NULL, *buf_recv2 = NULL;

  MPI_Request request;
  MPI_Status status;
  MPI_Comm *send_comm, *recv_comm;

  // DSM Multibox: Reference buffers/variables specific to this box.
  int& maxnlocal = boxBufs[box_id].maxnlocal;
  int& maxthreads = boxBufs[box_id].maxthreads;
  int*& nsend_thread = boxBufs[box_id].nsend_thread;
  int*& nrecv_thread = boxBufs[box_id].nrecv_thread;
  int*& nholes_thread = boxBufs[box_id].nholes_thread;
  int*& maxsend_thread = boxBufs[box_id].maxsend_thread;
  int**& exc_sendlist_thread = boxBufs[box_id].exc_sendlist_thread;
  AtomBuffer (&bufs_send)[3][8] = boxBufs[box_id].bufs_send[EXCHANGE_FUNCTION];
  AtomBuffer (&bufs_recv)[3][8] = boxBufs[box_id].bufs_recv[EXCHANGE_FUNCTION];
  AtomBuffer& internal_buf_recv_up = boxBufs[box_id].internal_buf_recv_up[EXCHANGE_FUNCTION];
  AtomBuffer& internal_buf_recv_down = boxBufs[box_id].internal_buf_recv_down[EXCHANGE_FUNCTION];

  // Loop over box layers: box_id, box_id-1, box_id+1
  // swapnum starts at 100 to ensure tag is unique across communicate, borders and exchange.
  // Messages from incorrect functions may match otherwise
  int swapnum = 100;
  for (int box_layer_index = 0; box_layer_index < 3; box_layer_index++) {
    // DSM: Additionally, now loop over 5 dimensions (sides + corners) rather than only 3 (sides)
    for(idim = 0; idim < 5; idim++) {

      // Skip y-dimension. This is internal swap between boxes on this process and handled outside this function
      if (idim == 1) { continue; }

      //#pragma omp barrier
      //#pragma omp master
      {
        int sendDirection, recvDirection, send_target_box_id, recv_target_box_id, sendtag, recvtag;
        layer_to_targets(&atom, box_layer_index, &send_target_box_id, &recv_target_box_id);

        // Determine which buffers we are using in first communication
        // Map from idim to UP, TOP_RIGHT, etc. neighbour ID
        sendDirection = dim_to_neigh(idim, 0);
        recvDirection = dim_to_neigh(idim, 1);
        buf_send = &bufs_send[box_layer_index][sendDirection];
        buf_recv1 = &bufs_recv[box_layer_index][recvDirection];
        send_comm = &boxBufs[send_target_box_id].comm[EXCHANGE_FUNCTION];
        recv_comm = &boxBufs[atom.box_id].comm[EXCHANGE_FUNCTION];
        sendtag = swapnum;
        recvtag = sendtag;
        ++swapnum;
        nsend = buf_send->natoms * 7; // *7 to represent how many MMD_floats are due to be sent
        buf_send->natoms = 0; // Reset atom counter for next use of this buffer

        /* send/recv atoms in both directions
           only if neighboring procs are different */

        nrecv1 = nrecv2 = 0;

        // DSM First send/recv number of elements we're expecting to/from neighbour, then irecv/send the actual data.
        // DSM: Fixed this to be Irecv/send pair. Had a chance of deadlock in the original code.
        MPI_Irecv(&nrecv1, 1, MPI_INT, procneigh[idim][1], recvtag, *recv_comm, &request);
        MPI_Send(&nsend, 1, MPI_INT, procneigh[idim][0], sendtag, *send_comm);
        MPI_Wait(&request, &status);

        // Increase size of receive buffer if receiving more than it can hold
        if(nrecv1 > buf_recv1->maxsize) {
          buf_recv1->growrecv(nrecv1);
        }

        // Perform exchange with first neighbour
        // DSM TODO: Could we not make this a static check? #define MMD_MPI_FLOAT as MPI_FLOAT or MPI_DOUBLE?
        if(sizeof(MMD_float) == 4) {
          if (nrecv1 != 0) {
            MPI_Irecv(buf_recv1->buf, nrecv1, MPI_FLOAT, procneigh[idim][1], recvtag, *recv_comm, &request);
          }

          if (nsend != 0) {
            MPI_Send(buf_send->buf, nsend, MPI_FLOAT, procneigh[idim][0], sendtag, *send_comm);
          }
        } else {
          if (nrecv1 != 0) {
            MPI_Irecv(buf_recv1->buf, nrecv1, MPI_DOUBLE, procneigh[idim][1], recvtag, *recv_comm, &request);
          }

          if (nsend != 0) {
            MPI_Send(buf_send->buf, nsend, MPI_DOUBLE, procneigh[idim][0], sendtag, *send_comm);
          }
        }

        if (nrecv1 != 0) {
          MPI_Wait(&request, &status); // End of first exchange
        }

        // Determine which buffers we are using in second communication
        sendDirection = dim_to_neigh(idim, 1);
        recvDirection = dim_to_neigh(idim, 0);
        buf_send = &bufs_send[box_layer_index][sendDirection];
        buf_recv2 = &bufs_recv[box_layer_index][recvDirection];
        sendtag = swapnum;
        recvtag = sendtag;
        ++swapnum;
        // TODO: Do we need to reverse send/recv target IDs and comms here?
        nsend = buf_send->natoms * 7; // *7 to represent how many MMD_floats are due to be sent
        buf_send->natoms = 0; // Reset atom counter for next use of this buffer

        // DSM: This branch is checking procgrid[idim]!=2. Can't ==1 as that is checked at start of loop (continue if
        // true) and can't <=0 as all dimensions need at least 1 processor.
        // DSM: i.e. We don't need to do a second send/recv exchange in the case where there are only 2 processors in this
        // dimension. Do not need to do exchanges in both directions every time. This communication is for exchanging
        // ownership of atoms, not updating atom coordinates - PBC do not apply.
        //if(procgrid[idim] > 2) {
        // Update: Unconditionally do second swap now. Required by 26-way communication
        MPI_Irecv(&nrecv2, 1, MPI_INT, procneigh[idim][0], recvtag, *recv_comm, &request);
        MPI_Send(&nsend, 1, MPI_INT, procneigh[idim][1], sendtag, *send_comm);
        MPI_Wait(&request, &status);

        // Increase size of receive buffer if receiving more than it can hold
        if(nrecv2 > buf_recv2->maxsize) {
          buf_recv2->growrecv(nrecv2);
        }

        // DSM: The second exchange,
        if(sizeof(MMD_float) == 4) {
          if (nrecv2 != 0) {
            MPI_Irecv(buf_recv2->buf, nrecv2, MPI_FLOAT, procneigh[idim][0], recvtag, *recv_comm, &request);
          }

          if (nsend != 0) {
            MPI_Send(buf_send->buf, nsend, MPI_FLOAT, procneigh[idim][1], sendtag, *send_comm);
          }
        } else {
          if (nrecv2 != 0) {
            MPI_Irecv(buf_recv2->buf, nrecv2, MPI_DOUBLE, procneigh[idim][0], recvtag, *recv_comm, &request);
          }

          if (nsend != 0) {
            MPI_Send(buf_send->buf, nsend, MPI_DOUBLE, procneigh[idim][1], sendtag, *send_comm);
          }
        }

        if (nrecv2 != 0) {
          MPI_Wait(&request, &status); // End of second exchange
        }
      }

      /* check incoming atoms to see if they are in my box
         if they are, add to my list */
      // DSM: Removed this check as all atoms received now must be in my box (otherwise they would not have been sent)

      // Again, /7 as 7 MMD_floats per atom received (position, velocity and type)
      nrecv1 /= 7;
      nrecv2 /= 7;

      // Add new atoms from each buffer to end of list. Update local count of atoms after each unpack loop.
      for (i = 0; i < nrecv1; ++i) {
        atom.unpack_exchange(atom.nlocal + i, &buf_recv1->buf[i * 7]);
      }
      atom.nlocal += nrecv1;

      for (i = 0; i < nrecv2; ++i) {
        atom.unpack_exchange(atom.nlocal + i, &buf_recv2->buf[i * 7]);
      }
      atom.nlocal += nrecv2;

      // #pragma omp barrier

    } // end of loop over dims
  } // end of loop over box layers
}

void Comm::exchange_blocking_neighbourtasks(Atom* atom, int box_id)
{
  // Skip exchange entirely if running with only one process and one box
  // Not just an optimisation; algorithm breaks otherwise.
  if (threads->mpi_num_threads == 1 && atom->boxes_per_process == 1) {
    return;
  }

  if(do_safeexchange)
    return exchange_all(*atom);

  // No outer loop over box layers this implementation. All 3 layers captured in each task

  // Determine box IDs to send to/receive from on all layers
  int send_target_box_id[3], recv_target_box_id[3];
  layer_to_targets(atom, 0, &send_target_box_id[0], &recv_target_box_id[0]);
  layer_to_targets(atom, 1, &send_target_box_id[1], &recv_target_box_id[1]);
  layer_to_targets(atom, 2, &send_target_box_id[2], &recv_target_box_id[2]);

  // swapnum starts at 100 to ensure tag is unique across communicate, borders and exchange.
  // Messages from incorrect functions may match otherwise
  int swapnum = 100; // to ensure tag is unique across communicate, borders and exchange
  // Loop over 5 dimensions (sides + corners) rather than original 3 (sides)
  for(int idim = 0; idim < 5; idim++) {

    // Skip y-dimension. This is internal swap between boxes on this process and handled outside this function
    if (idim == 1) { continue; }

    // Perform send with first neighbour
    #pragma oss task label (exchange_send_1) \
                     in(exchangePackSentinels[box_id]) \
                     out(exchangeSend1Sentinels[box_id][idim]) \
                     firstprivate(atom, box_id, idim, swapnum, send_target_box_id)
    {

      for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {
        // Determine which buffers we are using in first communication
        // Map from idim to UP, TOP_RIGHT, etc. neighbour ID
        int sendDirection = dim_to_neigh(idim, 0);
        AtomBuffer *buf_send = &boxBufs[box_id].bufs_send[EXCHANGE_FUNCTION][box_layer_index][sendDirection];
        MPI_Comm *send_comm = &boxBufs[send_target_box_id[box_layer_index]].comm[EXCHANGE_FUNCTION];
        int sendtag = swapnum;
        ++swapnum;
        int nsend = buf_send->natoms * 7; // *7 to represent how many MMD_floats are due to be sent
        buf_send->natoms = 0; // Reset atom counter for next use of this buffer

        MPI_Send(&nsend, 1, MPI_INT, procneigh[idim][0], sendtag, *send_comm);

        // DSM TODO: Could we not make this a static check? #define MMD_MPI_FLOAT as MPI_FLOAT or MPI_DOUBLE?
        if (nsend != 0) {
          if (sizeof(MMD_float) == 4) {
            MPI_Send(buf_send->buf, nsend, MPI_FLOAT, procneigh[idim][0], sendtag, *send_comm);
          } else {
            MPI_Send(buf_send->buf, nsend, MPI_DOUBLE, procneigh[idim][0], sendtag, *send_comm);
          }
        }
      }

    }

    // Perform recv with first neighbour
    // Recvs can be started even before packing tasks are finished
    #pragma oss task label (exchange_recv_1) \
                     in(initialIntegrateSentinels[box_id]) \
                     out(exchangeRecv1Sentinels[box_id][idim]) \
                     firstprivate(atom, box_id, idim, swapnum)
    {

      for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {
        // Determine which buffers we are using in first communication
        // Map from idim to UP, TOP_RIGHT, etc. neighbour ID
        int recvDirection = dim_to_neigh(idim, 1);
        AtomBuffer *buf_recv1 = &boxBufs[box_id].bufs_recv[EXCHANGE_FUNCTION][box_layer_index][recvDirection];
        MPI_Comm *recv_comm = &boxBufs[atom->box_id].comm[EXCHANGE_FUNCTION];
        int recvtag = swapnum;
        ++swapnum;

        boxBufs[box_id].nrecvs1[box_layer_index][idim] = 0;
        MPI_Recv(&boxBufs[box_id].nrecvs1[box_layer_index][idim], 1, MPI_INT, procneigh[idim][1], recvtag, *recv_comm, MPI_STATUS_IGNORE);

        // Increase size of receive buffer if receiving more than it can hold
        if (boxBufs[box_id].nrecvs1[box_layer_index][idim] > buf_recv1->maxsize) {
          buf_recv1->growrecv(boxBufs[box_id].nrecvs1[box_layer_index][idim]);
        }

        // Perform exchange with first neighbour
        // DSM TODO: Could we not make this a static check? #define MMD_MPI_FLOAT as MPI_FLOAT or MPI_DOUBLE?
        if (boxBufs[box_id].nrecvs1[box_layer_index][idim] != 0) {
          if (sizeof(MMD_float) == 4) {
            MPI_Recv(buf_recv1->buf, boxBufs[box_id].nrecvs1[box_layer_index][idim], MPI_FLOAT, procneigh[idim][1], recvtag, *recv_comm, MPI_STATUS_IGNORE);
          } else {
            MPI_Recv(buf_recv1->buf, boxBufs[box_id].nrecvs1[box_layer_index][idim], MPI_DOUBLE, procneigh[idim][1], recvtag, *recv_comm, MPI_STATUS_IGNORE);
          }
        }

        // Save recv buffer pointer for use in unpack task
        boxBufs[box_id].buf_recvs1[box_layer_index][idim] = buf_recv1;
      }

    }

    // 3 send/recv pair completes
    swapnum+=3;

    // Perform second set of exchanges
    #pragma oss task label (exchange_send_2) \
                     in(exchangePackSentinels[box_id]) \
                     out(exchangeSend2Sentinels[box_id][idim]) \
                     firstprivate(atom, box_id, idim, swapnum, send_target_box_id)
    {

      for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {
        // Determine which buffers we are using in second communication
        int sendDirection = dim_to_neigh(idim, 1);
        AtomBuffer *buf_send = &boxBufs[box_id].bufs_send[EXCHANGE_FUNCTION][box_layer_index][sendDirection];
        MPI_Comm *send_comm = &boxBufs[send_target_box_id[box_layer_index]].comm[EXCHANGE_FUNCTION];
        int sendtag = swapnum;
        ++swapnum;
        // TODO: Do we need to reverse send/recv target IDs and comms here?
        int nsend = buf_send->natoms * 7; // *7 to represent how many MMD_floats are due to be sent
        buf_send->natoms = 0; // Reset atom counter for next use of this buffer

        MPI_Send(&nsend, 1, MPI_INT, procneigh[idim][1], sendtag, *send_comm);

        // DSM: The second exchange,
        if (nsend != 0) {
          if (sizeof(MMD_float) == 4) {
            MPI_Send(buf_send->buf, nsend, MPI_FLOAT, procneigh[idim][1], sendtag, *send_comm);
          } else {
            MPI_Send(buf_send->buf, nsend, MPI_DOUBLE, procneigh[idim][1], sendtag, *send_comm);
          }
        }
      }

    }

    #pragma oss task label (exchange_recv_2) \
                     in(initialIntegrateSentinels[box_id]) \
                     out(exchangeRecv2Sentinels[box_id][idim]) \
                     firstprivate(atom, box_id, idim, swapnum)
    {

      for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {
        // Determine which buffers we are using in second communication
        int recvDirection = dim_to_neigh(idim, 0);
        AtomBuffer *buf_recv2 = &boxBufs[box_id].bufs_recv[EXCHANGE_FUNCTION][box_layer_index][recvDirection];
        MPI_Comm *recv_comm = &boxBufs[atom->box_id].comm[EXCHANGE_FUNCTION];
        int recvtag = swapnum;
        ++swapnum;

        boxBufs[box_id].nrecvs2[box_layer_index][idim] = 0;
        MPI_Recv(&boxBufs[box_id].nrecvs2[box_layer_index][idim], 1, MPI_INT, procneigh[idim][0], recvtag, *recv_comm, MPI_STATUS_IGNORE);

        // Increase size of receive buffer if receiving more than it can hold
        if (boxBufs[box_id].nrecvs2[box_layer_index][idim] > buf_recv2->maxsize) {
          buf_recv2->growrecv(boxBufs[box_id].nrecvs2[box_layer_index][idim]);
        }

        if (boxBufs[box_id].nrecvs2[box_layer_index][idim] != 0) {
          if (sizeof(MMD_float) == 4) {
            MPI_Recv(buf_recv2->buf, boxBufs[box_id].nrecvs2[box_layer_index][idim], MPI_FLOAT, procneigh[idim][0], recvtag, *recv_comm, MPI_STATUS_IGNORE);
          } else {
            MPI_Recv(buf_recv2->buf, boxBufs[box_id].nrecvs2[box_layer_index][idim], MPI_DOUBLE, procneigh[idim][0], recvtag, *recv_comm, MPI_STATUS_IGNORE);
          }
        }

        boxBufs[box_id].buf_recvs2[box_layer_index][idim] = buf_recv2;
      }

    }

    // Second set of 3 send/recv pairs completed
    swapnum+=3;

    // Unpack both neighbours' buffers. Commutative so multiple unpack tasks for same box do not
    // overwrite the x array concurrently
    // Is dependency on pack needed? Yes, what if unpack tries to grow array while a pack task is running
    #pragma oss task label (exchange_unpack) commutative(atom->nlocal) \
                     in(exchangePackSentinels[box_id]) \
                     in(exchangeRecv1Sentinels[box_id][idim]) \
                     in(exchangeRecv2Sentinels[box_id][idim]) \
                     out(exchangeSentinels[box_id]) \
                     firstprivate(atom, idim)
    {
      for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {
        AtomBuffer *buf_recv1 = boxBufs[box_id].buf_recvs1[box_layer_index][idim];
        int nrecv1 = boxBufs[box_id].nrecvs1[box_layer_index][idim];
        AtomBuffer *buf_recv2 = boxBufs[box_id].buf_recvs2[box_layer_index][idim];
        int nrecv2 = boxBufs[box_id].nrecvs2[box_layer_index][idim];

        // Again, /7 as 7 MMD_floats per atom received (position, velocity and type)
        nrecv1 /= 7;
        nrecv2 /= 7;

        // Add new atoms from each buffer to end of list. Update local count of atoms after each unpack loop.
        for (int i = 0; i < nrecv1; ++i) {
          atom->unpack_exchange(atom->nlocal + i, &buf_recv1->buf[i * 7]);
        }
        atom->nlocal += nrecv1;

        for (int i = 0; i < nrecv2; ++i) {
          atom->unpack_exchange(atom->nlocal + i, &buf_recv2->buf[i * 7]);
        }
        atom->nlocal += nrecv2;
      }
    }
  } // end of loop over dims
}

void Comm::exchange_nonblocking_neighbourtasks_tampi_iwaitall(Atom* atom, int box_id)
{
  // Skip exchange entirely if running with only one process and one box
  // Not just an optimisation; algorithm breaks otherwise.
  if (threads->mpi_num_threads == 1 && atom->boxes_per_process == 1) {
    return;
  }

  if(do_safeexchange)
    return exchange_all(*atom);

  // No outer loop over box layers this implementation. All 3 layers captured in each task

  // Determine box IDs to send to/receive from on all layers
  int send_target_box_id[3], recv_target_box_id[3];
  layer_to_targets(atom, 0, &send_target_box_id[0], &recv_target_box_id[0]);
  layer_to_targets(atom, 1, &send_target_box_id[1], &recv_target_box_id[1]);
  layer_to_targets(atom, 2, &send_target_box_id[2], &recv_target_box_id[2]);

  // swapnum starts at 100 to ensure tag is unique across communicate, borders and exchange.
  // Messages from incorrect functions may match otherwise
  int swapnum = 100; // to ensure tag is unique across communicate, borders and exchange
  // Loop over 5 dimensions (sides + corners) rather than original 3 (sides)
  for(int idim = 0; idim < 5; idim++) {

    // Skip y-dimension. This is internal swap between boxes on this process and handled outside this function
    if (idim == 1) { continue; }

    // Perform send with first neighbour
    #pragma oss task label (exchange_send_1) \
                     in(exchangePackSentinels[box_id]) \
                     out(exchangeSend1Sentinels[box_id][idim]) \
                     firstprivate(atom, box_id, idim, swapnum, send_target_box_id)
    {
      MPI_Request requests[] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL};
      int req_counter = 0;
      for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {
        // Determine which buffers we are using in first communication
        // Map from idim to UP, TOP_RIGHT, etc. neighbour ID
        int sendDirection = dim_to_neigh(idim, 0);
        AtomBuffer *buf_send = &boxBufs[box_id].bufs_send[EXCHANGE_FUNCTION][box_layer_index][sendDirection];
        MPI_Comm *send_comm = &boxBufs[send_target_box_id[box_layer_index]].comm[EXCHANGE_FUNCTION];
        int sendtag = swapnum;
        ++swapnum;
        int nsend = buf_send->natoms * 7; // *7 to represent how many MMD_floats are due to be sent
        buf_send->natoms = 0; // Reset atom counter for next use of this buffer

        MPI_Send(&nsend, 1, MPI_INT, procneigh[idim][0], sendtag, *send_comm);

        // DSM TODO: Could we not make this a static check? #define MMD_MPI_FLOAT as MPI_FLOAT or MPI_DOUBLE?
        if (nsend != 0) {
          if (sizeof(MMD_float) == 4) {
            MPI_Isend(buf_send->buf, nsend, MPI_FLOAT, procneigh[idim][0], sendtag, *send_comm, &requests[req_counter]);
          } else {
            MPI_Isend(buf_send->buf, nsend, MPI_DOUBLE, procneigh[idim][0], sendtag, *send_comm, &requests[req_counter]);
          }
          ++req_counter;
        }
      }

      if (req_counter > 0) {
        TAMPI_Iwaitall(req_counter, requests, MPI_STATUSES_IGNORE);
      }
    }

    // Perform recv with first neighbour
    // Recvs can be started even before packing tasks are finished
    #pragma oss task label (exchange_recv_1) \
                     in(initialIntegrateSentinels[box_id]) \
                     out(exchangeRecv1Sentinels[box_id][idim]) \
                     firstprivate(atom, box_id, idim, swapnum)
    {
      MPI_Request requests[] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL};
      int req_counter = 0;
      for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {
        // Determine which buffers we are using in first communication
        // Map from idim to UP, TOP_RIGHT, etc. neighbour ID
        int recvDirection = dim_to_neigh(idim, 1);
        AtomBuffer *buf_recv1 = &boxBufs[box_id].bufs_recv[EXCHANGE_FUNCTION][box_layer_index][recvDirection];
        MPI_Comm *recv_comm = &boxBufs[atom->box_id].comm[EXCHANGE_FUNCTION];
        int recvtag = swapnum;
        ++swapnum;

        boxBufs[box_id].nrecvs1[box_layer_index][idim] = 0;
        MPI_Recv(&boxBufs[box_id].nrecvs1[box_layer_index][idim], 1, MPI_INT, procneigh[idim][1], recvtag, *recv_comm, MPI_STATUS_IGNORE);

        // Increase size of receive buffer if receiving more than it can hold
        if (boxBufs[box_id].nrecvs1[box_layer_index][idim] > buf_recv1->maxsize) {
          buf_recv1->growrecv(boxBufs[box_id].nrecvs1[box_layer_index][idim]);
        }

        // Perform exchange with first neighbour
        // DSM TODO: Could we not make this a static check? #define MMD_MPI_FLOAT as MPI_FLOAT or MPI_DOUBLE?
        if (boxBufs[box_id].nrecvs1[box_layer_index][idim] != 0) {
          if (sizeof(MMD_float) == 4) {
            MPI_Irecv(buf_recv1->buf, boxBufs[box_id].nrecvs1[box_layer_index][idim], MPI_FLOAT, procneigh[idim][1], recvtag, *recv_comm, &requests[req_counter]);
          } else {
            MPI_Irecv(buf_recv1->buf, boxBufs[box_id].nrecvs1[box_layer_index][idim], MPI_DOUBLE, procneigh[idim][1], recvtag, *recv_comm, &requests[req_counter]);
          }
          ++req_counter;
        }

        // Save recv buffer pointer for use in unpack task
        boxBufs[box_id].buf_recvs1[box_layer_index][idim] = buf_recv1;
      }

      if (req_counter > 0) {
        TAMPI_Iwaitall(req_counter, requests, MPI_STATUSES_IGNORE);
      }
    }

    // 3 send/recv pair completes
    swapnum+=3;

    // Perform second set of exchanges
    #pragma oss task label (exchange_send_2) \
                     in(exchangePackSentinels[box_id]) \
                     out(exchangeSend2Sentinels[box_id][idim]) \
                     firstprivate(atom, box_id, idim, swapnum, send_target_box_id)
    {
      MPI_Request requests[] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL};
      int req_counter = 0;
      for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {
        // Determine which buffers we are using in second communication
        int sendDirection = dim_to_neigh(idim, 1);
        AtomBuffer *buf_send = &boxBufs[box_id].bufs_send[EXCHANGE_FUNCTION][box_layer_index][sendDirection];
        MPI_Comm *send_comm = &boxBufs[send_target_box_id[box_layer_index]].comm[EXCHANGE_FUNCTION];
        int sendtag = swapnum;
        ++swapnum;
        // TODO: Do we need to reverse send/recv target IDs and comms here?
        int nsend = buf_send->natoms * 7; // *7 to represent how many MMD_floats are due to be sent
        buf_send->natoms = 0; // Reset atom counter for next use of this buffer

        MPI_Send(&nsend, 1, MPI_INT, procneigh[idim][1], sendtag, *send_comm);

        // DSM: The second exchange,
        if (nsend != 0) {
          if (sizeof(MMD_float) == 4) {
            MPI_Isend(buf_send->buf, nsend, MPI_FLOAT, procneigh[idim][1], sendtag, *send_comm, &requests[req_counter]);
          } else {
            MPI_Isend(buf_send->buf, nsend, MPI_DOUBLE, procneigh[idim][1], sendtag, *send_comm, &requests[req_counter]);
          }
          ++req_counter;
        }
      }

      if (req_counter > 0) {
        TAMPI_Iwaitall(req_counter, requests, MPI_STATUSES_IGNORE);
      }
    }

    #pragma oss task label (exchange_recv_2) \
                     in(initialIntegrateSentinels[box_id]) \
                     out(exchangeRecv2Sentinels[box_id][idim]) \
                     firstprivate(atom, box_id, idim, swapnum)
    {
      MPI_Request requests[] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL};
      int req_counter = 0;
      for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {
        // Determine which buffers we are using in second communication
        int recvDirection = dim_to_neigh(idim, 0);
        AtomBuffer *buf_recv2 = &boxBufs[box_id].bufs_recv[EXCHANGE_FUNCTION][box_layer_index][recvDirection];
        MPI_Comm *recv_comm = &boxBufs[atom->box_id].comm[EXCHANGE_FUNCTION];
        int recvtag = swapnum;
        ++swapnum;

        boxBufs[box_id].nrecvs2[box_layer_index][idim] = 0;
        MPI_Recv(&boxBufs[box_id].nrecvs2[box_layer_index][idim], 1, MPI_INT, procneigh[idim][0], recvtag, *recv_comm, MPI_STATUS_IGNORE);

        // Increase size of receive buffer if receiving more than it can hold
        if (boxBufs[box_id].nrecvs2[box_layer_index][idim] > buf_recv2->maxsize) {
          buf_recv2->growrecv(boxBufs[box_id].nrecvs2[box_layer_index][idim]);
        }

        if (boxBufs[box_id].nrecvs2[box_layer_index][idim] != 0) {
          if (sizeof(MMD_float) == 4) {
            MPI_Irecv(buf_recv2->buf, boxBufs[box_id].nrecvs2[box_layer_index][idim], MPI_FLOAT, procneigh[idim][0], recvtag, *recv_comm, &requests[req_counter]);
          } else {
            MPI_Irecv(buf_recv2->buf, boxBufs[box_id].nrecvs2[box_layer_index][idim], MPI_DOUBLE, procneigh[idim][0], recvtag, *recv_comm, &requests[req_counter]);
          }
          ++req_counter;
        }

        boxBufs[box_id].buf_recvs2[box_layer_index][idim] = buf_recv2;
      }

      if (req_counter > 0) {
        TAMPI_Iwaitall(req_counter, requests, MPI_STATUSES_IGNORE);
      }
    }

    // Second set of 3 send/recv pairs completed
    swapnum+=3;

    // Unpack both neighbours' buffers. Commutative so multiple unpack tasks for same box do not
    // overwrite the x array concurrently
    // Is dependency on pack needed? Yes, what if unpack tries to grow array while a pack task is running
    #pragma oss task label (exchange_unpack) commutative(atom->nlocal) \
                     in(exchangePackSentinels[box_id]) \
                     in(exchangeRecv1Sentinels[box_id][idim]) \
                     in(exchangeRecv2Sentinels[box_id][idim]) \
                     out(exchangeSentinels[box_id]) \
                     firstprivate(atom, idim)
    {
      for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {
        AtomBuffer *buf_recv1 = boxBufs[box_id].buf_recvs1[box_layer_index][idim];
        int nrecv1 = boxBufs[box_id].nrecvs1[box_layer_index][idim];
        AtomBuffer *buf_recv2 = boxBufs[box_id].buf_recvs2[box_layer_index][idim];
        int nrecv2 = boxBufs[box_id].nrecvs2[box_layer_index][idim];

        // Again, /7 as 7 MMD_floats per atom received (position, velocity and type)
        nrecv1 /= 7;
        nrecv2 /= 7;

        // Add new atoms from each buffer to end of list. Update local count of atoms after each unpack loop.
        for (int i = 0; i < nrecv1; ++i) {
          atom->unpack_exchange(atom->nlocal + i, &buf_recv1->buf[i * 7]);
        }
        atom->nlocal += nrecv1;

        for (int i = 0; i < nrecv2; ++i) {
          atom->unpack_exchange(atom->nlocal + i, &buf_recv2->buf[i * 7]);
        }
        atom->nlocal += nrecv2;
      }
    }
  } // end of loop over dims
}

// DSM New functions that perform internal exchanges between boxes on this process
void Comm::exchange_internal(Atom &atom, int box_id)
{
  // Skip exchange entirely if running with only one process and one box
  // Not just an optimisation; algorithm breaks otherwise.
  if (threads->mpi_num_threads == 1 && atom.boxes_per_process == 1) {
    return;
  }

  // References for this box
  AtomBuffer& internal_buf_recv_down = boxBufs[box_id].internal_buf_recv_down[EXCHANGE_FUNCTION];
  AtomBuffer& internal_buf_recv_up = boxBufs[box_id].internal_buf_recv_up[EXCHANGE_FUNCTION];

  // Data has been pushed into my internal buffers. Unpack into end of atom list
  for (int i = 0; i < internal_buf_recv_down.natoms; ++i) {
    atom.unpack_exchange(atom.nlocal + i, &internal_buf_recv_down.buf[i * 7]);
  }
  atom.nlocal += internal_buf_recv_down.natoms;

  for (int i = 0; i < internal_buf_recv_up.natoms; ++i) {
    atom.unpack_exchange(atom.nlocal + i, &internal_buf_recv_up.buf[i * 7]);
  }
  atom.nlocal += internal_buf_recv_up.natoms;
}

void Comm::exchange_internal_send(AtomBuffer* buf_send, AtomBuffer* buf_recv) {
  // Skip exchange entirely if running with only one process and one box
  // Not just an optimisation; algorithm breaks otherwise.
  if (threads->mpi_num_threads == 1 && boxes_per_process == 1) {
    return;
  }

  // Perform memory copy into neighbouring box on this proc's recv buffer
  // Send buffers have already been packed as part of main exchange() function
  int nsend = buf_send->natoms*7; // *7 as 7 MMD_floats per atom
  // Increase size of receive buffers if receiving more than can be held
  if(nsend > buf_recv->maxsize) {
    buf_recv->growrecv(nsend);
  }
  // Push to neighbour.
  memcpy(buf_recv->buf, buf_send->buf, nsend * sizeof(MMD_float));
  // Set receiving buffer sizes
  buf_recv->natoms = buf_send->natoms;
  // Reset atom counter for next use of this buffer
  buf_send->natoms = 0;
}

void Comm::exchange_internal_recv(Atom* atom, AtomBuffer* recv_buf) {
  // Skip exchange entirely if running with only one process and one box
  // Not just an optimisation; algorithm breaks otherwise.
  if (threads->mpi_num_threads == 1 && boxes_per_process == 1) {
    return;
  }

  // Data has been pushed into my internal buffers. Unpack into end of atom list
  for (int i = 0; i < recv_buf->natoms; ++i) {
    atom->unpack_exchange(atom->nlocal + i, &recv_buf->buf[i * 7]);
  }
  atom->nlocal += recv_buf->natoms;
}

void Comm::borders_internal(Atom &atom, int box_id)
{
  // References for this box
  AtomBuffer& internal_buf_recv_down = boxBufs[box_id].internal_buf_recv_down[BORDERS_FUNCTION];
  AtomBuffer& internal_buf_recv_up = boxBufs[box_id].internal_buf_recv_up[BORDERS_FUNCTION];
  int* (&recvnum)[3] = boxBufs[box_id].recvnum;
  int* (&comm_recv_size)[3] = boxBufs[box_id].comm_recv_size;
  int* (&firstrecv)[3] = boxBufs[box_id].firstrecv;

  // Number of local and ghost atoms in this box
  int natoms = atom.nlocal + atom.nghost;

  // Data has been pushed into my internal buffers. Unpack into end of atom list
  for (int i = 0; i < internal_buf_recv_down.natoms; ++i) {
    atom.unpack_border(natoms + i, &internal_buf_recv_down.buf[i * atom.border_size]);
  }
  // Fix the missing pointers/values skipped in borders(). [2] is always y-ve swap and [3] always y+ve
  // box layer is always [0], other layers skip the y dimension
  recvnum[0][2] = internal_buf_recv_down.natoms;
  comm_recv_size[0][2] = internal_buf_recv_down.natoms * atom.comm_size;
  firstrecv[0][2] = natoms;
  atom.nghost += internal_buf_recv_down.natoms;

  natoms = atom.nlocal + atom.nghost;
  for (int i = 0; i < internal_buf_recv_up.natoms; ++i) {
    atom.unpack_border(natoms + i, &internal_buf_recv_up.buf[i * atom.border_size]);
  }
  recvnum[0][3] = internal_buf_recv_up.natoms;
  comm_recv_size[0][3] = internal_buf_recv_up.natoms * atom.comm_size;
  firstrecv[0][3] = natoms;
  atom.nghost += internal_buf_recv_up.natoms;
}

void Comm::communicate_internal(Atom &atom, int box_id)
{
  // References for this box
  AtomBuffer& internal_buf_recv_down = boxBufs[box_id].internal_buf_recv_down[COMMUNICATE_FUNCTION];
  AtomBuffer& internal_buf_recv_up = boxBufs[box_id].internal_buf_recv_up[COMMUNICATE_FUNCTION];
  int* (&recvnum)[3] = boxBufs[box_id].recvnum;
  // Values of firstrecv are set in borders_internal(): [2] is always y-ve swap and [3] always y+ve.
  int* (&firstrecv)[3] = boxBufs[box_id].firstrecv;

  // Data has been pushed into my internal buffers. Unpack into end of atom list
  // Unlike other unpack routines, unpack_comm has an internal loop over atoms
  atom.unpack_comm(recvnum[0][2], firstrecv[0][2], internal_buf_recv_down.buf);
  atom.unpack_comm(recvnum[0][3], firstrecv[0][3], internal_buf_recv_up.buf);
}

void Comm::communicate_internal_send(Atom* atom, AtomBuffer* buf_send, AtomBuffer* buf_recv, int* sendlist, int nsend, int sendnum) {
  // Increase size of receive buffer, if needed, then pack data directly into recv buffer using periodic boundary
  // conditions from send buffer
  if (nsend > buf_recv->maxsize) {
    buf_recv->growrecv(nsend);
  }
  atom->pack_comm(sendnum, sendlist, buf_recv->buf, buf_send->pbc_any, buf_send->pbc_x,
                  buf_send->pbc_y, buf_send->pbc_z);
}

void Comm::exchange_all(Atom &atom)
{
  // DSM Removed for multibox (doesn't build)
//  int i, m, n, idim, nsend, nrecv, nrecv1, nrecv2, nlocal;
//  MMD_float lo, hi, value;
//  MMD_float* x;
//
//  MPI_Request request;
//  MPI_Status status;
//
//  /* enforce PBC */
//
//  atom.pbc();
//
//  /* loop over dimensions */
//  int iswap = 0;
//
//  for(idim = 0; idim < 3; idim++) {
//
//    /* only exchange if more than one proc in this dimension */
//
//    if(procgrid[idim] == 1) {
//      iswap += 2 * need[idim];
//      continue;
//    }
//
//    /* fill buffer with atoms leaving my box
//    *        when atom is deleted, fill it in with last atom */
//
//    i = nsend = 0;
//
//    if(idim == 0) {
//      lo = atom.box.xlo;
//      hi = atom.box.xhi;
//    } else if(idim == 1) {
//      lo = atom.box.ylo;
//      hi = atom.box.yhi;
//    } else {
//      lo = atom.box.zlo;
//      hi = atom.box.zhi;
//    }
//
//    x = atom.x;
//
//    nlocal = atom.nlocal;
//
//    while(i < nlocal) {
//      if(x[i * PAD + idim] < lo || x[i * PAD + idim] >= hi) {
//        if(nsend > maxsend) growsend(nsend);
//
//        nsend += atom.pack_exchange(i, &buf_send[nsend]);
//        atom.copy(nlocal - 1, i);
//        nlocal--;
//      } else i++;
//    }
//
//    atom.nlocal = nlocal;
//
//    /* send/recv atoms in both directions
//    *        only if neighboring procs are different */
//    for(int ineed = 0; ineed < 2 * need[idim]; ineed += 1) {
//      if(ineed < procgrid[idim] - 1) {
//        MPI_Send(&nsend, 1, MPI_INT, sendproc_exc[iswap], 0, MPI_COMM_WORLD);
//        MPI_Recv(&nrecv, 1, MPI_INT, recvproc_exc[iswap], 0, MPI_COMM_WORLD, &status);
//
//        if(nrecv > maxrecv) growrecv(nrecv);
//
//        if(sizeof(MMD_float) == 4) {
//          MPI_Irecv(buf_recv, nrecv, MPI_FLOAT, recvproc_exc[iswap], 0,
//                    MPI_COMM_WORLD, &request);
//          MPI_Send(buf_send, nsend, MPI_FLOAT, sendproc_exc[iswap], 0, MPI_COMM_WORLD);
//        } else {
//          MPI_Irecv(buf_recv, nrecv, MPI_DOUBLE, recvproc_exc[iswap], 0,
//                    MPI_COMM_WORLD, &request);
//          MPI_Send(buf_send, nsend, MPI_DOUBLE, sendproc_exc[iswap], 0, MPI_COMM_WORLD);
//        }
//
//        MPI_Wait(&request, &status);
//
//        /* check incoming atoms to see if they are in my box
//        *        if they are, add to my list */
//
//        n = atom.nlocal;
//        m = 0;
//
//        while(m < nrecv) {
//          value = buf_recv[m + idim];
//
//          if(value >= lo && value < hi)
//            m += atom.unpack_exchange(n++, &buf_recv[m]);
//          else m += atom.skip_exchange(&buf_recv[m]);
//        }
//
//        atom.nlocal = n;
//      }
//
//      iswap += 1;
//
//    }
//  }
}

/* borders:
   make lists of nearby atoms to send to neighboring procs at every timestep
   one list is created for every swap that will be made
   as list is made, actually do swaps
   this does equivalent of a communicate (so don't need to explicitly
     call communicate routine on reneighboring timestep)
   this routine is called before every reneighboring
*/

// DSM borders() is always called after exchange()
// DSM: exchange() is always called before borders()
void Comm::borders(Atom* atoms[]) {
  // Select between blocking and non-blocking communication modes
  if (nonblocking_enabled) {
    borders_nonblocking(atoms);
  } else {
    for(int box_index = 0; box_index < atoms[0]->boxes_per_process; ++box_index) {
      // Dependencies on both sort and exchange completing as sort operation may be skipped
      #pragma oss task label(borders_pack) \
                       in(sortSentinels[box_index]) \
                       in(exchangePBCSentinels[box_index]) \
                       out(bordersPackSentinels[box_index]) \
                       firstprivate(box_index)
      borders_pack(atoms[box_index]);
    }

    for(int box_index = 0; box_index < atoms[0]->boxes_per_process; ++box_index) {
      //#pragma oss task label(borders_blocking) in(bordersPackSentinels[box_index]) out(bordersSentinels[box_index]) firstprivate(box_index)
      //borders_blocking(*atoms[box_index], box_index);

      // Tasks created within this function
      //borders_blocking_neighbourtasks(atoms[box_index], box_index);
      blocking_nonblocking_neighbourtasks_tampi_iwaitall(atoms[box_index], box_index);
    }

    for(int box_index = 0; box_index < atoms[0]->boxes_per_process; ++box_index) {
      // Internal swaps between boxes
      #pragma oss task label(borders_internal_send) in(bordersPackSentinels[box_index]) out(bordersInternalSendSentinels[box_index]) firstprivate(box_index)
      borders_internal_send(atoms[box_index]);
    }
  } // End of blocking/nonblocking branch

  for(int box_index = 0; box_index < atoms[0]->boxes_per_process; ++box_index) {
    // After swaps over all boxes completed with their off-process neighbours,
    // perform internal exchanges with all boxes on this process.
    // Commutative on other unpack tasks on this box to avoid multiple concurrent attempts to write into x array
    // TODO: Why does this need to depend on send? Crashes if only depends on pack
    #pragma oss task label(borders_internal_recv) commutative(atoms[box_index]->nghost) \
                       in(bordersInternalSendSentinels[atoms[box_index]->boxneigh_positive]) \
                       in(bordersInternalSendSentinels[atoms[box_index]->boxneigh_negative]) \
                       in(bordersSentinels[box_index]) \
                       in(bordersSendSentinels[box_index]) \
                       out(bordersInternalSentinels[box_index]) firstprivate(box_index)
    borders_internal(*atoms[box_index], box_index);
  }
}

// Determine which slabs, if any, all atoms are currently within. Copy atoms within a slab to the appropriate send
// buffers for those directions.
void Comm::borders_pack(Atom* atom) {

  int box_id = atom->box_id;
  AtomBuffer (&bufs_send)[3][8] = boxBufs[box_id].bufs_send[BORDERS_FUNCTION];
  AtomBuffer (&bufs_recv)[3][8] = boxBufs[box_id].bufs_recv[BORDERS_FUNCTION];
  int* (&sendlists)[3][8] = boxBufs[box_id].sendlists;
  int (&maxsendlists)[3][8] = boxBufs[box_id].maxsendlists;
  AtomBuffer& internal_buf_send_up = boxBufs[box_id].internal_buf_send_up[BORDERS_FUNCTION];
  AtomBuffer& internal_buf_send_down = boxBufs[box_id].internal_buf_send_down[BORDERS_FUNCTION];
  AtomBuffer& internal_buf_recv_up = boxBufs[box_id].internal_buf_recv_up[BORDERS_FUNCTION];
  AtomBuffer& internal_buf_recv_down = boxBufs[box_id].internal_buf_recv_down[BORDERS_FUNCTION];
  int*& internal_sendlist_up = boxBufs[box_id].internal_sendlist_up;
  int& internal_maxsendlist_up = boxBufs[box_id].internal_maxsendlist_up;
  int*& internal_sendlist_down = boxBufs[box_id].internal_sendlist_down;
  int& internal_maxsendlist_down = boxBufs[box_id].internal_maxsendlist_down;

  /* erase all ghost atoms */

  atom->nghost = 0;

  // Loop over atoms here instead of in iswaps loop. Only one loop over the atoms array rather than one per
  // dimension as in original implementation.
  // Dissimilar from exchange() calculations: per slab rather than per box and each atom can be packed into multiple
  // buffers. An atom in the top-left corner slab will also need to be sent directly left, for example.
  // 26 exchanges = 8 process neighbours * 3 box layers + 2 internal exchanges between boxes on same process
  int nfirst = 0;
  int nlast = atom->nlocal+atom->nghost; // DSM: nghost will always be 0 here, we sort all local atoms in a single pass now
  MMD_float xcoord, ycoord, zcoord;
  int type;
  for(int i = nfirst; i < nlast; i++) { // DSM: From 0 to atom.nlocal+atom.nghost
    // Get atom's current coordinates and type (needed for packing)
    xcoord = atom->x[i * PAD];
    ycoord = atom->x[i * PAD + 1];
    zcoord = atom->x[i * PAD + 2];
    type = atom->type[i];

    // First, check which layer this atom is on
    // box_id (same layer), box_id-1 (below), box_id+1 (above)
    // Index into array of send/recv buffers.
    // Possible for atom to be in all 3 layers if box is small enough that bottom slab boundary overlaps with top:
    // e.g. bottom slab 48-51 units and top 50-53 units.
    // Use these flags to control which layers we need to loop over
    int atom_layers[3] = {SAME_LAYER, NO_LAYER, NO_LAYER};

    // Check if atom is within a y slab and set layers flags accordingly.
    if (ycoord >= atom->box.yneg_slab_lo && ycoord <= atom->box.yneg_slab_hi) {
      // In -ve y slab: box_id-1 / "lower" layer.
      atom_layers[1] = LOWER_LAYER;
      // In y slab, necessarily pack into internal buffer to be sent down
      internal_buf_send_down.pack_border(i, xcoord, ycoord, zcoord, type, &internal_sendlist_down, internal_maxsendlist_down);
    }
    // NOT an else-if. Atom can be in all 3 layers!
    if (ycoord >= atom->box.ypos_slab_lo && ycoord <= atom->box.ypos_slab_hi) {
      // In +ve y slab: box_id+1 / "upper" layer
      atom_layers[2] = UPPER_LAYER;
      internal_buf_send_up.pack_border(i, xcoord, ycoord, zcoord, type, &internal_sendlist_up, internal_maxsendlist_up);
    }

    // This loop will be executed 1-3 times:
    //   * Always for neighbours on the same layer of boxes
    //   * Sometimes for neighbours on the above and below layers of boxes (can be both!), dependent on result of above checks
    for (int layer_index=0; layer_index < 3; ++layer_index) {
      int atom_layer = atom_layers[layer_index];

      // Skip if no second layer needing checked
      if (atom_layer == NO_LAYER) {
        continue;
      }

      if (xcoord >= atom->box.xneg_slab_lo && xcoord <= atom->box.xneg_slab_hi) { // Is atom in x -ve slab?
        // x-ve: send "left"
        bufs_send[atom_layer][LEFT].pack_border(i, xcoord, ycoord, zcoord, type,
                                                &sendlists[atom_layer][LEFT], maxsendlists[atom_layer][LEFT]);

        // Is atom also in a z dimension slab?
        if (zcoord >= atom->box.zneg_slab_lo && zcoord <= atom->box.zneg_slab_hi) {
          // x-ve and z-ve: send to "bottom-left"
          bufs_send[atom_layer][BOTTOM_LEFT].pack_border(i, xcoord, ycoord, zcoord, type,
                                                         &sendlists[atom_layer][BOTTOM_LEFT], maxsendlists[atom_layer][BOTTOM_LEFT]);
        }

        // NOT an else-if. Possible for atom to be in -ve and +ve slabs given certain decompositions
        if (zcoord >= atom->box.zpos_slab_lo && zcoord <= atom->box.zpos_slab_hi) {
          // x-ve and z+ve: send to "top-left"
          bufs_send[atom_layer][TOP_LEFT].pack_border(i, xcoord, ycoord, zcoord, type,
                                                      &sendlists[atom_layer][TOP_LEFT], maxsendlists[atom_layer][TOP_LEFT]);
        }

      }

      // NOT an else-if. Possible for atom to be in -ve and +ve slabs given certain decompositions
      if (xcoord >= atom->box.xpos_slab_lo && xcoord <= atom->box.xpos_slab_hi) { // Is atom in x +ve slab?
        // x+ve: send "right"
        bufs_send[atom_layer][RIGHT].pack_border(i, xcoord, ycoord, zcoord, type,
                                                 &sendlists[atom_layer][RIGHT], maxsendlists[atom_layer][RIGHT]);

        // Is atom also in a z dimension slab?
        if (zcoord >= atom->box.zneg_slab_lo && zcoord <= atom->box.zneg_slab_hi) {
          // x+ve and z-ve: send to "bottom-right"
          bufs_send[atom_layer][BOTTOM_RIGHT].pack_border(i, xcoord, ycoord, zcoord, type,
                                                          &sendlists[atom_layer][BOTTOM_RIGHT], maxsendlists[atom_layer][BOTTOM_RIGHT]);
        }

        if (zcoord >= atom->box.zpos_slab_lo && zcoord <= atom->box.zpos_slab_hi) {
          // x+ve and z+ve: send to "top-right"
          bufs_send[atom_layer][TOP_RIGHT].pack_border(i, xcoord, ycoord, zcoord, type,
                                                       &sendlists[atom_layer][TOP_RIGHT], maxsendlists[atom_layer][TOP_RIGHT]);
        }

      }

      if (zcoord >= atom->box.zneg_slab_lo && zcoord <= atom->box.zneg_slab_hi) { // Is atom in z -ve slab?
        // z-ve: send "down"
        bufs_send[atom_layer][DOWN].pack_border(i, xcoord, ycoord, zcoord, type,
                                                &sendlists[atom_layer][DOWN], maxsendlists[atom_layer][DOWN]);

      }

      if (zcoord >= atom->box.zpos_slab_lo && zcoord <= atom->box.zpos_slab_hi) { // Is atom in z +ve slab?
        // z+ve: send "up"
        bufs_send[atom_layer][UP].pack_border(i, xcoord, ycoord, zcoord, type,
                                              &sendlists[atom_layer][UP], maxsendlists[atom_layer][UP]);

      }
    } // end of loop over layers

  } // end of loop over atoms
}

// Perform memory copy into neighbouring box on this proc's recv buffer.
void Comm::borders_internal_send(Atom* atom) {

  int box_id = atom->box_id;
  int* (&sendnum)[3] = boxBufs[box_id].sendnum;
  int* (&comm_send_size)[3] = boxBufs[box_id].comm_send_size;
  AtomBuffer& internal_buf_send_up = boxBufs[box_id].internal_buf_send_up[BORDERS_FUNCTION];
  AtomBuffer& internal_buf_send_down = boxBufs[box_id].internal_buf_send_down[BORDERS_FUNCTION];

  // First exchange to lower layer (usually box_id-1 but periodic)
  int box_neighbour = atom->boxneigh_negative;
  AtomBuffer* buf_send = &internal_buf_send_down;
  AtomBuffer* buf_recv = &(boxBufs[box_neighbour].internal_buf_recv_down[BORDERS_FUNCTION]);
  int nsend = buf_send->natoms * atom->border_size; // natoms * num MMD_floats per atom (e.g. 4: x,y,z,type)
  // Increase size of receive buffers if receiving more than can be held
  if(nsend > buf_recv->maxsize) {
    buf_recv->growrecv(nsend);
  }
  // Push to neighbour.
  memcpy(buf_recv->buf, buf_send->buf, nsend * sizeof(MMD_float));
  // Update relevant pointers and counters for this swap
  // Recv side done in seperate internal swaps function
  // Always box layer 0. box layers 1 and 2 are completely skipped for y dimension
  // Lower layer send is always swap [2]. Upper layer is always swap [3]
  sendnum[0][2] = buf_send->natoms;
  comm_send_size[0][2] = buf_send->natoms * atom->comm_size;
  // Set receiving buffer sizes
  buf_recv->natoms = buf_send->natoms;
  // Reset atom counter for next use of this buffer
  buf_send->natoms = 0;

  // Second exchange to upper layer (usually box_id+1)
  box_neighbour = atom->boxneigh_positive;
  buf_send = &internal_buf_send_up;
  buf_recv = &(boxBufs[box_neighbour].internal_buf_recv_up[BORDERS_FUNCTION]);
  nsend = buf_send->natoms * atom->border_size;
  if(nsend > buf_recv->maxsize) {
    buf_recv->growrecv(nsend);
  }
  memcpy(buf_recv->buf, buf_send->buf, nsend * sizeof(MMD_float));
  sendnum[0][3] = buf_send->natoms;
  comm_send_size[0][3] = buf_send->natoms * atom->comm_size;
  buf_recv->natoms = buf_send->natoms;
  buf_send->natoms = 0;
}

void Comm::borders_blocking(Atom &atom, int box_id)
{
  int i, m, n, iswap, idim, ineed, nsend, nrecv, nall, sendtag, recvtag;
  int swapnum = 0;
  AtomBuffer *buf = NULL, *buf_send = NULL, *buf_recv = NULL;
  MMD_float* x;
  MPI_Request request;
  MPI_Status status;
  MPI_Comm *send_comm, *recv_comm;

  // DSM Multibox: Reference buffers/variables specific to this box.
  int& nswap = boxBufs[box_id].nswap;
  int& maxnlocal = boxBufs[box_id].maxnlocal;
  int*& sendproc = boxBufs[box_id].sendproc;
  int*& recvproc = boxBufs[box_id].recvproc;
  int*& sendneigh = boxBufs[box_id].sendneigh;
  int*& recvneigh = boxBufs[box_id].recvneigh;
  int* (&comm_send_size)[3] = boxBufs[box_id].comm_send_size;
  int* (&comm_recv_size)[3] = boxBufs[box_id].comm_recv_size;
  int& nrecv_atoms = boxBufs[box_id].nrecv_atoms;
  int* (&sendnum)[3] = boxBufs[box_id].sendnum;
  int* (&recvnum)[3] = boxBufs[box_id].recvnum;
  int* (&firstrecv)[3] = boxBufs[box_id].firstrecv;
  int* (&reverse_send_size)[3] = boxBufs[box_id].reverse_send_size;
  int* (&reverse_recv_size)[3] = boxBufs[box_id].reverse_recv_size;
  AtomBuffer (&bufs_send)[3][8] = boxBufs[box_id].bufs_send[BORDERS_FUNCTION];
  AtomBuffer (&bufs_recv)[3][8] = boxBufs[box_id].bufs_recv[BORDERS_FUNCTION];
  int* (&sendlists)[3][8] = boxBufs[box_id].sendlists;
  int (&maxsendlists)[3][8] = boxBufs[box_id].maxsendlists;

  // Loop over box layers: box_id, box_id-1, box_id+1
  for (int box_layer_index = 0; box_layer_index < 3; box_layer_index++) {
    iswap = 0;
    // Loop over 5 dimensions x, y, z, diagonal-right, diagonal-left
    for(idim = 0; idim < 5; idim++) {

      // DSM: Do not perform MPI communication in the y / idim==1 dimension
      // Memory copy into neighbouring box on this proc's recv buffer will be done outside this process
      if (idim == 1) {
        iswap += 2; // Skip over both sides of the y dimension swap here
        continue;
      }

      for(ineed = 0; ineed < 2 * need[idim]; ineed++) { // DSM Multibox: Will always be ineed < 2 here, need guaranteed to be 1.

        /* swap atoms with other proc
        put incoming ghosts at end of my atom arrays
        if swapping with self, simply copy, no messages */

        //#pragma omp master
        {
          // DSM Two send/recv pairs per dimension
          int send_target_box_id, recv_target_box_id;
          layer_to_targets(&atom, box_layer_index, &send_target_box_id, &recv_target_box_id);

          // Determine which buffers we are using in communication
          int sendNeighbour = sendneigh[iswap];
          int recvNeighbour = recvneigh[iswap];
          buf_send = &bufs_send[box_layer_index][sendNeighbour];
          buf_recv = &bufs_recv[box_layer_index][recvNeighbour];
          send_comm = &boxBufs[send_target_box_id].comm[BORDERS_FUNCTION];
          recv_comm = &boxBufs[atom.box_id].comm[BORDERS_FUNCTION];
          sendtag = swapnum; // Not the same as iswap - swapnum counts over all box layers
          recvtag = sendtag;
          ++swapnum;
          nsend = buf_send->natoms; // Number of atoms to send (*not* MMD_Floats)
          buf_send->natoms = 0; // Reset atom counter for next use of this buffer

          if(sendproc[iswap] == me && send_target_box_id == atom.box_id) {
            // Simply copy if swapping with yourself
            nrecv = nsend;
            buf = buf_send;
          } else {
            // DSM: Changed to Irecv/Send pair from Send/Recv. No guarantee original code progresses (likely will due to buffered send)
            // DSM TODO: Double check the comm index makes sense - should the receive be box_comms[sender_id]? What happens if we have multiple sends/receives to the same box at the same time?
            // DSM TODO: Optimise this - we can skip a send/recv pair if our nsend/nrecv == 0
            MPI_Irecv(&nrecv, 1, MPI_INT, recvproc[iswap], recvtag, *recv_comm, &request);
            MPI_Send(&nsend, 1, MPI_INT, sendproc[iswap], sendtag,  *send_comm);
            MPI_Wait(&request, &status);

            // Increase receiving buffer size if needed
            if(nrecv * atom.border_size > buf_recv->maxsize) buf_recv->growrecv(nrecv * atom.border_size);

            if(sizeof(MMD_float) == 4) {
              MPI_Irecv(buf_recv->buf, nrecv * atom.border_size, MPI_FLOAT,
                        recvproc[iswap], recvtag, *recv_comm, &request);
              MPI_Send(buf_send->buf, nsend * atom.border_size, MPI_FLOAT,
                       sendproc[iswap], sendtag,  *send_comm);
            } else {
              MPI_Irecv(buf_recv->buf, nrecv * atom.border_size, MPI_DOUBLE,
                        recvproc[iswap], recvtag, *recv_comm, &request);
              MPI_Send(buf_send->buf, nsend * atom.border_size, MPI_DOUBLE,
                       sendproc[iswap], sendtag, *send_comm);
            }

            MPI_Wait(&request, &status);
            buf = buf_recv;
          }

          nrecv_atoms = nrecv;
        }
        /* unpack buffer */

        n = atom.nlocal + atom.nghost;
        nrecv = nrecv_atoms;

        //#pragma omp for
        for(int i = 0; i < nrecv; i++)
          atom.unpack_border(n + i, &buf->buf[i * 4]); // DSM TODO: Shouldn't this be atom.border_size instead of "4"?

        // #pragma omp barrier

        /* set all pointers & counters */

        //#pragma omp master
        {
          // DSM These now duplicated per layer
          sendnum[box_layer_index][iswap] = nsend;
          recvnum[box_layer_index][iswap] = nrecv;
          // DSM: * atom.comm_size (constant '3' defined in constructor) since need 3 doubles per atom (x,y,z coordinates)
          comm_send_size[box_layer_index][iswap] = nsend * atom.comm_size;
          comm_recv_size[box_layer_index][iswap] = nrecv * atom.comm_size;
          reverse_send_size[box_layer_index][iswap] = nrecv * atom.reverse_size;
          reverse_recv_size[box_layer_index][iswap] = nsend * atom.reverse_size;
          firstrecv[box_layer_index][iswap] = atom.nlocal + atom.nghost;
          atom.nghost += nrecv;
        }

        //#pragma omp barrier
        iswap++;
      } // DSM End of loop over ineed/number of swaps
    } // DSM End of loop over dimensions
  } // DSM End of loop over box layers

  /* insure buffers are large enough for reverse comm */
  // DSM Multibox: Don't need to worry about this path, reverse_communicate is only done in halfneigh calculation. We
  // are ignoring this method.

  /*
  int max1, max2;
  max1 = max2 = 0;

  for(iswap = 0; iswap < nswap; iswap++) {
    max1 = MAX(max1, reverse_send_size[iswap]);
    max2 = MAX(max2, reverse_recv_size[iswap]);
  }

  if(max1 > maxsend) growsend(max1, box_id);

  if(max2 > maxrecv) growrecv(max2, box_id);
  */
}

void Comm::borders_blocking_neighbourtasks(Atom* atom, int box_id)
{
  int iswap, idim, ineed;
  int swapnum = 200; // to ensure tag is unique accross communicate, borders and exchange

  // Determine box IDs to send to/receive from on all layers
  int send_target_box_id[3], recv_target_box_id[3];
  layer_to_targets(atom, 0, &send_target_box_id[0], &recv_target_box_id[0]);
  layer_to_targets(atom, 1, &send_target_box_id[1], &recv_target_box_id[1]);
  layer_to_targets(atom, 2, &send_target_box_id[2], &recv_target_box_id[2]);

  // No outer loop over box layers this implementation. All 3 layers captured in each task
  iswap = 0;
  // Loop over 5 dimensions x, y, z, diagonal-right, diagonal-left
  for(idim = 0; idim < 5; idim++) {

    // DSM: Do not perform MPI communication in the y / idim==1 dimension
    // Memory copy into neighbouring box on this proc's recv buffer will be done outside this process
    if (idim == 1) {
      iswap += 2; // Skip over both sides of the y dimension swap here
      continue;
    }

    for(ineed = 0; ineed < 2 * need[idim]; ineed++) { // DSM Multibox: Will always be ineed < 2 here, need guaranteed to be 1.

      /* swap atoms with other proc
      put incoming ghosts at end of my atom arrays
      if swapping with self, simply copy, no messages */

      #pragma oss task label(borders_send) \
                       in(bordersPackSentinels[box_id]) \
                       out(bordersSendSentinels[box_id]) \
                       firstprivate(atom, box_id, iswap, swapnum)
      {
        for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {
          // Determine which buffers we are using in communication
          int sendNeighbour = boxBufs[box_id].sendneigh[iswap];
          AtomBuffer* buf_send = &boxBufs[box_id].bufs_send[BORDERS_FUNCTION][box_layer_index][sendNeighbour];
          MPI_Comm* send_comm = &boxBufs[send_target_box_id[box_layer_index]].comm[BORDERS_FUNCTION];
          int sendtag = swapnum; // Not the same as iswap - swapnum counts over all box layers
          ++swapnum;
          int nsend = buf_send->natoms; // Number of atoms to send (*not* MMD_Floats)

          MPI_Send(&nsend, 1, MPI_INT, boxBufs[box_id].sendproc[iswap], sendtag, *send_comm);

          // Skip a send/recv pair if our nsend/nrecv == 0
          if (nsend > 0) {
            if (sizeof(MMD_float) == 4) {
              MPI_Send(buf_send->buf, nsend * atom->border_size, MPI_FLOAT,
                       boxBufs[box_id].sendproc[iswap], sendtag, *send_comm);
            } else {
              MPI_Send(buf_send->buf, nsend * atom->border_size, MPI_DOUBLE,
                       boxBufs[box_id].sendproc[iswap], sendtag, *send_comm);
            }
          }

          /* set all pointers & counters */

          boxBufs[box_id].sendnum[box_layer_index][iswap] = nsend;
          // DSM: * atom.comm_size (constant '3' defined in constructor) since need 3 doubles per atom (x,y,z coordinates)
          boxBufs[box_id].comm_send_size[box_layer_index][iswap] = nsend * atom->comm_size;
          // Reset atom counter for next use of this buffer
          buf_send->natoms = 0;
        }
      }

      // in-dependency on first operation in iteration (initialIntegrate).
      // Do not need to wait for packs or sends, recvs can be posted immediately.
      #pragma oss task label(borders_recv) \
                       in(initialIntegrateSentinels[box_id]) \
                       out(bordersRecvSentinels[box_id][idim][iswap]) \
                       firstprivate(atom, box_id, iswap, swapnum, idim, ineed)
      {
        for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {
          // Determine which buffers we are using in communication
          int recvNeighbour = boxBufs[box_id].recvneigh[iswap];
          AtomBuffer* buf_recv = &boxBufs[box_id].bufs_recv[BORDERS_FUNCTION][box_layer_index][recvNeighbour];
          MPI_Comm* recv_comm = &boxBufs[atom->box_id].comm[BORDERS_FUNCTION];
          int recvtag = swapnum; // Not the same as iswap - swapnum counts over all box layers
          ++swapnum;

          int nrecv = 0;
          // DSM TODO: Double check the comm index makes sense - should the receive be box_comms[sender_id]? What happens if we have multiple sends/receives to the same box at the same time?
          MPI_Recv(&nrecv, 1, MPI_INT, boxBufs[box_id].recvproc[iswap], recvtag, *recv_comm, MPI_STATUS_IGNORE);

          // Skip a send/recv pair if our nsend/nrecv == 0
          if (nrecv > 0) {
            // Increase receiving buffer size if needed
            if (nrecv * atom->border_size > buf_recv->maxsize) buf_recv->growrecv(nrecv * atom->border_size);

            if (sizeof(MMD_float) == 4) {
              MPI_Recv(buf_recv->buf, nrecv * atom->border_size, MPI_FLOAT,
                       boxBufs[box_id].recvproc[iswap], recvtag, *recv_comm, MPI_STATUS_IGNORE);
            } else {
              MPI_Recv(buf_recv->buf, nrecv * atom->border_size, MPI_DOUBLE,
                       boxBufs[box_id].recvproc[iswap], recvtag, *recv_comm, MPI_STATUS_IGNORE);
            }
          }

          // Save values for unpack task
          boxBufs[box_id].nrecvs_borders[box_layer_index][idim][ineed] = nrecv;
          boxBufs[box_id].buf_recvs_borders[box_layer_index][idim][ineed] = buf_recv;
        }
      }

      // Commutative to avoid multiple tasks unpacking into same box's x array concurrently
      // Similarly, depends on pack having completed TODO: Could this be removed? Unlikely, what if unpack_border reallocs array as pack is occurring?
      // Additionally depends on this iteration's recv task having finished, rather than all recvs for this box
      #pragma oss task label(borders_unpack) commutative(atom->nghost) \
                       in(bordersRecvSentinels[box_id][idim][iswap]) \
                       in(bordersPackSentinels[box_id]) \
                       out(bordersUnpackSentinels[box_id]) \
                       firstprivate(atom, box_id, iswap, swapnum, idim, ineed)
      {
        int n = 0, nrecv = 0;
        AtomBuffer *buf = NULL;
        for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {
          n = atom->nlocal + atom->nghost;
          nrecv = boxBufs[box_id].nrecvs_borders[box_layer_index][idim][ineed];
          buf = boxBufs[box_id].buf_recvs_borders[box_layer_index][idim][ineed];

          //#pragma omp for
          for (int i = 0; i < nrecv; i++)
            // DSM TODO: Shouldn't this be atom.border_size instead of "4"?
            atom->unpack_border(n + i, &buf->buf[i * 4]);

          /* set all pointers & counters */

          boxBufs[box_id].recvnum[box_layer_index][iswap] = nrecv;
          // DSM: * atom.comm_size (constant '3' defined in constructor) since need 3 doubles per atom (x,y,z coordinates)
          boxBufs[box_id].comm_recv_size[box_layer_index][iswap] = nrecv * atom->comm_size;
          boxBufs[box_id].firstrecv[box_layer_index][iswap] = atom->nlocal + atom->nghost;
          atom->nghost += nrecv;
        }
      }

      // 3 send-recv pairs (1 per layer) done per iswap iteration
      // i.e. per pair of send and recv tasks
      swapnum+=3;

      iswap++;
    } // DSM End of loop over ineed/number of swaps
  } // DSM End of loop over dimensions
}

void Comm::blocking_nonblocking_neighbourtasks_tampi_iwaitall(Atom* atom, int box_id)
{
  int iswap, idim, ineed;
  int swapnum = 200; // to ensure tag is unique accross communicate, borders and exchange

  // Determine box IDs to send to/receive from on all layers
  int send_target_box_id[3], recv_target_box_id[3];
  layer_to_targets(atom, 0, &send_target_box_id[0], &recv_target_box_id[0]);
  layer_to_targets(atom, 1, &send_target_box_id[1], &recv_target_box_id[1]);
  layer_to_targets(atom, 2, &send_target_box_id[2], &recv_target_box_id[2]);

  // No outer loop over box layers this implementation. All 3 layers captured in each task
  iswap = 0;
  // Loop over 5 dimensions x, y, z, diagonal-right, diagonal-left
  for(idim = 0; idim < 5; idim++) {

    // DSM: Do not perform MPI communication in the y / idim==1 dimension
    // Memory copy into neighbouring box on this proc's recv buffer will be done outside this process
    if (idim == 1) {
      iswap += 2; // Skip over both sides of the y dimension swap here
      continue;
    }

    for(ineed = 0; ineed < 2 * need[idim]; ineed++) { // DSM Multibox: Will always be ineed < 2 here, need guaranteed to be 1.

      /* swap atoms with other proc
      put incoming ghosts at end of my atom arrays
      if swapping with self, simply copy, no messages */

      #pragma oss task label(borders_send) \
                       in(bordersPackSentinels[box_id]) \
                       out(bordersSendSentinels[box_id]) \
                       firstprivate(atom, box_id, iswap, swapnum)
      {
        MPI_Request requests[] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL};
        int req_counter = 0;
        for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {
          // Determine which buffers we are using in communication
          int sendNeighbour = boxBufs[box_id].sendneigh[iswap];
          AtomBuffer* buf_send = &boxBufs[box_id].bufs_send[BORDERS_FUNCTION][box_layer_index][sendNeighbour];
          MPI_Comm* send_comm = &boxBufs[send_target_box_id[box_layer_index]].comm[BORDERS_FUNCTION];
          int sendtag = swapnum; // Not the same as iswap - swapnum counts over all box layers
          ++swapnum;
          int nsend = buf_send->natoms; // Number of atoms to send (*not* MMD_Floats)

          MPI_Send(&nsend, 1, MPI_INT, boxBufs[box_id].sendproc[iswap], sendtag, *send_comm);

          // Skip a send/recv pair if our nsend/nrecv == 0
          if (nsend > 0) {
            if (sizeof(MMD_float) == 4) {
              MPI_Isend(buf_send->buf, nsend * atom->border_size, MPI_FLOAT,
                       boxBufs[box_id].sendproc[iswap], sendtag, *send_comm, &requests[req_counter]);
            } else {
              MPI_Isend(buf_send->buf, nsend * atom->border_size, MPI_DOUBLE,
                       boxBufs[box_id].sendproc[iswap], sendtag, *send_comm, &requests[req_counter]);
            }
            ++req_counter;
          }

          /* set all pointers & counters */

          boxBufs[box_id].sendnum[box_layer_index][iswap] = nsend;
          // DSM: * atom.comm_size (constant '3' defined in constructor) since need 3 doubles per atom (x,y,z coordinates)
          boxBufs[box_id].comm_send_size[box_layer_index][iswap] = nsend * atom->comm_size;
          // Reset atom counter for next use of this buffer
          buf_send->natoms = 0;
        }

        if (req_counter > 0) {
          TAMPI_Iwaitall(req_counter, requests, MPI_STATUSES_IGNORE);
        }
      }

      // in-dependency on first operation in iteration (initialIntegrate).
      // Do not need to wait for packs or sends, recvs can be posted immediately.
      #pragma oss task label(borders_recv) \
                       in(initialIntegrateSentinels[box_id]) \
                       out(bordersRecvSentinels[box_id][idim][iswap]) \
                       firstprivate(atom, box_id, iswap, swapnum, idim, ineed)
      {
        MPI_Request requests[] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL};
        int req_counter = 0;
        for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {
          // Determine which buffers we are using in communication
          int recvNeighbour = boxBufs[box_id].recvneigh[iswap];
          AtomBuffer* buf_recv = &boxBufs[box_id].bufs_recv[BORDERS_FUNCTION][box_layer_index][recvNeighbour];
          MPI_Comm* recv_comm = &boxBufs[atom->box_id].comm[BORDERS_FUNCTION];
          int recvtag = swapnum; // Not the same as iswap - swapnum counts over all box layers
          ++swapnum;

          int nrecv = 0;
          // DSM TODO: Double check the comm index makes sense - should the receive be box_comms[sender_id]? What happens if we have multiple sends/receives to the same box at the same time?
          MPI_Recv(&nrecv, 1, MPI_INT, boxBufs[box_id].recvproc[iswap], recvtag, *recv_comm, MPI_STATUS_IGNORE);

          // Skip a send/recv pair if our nsend/nrecv == 0
          if (nrecv > 0) {
            // Increase receiving buffer size if needed
            if (nrecv * atom->border_size > buf_recv->maxsize) buf_recv->growrecv(nrecv * atom->border_size);

            if (sizeof(MMD_float) == 4) {
              MPI_Irecv(buf_recv->buf, nrecv * atom->border_size, MPI_FLOAT,
                       boxBufs[box_id].recvproc[iswap], recvtag, *recv_comm, &requests[req_counter]);
            } else {
              MPI_Irecv(buf_recv->buf, nrecv * atom->border_size, MPI_DOUBLE,
                       boxBufs[box_id].recvproc[iswap], recvtag, *recv_comm, &requests[req_counter]);
            }
            ++req_counter;
          }

          // Save values for unpack task
          boxBufs[box_id].nrecvs_borders[box_layer_index][idim][ineed] = nrecv;
          boxBufs[box_id].buf_recvs_borders[box_layer_index][idim][ineed] = buf_recv;
        }

        if (req_counter > 0) {
          TAMPI_Iwaitall(req_counter, requests, MPI_STATUSES_IGNORE);
        }
      }

      // Commutative to avoid multiple tasks unpacking into same box's x array concurrently
      // Similarly, depends on pack having completed TODO: Could this be removed? Unlikely, what if unpack_border reallocs array as pack is occurring?
      // Additionally depends on this iteration's recv task having finished, rather than all recvs for this box
      #pragma oss task label(borders_unpack) commutative(atom->nghost) \
                       in(bordersRecvSentinels[box_id][idim][iswap]) \
                       in(bordersPackSentinels[box_id]) \
                       out(bordersUnpackSentinels[box_id]) \
                       firstprivate(atom, box_id, iswap, swapnum, idim, ineed)
      {
        int n = 0, nrecv = 0;
        AtomBuffer *buf = NULL;
        for (int box_layer_index = 0; box_layer_index < 3; ++box_layer_index) {
          n = atom->nlocal + atom->nghost;
          nrecv = boxBufs[box_id].nrecvs_borders[box_layer_index][idim][ineed];
          buf = boxBufs[box_id].buf_recvs_borders[box_layer_index][idim][ineed];

          //#pragma omp for
          for (int i = 0; i < nrecv; i++)
            // DSM TODO: Shouldn't this be atom.border_size instead of "4"?
            atom->unpack_border(n + i, &buf->buf[i * 4]);

          /* set all pointers & counters */

          boxBufs[box_id].recvnum[box_layer_index][iswap] = nrecv;
          // DSM: * atom.comm_size (constant '3' defined in constructor) since need 3 doubles per atom (x,y,z coordinates)
          boxBufs[box_id].comm_recv_size[box_layer_index][iswap] = nrecv * atom->comm_size;
          boxBufs[box_id].firstrecv[box_layer_index][iswap] = atom->nlocal + atom->nghost;
          atom->nghost += nrecv;
        }
      }

      // 3 send-recv pairs (1 per layer) done per iswap iteration
      // i.e. per pair of send and recv tasks
      swapnum+=3;

      iswap++;
    } // DSM End of loop over ineed/number of swaps
  } // DSM End of loop over dimensions
}

void Comm::borders_nonblocking(Atom* atoms[])
{
  // Arrays for isend/irecv requests
  // nswaps communications per box layer, 3 box layers per box on this process.
  // Internal swaps in y-dimension are skipped => subtract 6 per box (2 y-swaps per box layer)
  // nswaps and boxes_per_process constant per box
  int nrequests = boxBufs[0].nswap * 3 * atoms[0]->boxes_per_process - 6 * atoms[0]->boxes_per_process;
  int req_i = 0;
  MPI_Request recv_requests[nrequests];
  MPI_Request send_requests[nrequests];
  MPI_Comm *send_comm, *recv_comm;

  // Pack all atoms in slab boundaries into buffers
  for(int box_id = 0; box_id < atoms[0]->boxes_per_process; ++box_id) {
    borders_pack(atoms[box_id]);
  }

  for(int box_id = 0; box_id < atoms[0]->boxes_per_process; ++box_id) {

    Atom& atom = *atoms[box_id];

    int i, m, n, iswap, idim, ineed, nall, nfirst, nlast, type;
    AtomBuffer *buf_send = NULL, *buf_recv = NULL;
    MMD_float* x;
    MMD_float xcoord, ycoord, zcoord;
    // Swap counter for this box
    int swapnum = 200; // to ensure tag is unique accross communicate, borders and exchange

    // DSM Multibox: Reference buffers/variables specific to this box.
    int& nswap = boxBufs[box_id].nswap;
    int& maxnlocal = boxBufs[box_id].maxnlocal;
    int*& sendproc = boxBufs[box_id].sendproc;
    int*& recvproc = boxBufs[box_id].recvproc;
    int*& sendneigh = boxBufs[box_id].sendneigh;
    int*& recvneigh = boxBufs[box_id].recvneigh;
    int& nrecv_atoms = boxBufs[box_id].nrecv_atoms;
    int* (&sendnum)[3] = boxBufs[box_id].sendnum;
    int* (&recvnum)[3] = boxBufs[box_id].recvnum;
    int* (&firstrecv)[3] = boxBufs[box_id].firstrecv;
    int* (&comm_send_size)[3] = boxBufs[box_id].comm_send_size;
    int* (&comm_recv_size)[3] = boxBufs[box_id].comm_recv_size;
    int* (&reverse_send_size)[3] = boxBufs[box_id].reverse_send_size;
    int* (&reverse_recv_size)[3] = boxBufs[box_id].reverse_recv_size;
    AtomBuffer (&bufs_send)[3][8] = boxBufs[box_id].bufs_send[BORDERS_FUNCTION];
    AtomBuffer (&bufs_recv)[3][8] = boxBufs[box_id].bufs_recv[BORDERS_FUNCTION];
    int* (&sendlists)[3][8] = boxBufs[box_id].sendlists;
    int (&maxsendlists)[3][8] = boxBufs[box_id].maxsendlists;
    AtomBuffer& internal_buf_send_up = boxBufs[box_id].internal_buf_send_up[BORDERS_FUNCTION];
    AtomBuffer& internal_buf_send_down = boxBufs[box_id].internal_buf_send_down[BORDERS_FUNCTION];
    AtomBuffer& internal_buf_recv_up = boxBufs[box_id].internal_buf_recv_up[BORDERS_FUNCTION];
    AtomBuffer& internal_buf_recv_down = boxBufs[box_id].internal_buf_recv_down[BORDERS_FUNCTION];
    int*& internal_sendlist_up = boxBufs[box_id].internal_sendlist_up;
    int& internal_maxsendlist_up = boxBufs[box_id].internal_maxsendlist_up;
    int*& internal_sendlist_down = boxBufs[box_id].internal_sendlist_down;
    int& internal_maxsendlist_down = boxBufs[box_id].internal_maxsendlist_down;

    // DSM: This is the atom position array
    x = atom.x;

    // Loop over box layers: box_id, box_id-1, box_id+1
    for (int box_layer_index = 0; box_layer_index < 3; box_layer_index++) {
      iswap = 0;
      // Loop over 5 dimensions x, y, z, diagonal-right, diagonal-left
      for(idim = 0; idim < 5; idim++) {

        // DSM: Do not perform MPI communication in the y / idim==1 dimension
        // Instead, perform memory copy into neighbouring box on this proc's recv buffer.
        // The borders_internal function, called following completion of all borders operations, will then pull this
        // data into each box's atoms list.
        if (idim == 1) {
          // Only perform this copy step once
          if (box_layer_index == 0) {
            // First exchange to lower layer (usually box_id-1 but periodic)
            int box_neighbour = atom.boxneigh_negative;
            buf_send = &internal_buf_send_down;
            buf_recv = &(boxBufs[box_neighbour].internal_buf_recv_down[BORDERS_FUNCTION]);
            int nsend = buf_send->natoms * atom.border_size; // natoms * num MMD_floats per atom (e.g. 4: x,y,z,type)
            // Increase size of receive buffers if receiving more than can be held
            if(nsend > buf_recv->maxsize) {
              buf_recv->growrecv(nsend);
            }
            // Push to neighbour.
            memcpy(buf_recv->buf, buf_send->buf, nsend * sizeof(MMD_float));
            // Update relevant pointers and counters for this swap
            // Recv side done in seperate internal swaps function
            sendnum[box_layer_index][iswap] = buf_send->natoms;
            comm_send_size[box_layer_index][iswap] = buf_send->natoms * atom.comm_size;
            // Set receiving buffer sizes
            buf_recv->natoms = buf_send->natoms;
            // Reset atom counter for next use of this buffer
            buf_send->natoms = 0;
            iswap++;

            // Second exchange to upper layer (usually box_id+1)
            box_neighbour = atom.boxneigh_positive;
            buf_send = &internal_buf_send_up;
            buf_recv = &(boxBufs[box_neighbour].internal_buf_recv_up[BORDERS_FUNCTION]);
            nsend = buf_send->natoms * atom.border_size;
            if(nsend > buf_recv->maxsize) {
              buf_recv->growrecv(nsend);
            }
            memcpy(buf_recv->buf, buf_send->buf, nsend * sizeof(MMD_float));
            sendnum[box_layer_index][iswap] = buf_send->natoms;
            comm_send_size[box_layer_index][iswap] = buf_send->natoms * atom.comm_size;
            buf_recv->natoms = buf_send->natoms;
            buf_send->natoms = 0;
            iswap++;

          } else {
            iswap += 2; // Skip over both sides of the y dimension swap here
          }
          continue; // Skip to next dimension
        }

        for(ineed = 0; ineed < 2 * need[idim]; ineed++) { // DSM Multibox: Will always be ineed < 2 here, need guaranteed to be 1.

          /* swap atoms with other proc
          put incoming ghosts at end of my atom arrays
          if swapping with self, simply copy, no messages */

          //#pragma omp master
          {
            // DSM Two send/recv pairs per dimension
            int send_target_box_id, recv_target_box_id, sendtag, recvtag;
            layer_to_targets(&atom, box_layer_index, &send_target_box_id, &recv_target_box_id);

            // Determine which buffers we are using in communication
            int sendNeighbour = sendneigh[iswap];
            buf_send = &bufs_send[box_layer_index][sendNeighbour];
            sendnum[box_layer_index][iswap] = buf_send->natoms; // Number of atoms to send (*not* MMD_Floats)
            buf_send->natoms = 0; // Reset atom counter for next use of this buffer
            send_comm = &boxBufs[send_target_box_id].comm[BORDERS_FUNCTION];
            recv_comm = &boxBufs[atom.box_id].comm[BORDERS_FUNCTION];
            sendtag = swapnum;
            recvtag = sendtag;
            ++swapnum;

            if(sendproc[iswap] == me && send_target_box_id == atom.box_id) {
              // Simply copy if swapping with yourself
              recvnum[box_layer_index][iswap] = sendnum[box_layer_index][iswap];
            } else {
              // DSM: Changed to Irecv/Send pair from Send/Recv. No guarantee original code progresses (likely will due to buffered send)
              // DSM TODO: Double check the comm index makes sense - should the receive be box_comms[sender_id]? What happens if we have multiple sends/receives to the same box at the same time?
              // DSM TODO: Optimise this - we can skip a send/recv pair if our nsend/nrecv == 0
              MPI_Irecv(&recvnum[box_layer_index][iswap], 1, MPI_INT, recvproc[iswap], recvtag, *recv_comm, &recv_requests[req_i]);
              MPI_Isend(&sendnum[box_layer_index][iswap], 1, MPI_INT, sendproc[iswap], sendtag, *send_comm, &send_requests[req_i]);
              req_i++;
            } // End of branch
            iswap++;
          } // End of omp block
        } // End of loop over swaps (ineed)
      } // End of loop over dimensions (idim)
    } // End of loop over box layers
  } // End of loop over boxes per process

  // All size messages posted. Wait for completion
  // req_i rather than nrequests as may have skipped requests if sending to self
  MPI_Waitall(req_i, recv_requests, MPI_STATUSES_IGNORE);
  MPI_Waitall(req_i, send_requests, MPI_STATUSES_IGNORE);

  // Now post actual data
  req_i = 0;
  for(int box_id = 0; box_id < atoms[0]->boxes_per_process; ++box_id) {

    Atom &atom = *atoms[box_id];
    // Swap counter for this box
    int swapnum = 200; // to ensure tag is unique accross communicate, borders and exchange

    for (int box_layer_index = 0; box_layer_index < 3; box_layer_index++) {

      int iswap = 0;
      int*& sendproc = boxBufs[box_id].sendproc;
      int*& recvproc = boxBufs[box_id].recvproc;
      int*& sendneigh = boxBufs[box_id].sendneigh;
      int*& recvneigh = boxBufs[box_id].recvneigh;
      int* (&sendnum)[3] = boxBufs[box_id].sendnum;
      int* (&recvnum)[3] = boxBufs[box_id].recvnum;
      AtomBuffer (&bufs_send)[3][8] = boxBufs[box_id].bufs_send[BORDERS_FUNCTION];
      AtomBuffer (&bufs_recv)[3][8] = boxBufs[box_id].bufs_recv[BORDERS_FUNCTION];
      AtomBuffer *buf_send = NULL, *buf_recv = NULL;

      for (int idim = 0; idim < 5; idim++) {

        // Skip y-dimension
        if (idim == 1) {
          iswap += 2;
          continue;
        }

        for (int ineed = 0; ineed < 2 * need[idim]; ineed++) {

          int nsend = sendnum[box_layer_index][iswap];
          int nrecv = recvnum[box_layer_index][iswap];
          int send_target_box_id, recv_target_box_id, sendtag, recvtag;
          layer_to_targets(&atom, box_layer_index, &send_target_box_id, &recv_target_box_id);

          // Determine which buffers we are using in communication
          int sendNeighbour = sendneigh[iswap];
          int recvNeighbour = recvneigh[iswap];
          buf_send = &bufs_send[box_layer_index][sendNeighbour];
          buf_recv = &bufs_recv[box_layer_index][recvNeighbour];
          send_comm = &boxBufs[send_target_box_id].comm[BORDERS_FUNCTION];
          recv_comm = &boxBufs[atom.box_id].comm[BORDERS_FUNCTION];
          sendtag = swapnum;
          recvtag = sendtag;
          ++swapnum;

          // Skip sending to self (data will be unpacked directly from send buffer in next step)
          if(sendproc[iswap] == me && send_target_box_id == atom.box_id) {
            iswap++;
            continue;
          }

          // Increase receiving buffer size if needed
          if (nrecv * atom.border_size > buf_recv->maxsize) buf_recv->growrecv(nrecv * atom.border_size);

          if (sizeof(MMD_float) == 4) {
            MPI_Irecv(buf_recv->buf, nrecv * atom.border_size, MPI_FLOAT,
                      recvproc[iswap], recvtag, *recv_comm, &recv_requests[req_i]);
            MPI_Isend(buf_send->buf, nsend * atom.border_size, MPI_FLOAT,
                     sendproc[iswap], sendtag, *send_comm, &send_requests[req_i]);
          } else {
            MPI_Irecv(buf_recv->buf, nrecv * atom.border_size, MPI_DOUBLE,
                      recvproc[iswap], recvtag, *recv_comm, &recv_requests[req_i]);
            MPI_Isend(buf_send->buf, nsend * atom.border_size, MPI_DOUBLE,
                     sendproc[iswap], sendtag, *send_comm, &send_requests[req_i]);
          }
          iswap++;
          req_i++;
        } // End of loop over ineed (swaps)
      } // End of loop over idim
    } // End of loop over box layers
  } // End of loop over boxes per process

  // Wait for all data communication to complete
  // req_i rather than nrequests as may have skipped requests if sending to self
  MPI_Waitall(req_i, recv_requests, MPI_STATUSES_IGNORE);
  MPI_Waitall(req_i, send_requests, MPI_STATUSES_IGNORE);

  // Now unpack data
  for(int box_id = 0; box_id < atoms[0]->boxes_per_process; ++box_id) {
    Atom &atom = *atoms[box_id];
    for (int box_layer_index = 0; box_layer_index < 3; box_layer_index++) {

      int iswap = 0;
      int*& sendproc = boxBufs[box_id].sendproc;
      int*& recvproc = boxBufs[box_id].recvproc;
      int*& sendneigh = boxBufs[box_id].sendneigh;
      int*& recvneigh = boxBufs[box_id].recvneigh;
      int* (&sendnum)[3] = boxBufs[box_id].sendnum;
      int* (&recvnum)[3] = boxBufs[box_id].recvnum;
      AtomBuffer (&bufs_send)[3][8] = boxBufs[box_id].bufs_send[BORDERS_FUNCTION];
      AtomBuffer (&bufs_recv)[3][8] = boxBufs[box_id].bufs_recv[BORDERS_FUNCTION];
      int* (&firstrecv)[3] = boxBufs[box_id].firstrecv;
      int* (&comm_send_size)[3] = boxBufs[box_id].comm_send_size;
      int* (&comm_recv_size)[3] = boxBufs[box_id].comm_recv_size;
      int* (&reverse_send_size)[3] = boxBufs[box_id].reverse_send_size;
      int* (&reverse_recv_size)[3] = boxBufs[box_id].reverse_recv_size;

      AtomBuffer *buf = NULL, *buf_send = NULL;

      for (int idim = 0; idim < 5; idim++) {

        // Skip y-dimension
        if (idim == 1) {
          iswap += 2;
          continue;
        }

        for (int ineed = 0; ineed < 2 * need[idim]; ineed++) {

          int nsend = sendnum[box_layer_index][iswap];
          int nrecv = recvnum[box_layer_index][iswap];
          int send_target_box_id, recv_target_box_id;
          layer_to_targets(&atom, box_layer_index, &send_target_box_id, &recv_target_box_id);

          // Determine which buffers we are using in communication
          int sendNeighbour = sendneigh[iswap];
          int recvNeighbour = recvneigh[iswap];
          buf_send = &bufs_send[box_layer_index][sendNeighbour];
          buf = &bufs_recv[box_layer_index][recvNeighbour];

          // Sending to self. Just unpack the send buffer
          if(sendproc[iswap] == me && send_target_box_id == atom.box_id) {
            buf = buf_send;
          }

          int n = atom.nlocal + atom.nghost;

          //#pragma omp for
          for(int i = 0; i < nrecv; i++) {
            atom.unpack_border(n + i, &buf->buf[i * 4]); // DSM TODO: Shouldn't this be atom.border_size instead of "4"?
          }

          /* set all pointers & counters */

          //#pragma omp master
          {
            // DSM: * atom.comm_size (constant '3' defined in constructor) since need 3 doubles per atom (x,y,z coordinates)
            comm_send_size[box_layer_index][iswap] = nsend * atom.comm_size;
            comm_recv_size[box_layer_index][iswap] = nrecv * atom.comm_size;
            reverse_send_size[box_layer_index][iswap] = nrecv * atom.reverse_size;
            reverse_recv_size[box_layer_index][iswap] = nsend * atom.reverse_size;
            firstrecv[box_layer_index][iswap] = atom.nlocal + atom.nghost;
            atom.nghost += nrecv;
          }

          iswap++;
        } // End of loop over ineed (swaps)
      } // End of loop over idim
    } // End of loop over box layers
  } // End of loop over boxes per process
} // End of function

// DSM Multibox: changed below buffer reallocation functions to work per box

/* realloc the size of the send buffer as needed with BUFFACTOR & BUFEXTRA */

void AtomBuffer::growsend(int n) {
  maxsize = static_cast<int>(BUFFACTOR * n);
  buf = (MMD_float*) realloc(buf, (maxsize + BUFEXTRA) * sizeof(MMD_float));

  /*
  if (!buf) {
    printf("Failed to grow send buffer by %d! Out of memory? Aborting!\n", n);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  */
}

/* free/malloc the size of the recv buffer as needed with BUFFACTOR */

void AtomBuffer::growrecv(int n) {
  maxsize = static_cast<int>(BUFFACTOR * n);
  free(buf);
  buf = (MMD_float*) malloc(maxsize * sizeof(MMD_float));
}

// Replaces equivalent in atom.pack_border()
// Moved to an AtomBuffer method to facilitate PBC correction and updating this buffer's sendlist
// Takes atom x, y, z coordinates and type as arguments rather then indexing into atom.x[] and atom.type[]
// Still takes index i to record in the sendlist
int AtomBuffer::pack_border(int i, MMD_float x, MMD_float y, MMD_float z, int type, int** sendlist, int& maxsendlist)
{
  int m = natoms*4; // Add new atom to end of buffer

  // Moved buffer size checks into this function rather than performing them before every call to this function
  if ((natoms+1)*4 > maxsize) {
    growsend(maxsize);
  }

  // DSM: Is this branch necessary? pbc_flags are set to 0 if not periodic in a dimension - are the branches not equivalent?
  // Does adding 0.0 cause rounding errors?
  if(pbc_any == 0) {
    buf[m++] = x;
    buf[m++] = y;
    buf[m++] = z;
    buf[m++] = type;
  } else {
    buf[m++] = x + pbc_x;
    buf[m++] = y + pbc_y;
    buf[m++] = z + pbc_z;
    buf[m++] = type;
  }

  // Mark this atom for sending in future communicate() calls
  // realloc the size of the sendlist as needed with BUFFACTOR
  // Comm::growlist made obsolete by this method
  if (natoms+1 >= maxsendlist) {
    // Replaces old "growlist" function
    // No longer scaling at same rate as original code. Original code knew total number of atoms to send at pack time
    // and would base buffer growth on this (BUFFACTOR * natoms). This does BUFFACTOR * old maximum size - likely to
    // cause more reallocs
    maxsendlist = static_cast<int>(BUFFACTOR * maxsendlist);
    *sendlist = (int*) realloc(*sendlist, maxsendlist * sizeof(int));
  }
  (*sendlist)[natoms] = i;

  // Increment count of atoms held in this buffer
  natoms++;

  return m;
}
