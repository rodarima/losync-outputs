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

#ifndef COMM_H
#define COMM_H

#include "atom.h"
#include "threadData.h"
#include "timer.h"
#include <mpi.h> // DSM Multibox. Needed for communicator arrays

#ifdef USE_TAMPI
#include <TAMPI.h>
#endif

// DSM: Holds pointer to memory buffer of MMD_floats and its maximum size
class AtomBuffer {
  public:
    int maxsize;        // maximum number of MMD_floats that can be sent by/received into corresponding buffer
    int natoms;         // number of atoms within this buffer
    MMD_float* buf;     // buffer for all comms

    // PBC flags moved here - one set per (send) buffer rather than one per swap
    // Remove separate pbc_flagx, pbc_flagy, pbc_flagz variables and replaced with value to correct by in each
    // dimension. e.g. pbc_x will be one of -box.xprd, 0, box.xprd. Saves having to compute the product every time in
    // packing routines. pbc_any kept TODO: is pbc_any necessary?
    int pbc_any;        // whether any PBC on this buffer
    MMD_float pbc_x;    // PBC correction in x for this swap
    MMD_float pbc_y;    // same in y
    MMD_float pbc_z;    // same in z

    // Used for internal buffers only. Required to know where to unpack 1st atom to outside the communicate() function
    int internal_firstrecv;

    void growsend(int); // Increase size of buffer to be able to store at least the given number of items
    void growrecv(int); // As above but different factor for size increase. More suited to recieve buffer pattern

    // Pack this buffer for border communication. Replaces equivalent atom.pack_border()
    int pack_border(int, MMD_float, MMD_float, MMD_float, int, int**, int&);
};

// DSM Multibox: Move majority of Comm attributes into a struct to duplicate per box. Better solution would be to
// rewrite program to use multiple Comm object instances but this would require greater structural changes.
struct BoxBufs {
    int nswap;                        // # of swaps to perform

    // DSM These now duplicated per box layer
    int* sendnum[3], *recvnum[3];     // # of atoms to send/recv in each swap
    int* firstrecv[3];                // where to put 1st recv atom in each swap
    int* comm_send_size[3];           // # of values to send in each comm
    int* comm_recv_size[3];           // # of values to recv in each comm
    int* reverse_send_size[3];        // # of values to send in each reverse
    int* reverse_recv_size[3];        // # of values to recv in each reverse

    int* sendproc, *recvproc;         // proc to send/recv with at each swap
    int* sendproc_exc, *recvproc_exc; // proc to send/recv with at each swap for safe exchange
    // DSM Multibox change
    int* sendneigh, *recvneigh;       // neighbour directions (UP, RIGHT, etc.) to send/recv with at each swap
    MPI_Request request;              // outstanding MPI request for this box
    MPI_Comm comm[3];                 // communicator for this box, 3 as 1 per communication function (as bufs below)

    //MMD_float* buf;
    // DSM Replaced with 3*24 buffers,
    //       one for each exchange between neighbours
    //   and one for each function that performs communication: communicate(), exchange(), borders()
    //     = 3 functions * 3 layers * 8 neighbours
    // Indexed by [function][box_layer][neighbour].
    // 3 functions: COMMUNICATE_FUNCTION, EXCHANGE_FUNCTION, BORDERS_FUNCTION
    // 3 box layers: box_id (same layer), box_id-1 (below), box_id+1 (above)
    // 8 neighbours: UP, TOP_RIGHT, RIGHT, etc.

    AtomBuffer bufs_send[3][3][8];                 // send buffers for all comm
    AtomBuffer bufs_recv[3][3][8];                 // recv buffers for all comm
    int* sendlists[3][8];                          // list of atoms for each buffer to send in each swap. Used in borders() and communicate()
    int maxsendlists[3][8];                        // size of sendlist

    // DSM buffers for special case: internal memory transfer in y-direction
    // Again 3 buffers, 1 per function that performs communication
    AtomBuffer internal_buf_send_up[3];
    AtomBuffer internal_buf_send_down[3];
    AtomBuffer internal_buf_recv_up[3];
    AtomBuffer internal_buf_recv_down[3];
    int* internal_sendlist_up;
    int internal_maxsendlist_up;
    int* internal_sendlist_down;
    int internal_maxsendlist_down;

    // Used in exchange. Saves recv buffer pointers and sizes from recv tasks to be used in unpack tasks
    // Attributes to avoid going out of scope before tasks are complete
    // 3 box layers, 5 dimensions per swap
    AtomBuffer* buf_recvs1[3][5];
    int nrecvs1[3][5];
    AtomBuffer* buf_recvs2[3][5];
    int nrecvs2[3][5];
    // As above but for borders
    // 3 box layers, 5 dimensions, 2 swaps
    AtomBuffer* buf_recvs_borders[3][5][2];
    int nrecvs_borders[3][5][2];

    // Buffers used exclusively in communicate meighbour tasks implmentation using a single message per neighbour
    AtomBuffer buf_send_single;
    AtomBuffer buf_recv_single;

    // DSM Replaced these with xneg_slab_lo, etc. constants in atom.box struct
    //MMD_float* slablo, *slabhi;          // bounds of slabs to send to other procs

    int copy_size;
    int* nsend_thread;
    int* nrecv_thread;
    int* nholes_thread;
    int** exc_sendlist_thread;
    int* send_flag;
    int* maxsend_thread;
    int maxthreads;
    int maxnlocal;
    int nrecv_atoms;
};

class Comm
{
  public:
    Comm(int); // DSM Multibox: constructor now takes boxes_per_process as an argument
    ~Comm();
    //void free_box_comms(); // DSM Multibox: for freeing box communicators
    int setup(MMD_float, int, int, Atom* [], int); // DSM Multibox: Array, two process grid int and nonblocking flag arguments
    // Communicate
    void communicate(Atom**); // DSM Multibox: Now takes array of Atom** structs, one per box
    void communicate_nonblocking(Atom**); // Nonblocking implementation
    void communicate_blocking(Atom &, int); // DSM Multibox: Blocking implementation. Now also takes int box_id
    void communicate_blocking_isend(Atom &, int);
    void communicate_blocking_alltasks(Atom *, int);
    void communicate_blocking_alltasks_recvfirst(Atom *, int);
    void communicate_nonblocking_neighbourtasks(Atom *, int);
    void communicate_blocking_neighbourtasks(Atom *, int);
    void communicate_blocking_single_message_neighbourtasks(Atom *, int);
    void communicate_nonblocking_alltasks_tampi_iwait(Atom *, int);
    void communicate_nonblocking_neighbourtasks_tampi_iwaitall(Atom *, int);
    void communicate_internal(Atom &, int); // communicate() equivalent of exchange_internal()
    void communicate_internal_send(Atom*, AtomBuffer*, AtomBuffer*, int*, int, int);
    void reverse_communicate(Atom &);
    // Exchange
    void exchange(Atom**);
    void exchange_nonblocking(Atom**);
    void exchange_blocking(Atom &, int);
    void exchange_blocking_neighbourtasks(Atom*, int);
    void exchange_nonblocking_neighbourtasks_tampi_iwaitall(Atom*, int);
    void exchange_pack(Atom*, AtomBuffer bufs_send[3][8], AtomBuffer*, AtomBuffer*);
    void exchange_internal(Atom &, int); // DSM Multibox: New function that performs exchange between boxes on same proc
    void exchange_internal_send(AtomBuffer*, AtomBuffer*);
    void exchange_internal_recv(Atom*, AtomBuffer*);
    void exchange_all(Atom &);
    // Borders
    void borders(Atom**);
    void borders_nonblocking(Atom**);
    void borders_blocking(Atom &, int); // DSM Multibox: Now also takes int box_id
    void borders_blocking_neighbourtasks(Atom *, int);
    void blocking_nonblocking_neighbourtasks_tampi_iwaitall(Atom *, int);
    void borders_internal(Atom &, int); // borders() equivalent of exchange_internal()
    void borders_internal_send(Atom*);
    void borders_pack(Atom* atom);
    // DSM Multibox: changed below buffer reallocation functions to work per box and per buffer
    void growsend(int, int);
    void growrecv(int, int);
    void growlist(int, int, int);
    // Sentinels for tasking
    char *initialIntegrateSentinels;
    char *sortSentinels;
    char *forceComputeSentinels;
    char *finalIntegrateSentinels;
    char *communicateSentinels;
    char *communicateInternalPackSentinels;
    char *communicateInternalUnpackSentinels;
    char **communicatePackSentinels;
    char **communicateSendSentinel;
    char **communicateRecvSentinels;
    char *exchangeSentinels;
    char *exchangePackSentinels;
    char *exchangeInternalSendSentinels;
    char *exchangeInternalRecvSentinels;
    char *exchangePBCSentinels;
    char **exchangeSend1Sentinels;
    char **exchangeRecv1Sentinels;
    char **exchangeSend2Sentinels;
    char **exchangeRecv2Sentinels;
    char *bordersSentinels;
    char *bordersPackSentinels;
    char *bordersUnpackSentinels;
    char *bordersSendSentinels;
    char (*bordersRecvSentinels)[5][2]; // 5 dimensions, 2 swaps per dimension
    char *bordersInternalSentinels;
    char *bordersInternalSendSentinels;
    char *neighbourBuildSentinels;

  public:
    int me;                           // my proc ID
    // DSM: Multibox change.
    int boxes_per_process;            // number of boxes on each MPI process
    BoxBufs* boxBufs;                 // duplicated buffers and counters required per box

    // DSM: Remaining attributes not duplicated in boxBufs, i.e. constants per box (or values relating to broken code)
    // DSM To switch to 26-way communication, procneigh is now 5 dimensional array:
    //   - Sides (original 2 dimensions): x, z
    //   - Corners (2 new dimensions): "top-right" & "bottom-right", "top-left" & "bottom-left"
    //   - Internal communication between boxes: y
    // (Sides+Corners)*3 for box_id+1, box_id, and box_id-1 layers = 4*3 = 12 dimensions. NB: Do not need to have
    // separate elements in procneigh for each box layer. Process ranks will be the same on each layer, only target box
    // IDs differ.
    // + 1 for internal communication in y dimension = 13 dimensions total.
    // Exchange in +ve and -ve direction in each dimension => 13*2 = 26-way communication between process neighbours
    int procneigh[5][2];              // my 26 proc neighs (procneigh[1][0:1] will always be my rank as boxes split in y dimension)
    int procgrid[3];                  // # of procs in each dim
    int need[5];                      // how many procs away needed in each dim
    int boxgrid[3];                   // # of boxes in each dim (analogous to procgrid)

    ThreadData* threads;		    //

    int check_safeexchange;           // if sets give warnings if an atom moves further than subdomain size
    int do_safeexchange;		          // exchange atoms with all subdomains within neighbor cutoff
    Timer* timer;

    // DSM Flag enabling/disabling non-blocking communication mode
    int nonblocking_enabled;
};

#endif
