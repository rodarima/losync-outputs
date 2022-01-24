/// \file
/// Functions to maintain link cell structures for fast pair finding.

#ifndef __LINK_CELLS_H_
#define __LINK_CELLS_H_

#include "mytype.h"
#include "initAtoms.h"

#include <mpi.h>

/// The maximum number of atoms that can be stored in a link cell.
#define MAXATOMS 64 

struct DomainSt;
struct AtomsSt;

/// Link cell data.  For convenience, we keep a copy of the localMin and
/// localMax coordinates that are also found in the DomainsSt.
typedef struct LinkCellSt
{
    int gridSize[3];           //!< number of boxes in each dimension on processor
    int nLocalBoxes;           //!< total number of local boxes on processor
    int nHaloBoxes;            //!< total number of remote halo/ghost boxes on processor
    int nTotalBoxes;           //!< total number of boxes on processor
    //!< nLocalBoxes + nHaloBoxes
    real3 localMin;            //!< minimum local bounds on processor
    real3 localMax;            //!< maximum local bounds on processor
    real3 boxSize;             //!< size of box in each dimension
    real3 invBoxSize;          //!< inverse size of box in each dimension

    int* nAtoms;               //!< total number of atoms in each box
    int** nbrBoxes;            //!< neighbor boxes for each box

    // New additions
    int* faces;                //!< index of the face each box is located on
    int numColumns;            //!< number of columns
    int* columnNumNeighbours;  //!< number of neighbours of each column
    int*** columnCommNeighbours;   //!< neighbours of each column
    Atoms** pendingAtoms;         //!< atoms that are due to be transferred to the given column
    int* numInternalNeighbours; //!< number of neighbouring columns each column has in the local domain (including halos)
    int** internalNeighbours;   //!< list of column ids that neighbour each column in the local domain (including halos)

    char*** sendBufs;
    // Separate set of send buffers for each box
    char*** recvBufs;
    // Separate set of send requests for each box
    MPI_Request** sendRequests;
    // Separate set of receive requests for each box
    MPI_Request** recvRequests;
    // Number of bytes to be sent for each box
    int** nSends;
    // Number of bytes to be received from each neighbour
    int** nRecvs;
    // Need to store receive statuses so that message info can be accessed
    // outside the current task
    MPI_Status** recvStatus; 

    // Sentinels
    char* firstMallocSentinel;
    char* secondMallocSentinel;
    char* loadSentinel;
    char* sendSentinel;
    char* recvSentinel;
    char* unloadSentinel;
    char* firstVelocitySentinel;
    char* secondVelocitySentinel;
    char* positionSentinel;
    char* forceSentinel;
    char* zeroSentinel;
    char* firstFreeSentinel;
    char* secondFreeSentinel;
    char* firstUpdateSentinel;
    char* secondUpdateSentinel;

    // Various buffers
    // Separate set of send buffers for each box
    char*** sendBufs;
    // Separate set of send buffers for each box
    char*** recvBufs;
    // Separate set of send requests for each box
    MPI_Request** sendRequests;
    // Separate set of receive requests for each box
    MPI_Request** recvRequests;
    // Number of bytes to be sent for each box
    int** nSends;
    // Number of bytes to be received from each neighbour
    int** nRecvs;
    // Need to store receive statuses so that message info can be accessed
    // outside the current task
    MPI_Status** recvStatus; 
    // Number of neighbours each box has
    int* nColumnNeighbours;


} LinkCell;

LinkCell* initLinkCells(const struct DomainSt* domain, real_t cutoff);
void destroyLinkCells(LinkCell** boxes);

int getNeighborBoxes(LinkCell* boxes, int iBox, int* nbrBoxes);
void putAtomInBox(LinkCell* boxes, struct AtomsSt* atoms,
        const int gid, const int iType,
        const real_t x,  const real_t y,  const real_t z,
        const real_t px, const real_t py, const real_t pz);
int getBoxFromTuple(LinkCell* boxes, int x, int y, int z);
int getFaceFromBox(int iBox);

void moveAtom(LinkCell* boxes, struct AtomsSt* atoms, int iId, int iBox, int jBox);

/// Update link cell data structures when the atoms have moved.
void updateLinkCells(LinkCell* boxes, struct AtomsSt* atoms);

int maxOccupancy(LinkCell* boxes);

// New addition
void addNeighbour(int* neighbourList, int *nNeighbours, int newNeighbour);
void getTuple(LinkCell* boxes, int iBox, int* ixp, int* iyp, int* izp);


#endif
