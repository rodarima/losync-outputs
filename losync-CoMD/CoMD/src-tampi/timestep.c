/// \file
/// Leapfrog time integrator

#include "timestep.h"

#include <omp.h>
#include "mpi.h"

#include "CoMDTypes.h"
#include "linkCells.h"
#include "parallel.h"
#include "performanceTimers.h"

static void firstAdvanceVelocity(SimFlat* s, int nBoxes, real_t dt);
static void secondAdvanceVelocity(SimFlat* s, int nBoxes, real_t dt);
static void advancePosition(SimFlat* s, int nBoxes, real_t dt);


/// Advance the simulation time to t+dt using a leap frog method
/// (equivalent to velocity verlet).
///
/// Forces must be computed before calling the integrator the first time.
///
///  - Advance velocities half time step using forces
///  - Advance positions full time step using velocities
///  - Update link cells and exchange remote particles
///  - Compute forces
///  - Update velocities half time step using forces
///
/// This leaves positions, velocities, and forces at t+dt, with the
/// forces ready to perform the half step velocity update at the top of
/// the next call.
///
/// After nSteps the kinetic energy is computed for diagnostic output.
double timestep(SimFlat* s, int nSteps, real_t dt)
{
    computeForce(s,1);
#pragma oss taskwait
    for (int ii=0; ii<nSteps; ++ii)
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        startTimer(velocityTimer);
        firstAdvanceVelocity(s, s->boxes->nLocalBoxes, 0.5*dt); 
        stopTimer(velocityTimer);

        startTimer(positionTimer);
        advancePosition(s, s->boxes->nLocalBoxes, dt);
        stopTimer(positionTimer);

        startTimer(redistributeTimer);
        redistributeAtoms(s);
        stopTimer(redistributeTimer);

        startTimer(computeForceTimer);
        computeForce(s,0);
        stopTimer(computeForceTimer);

        startTimer(velocityTimer);
        secondAdvanceVelocity(s, s->boxes->nLocalBoxes, 0.5*dt); 

        stopTimer(velocityTimer);
    }

#pragma oss taskwait
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    computeForce(s,1);
    kineticEnergy(s);

    return s->ePotential;
}

void computeForce(SimFlat* s, int withEnergy)
{
    s->pot->force(s, withEnergy);
}

void firstAdvanceVelocity(SimFlat* s, int nBoxes, real_t dt)
{
    LinkCell* ll = ((SimFlat*) s)->boxes;
    int nCells = ll->gridSize[2];

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    for (int x=-1;x<ll->gridSize[0]+1;++x)
    {
        for (int y=-1;y<ll->gridSize[1]+1;++y)
        {
            int iColumn = y+1+(x+1)*(ll->gridSize[1]+2);
#pragma oss task firstprivate(iColumn) \
            in(ll->secondVelocitySentinel[iColumn]) \
            out(ll->firstVelocitySentinel[iColumn])
            {
                for (int iCell=-1; iCell<nCells+1; ++iCell)
                {
                    int iBox = getBoxFromTuple(s->boxes,x,y,iCell);

                    for (int iOff=MAXATOMS*iBox,ii=0; ii<s->boxes->nAtoms[iBox]; ii++,iOff++)
                    {
                        s->atoms->p[iOff][0] += dt*s->atoms->f[iOff][0];
                        s->atoms->p[iOff][1] += dt*s->atoms->f[iOff][1];
                        s->atoms->p[iOff][2] += dt*s->atoms->f[iOff][2];
                    }
                }
            }
        }
    }
}

void secondAdvanceVelocity(SimFlat* s, int nBoxes, real_t dt)
{
    LinkCell* ll = ((SimFlat*) s)->boxes;
    int nCells = ll->gridSize[2];

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    for (int x=-1;x<ll->gridSize[0]+1;++x)
    {
        for (int y=-1;y<ll->gridSize[1]+1;++y)
        {
            int iColumn = y+1+(x+1)*(ll->gridSize[1]+2);
#pragma oss task firstprivate(iColumn) \
            in(ll->forceSentinel[iColumn]) \
            out(ll->secondVelocitySentinel[iColumn])
            {
                for (int iCell=-1; iCell<nCells+1; ++iCell)
                {
                    int iBox = getBoxFromTuple(s->boxes,x,y,iCell);

                    for (int iOff=MAXATOMS*iBox,ii=0; ii<s->boxes->nAtoms[iBox]; ii++,iOff++)
                    {
                        s->atoms->p[iOff][0] += dt*s->atoms->f[iOff][0];
                        s->atoms->p[iOff][1] += dt*s->atoms->f[iOff][1];
                        s->atoms->p[iOff][2] += dt*s->atoms->f[iOff][2];
                    }
                }
            }
        }
    }
}

void advancePosition(SimFlat* s, int nBoxes, real_t dt)
{
    LinkCell* ll = ((SimFlat*) s)->boxes;
    int nCells = ll->gridSize[2];

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    for (int x=-1;x<ll->gridSize[0]+1;++x)
    {
        for (int y=-1;y<ll->gridSize[1]+1;++y)
        {
            int iColumn = y+1+(x+1)*(ll->gridSize[1]+2);
#pragma oss task firstprivate(iColumn) \
            in(ll->firstVelocitySentinel[iColumn]) \
            out(ll->positionSentinel[iColumn])
            {
                for (int iCell=-1; iCell<nCells+1; ++iCell)
                {
                    int iBox = getBoxFromTuple(s->boxes,x,y,iCell);

                    for (int iOff=MAXATOMS*iBox,ii=0; ii<ll->nAtoms[iBox]; ii++,iOff++)
                    {
                        int iSpecies = s->atoms->iSpecies[iOff];
                        real_t invMass = 1.0/s->species[iSpecies].mass;
                        s->atoms->r[iOff][0] += dt*s->atoms->p[iOff][0]*invMass;
                        s->atoms->r[iOff][1] += dt*s->atoms->p[iOff][1]*invMass;
                        s->atoms->r[iOff][2] += dt*s->atoms->p[iOff][2]*invMass;
                    }
                }
            }
        }
    }
}

/// Calculates total kinetic and potential energy across all tasks.  The
/// local potential energy is a by-product of the force routine.
void kineticEnergy(SimFlat* s)
{
    real_t eLocal[2];
    real_t kenergy = 0.0;
    eLocal[0] = s->ePotential;
    eLocal[1] = 0;
    for (int iBox=0; iBox<s->boxes->nLocalBoxes; iBox++)
    {
        for (int iOff=MAXATOMS*iBox,ii=0; ii<s->boxes->nAtoms[iBox]; ii++,iOff++)
        {
            int iSpecies = s->atoms->iSpecies[iOff];
            real_t invMass = 0.5/s->species[iSpecies].mass;
            kenergy += ( s->atoms->p[iOff][0] * s->atoms->p[iOff][0] +
                    s->atoms->p[iOff][1] * s->atoms->p[iOff][1] +
                    s->atoms->p[iOff][2] * s->atoms->p[iOff][2] )*invMass;
        }
    }

    eLocal[1] = kenergy;

    real_t eSum[2];
    startTimer(commReduceTimer);
    addRealParallel(eLocal, eSum, 2);
    stopTimer(commReduceTimer);

    s->ePotential = eSum[0];
    s->eKinetic = eSum[1];
}

/// \details
/// This function provides one-stop shopping for the sequence of events
/// that must occur for a proper exchange of halo atoms after the atom
/// positions have been updated by the integrator.
///
/// - updateLinkCells: Since atoms have moved, some may be in the wrong
///   link cells.
/// - haloExchange (atom version): Sends atom data to remote tasks. 
/// - sort: Sort the atoms.
///
/// \see updateLinkCells
/// \see initAtomHaloExchange
/// \see sortAtomsInCell
void redistributeAtoms(SimFlat* sim)
{
    updateLinkCells(sim->boxes, sim->atoms);

    startTimer(atomHaloTimer);
    haloExchange(sim->atomExchange, sim);
    stopTimer(atomHaloTimer);

}
