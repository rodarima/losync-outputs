/*
 * Mesh.cc
 *
 *  Created on: Jan 5, 2012
 *      Author: cferenba
 *
 * Copyright (c) 2012, Triad National Security, LLC.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style open-source
 * license; see top-level LICENSE file for full license text.
 */

#include "Mesh.hh"

#include <stdint.h>
#include <cmath>
#include <iostream>
#include <algorithm>

#include "Vec2.hh"
#include "Memory.hh"
#include "Parallel.hh"
#include "InputFile.hh"
#include "GenMesh.hh"
#include "WriteXY.hh"
#include "ExportGold.hh"

using namespace std;


Mesh::Mesh(const InputFile* inp) :
    gmesh(NULL), egold(NULL), wxy(NULL) {

    using Parallel::mype;

    chunksize = inp->getInt("chunksize", 0);
    if (chunksize < 0) {
        if (mype == 0)
            cerr << "Error: bad chunksize " << chunksize << endl;
        exit(1);
    }

    subregion = inp->getDoubleList("subregion", vector<double>());
    if (subregion.size() != 0 && subregion.size() != 4) {
        if (mype == 0)
            cerr << "Error:  subregion must have 4 entries" << endl;
        exit(1);
    }

    writexy = inp->getInt("writexy", 0);
    writegold = inp->getInt("writegold", 0);

    gmesh = new GenMesh(inp);
    wxy = new WriteXY(this);
    egold = new ExportGold(this);

    init();
}


Mesh::~Mesh() {
    delete gmesh;
    delete wxy;
    delete egold;
}


void Mesh::init() {

    // generate mesh
    vector<double2> nodepos;
    vector<int> cellstart, cellsize, cellnodes;
    vector<int> slavemstrpes, slavemstrcounts, slavepoints;
    vector<int> masterslvpes, masterslvcounts, masterpoints;
    gmesh->generate(nodepos, cellstart, cellsize, cellnodes,
            slavemstrpes, slavemstrcounts, slavepoints,
            masterslvpes, masterslvcounts, masterpoints);

    nump = nodepos.size();
    numz = cellstart.size();
    nums = cellnodes.size();
    numc = nums;

    // copy cell sizes to mesh
    znump = Memory::alloc<int>(numz);
    copy(cellsize.begin(), cellsize.end(), znump);

    // populate maps:
    // use the cell* arrays to populate the side maps
    initSides(cellstart, cellsize, cellnodes);
    // release memory from cell* arrays
    cellstart.resize(0);
    cellsize.resize(0);
    cellnodes.resize(0);
    // now populate edge maps using side maps
    initEdges();

    // populate chunk information
    initChunks();

    // create inverse map for corner-to-point gathers
    initInvMap();

    // calculate parallel data structures
    initParallel(slavemstrpes, slavemstrcounts, slavepoints,
            masterslvpes, masterslvcounts, masterpoints);
    // release memory from parallel-related arrays
    slavemstrpes.resize(0);
    slavemstrcounts.resize(0);
    slavepoints.resize(0);
    masterslvpes.resize(0);
    masterslvcounts.resize(0);
    masterpoints.resize(0);

    // write mesh statistics
    writeStats();

    // allocate remaining arrays
    px = Memory::alloc<double2>(nump);
    ex = Memory::alloc<double2>(nume);
    zx = Memory::alloc<double2>(numz);
    px0 = Memory::alloc<double2>(nump);
    pxp = Memory::alloc<double2>(nump);
    exp = Memory::alloc<double2>(nume);
    zxp = Memory::alloc<double2>(numz);
    sarea = Memory::alloc<double>(nums);
    svol = Memory::alloc<double>(nums);
    zarea = Memory::alloc<double>(numz);
    zvol = Memory::alloc<double>(numz);
    sareap = Memory::alloc<double>(nums);
    svolp = Memory::alloc<double>(nums);
    zareap = Memory::alloc<double>(numz);
    zvolp = Memory::alloc<double>(numz);
    zvol0 = Memory::alloc<double>(numz);
    ssurfp = Memory::alloc<double2>(nums);
    elen = Memory::alloc<double>(nume);
    zdl = Memory::alloc<double>(numz);
    smf = Memory::alloc<double>(nums);

    // do a few initial calculations
    //#pragma omp parallel for schedule(static)
    for (int pch = 0; pch < numpch; ++pch) {
        int pfirst = pchpfirst[pch];
        int plast = pchplast[pch];
        // copy nodepos into px, distributed across threads
        for (int p = pfirst; p < plast; ++p)
            px[p] = nodepos[p];

    }

    numsbad = 0;
    //#pragma omp parallel for schedule(static)
    for (int sch = 0; sch < numsch; ++sch) {
        int sfirst = schsfirst[sch];
        int slast = schslast[sch];
        calcCtrs(px, ex, zx, sfirst, slast);
        calcVols(px, zx, sarea, svol, zarea, zvol, sfirst, slast);
        calcSideFracs(sarea, zarea, smf, sfirst, slast);
    }
    checkBadSides();

    // Used to simplify task dependencies in sumOnProc
    sumOnProctaskindices = Memory::alloc<int*>(numpch);
    sumOnProctasksizes = Memory::alloc<int>(numpch);
    for (int pch = 0; pch < numpch; ++pch) {
        int pfirst = pchpfirst[pch];
        int plast = pchplast[pch];
        // Indices for only this chunk. Vector as mappc stopping point unknown
        std::vector<int> pchindices;

        for (int p = pfirst; p < plast; ++p) {
            for (int c = mappcfirst[p]; c >= 0; c = mapccnext[c]) {
                pchindices.push_back(c);
            }
        }

        // Copy data into array for OmpSs-2. Fails to compile with vector
        // Can't just store int pointer as data() is freed when vector goes out
        // of scope.
        sumOnProctasksizes[pch] = pchindices.size();
        sumOnProctaskindices[pch] = Memory::alloc<int>(pchindices.size());
        std::copy(pchindices.begin(), pchindices.end(), sumOnProctaskindices[pch]);
    }
}


void Mesh::initSides(
        const vector<int>& cellstart,
        const vector<int>& cellsize,
        const vector<int>& cellnodes) {

    mapsp1 = Memory::alloc<int>(nums);
    mapsp2 = Memory::alloc<int>(nums);
    mapsz  = Memory::alloc<int>(nums);
    mapss3 = Memory::alloc<int>(nums);
    mapss4 = Memory::alloc<int>(nums);

    for (int z = 0; z < numz; ++z) {
        int sbase = cellstart[z];
        int size = cellsize[z];
        for (int n = 0; n < size; ++n) {
            int s = sbase + n;
            int snext = sbase + (n + 1 == size ? 0 : n + 1);
            int slast = sbase + (n == 0 ? size : n) - 1;
            mapsz[s] = z;
            mapsp1[s] = cellnodes[s];
            mapsp2[s] = cellnodes[snext];
            mapss3[s] = slast;
            mapss4[s] = snext;
        } // for n
    } // for z

}


void Mesh::initEdges() {

    vector<vector<int> > edgepp(nump), edgepe(nump);

    mapse = Memory::alloc<int>(nums);

    int e = 0;
    for (int s = 0; s < nums; ++s) {
        int p1 = min(mapsp1[s], mapsp2[s]);
        int p2 = max(mapsp1[s], mapsp2[s]);

        vector<int>& vpp = edgepp[p1];
        vector<int>& vpe = edgepe[p1];
        int i = find(vpp.begin(), vpp.end(), p2) - vpp.begin();
        if (i == vpp.size()) {
            // (p, p2) isn't in the edge list - add it
            vpp.push_back(p2);
            vpe.push_back(e);
            ++e;
        }
        mapse[s] = vpe[i];
    }  // for s

    nume = e;

}


void Mesh::initChunks() {

    if (chunksize == 0) chunksize = max(nump, nums);

    // compute side chunks
    // use 'chunksize' for maximum chunksize; decrease as needed
    // to ensure that no zone has its sides split across chunk
    // boundaries
    int s1, s2 = 0;
    while (s2 < nums) {
        s1 = s2;
        s2 = min(s2 + chunksize, nums);
        while (s2 < nums && mapsz[s2] == mapsz[s2-1])
            --s2;
        schsfirst.push_back(s1);
        schslast.push_back(s2);
        schzfirst.push_back(mapsz[s1]);
        schzlast.push_back(mapsz[s2-1] + 1);
    }
    numsch = schsfirst.size();

    // compute point chunks
    int p1, p2 = 0;
    while (p2 < nump) {
        p1 = p2;
        p2 = min(p2 + chunksize, nump);
        pchpfirst.push_back(p1);
        pchplast.push_back(p2);
    }
    numpch = pchpfirst.size();

    // compute zone chunks
    int z1, z2 = 0;
    while (z2 < numz) {
        z1 = z2;
        z2 = min(z2 + chunksize, numz);
        zchzfirst.push_back(z1);
        zchzlast.push_back(z2);
    }
    numzch = zchzfirst.size();

}


void Mesh::initInvMap() {
    mappcfirst = Memory::alloc<int>(nump);
    mapccnext = Memory::alloc<int>(nums);

    vector<pair<int, int> > pcpair(nums);
    for (int c = 0; c < numc; ++c)
        pcpair[c] = make_pair(mapsp1[c], c);
    sort(pcpair.begin(), pcpair.end());
    for (int i = 0; i < numc; ++i) {
        int p = pcpair[i].first;
        int c = pcpair[i].second;

        if (i == 0) mappcfirst[p] = c;
        else {
            int pm = pcpair[i-1].first;
            if (p != pm) mappcfirst[p] = c;
        }
        mapccnext[c] = -1;
        if (i+1 != numc) {
            int pp = pcpair[i+1].first;
            int cp = pcpair[i+1].second;
            if (p == pp) mapccnext[c] = cp;
        }
    }
}


void Mesh::initParallel(
        const vector<int>& slavemstrpes,
        const vector<int>& slavemstrcounts,
        const vector<int>& slavepoints,
        const vector<int>& masterslvpes,
        const vector<int>& masterslvcounts,
        const vector<int>& masterpoints) {
    if (Parallel::numpe == 1) return;

    nummstrpe = slavemstrpes.size();
    mapmstrpepe = Memory::alloc<int>(nummstrpe);
    copy(slavemstrpes.begin(), slavemstrpes.end(), mapmstrpepe);
    mstrpenumslv = Memory::alloc<int>(nummstrpe);
    copy(slavemstrcounts.begin(), slavemstrcounts.end(), mstrpenumslv);
    mapmstrpeslv1 = Memory::alloc<int>(nummstrpe);
    int count = 0;
    for (int mstrpe = 0; mstrpe < nummstrpe; ++mstrpe) {
        mapmstrpeslv1[mstrpe] = count;
        count += mstrpenumslv[mstrpe];
    }
    numslv = slavepoints.size();
    mapslvp = Memory::alloc<int>(numslv);
    copy(slavepoints.begin(), slavepoints.end(), mapslvp);

    numslvpe = masterslvpes.size();
    mapslvpepe = Memory::alloc<int>(numslvpe);
    copy(masterslvpes.begin(), masterslvpes.end(), mapslvpepe);
    slvpenumprx = Memory::alloc<int>(numslvpe);
    copy(masterslvcounts.begin(), masterslvcounts.end(), slvpenumprx);
    mapslvpeprx1 = Memory::alloc<int>(numslvpe);
    count = 0;
    for (int slvpe = 0; slvpe < numslvpe; ++slvpe) {
        mapslvpeprx1[slvpe] = count;
        count += slvpenumprx[slvpe];
    }
    numprx = masterpoints.size();
    mapprxp = Memory::alloc<int>(numprx);
    copy(masterpoints.begin(), masterpoints.end(), mapprxp);

    // Send buffer allocation moved here to avoid scoping issues with tasking
    // Communication functions are generics but limited to either double or
    // double2 types in practice.
    // sumToPoints() called twice in main loop but first call is always T<double>
    // and second is always T<double2> => No concerns of buffer reuse before ready
    prxvar_double = Memory::alloc<double>(numprx);
    prxvar_double2 = Memory::alloc<double2>(numprx);
    slvvar_gather_double = Memory::alloc<double>(numslv);
    slvvar_gather_double2 = Memory::alloc<double2>(numslv);
    slvvar_scatter_double = Memory::alloc<double>(numslv);
    slvvar_scatter_double2 = Memory::alloc<double2>(numslv);
}


void Mesh::writeStats() {

    int64_t gnump = nump;
    // make sure that boundary points aren't double-counted;
    // only count them if they are masters
    if (Parallel::numpe > 1) gnump -= numslv;
    int64_t gnumz = numz;
    int64_t gnums = nums;
    int64_t gnume = nume;
    int gnumpch = numpch;
    int gnumzch = numzch;
    int gnumsch = numsch;

    Parallel::globalSum(gnump);
    Parallel::globalSum(gnumz);
    Parallel::globalSum(gnums);
    Parallel::globalSum(gnume);
    Parallel::globalSum(gnumpch);
    Parallel::globalSum(gnumzch);
    Parallel::globalSum(gnumsch);

    if (Parallel::mype > 0) return;

    cout << "--- Mesh Information ---" << endl;
    cout << "Points:  " << gnump << endl;
    cout << "Zones:  "  << gnumz << endl;
    cout << "Sides:  "  << gnums << endl;
    cout << "Edges:  "  << gnume << endl;
    cout << "Side chunks:  " << gnumsch << endl;
    cout << "Point chunks:  " << gnumpch << endl;
    cout << "Zone chunks:  " << gnumzch << endl;
    cout << "Chunk size:  " << chunksize << endl;
    cout << "------------------------" << endl;

}


void Mesh::write(
        const string& probname,
        const int cycle,
        const double time,
        const double* zr,
        const double* ze,
        const double* zp) {

    if (writexy) {
        if (Parallel::mype == 0)
            cout << "Writing .xy file..." << endl;
        wxy->write(probname, zr, ze, zp);
    }
    if (writegold) {
        if (Parallel::mype == 0)
            cout << "Writing gold file..." << endl;
        egold->write(probname, cycle, time, zr, ze, zp);
    }

}


vector<int> Mesh::getXPlane(const double c) {

    vector<int> mapbp;
    const double eps = 1.e-12;

    for (int p = 0; p < nump; ++p) {
        if (fabs(px[p].x - c) < eps) {
            mapbp.push_back(p);
        }
    }
    return mapbp;

}


vector<int> Mesh::getYPlane(const double c) {

    vector<int> mapbp;
    const double eps = 1.e-12;

    for (int p = 0; p < nump; ++p) {
        if (fabs(px[p].y - c) < eps) {
            mapbp.push_back(p);
        }
    }
    return mapbp;

}


void Mesh::getPlaneChunks(
        const int numb,
        const int* mapbp,
        vector<int>& pchbfirst,
        vector<int>& pchblast) {

    pchbfirst.resize(0);
    pchblast.resize(0);

    // compute boundary point chunks
    // (boundary points contained in each point chunk)
    int bf, bl = 0;
    for (int pch = 0; pch < numpch; ++pch) {
         int pl = pchplast[pch];
         bf = bl;
         bl = lower_bound(&mapbp[bf], &mapbp[numb], pl) - &mapbp[0];
         pchbfirst.push_back(bf);
         pchblast.push_back(bl);
    }

}


void Mesh::calcCtrs(
        const double2* px,
        double2* ex,
        double2* zx,
        const int sfirst,
        const int slast) {

    int zfirst = mapsz[sfirst];
    int zlast = (slast < nums ? mapsz[slast] : numz);
    fill(&zx[zfirst], &zx[zlast], double2(0., 0.));

    for (int s = sfirst; s < slast; ++s) {
        int p1 = mapsp1[s];
        int p2 = mapsp2[s];
        int e = mapse[s];
        int z = mapsz[s];
        ex[e] = 0.5 * (px[p1] + px[p2]);
        zx[z] += px[p1];
    }

    for (int z = zfirst; z < zlast; ++z) {
        zx[z] /= (double) znump[z];
    }

}


void Mesh::calcVols(
        const double2* px,
        const double2* zx,
        double* sarea,
        double* svol,
        double* zarea,
        double* zvol,
        const int sfirst,
        const int slast) {

    int zfirst = mapsz[sfirst];
    int zlast = (slast < nums ? mapsz[slast] : numz);
    fill(&zvol[zfirst], &zvol[zlast], 0.);
    fill(&zarea[zfirst], &zarea[zlast], 0.);

    const double third = 1. / 3.;
    int count = 0;
    for (int s = sfirst; s < slast; ++s) {
        int p1 = mapsp1[s];
        int p2 = mapsp2[s];
        int z = mapsz[s];

        // compute side volumes, sum to zone
        double sa = 0.5 * cross(px[p2] - px[p1], zx[z] - px[p1]);
        double sv = third * sa * (px[p1].x + px[p2].x + zx[z].x);
        sarea[s] = sa;
        svol[s] = sv;
        zarea[z] += sa;
        zvol[z] += sv;

        // check for negative side volumes
        if (sv <= 0.) count += 1;

    } // for s

    if (count > 0) {
        //#pragma omp atomic
        //numsbad += count;
        // Replaces checkBadSides() to avoid need for atomic
        cerr << "Error: negative side volume" << endl;
        cerr << "Exiting..." << endl;
        exit(1);
    }

}


void Mesh::checkBadSides() {

    // if there were negative side volumes, error exit
    if (numsbad > 0) {
        cerr << "Error: " << numsbad << " negative side volumes" << endl;
        cerr << "Exiting..." << endl;
        exit(1);
    }

}


void Mesh::calcSideFracs(
        const double* sarea,
        const double* zarea,
        double* smf,
        const int sfirst,
        const int slast) {

    #pragma ivdep
    for (int s = sfirst; s < slast; ++s) {
        int z = mapsz[s];
        smf[s] = sarea[s] / zarea[z];
    }
}


void Mesh::calcSurfVecs(
        const double2* zx,
        const double2* ex,
        double2* ssurf,
        const int sfirst,
        const int slast) {

    #pragma ivdep
    for (int s = sfirst; s < slast; ++s) {
        int z = mapsz[s];
        int e = mapse[s];

        ssurf[s] = rotateCCW(ex[e] - zx[z]);

    }

}


void Mesh::calcEdgeLen(
        const double2* px,
        double* elen,
        const int sfirst,
        const int slast) {

    for (int s = sfirst; s < slast; ++s) {
        const int p1 = mapsp1[s];
        const int p2 = mapsp2[s];
        const int e = mapse[s];

        elen[e] = length(px[p2] - px[p1]);

    }
}


void Mesh::calcCharLen(
        const double* sarea,
        double* zdl,
        const int sfirst,
        const int slast) {

    int zfirst = mapsz[sfirst];
    int zlast = (slast < nums ? mapsz[slast] : numz);
    fill(&zdl[zfirst], &zdl[zlast], 1.e99);

    for (int s = sfirst; s < slast; ++s) {
        int z = mapsz[s];
        int e = mapse[s];

        double area = sarea[s];
        double base = elen[e];
        double fac = (znump[z] == 3 ? 3. : 4.);
        double sdl = fac * area / base;
        zdl[z] = min(zdl[z], sdl);
    }
}

double* Mesh::giveGatherBuf(double* array) { return slvvar_gather_double; }
double2* Mesh::giveGatherBuf(double2* array) { return slvvar_gather_double2; }

template <typename T>
void Mesh::parallelGather(
        const T* pvar,
        T* prxvar,
        const int tagmod) {
#ifdef USE_MPI
    // This routine gathers slave values for which MYPE owns the masters.
    const int tagmpi = 100 + tagmod;
    const int type_size = sizeof(T);

    // TODO: This entirely defeats the point of using generics
    T* slvvar = giveGatherBuf(prxvar);

    // Post receives for incoming messages from slaves.
    // Store results in proxy buffer.
    for (int slvpe = 0; slvpe < numslvpe; ++slvpe) {
        int pe = mapslvpepe[slvpe];
        int nprx = slvpenumprx[slvpe];
        int prx1 = mapslvpeprx1[slvpe];

        #pragma oss task label("Parallel Gather: Post Irecvs") \
                         out(prxvar[prx1:nprx*type_size - 1]) \
                         firstprivate(prx1, nprx, type_size, pe, tagmpi)
        {
            MPI_Request request;
            MPI_Irecv(&prxvar[prx1], nprx * type_size, MPI_BYTE,
                    pe, tagmpi, MPI_COMM_WORLD, &request);
            TAMPI_Iwait(&request, MPI_STATUS_IGNORE);
        }
    }

    // Load slave data buffer from points.
    // TODO Merge with send task?
    #pragma oss task label("Parallel Gather: Pack send buffers") \
                     in({pvar[mapslvp[p]], p=0;numslv}) \
                     out(slvvar[0:numslv-1])
    {
        for (int slv = 0; slv < numslv; ++slv) {
            int p = mapslvp[slv];
            slvvar[slv] = pvar[p];
        }
    }

    // Send slave data to master PEs.
    for (int mstrpe = 0; mstrpe < nummstrpe; ++mstrpe) {
        int pe = mapmstrpepe[mstrpe];
        int nslv = mstrpenumslv[mstrpe];
        int slv1 = mapmstrpeslv1[mstrpe];

        #pragma oss task label("Parallel Gather: Sends") \
                         in(slvvar[slv1:nslv-1]) \
                         firstprivate(slv1, nslv, type_size, pe, tagmpi)
        {
            MPI_Send(&slvvar[slv1], nslv * type_size, MPI_BYTE,
                    pe, tagmpi, MPI_COMM_WORLD);
        }

    }
#endif
}


template <typename T>
void Mesh::parallelSum(
        T* pvar,
        T* prxvar) {
#ifdef USE_MPI
    #pragma oss task label("Parallel sum") \
                     inout(prxvar[0:numprx-1], {pvar[mapprxp[p]], p=0;numprx})
    {
        // Compute sum of all (proxy/master) sets.
        // Store results in master.
        for (int prx = 0; prx < numprx; ++prx) {
            int p = mapprxp[prx];
            pvar[p] += prxvar[prx];
        }

        // Copy updated master data back to proxies.
        for (int prx = 0; prx < numprx; ++prx) {
            int p = mapprxp[prx];
            prxvar[prx] = pvar[p];
        }
    }
#endif
}

double* Mesh::giveScatterBuf(double* array) { return slvvar_scatter_double; }
double2* Mesh::giveScatterBuf(double2* array) { return slvvar_scatter_double2; }

template <typename T>
void Mesh::parallelScatter(
        T* pvar,
        const T* prxvar,
        const int tagmod) {
#ifdef USE_MPI
    // This routine scatters master values on MYPE to all slave copies
    // owned by other PEs.
    const int tagmpi = 200 + tagmod;
    const int type_size = sizeof(T);

    // TODO: This entirely defeats the point of using generics
    T* slvvar = giveScatterBuf(pvar);

    // Post receives for incoming messages from masters.
    // Store results in slave buffer.
    for (int mstrpe = 0; mstrpe < nummstrpe; ++mstrpe) {
        int pe = mapmstrpepe[mstrpe];
        int nslv = mstrpenumslv[mstrpe];
        int slv1 = mapmstrpeslv1[mstrpe];

        #pragma oss task label("Parallel Scatter: Post Irecvs") \
                         out(slvvar[slv1:nslv-1]) \
                         firstprivate(slv1, nslv, type_size, pe, tagmpi)
        {
            MPI_Request request;
            MPI_Irecv(&slvvar[slv1], nslv * type_size, MPI_BYTE,
                    pe, tagmpi, MPI_COMM_WORLD, &request);
            TAMPI_Iwait(&request, MPI_STATUS_IGNORE);
        }
    }

    // Send updated slave data from proxy buffer back to slave PEs.
    for (int slvpe = 0; slvpe < numslvpe; ++slvpe) {
        int pe = mapslvpepe[slvpe];
        int nprx = slvpenumprx[slvpe];
        int prx1 = mapslvpeprx1[slvpe];

        #pragma oss task label("Parallel Scatter: Sends") \
                         in(prxvar[prx1:nprx-1]) \
                         firstprivate(prx1, nprx, type_size, pe, tagmpi)
        {
            MPI_Send((void*)&prxvar[prx1], nprx * type_size, MPI_BYTE,
                    pe, tagmpi, MPI_COMM_WORLD);
        }
    }

    // Store slave data from buffer back to points.
    // TODO Merge with recv task?
    #pragma oss task label("Parallel Scatter: Unpack recv buffers") \
                     in(slvvar[0:numslv-1]) \
                     out({pvar[mapslvp[p]], p=0;numslv})
    {
        for (int slv = 0; slv < numslv; ++slv) {
            int p = mapslvp[slv];
            pvar[p] = slvvar[slv];
        }
    }
#endif
}

double* Mesh::giveprxvar(double* arrayt) { return prxvar_double; }
double2* Mesh::giveprxvar(double2* arrayt) { return prxvar_double2; }

template <typename T>
void Mesh::sumAcrossProcs(T* pvar, const int tagmod) {
    if (Parallel::numpe == 1) return;

    T* prxvar = giveprxvar(pvar);
    parallelGather(pvar, &prxvar[0], tagmod);
    parallelSum(pvar, &prxvar[0]);
    parallelScatter(pvar, &prxvar[0], tagmod);
}


template <typename T>
void Mesh::sumOnProc(
        const T* cvar,
        T* pvar) {

    //#pragma omp parallel for schedule(static)
    for (int pch = 0; pch < numpch; ++pch) {
        int pfirst = pchpfirst[pch];
        int plast = pchplast[pch];

        #pragma oss task label("Communicate: sum on process") \
                         in({cvar[sumOnProctaskindices[pch][i]], i=0;sumOnProctasksizes[pch]}) \
                         out(pvar[pfirst:plast-1]) \
                         firstprivate(pfirst, plast)
        {
            for (int p = pfirst; p < plast; ++p) {
                T x = T();
                for (int c = mappcfirst[p]; c >= 0; c = mapccnext[c]) {
                    x += cvar[c];
                }
                pvar[p] = x;
            }  // for p
        } // task
    }  // for pch
}


template <>
void Mesh::sumToPoints(
        const double* cvar,
        double* pvar,
        int tagmod) {

    #pragma oss taskwait
    sumOnProc(cvar, pvar);
    #pragma oss taskwait
    if (Parallel::numpe > 1)
        sumAcrossProcs(pvar, tagmod);
    #pragma oss taskwait

}


template <>
void Mesh::sumToPoints(
        const double2* cvar,
        double2* pvar,
        int tagmod) {

    #pragma oss taskwait
    sumOnProc(cvar, pvar);
    #pragma oss taskwait
    if (Parallel::numpe > 1)
        sumAcrossProcs(pvar, tagmod);
    #pragma oss taskwait

}

