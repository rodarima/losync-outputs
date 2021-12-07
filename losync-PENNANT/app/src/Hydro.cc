/*
 * Hydro.cc
 *
 *  Created on: Dec 22, 2011
 *      Author: cferenba
 *
 * Copyright (c) 2012, Triad National Security, LLC.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style open-source
 * license; see top-level LICENSE file for full license text.
 */

#include "Hydro.hh"

#include <string>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <iostream>
#include <iomanip>

#include "Parallel.hh"
#include "Memory.hh"
#include "InputFile.hh"
#include "Mesh.hh"
#include "PolyGas.hh"
#include "TTS.hh"
#include "QCS.hh"
#include "HydroBC.hh"

using namespace std;


Hydro::Hydro(const InputFile* inp, Mesh* m) : mesh(m) {
    cfl = inp->getDouble("cfl", 0.6);
    cflv = inp->getDouble("cflv", 0.1);
    rinit = inp->getDouble("rinit", 1.);
    einit = inp->getDouble("einit", 0.);
    rinitsub = inp->getDouble("rinitsub", 1.);
    einitsub = inp->getDouble("einitsub", 0.);
    uinitradial = inp->getDouble("uinitradial", 0.);
    bcx = inp->getDoubleList("bcx", vector<double>());
    bcy = inp->getDoubleList("bcy", vector<double>());

    pgas = new PolyGas(inp, this);
    tts = new TTS(inp, this);
    qcs = new QCS(inp, this);

    const double2 vfixx = double2(1., 0.);
    const double2 vfixy = double2(0., 1.);
    for (int i = 0; i < bcx.size(); ++i)
        bcs.push_back(new HydroBC(mesh, vfixx, mesh->getXPlane(bcx[i])));
    for (int i = 0; i < bcy.size(); ++i)
        bcs.push_back(new HydroBC(mesh, vfixy, mesh->getYPlane(bcy[i])));

    // Array of timestep sizes for all cycles
    int cstop = inp->getInt("cstop", 999999);
    dt = Memory::alloc<double>(cstop);

    init();
}


Hydro::~Hydro() {

    delete tts;
    delete qcs;
    for (int i = 0; i < bcs.size(); ++i) {
        delete bcs[i];
    }
}


void Hydro::init() {

    const int numpch = mesh->numpch;
    const int numzch = mesh->numzch;
    const int nump = mesh->nump;
    const int numz = mesh->numz;
    const int nums = mesh->nums;

    const double2* zx = mesh->zx;
    const double* zvol = mesh->zvol;

    // allocate arrays
    pu = Memory::alloc<double2>(nump);
    pu0 = Memory::alloc<double2>(nump);
    pap = Memory::alloc<double2>(nump);
    pf = Memory::alloc<double2>(nump);
    pmaswt = Memory::alloc<double>(nump);
    cmaswt = Memory::alloc<double>(nums);
    zm = Memory::alloc<double>(numz);
    zr = Memory::alloc<double>(numz);
    zrp = Memory::alloc<double>(numz);
    ze = Memory::alloc<double>(numz);
    zetot = Memory::alloc<double>(numz);
    zw = Memory::alloc<double>(numz);
    zwrate = Memory::alloc<double>(numz);
    zp = Memory::alloc<double>(numz);
    zss = Memory::alloc<double>(numz);
    zdu = Memory::alloc<double>(numz);
    sfp = Memory::alloc<double2>(nums);
    sfq = Memory::alloc<double2>(nums);
    sft = Memory::alloc<double2>(nums);
    cftot = Memory::alloc<double2>(nums);
    dtrec = Memory::alloc<double>(numzch);
    msgdtrec = Memory::alloc<char*>(numzch);

    // initialize hydro vars
    //#pragma omp parallel for schedule(static)
    for (int zch = 0; zch < numzch; ++zch) {
        int zfirst = mesh->zchzfirst[zch];
        int zlast = mesh->zchzlast[zch];

        fill(&zr[zfirst], &zr[zlast], rinit);
        fill(&ze[zfirst], &ze[zlast], einit);
        fill(&zwrate[zfirst], &zwrate[zlast], 0.);

        const vector<double>& subrgn = mesh->subregion;
        if (!subrgn.empty()) {
            const double eps = 1.e-12;
            #pragma ivdep
            for (int z = zfirst; z < zlast; ++z) {
                if (zx[z].x > (subrgn[0] - eps) &&
                    zx[z].x < (subrgn[1] + eps) &&
                    zx[z].y > (subrgn[2] - eps) &&
                    zx[z].y < (subrgn[3] + eps)) {
                    zr[z] = rinitsub;
                    ze[z] = einitsub;
                }
            }
        }

        #pragma ivdep
        for (int z = zfirst; z < zlast; ++z) {
            zm[z] = zr[z] * zvol[z];
            zetot[z] = ze[z] * zm[z];
        }
    }  // for sch

    //#pragma omp parallel for schedule(static)
    for (int pch = 0; pch < numpch; ++pch) {
        int pfirst = mesh->pchpfirst[pch];
        int plast = mesh->pchplast[pch];
        if (uinitradial != 0.)
            initRadialVel(uinitradial, pfirst, plast);
        else
            fill(&pu[pfirst], &pu[plast], double2(0., 0.));
    }  // for pch

    // Replace single dt hydro reset with initialisation of new dt array
    //resetDtHydro();
    for (int zch = 0; zch < numzch; ++zch) {
        dtrec[zch] = 1.e99;
        msgdtrec[zch] = Memory::alloc<char>(80);
        strcpy(msgdtrec[zch], "Hydro default");
    }

    // Used to simplify task dependencies in step 4a (boundary conditions)
    // All relevant data structures have been filled by this point.
    //   pchbfirst and pchblast set up by getPlaneChunks in HydroBC constructor
    // Does a run through of the nested loops in step 4a but saves all indices
    // in a flat 2d array to avoid the need for mapping arrays in task pragma.
    bcstaskindices = Memory::alloc<int*>(numpch);
    bcstasksizes = Memory::alloc<int>(numpch);
    for (int pch = 0; pch < numpch; ++pch) {
        // Indices for only this chunk. Vector as bfirst-blast span unknown
        std::vector<int> pchindices;

        for (int i = 0; i < bcs.size(); ++i) {
            int bfirst = bcs[i]->pchbfirst[pch];
            int blast = bcs[i]->pchblast[pch];
            for (int b = bfirst; b < blast; ++b) {
                pchindices.push_back(bcs[i]->mapbp[b]);
            }
        }

        // Copy data into array for OmpSs-2. Fails to compile with vector
        // Can't just store int pointer as data() is freed when vector goes out
        // of scope.
        bcstasksizes[pch] = pchindices.size();
        bcstaskindices[pch] = Memory::alloc<int>(pchindices.size());
        std::copy(pchindices.begin(), pchindices.end(), bcstaskindices[pch]);
    }
}


void Hydro::initRadialVel(
        const double vel,
        const int pfirst,
        const int plast) {
    const double2* px = mesh->px;
    const double eps = 1.e-12;

    #pragma ivdep
    for (int p = pfirst; p < plast; ++p) {
        double pmag = length(px[p]);
        if (pmag > eps)
            pu[p] = vel * px[p] / pmag;
        else
            pu[p] = double2(0., 0.);
    }
}


void Hydro::doCycle(const int cycle) {

    const int numpch = mesh->numpch;
    const int numsch = mesh->numsch;
    double2* px = mesh->px;
    double2* ex = mesh->ex;
    double2* zx = mesh->zx;
    double* sarea = mesh->sarea;
    double* svol = mesh->svol;
    double* zarea = mesh->zarea;
    double* zvol = mesh->zvol;
    double* sareap = mesh->sareap;
    double* svolp = mesh->svolp;
    double* zareap = mesh->zareap;
    double* zvolp = mesh->zvolp;
    double* zvol0 = mesh->zvol0;
    double2* ssurfp = mesh->ssurfp;
    double* elen = mesh->elen;
    double2* px0 = mesh->px0;
    double2* pxp = mesh->pxp;
    double2* exp = mesh->exp;
    double2* zxp = mesh->zxp;
    double* smf = mesh->smf;
    double* zdl = mesh->zdl;

    // Begin hydro cycle
    for (int pch = 0; pch < numpch; ++pch) {
        int pfirst = mesh->pchpfirst[pch];
        int plast = mesh->pchplast[pch];

        // Dependency blocks order:
        //   copies
        //   advPosHalf
        //
        // advPosHalf in(px0[pfirst:plast-1], pu0[pfirst:plast-1]) not included
        // in pragma as is output of previous step in task.
        #pragma oss task label("Step 1: advance mesh") \
                         in(px[pfirst:plast-1], pu[pfirst:plast-1]) \
                         out(px0[pfirst:plast-1], pu0[pfirst:plast-1]) \
                         \
                         in (dt[cycle-1]) \
                         out(pxp[pfirst:plast-1]) \
                         \
                         firstprivate(pfirst, plast, cycle)
        {
          // save off point variable values from previous cycle
          copy(&px[pfirst], &px[plast], &px0[pfirst]);
          copy(&pu[pfirst], &pu[plast], &pu0[pfirst]);

          // ===== Predictor step =====
          // 1. advance mesh to center of time step
          advPosHalf(px0, pu0, dt[cycle-1], pxp, pfirst, plast);
        }
    } // for pch

    for (int sch = 0; sch < numsch; ++sch) {
        int sfirst = mesh->schsfirst[sch];
        int slast = mesh->schslast[sch];
        int zfirst = mesh->schzfirst[sch];
        int zlast = mesh->schzlast[sch];

        // Dependency blocks order:
        //   copy
        //   calcCtrls
        //   calcVols
        //   calcSurfVecs
        //   calcEdgeLen
        //   calcCharLen
        //   calcRho
        //   calcCrnrMass
        //   calcStateAtHalf
        //   pgas->calcForce
        //   tts->calcForce
        //   qcs->calcForce
        //   sumCrnrForce
        //
        // "zlast = (slast < nums ? mapsz[slast] : numz);" in calcCtrls, calcVols and calcCharLen
        // Two separate types of task needed
        //
        // "zlast = (slast < nums ? mesh->mapsz[slast] : numz);" in qcs->calcforce not relevant.
        // Only applies to z0uc temporary array in setCornerDiv() which is freed at end of function.
        if (slast < mesh->nums) {
            #pragma oss task label("Steps 1a-4 standard: mesh & masses & material & forces") \
                             in(zvol[zfirst:zlast-1]) \
                             out(zvol0[zfirst:zlast-1]) \
                             \
                             in({pxp[mesh->mapsp1[p1]], p1=sfirst:slast-1}, {pxp[mesh->mapsp2[p2]], p2=sfirst:slast-1}) \
                             out({zxp[mesh->mapsz[z]], z=sfirst:slast-1}, {exp[mesh->mapse[e]], e=sfirst:slast-1}) \
                             \
                             in({pxp[mesh->mapsp1[p1]], p1=sfirst:slast-1}, {pxp[mesh->mapsp2[p2]], p2=sfirst:slast-1}, \
                                {zxp[mesh->mapsz[z]], z=sfirst:slast-1}) \
                             out({zvolp[mesh->mapsz[z]], z=sfirst:slast-1}, {zareap[mesh->mapsz[z]], z=sfirst:slast-1}, \
                                 svolp[sfirst:slast-1], sareap[sfirst:slast-1]) \
                             \
                             in({zxp[mesh->mapsz[z]], z=sfirst:slast-1}, {exp[mesh->mapse[e]], e=sfirst:slast-1}) \
                             out(ssurfp[sfirst:slast-1]) \
                             \
                             in({pxp[mesh->mapsp1[p1]], p1=sfirst:slast-1}, {pxp[mesh->mapsp2[p2]], p2=sfirst:slast-1}) \
                             out({elen[mesh->mapse[e]], e=sfirst:slast-1}) \
                             \
                             in(sareap[sfirst:slast-1], {elen[mesh->mapse[e]], e=sfirst:slast-1}, \
                                {mesh->znump[mesh->mapsz[z]], z=sfirst:slast-1}) \
                             out({zdl[mesh->mapsz[z]], z=sfirst:slast-1}) \
                             \
                             in(zm[zfirst:zlast-1], zvolp[zfirst:zlast-1]) \
                             out(zrp[zfirst:zlast-1]) \
                             \
                             in({zrp[mesh->mapsz[z]], z=sfirst:slast-1}, {zareap[mesh->mapsz[z]], z=sfirst:slast-1}, \
                                 smf[sfirst:slast-1], {smf[mesh->mapss3[s3]], s3=sfirst:slast-1}) \
                             out(cmaswt[sfirst:slast-1]) \
                             \
                             in(dt[cycle-1], zr[zfirst:zlast-1], ze[zfirst:zlast-1], zm[zfirst:zlast-1], \
                                zvolp[zfirst:zlast-1], zvol0[zfirst:zlast-1], zss[zfirst:zlast-1], \
                                zwrate[zfirst:zlast-1]) \
                             out(zp[zfirst:zlast-1], zss[zfirst:zlast-1]) \
                             \
                             in({zp[mesh->mapsz[s]], s=sfirst:slast-1}, ssurfp[sfirst:slast-1]) \
                             out(sfp[sfirst:slast-1]) \
                             \
                             in({zareap[mesh->mapsz[s]], s=sfirst:slast-1}, {zrp[mesh->mapsz[s]], s=sfirst:slast-1}, \
                                {zss[mesh->mapsz[s]], s=sfirst:slast-1}, \
                                sareap[sfirst:slast-1], smf[sfirst:slast-1], ssurfp[sfirst:slast-1]) \
                             out(sft[sfirst:slast-1]) \
                             \
                             in({mesh->znump[mesh->mapsz[z]], z=sfirst:slast-1}, \
                                {pu[mesh->mapsp2[mesh->mapss3[p]]], p=sfirst:slast-1}, \
                                {pxp[mesh->mapsp2[mesh->mapss3[p]]], p=sfirst:slast-1}, \
                                {pu[mesh->mapsp2[p2]], p2=sfirst:slast-1}, \
                                {exp[mesh->mapse[e2]], e2=sfirst:slast-1}, \
                                {zxp[mesh->mapsz[mesh->mapss3[z]]], z=sfirst:slast-1}, \
                                {pu[mesh->mapsp1[mesh->mapss3[p1]]], p1=sfirst:slast-1}, \
                                {exp[mesh->mapse[mesh->mapss3[e1]]], e1=sfirst:slast-1}, \
                                {elen[mesh->mapse[mesh->mapss3[e1]]], e1=sfirst:slast-1}, \
                                {elen[mesh->mapse[e2]], e2=sfirst:slast-1}, \
                                {zss[mesh->mapsz[z]], z=sfirst:slast-1}, \
                                {zrp[mesh->mapsz[z]], z=sfirst:slast-1}, \
                                {pu[mesh->mapsp2[mesh->mapss3[p]]], p=sfirst:slast-1}, \
                                {pu[mesh->mapsp1[mesh->mapss3[p1]]], p1=sfirst:slast-1}, \
                                {pu[mesh->mapsp2[p2]], p2=sfirst:slast-1}, \
                                {elen[mesh->mapse[mesh->mapss3[e1]]], e1=sfirst:slast-1}, \
                                {elen[mesh->mapse[e2]], e2=sfirst:slast-1}, \
                                {pu[mesh->mapsp1[p1]], p1=sfirst:slast-1}, \
                                {pu[mesh->mapsp2[p2]], p2=sfirst:slast-1}, \
                                {pxp[mesh->mapsp1[p1]], p1=sfirst:slast-1}, \
                                {pxp[mesh->mapsp2[p2]], p2=sfirst:slast-1}, \
                                {elen[mesh->mapse[e]], e=sfirst:slast-1}, \
                                {zss[mesh->mapsz[z]], z=sfirst:slast-1}) \
                             out(sfq[sfirst:slast-1], {zdu[mesh->mapsz[z]], z=sfirst:slast-1}) \
                             \
                             in(sfp[sfirst:slast-1], sfq[sfirst:slast-1], sft[sfirst:slast-1], \
                                {sfp[mesh->mapss3[s3]], s3=sfirst:slast-1}, {sfq[mesh->mapss3[s3]], s3=sfirst:slast-1}, \
                                {sft[mesh->mapss3[s3]], s3=sfirst:slast-1}) \
                             \
                             firstprivate(sfirst, slast, zfirst, zlast, cycle)
            {
                // save off zone variable values from previous cycle
                copy(&zvol[zfirst], &zvol[zlast], &zvol0[zfirst]);

                // 1a. compute new mesh geometry
                mesh->calcCtrs(pxp, exp, zxp, sfirst, slast);
                mesh->calcVols(pxp, zxp, sareap, svolp, zareap, zvolp,
                        sfirst, slast);
                mesh->calcSurfVecs(zxp, exp, ssurfp, sfirst, slast);
                mesh->calcEdgeLen(pxp, elen, sfirst, slast);
                mesh->calcCharLen(sareap, zdl, sfirst, slast);

                // 2. compute point masses
                calcRho(zm, zvolp, zrp, zfirst, zlast);
                calcCrnrMass(zrp, zareap, smf, cmaswt, sfirst, slast);

                // 3. compute material state (half-advanced)
                pgas->calcStateAtHalf(zr, zvolp, zvol0, ze, zwrate, zm, dt[cycle-1],
                        zp, zss, zfirst, zlast);

                // 4. compute forces
                pgas->calcForce(zp, ssurfp, sfp, sfirst, slast);
                tts->calcForce(zareap, zrp, zss, sareap, smf, ssurfp, sft,
                        sfirst, slast);
                qcs->calcForce(sfq, sfirst, slast);
                sumCrnrForce(sfp, sfq, sft, cftot, sfirst, slast);
            }
        } else {
            #pragma oss task label("Steps 1a-4 numz: mesh & masses & material & forces") \
                             in(zvol[zfirst:zlast-1]) \
                             out(zvol0[zfirst:zlast-1]) \
                             \
                             in({pxp[mesh->mapsp1[p1]], p1=sfirst:slast-1}, {pxp[mesh->mapsp2[p2]], p2=sfirst:slast-1}) \
                             out(zxp[mesh->mapsz[sfirst]:mesh->numz-1], {zxp[mesh->mapsz[z]], z=sfirst:slast-1}, \
                                 {exp[mesh->mapse[e]], e=sfirst:slast-1}) \
                             \
                             in({pxp[mesh->mapsp1[p1]], p1=sfirst:slast-1}, {pxp[mesh->mapsp2[p2]], p2=sfirst:slast-1}, \
                                {zxp[mesh->mapsz[z]], z=sfirst:slast-1}) \
                             out(zvolp[mesh->mapsz[sfirst]:mesh->numz-1], {zvolp[mesh->mapsz[z]], z=sfirst:slast-1}, \
                                 zareap[mesh->mapsz[sfirst]:mesh->numz-1], {zareap[mesh->mapsz[z]], z=sfirst:slast-1}, \
                                 svolp[sfirst:slast-1], sareap[sfirst:slast-1]) \
                             \
                             in({zxp[mesh->mapsz[z]], z=sfirst:slast-1}, {exp[mesh->mapse[e]], e=sfirst:slast-1}) \
                             out(ssurfp[sfirst:slast-1]) \
                             \
                             in({pxp[mesh->mapsp1[p1]], p1=sfirst:slast-1}, {pxp[mesh->mapsp2[p2]], p2=sfirst:slast-1}) \
                             out({elen[mesh->mapse[e]], e=sfirst:slast-1}) \
                             \
                             in(sareap[sfirst:slast-1], {elen[mesh->mapse[e]], e=sfirst:slast-1}, \
                                {mesh->znump[mesh->mapsz[z]], z=sfirst:slast-1}) \
                             out(zdl[mesh->mapsz[sfirst]:mesh->numz-1], {zdl[mesh->mapsz[z]], z=sfirst:slast-1}) \
                             \
                             in(zm[zfirst:zlast-1], zvolp[zfirst:zlast-1]) \
                             out(zrp[zfirst:zlast-1]) \
                             \
                             in({zrp[mesh->mapsz[z]], z=sfirst:slast-1}, {zareap[mesh->mapsz[z]], z=sfirst:slast-1}, \
                                 smf[sfirst:slast-1], {smf[mesh->mapss3[s3]], s3=sfirst:slast-1}) \
                             out(cmaswt[sfirst:slast-1]) \
                             \
                             in(dt[cycle-1], zr[zfirst:zlast-1], ze[zfirst:zlast-1], zm[zfirst:zlast-1], \
                                zvolp[zfirst:zlast-1], zvol0[zfirst:zlast-1], zss[zfirst:zlast-1], \
                                zwrate[zfirst:zlast-1]) \
                             out(zp[zfirst:zlast-1], zss[zfirst:zlast-1]) \
                             \
                             in({zp[mesh->mapsz[s]], s=sfirst:slast-1}, ssurfp[sfirst:slast-1]) \
                             out(sfp[sfirst:slast-1]) \
                             \
                             in({zareap[mesh->mapsz[s]], s=sfirst:slast-1}, {zrp[mesh->mapsz[s]], s=sfirst:slast-1}, \
                                {zss[mesh->mapsz[s]], s=sfirst:slast-1}, \
                                sareap[sfirst:slast-1], smf[sfirst:slast-1], ssurfp[sfirst:slast-1]) \
                             out(sft[sfirst:slast-1]) \
                             \
                             in({mesh->znump[mesh->mapsz[z]], z=sfirst:slast-1}, \
                                {pu[mesh->mapsp2[mesh->mapss3[p]]], p=sfirst:slast-1}, \
                                {pxp[mesh->mapsp2[mesh->mapss3[p]]], p=sfirst:slast-1}, \
                                {pu[mesh->mapsp2[p2]], p2=sfirst:slast-1}, \
                                {exp[mesh->mapse[e2]], e2=sfirst:slast-1}, \
                                {zxp[mesh->mapsz[mesh->mapss3[z]]], z=sfirst:slast-1}, \
                                {pu[mesh->mapsp1[mesh->mapss3[p1]]], p1=sfirst:slast-1}, \
                                {exp[mesh->mapse[mesh->mapss3[e1]]], e1=sfirst:slast-1}, \
                                {elen[mesh->mapse[mesh->mapss3[e1]]], e1=sfirst:slast-1}, \
                                {elen[mesh->mapse[e2]], e2=sfirst:slast-1}, \
                                {zss[mesh->mapsz[z]], z=sfirst:slast-1}, \
                                {zrp[mesh->mapsz[z]], z=sfirst:slast-1}, \
                                {pu[mesh->mapsp2[mesh->mapss3[p]]], p=sfirst:slast-1}, \
                                {pu[mesh->mapsp1[mesh->mapss3[p1]]], p1=sfirst:slast-1}, \
                                {pu[mesh->mapsp2[p2]], p2=sfirst:slast-1}, \
                                {elen[mesh->mapse[mesh->mapss3[e1]]], e1=sfirst:slast-1}, \
                                {elen[mesh->mapse[e2]], e2=sfirst:slast-1}, \
                                {pu[mesh->mapsp1[p1]], p1=sfirst:slast-1}, \
                                {pu[mesh->mapsp2[p2]], p2=sfirst:slast-1}, \
                                {pxp[mesh->mapsp1[p1]], p1=sfirst:slast-1}, \
                                {pxp[mesh->mapsp2[p2]], p2=sfirst:slast-1}, \
                                {elen[mesh->mapse[e]], e=sfirst:slast-1}, \
                                {zss[mesh->mapsz[z]], z=sfirst:slast-1}) \
                             out(sfq[sfirst:slast-1], {zdu[mesh->mapsz[z]], z=sfirst:slast-1}) \
                             \
                             in(sfp[sfirst:slast-1], sfq[sfirst:slast-1], sft[sfirst:slast-1], \
                                {sfp[mesh->mapss3[s3]], s3=sfirst:slast-1}, {sfq[mesh->mapss3[s3]], s3=sfirst:slast-1}, \
                                {sft[mesh->mapss3[s3]], s3=sfirst:slast-1}) \
                             \
                             firstprivate(sfirst, slast, zfirst, zlast, cycle)
            {
                // save off zone variable values from previous cycle
                copy(&zvol[zfirst], &zvol[zlast], &zvol0[zfirst]);

                // 1a. compute new mesh geometry
                mesh->calcCtrs(pxp, exp, zxp, sfirst, slast);
                mesh->calcVols(pxp, zxp, sareap, svolp, zareap, zvolp,
                        sfirst, slast);
                mesh->calcSurfVecs(zxp, exp, ssurfp, sfirst, slast);
                mesh->calcEdgeLen(pxp, elen, sfirst, slast);
                mesh->calcCharLen(sareap, zdl, sfirst, slast);

                // 2. compute point masses
                calcRho(zm, zvolp, zrp, zfirst, zlast);
                calcCrnrMass(zrp, zareap, smf, cmaswt, sfirst, slast);

                // 3. compute material state (half-advanced)
                pgas->calcStateAtHalf(zr, zvolp, zvol0, ze, zwrate, zm, dt[cycle-1],
                        zp, zss, zfirst, zlast);

                // 4. compute forces
                pgas->calcForce(zp, ssurfp, sfp, sfirst, slast);
                tts->calcForce(zareap, zrp, zss, sareap, smf, ssurfp, sft,
                        sfirst, slast);
                qcs->calcForce(sfq, sfirst, slast);
                sumCrnrForce(sfp, sfq, sft, cftot, sfirst, slast);
            }
        }
    }  // for sch

    // Not needed anymore, check and abort done in calcVols()
    //mesh->checkBadSides();

    // sum corner masses, forces to points
    mesh->sumToPoints(cmaswt, pmaswt, 0);
    mesh->sumToPoints(cftot, pf, 1);

    //#pragma omp parallel for schedule(static)
    for (int pch = 0; pch < numpch; ++pch) {
        int pfirst = mesh->pchpfirst[pch];
        int plast = mesh->pchplast[pch];

        // Dependency blocks order:
        //   applyFixedBC // Uses helper arrays set up in Hydro init()
        //   calcAccel
        //   advPosFull
        #pragma oss task label("Steps 4a-6: BC & accel & adv") \
                         inout({pu0[bcstaskindices[pch][i]], i=0:bcstasksizes[pch]-1}, \
                               {pf[bcstaskindices[pch][i]], i=0:bcstasksizes[pch]-1}) \
                         \
                         in(pf[pfirst:plast-1], pmaswt[pfirst:plast-1]) \
                         out(pap[pfirst:plast-1]) \
                         \
                         in(dt[cycle-1], px0[pfirst:plast-1], pu0[pfirst:plast-1], pap[pfirst:plast-1]) \
                         out(px[pfirst:plast-1], pu[pfirst:plast-1]) \
                         \
                         firstprivate(pfirst, plast, cycle)
        {
            // 4a. apply boundary conditions
            for (int i = 0; i < bcs.size(); ++i) {
                int bfirst = bcs[i]->pchbfirst[pch];
                int blast = bcs[i]->pchblast[pch];
                bcs[i]->applyFixedBC(pu0, pf, bfirst, blast);
            }

            // 5. compute accelerations
            calcAccel(pf, pmaswt, pap, pfirst, plast);

            // ===== Corrector step =====
            // 6. advance mesh to end of time step
            advPosFull(px0, pu0, pap, dt[cycle-1], px, pu, pfirst, plast);
        }
    }  // for pch

    // dtrec reset moved into calcDtHydro task
    //resetDtHydro();

    //#pragma omp parallel for schedule(static)
    for (int sch = 0; sch < numsch; ++sch) {
        int sfirst = mesh->schsfirst[sch];
        int slast = mesh->schslast[sch];
        int zfirst = mesh->schzfirst[sch];
        int zlast = mesh->schzlast[sch];

        // Dependency blocks order:
        //   calcCtrls
        //   calcVols
        //
        // "zlast = (slast < nums ? mapsz[slast] : numz);" in calcCtrls and calcVols
        // Two separate types of task needed
        if (slast < mesh->nums) {
            #pragma oss task label("Step 6a standard: new mesh geometry") \
                             in({px[mesh->mapsp1[p1]], p1=sfirst:slast-1}, {px[mesh->mapsp2[p2]], p2=sfirst:slast-1}) \
                             out({zx[mesh->mapsz[z]], z=sfirst:slast-1}, {ex[mesh->mapse[e]], e=sfirst:slast-1}) \
                             \
                             in({px[mesh->mapsp1[p1]], p1=sfirst:slast-1}, {px[mesh->mapsp2[p2]], p2=sfirst:slast-1}, \
                                {zx[mesh->mapsz[z]], z=sfirst:slast-1}) \
                             out({zvol[mesh->mapsz[z]], z=sfirst:slast-1}, {zarea[mesh->mapsz[z]], z=sfirst:slast-1}, \
                                 svol[sfirst:slast-1], sarea[sfirst:slast-1]) \
                             \
                             firstprivate(sfirst, slast, zfirst, zlast)
            {
                // 6a. compute new mesh geometry
                mesh->calcCtrs(px, ex, zx, sfirst, slast);
                mesh->calcVols(px, zx, sarea, svol, zarea, zvol,
                        sfirst, slast);
            }
        } else {
            #pragma oss task label("Step 6a numz: new mesh geometry") \
                             in({px[mesh->mapsp1[p1]], p1=sfirst:slast-1}, {px[mesh->mapsp2[p2]], p2=sfirst:slast-1}) \
                             out(zx[mesh->mapsz[sfirst]:mesh->numz-1], {zx[mesh->mapsz[z]], z=sfirst:slast-1}, \
                                 {ex[mesh->mapse[e]], e=sfirst:slast-1}) \
                             \
                             in({px[mesh->mapsp1[p1]], p1=sfirst:slast-1}, {px[mesh->mapsp2[p2]], p2=sfirst:slast-1}, \
                                {zx[mesh->mapsz[z]], z=sfirst:slast-1}) \
                             out(zvol[mesh->mapsz[sfirst]:mesh->numz-1], {zvol[mesh->mapsz[z]], z=sfirst:slast-1}, \
                                 zarea[mesh->mapsz[sfirst]:mesh->numz-1], {zarea[mesh->mapsz[z]], z=sfirst:slast-1}, \
                                 svol[sfirst:slast-1], sarea[sfirst:slast-1]) \
                             \
                             firstprivate(sfirst, slast, zfirst, zlast)
            {
                // 6a. compute new mesh geometry
                mesh->calcCtrs(px, ex, zx, sfirst, slast);
                mesh->calcVols(px, zx, sarea, svol, zarea, zvol,
                        sfirst, slast);
            }
        }
    }  // for sch
    #pragma oss taskwait

    for (int zch = 0; zch < mesh->numzch; ++zch) {
        int zfirst = mesh->zchzfirst[zch];
        int zlast = mesh->zchzlast[zch];

        // Dependency blocks order:
        //   calcDtHydro (i.e. calcDtCourant + calcDtVolume)
        //
        // Formerly step 9.
        #pragma oss task label("Step 6.5: compute timestep for next cycle") \
                    in(dt[cycle-1], zdu[zfirst:zlast-1], zss[zfirst:zlast-1], zdl[zfirst:zlast-1], \
                       zvol[zfirst:zlast-1], zvol0[zfirst:zlast-1])\
                    out (dtrec) \
                    firstprivate(zfirst, zlast, zch, cycle)
        {
            // 6.5.  compute timestep for next cycle
            dtrec[zch] = 1.e99;
            strcpy(msgdtrec[zch], "Hydro default");
            calcDtHydro(zdl, zvol, zvol0, dt[cycle-1], zfirst, zlast, zch);
        }
    } // for zch

    for (int sch = 0; sch < numsch; ++sch) {
        int sfirst = mesh->schsfirst[sch];
        int slast = mesh->schslast[sch];
        int zfirst = mesh->schzfirst[sch];
        int zlast = mesh->schzlast[sch];

        // Dependency blocks order:
        //   fill & calcWork
        #pragma oss task label("Step 7: compute work") \
                         in(dt[cycle-1], sfp[sfirst:slast-1], sfq[sfirst:slast-1], \
                            {pu0[mesh->mapsp1[p1]], p1=sfirst:slast-1}, {pu0[mesh->mapsp2[p2]], p2=sfirst:slast-1}, \
                            {pu[mesh->mapsp1[p1]], p1=sfirst:slast-1}, {pu[mesh->mapsp2[p2]], p2=sfirst:slast-1}, \
                            {pxp[mesh->mapsp1[p1]], p1=sfirst:slast-1}, {pxp[mesh->mapsp2[p2]], p2=sfirst:slast-1}) \
                         out(zw[zfirst:zlast-1], \
                             {zw[mesh->mapsz[z]], z=sfirst:slast-1}, {zetot[mesh->mapsz[z]], z=sfirst:slast-1}) \
                         \
                         firstprivate(sfirst, slast, zfirst, zlast, cycle)
        {
            // 7. compute work
            fill(&zw[zfirst], &zw[zlast], 0.);
            calcWork(sfp, sfq, pu0, pu, pxp, dt[cycle-1], zw, zetot,
                    sfirst, slast);
        }
    }  // for sch

    // Not needed anymore, check and abort done in calcVols()
    //mesh->checkBadSides();

    //#pragma omp parallel for schedule(static)
    for (int zch = 0; zch < mesh->numzch; ++zch) {
        int zfirst = mesh->zchzfirst[zch];
        int zlast = mesh->zchzlast[zch];

        // Dependency blocks order:
        //   calcWorkRate
        //   calcEnergy
        //   calcRho
        #pragma oss task label("Steps 7a-8: work rate & energy & rho") \
                    in(dt[cycle-1], zvol[zfirst:zlast-1], zvol0[zfirst:zlast-1], zw[zfirst:zlast-1], zp[zfirst:zlast-1]) \
                    out(zwrate[zfirst:zlast-1]) \
                    \
                    in(zetot[zfirst:zlast-1], zm[zfirst:zlast-1]) \
                    out(ze[zfirst:zlast-1]) \
                    \
                    in(zm[zfirst:zlast-1], zvol[zfirst:zlast-1]) \
                    out(zr[zfirst:zlast-1]) \
                    \
                    firstprivate(zfirst, zlast, cycle)
        {
            // 7a. compute work rate
            calcWorkRate(zvol0, zvol, zw, zp, dt[cycle-1], zwrate, zfirst, zlast);

            // 8. update state variables
            calcEnergy(zetot, zm, ze, zfirst, zlast);
            calcRho(zm, zvol, zr, zfirst, zlast);
        }
    }  // for zch
    #pragma oss taskwait
}


void Hydro::advPosHalf(
        const double2* px0,
        const double2* pu0,
        const double dt,
        double2* pxp,
        const int pfirst,
        const int plast) {

    double dth = 0.5 * dt;

    #pragma ivdep
    for (int p = pfirst; p < plast; ++p) {
        pxp[p] = px0[p] + pu0[p] * dth;
    }
}


void Hydro::advPosFull(
        const double2* px0,
        const double2* pu0,
        const double2* pa,
        const double dt,
        double2* px,
        double2* pu,
        const int pfirst,
        const int plast) {

    #pragma ivdep
    for (int p = pfirst; p < plast; ++p) {
        pu[p] = pu0[p] + pa[p] * dt;
        px[p] = px0[p] + 0.5 * (pu[p] + pu0[p]) * dt;
    }

}


void Hydro::calcCrnrMass(
        const double* zr,
        const double* zarea,
        const double* smf,
        double* cmaswt,
        const int sfirst,
        const int slast) {

    #pragma ivdep
    for (int s = sfirst; s < slast; ++s) {
        int s3 = mesh->mapss3[s];
        int z = mesh->mapsz[s];

        double m = zr[z] * zarea[z] * 0.5 * (smf[s] + smf[s3]);
        cmaswt[s] = m;
    }
}


void Hydro::sumCrnrForce(
        const double2* sf,
        const double2* sf2,
        const double2* sf3,
        double2* cftot,
        const int sfirst,
        const int slast) {

    #pragma ivdep
    for (int s = sfirst; s < slast; ++s) {
        int s3 = mesh->mapss3[s];

        double2 f = (sf[s] + sf2[s] + sf3[s]) -
                    (sf[s3] + sf2[s3] + sf3[s3]);
        cftot[s] = f;
    }
}


void Hydro::calcAccel(
        const double2* pf,
        const double* pmass,
        double2* pa,
        const int pfirst,
        const int plast) {

    const double fuzz = 1.e-99;

    #pragma ivdep
    for (int p = pfirst; p < plast; ++p) {
        pa[p] = pf[p] / max(pmass[p], fuzz);
    }

}


void Hydro::calcRho(
        const double* zm,
        const double* zvol,
        double* zr,
        const int zfirst,
        const int zlast) {

    #pragma ivdep
    for (int z = zfirst; z < zlast; ++z) {
        zr[z] = zm[z] / zvol[z];
    }

}


void Hydro::calcWork(
        const double2* sf,
        const double2* sf2,
        const double2* pu0,
        const double2* pu,
        const double2* px,
        const double dt,
        double* zw,
        double* zetot,
        const int sfirst,
        const int slast) {

    // Compute the work done by finding, for each element/node pair,
    //   dwork= force * vavg
    // where force is the force of the element on the node
    // and vavg is the average velocity of the node over the time period

    const double dth = 0.5 * dt;

    for (int s = sfirst; s < slast; ++s) {
        int p1 = mesh->mapsp1[s];
        int p2 = mesh->mapsp2[s];
        int z = mesh->mapsz[s];

        double2 sftot = sf[s] + sf2[s];
        double sd1 = dot( sftot, (pu0[p1] + pu[p1]));
        double sd2 = dot(-sftot, (pu0[p2] + pu[p2]));
        double dwork = -dth * (sd1 * px[p1].x + sd2 * px[p2].x);

        zetot[z] += dwork;
        zw[z] += dwork;

    }

}


void Hydro::calcWorkRate(
        const double* zvol0,
        const double* zvol,
        const double* zw,
        const double* zp,
        const double dt,
        double* zwrate,
        const int zfirst,
        const int zlast) {
    double dtinv = 1. / dt;
    #pragma ivdep
    for (int z = zfirst; z < zlast; ++z) {
        double dvol = zvol[z] - zvol0[z];
        zwrate[z] = (zw[z] + zp[z] * dvol) * dtinv;
    }

}


void Hydro::calcEnergy(
        const double* zetot,
        const double* zm,
        double* ze,
        const int zfirst,
        const int zlast) {

    const double fuzz = 1.e-99;
    #pragma ivdep
    for (int z = zfirst; z < zlast; ++z) {
        ze[z] = zetot[z] / (zm[z] + fuzz);
    }

}


void Hydro::sumEnergy(
        const double* zetot,
        const double* zarea,
        const double* zvol,
        const double* zm,
        const double* smf,
        const double2* px,
        const double2* pu,
        double& ei,
        double& ek,
        const int zfirst,
        const int zlast,
        const int sfirst,
        const int slast) {

    // compute internal energy
    double sumi = 0.;
    for (int z = zfirst; z < zlast; ++z) {
        sumi += zetot[z];
    }
    // multiply by 2\pi for cylindrical geometry
    ei += sumi * 2 * M_PI;

    // compute kinetic energy
    // in each individual zone:
    // zone ke = zone mass * (volume-weighted average of .5 * u ^ 2)
    //         = zm sum(c in z) [cvol / zvol * .5 * u ^ 2]
    //         = sum(c in z) [zm * cvol / zvol * .5 * u ^ 2]
    double sumk = 0.;
    for (int s = sfirst; s < slast; ++s) {
        int s3 = mesh->mapss3[s];
        int p1 = mesh->mapsp1[s];
        int z = mesh->mapsz[s];

        double cvol = zarea[z] * px[p1].x * 0.5 * (smf[s] + smf[s3]);
        double cke = zm[z] * cvol / zvol[z] * 0.5 * length2(pu[p1]);
        sumk += cke;
    }
    // multiply by 2\pi for cylindrical geometry
    ek += sumk * 2 * M_PI;

}


void Hydro::calcDtCourant(
        const double* zdl,
        double& dtrec,
        char* msgdtrec,
        const int zfirst,
        const int zlast) {

    const double fuzz = 1.e-99;
    double dtnew = 1.e99;
    int zmin = -1;
    for (int z = zfirst; z < zlast; ++z) {
        double cdu = max(zdu[z], max(zss[z], fuzz));
        double zdthyd = zdl[z] * cfl / cdu;
        zmin = (zdthyd < dtnew ? z : zmin);
        dtnew = (zdthyd < dtnew ? zdthyd : dtnew);
    }

    if (dtnew < dtrec) {
        dtrec = dtnew;
        snprintf(msgdtrec, 80, "Hydro Courant limit for z = %d", zmin);
    }

}


void Hydro::calcDtVolume(
        const double* zvol,
        const double* zvol0,
        const double dtlast,
        double& dtrec,
        char* msgdtrec,
        const int zfirst,
        const int zlast) {

    double dvovmax = 1.e-99;
    int zmax = -1;
    for (int z = zfirst; z < zlast; ++z) {
        double zdvov = abs((zvol[z] - zvol0[z]) / zvol0[z]);
        zmax = (zdvov > dvovmax ? z : zmax);
        dvovmax = (zdvov > dvovmax ? zdvov : dvovmax);
    }
    double dtnew = dtlast * cflv / dvovmax;
    if (dtnew < dtrec) {
        dtrec = dtnew;
        snprintf(msgdtrec, 80, "Hydro dV/V limit for z = %d", zmax);
    }

}


void Hydro::calcDtHydro(
        const double* zdl,
        const double* zvol,
        const double* zvol0,
        const double dtlast,
        const int zfirst,
        const int zlast,
        const int zch) {

    double dtchunk = 1.e99;
    char msgdtchunk[80];

    calcDtCourant(zdl, dtchunk, msgdtchunk, zfirst, zlast);
    calcDtVolume(zvol, zvol0, dtlast, dtchunk, msgdtchunk,
            zfirst, zlast);
    if (dtchunk < dtrec[zch]) {
        dtrec[zch] = dtchunk;
        strncpy(msgdtrec[zch], msgdtchunk, 80);
    }

}


void Hydro::getDtHydro(
        double& dtnew,
        string& msgdtnew) {

    // Calculate local minimum across chunks here.
    // MUST only be called once all entries in dtrec array are filled.
    // This is a synchronisation point for all chunks.
    for (int zch = 0; zch < mesh->numzch; ++zch) {
        if (dtrec[zch] < dtnew) {
            dtnew = dtrec[zch];
            msgdtnew = string(msgdtrec[zch]);
        }
    }
}


/*void Hydro::resetDtHydro() {

    dtrec = 1.e99;
    strcpy(msgdtrec, "Hydro default");

}*/


void Hydro::writeEnergyCheck() {

    using Parallel::mype;

    double ei = 0.;
    double ek = 0.;
    //#pragma omp parallel for schedule(static)
    for (int sch = 0; sch < mesh->numsch; ++sch) {
        int sfirst = mesh->schsfirst[sch];
        int slast = mesh->schslast[sch];
        int zfirst = mesh->schzfirst[sch];
        int zlast = mesh->schzlast[sch];

        double eichunk = 0.;
        double ekchunk = 0.;
        sumEnergy(zetot, mesh->zarea, mesh->zvol, zm, mesh->smf,
                mesh->px, pu, eichunk, ekchunk,
                zfirst, zlast, sfirst, slast);
        //#pragma omp critical
        {
            ei += eichunk;
            ek += ekchunk;
        }
    }

    Parallel::globalSum(ei);
    Parallel::globalSum(ek);

    if (mype == 0) {
        cout << scientific << setprecision(6);
        cout << "Energy check:  "
             << "total energy  = " << setw(14) << ei + ek << endl;
        cout << "(internal = " << setw(14) << ei
             << ", kinetic = " << setw(14) << ek << ")" << endl;
    }

 }
