/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.
  Copyright 2014 STATOIL ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <opm/polymer/fullyimplicit/FullyImplicitTwophasePolymerSolver.hpp>

#include <opm/core/pressure/tpfa/trans_tpfa.h>
#include <opm/polymer/fullyimplicit/AutoDiffBlock.hpp>
#include <opm/polymer/fullyimplicit/AutoDiffHelpers.hpp>
#include <opm/polymer/fullyimplicit/IncompPropsAdInterface.hpp>
#include <opm/polymer/PolymerProperties.hpp>
#include <opm/polymer/PolymerState.hpp>
#include <opm/polymer/fullyimplicit/PolymerPropsAd.hpp>

#include <opm/core/grid.h>
#include <opm/core/linalg/LinearSolverInterface.hpp>
#include <opm/core/props/rock/RockCompressibility.hpp>
#include <opm/core/simulator/TwophaseState.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/well_controls.h>

#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <Eigen/Eigen>
#include <algorithm>

namespace Opm {





typedef AutoDiffBlock<double> ADB;
typedef ADB::V V;
typedef ADB::M M;
typedef Eigen::Array<double,
                     Eigen::Dynamic,
                     Eigen::Dynamic,
                     Eigen::RowMajor> DataBlock;





namespace {
    
    std::vector<int>
    buildAllCells(const int nc)
    {
        std::vector<int> all_cells(nc);
        for (int c = 0; c < nc; ++c) { all_cells[c] = c; }

        return all_cells;
    }
    struct Chop01 {
        double operator()(double x) const { return std::max(std::min(x, 1.0), 0.0); }
    };





    V computePerfPress(const UnstructuredGrid& grid, const Wells& wells, const V& rho, const double grav)
    {
        const int nw = wells.number_of_wells;
        const int nperf = wells.well_connpos[nw];
        const int dim = grid.dimensions;
        V wdp = V::Zero(nperf,1);
        assert(wdp.size() == rho.size());

        // Main loop, iterate over all perforations,
        // using the following formula:
        //    wdp(perf) = g*(perf_z - well_ref_z)*rho(perf)
        // where the total density rho(perf) is taken to be
        //    sum_p (rho_p*saturation_p) in the perforation cell.
        // [although this is computed on the outside of this function].
        for (int w = 0; w < nw; ++w) {
            const double ref_depth = wells.depth_ref[w];
            for (int j = wells.well_connpos[w]; j < wells.well_connpos[w + 1]; ++j) {
                const int cell = wells.well_cells[j];
                const double cell_depth = grid.cell_centroids[dim * cell + dim - 1];
                wdp[j] = rho[j]*grav*(cell_depth - ref_depth);
            }
        }
        return wdp;
    }

} //anonymous namespace





    FullyImplicitTwophasePolymerSolver::
    FullyImplicitTwophasePolymerSolver(const UnstructuredGrid&         grid,
                                       const IncompPropsAdInterface&   fluid,
                                       const PolymerPropsAd&           polymer_props_ad,
                                       const LinearSolverInterface&    linsolver,
                                       const Wells&                    wells,
                                       const double*                   gravity)
        : grid_ (grid)
        , fluid_(fluid)
        , polymer_props_ad_ (polymer_props_ad)
        , linsolver_(linsolver)
        , wells_(wells)
        , gravity_(gravity)
        , cells_ (buildAllCells(grid.number_of_cells))
        , ops_(grid)
        , wops_(wells)
//        , mob_(std::vector<ADB>(fluid.numPhases() + 1, ADB::null()))
		, cmax_(V::Zero(grid.number_of_cells))
        , rq_(fluid.numPhases() + 1)
        , residual_( { std::vector<ADB>(fluid.numPhases() + 1, ADB::null()), ADB::null(), ADB::null()})
     {
     }





    FullyImplicitTwophasePolymerSolver::
    WellOps::WellOps(const Wells& wells)
        : w2p(wells.well_connpos[ wells.number_of_wells ],
              wells.number_of_wells)
        , p2w(wells.number_of_wells,
              wells.well_connpos[ wells.number_of_wells ])
    {
        const int        nw   = wells.number_of_wells;
        const int* const wpos = wells.well_connpos;

        typedef Eigen::Triplet<double> Tri;

        std::vector<Tri> scatter, gather;
        scatter.reserve(wpos[nw]);
        gather .reserve(wpos[nw]);

        for (int w = 0, i = 0; w < nw; ++w) {
            for (; i < wpos[ w + 1 ]; ++i) {
                scatter.push_back(Tri(i, w, 1.0));
                gather .push_back(Tri(w, i, 1.0));
            }
        }

        w2p.setFromTriplets(scatter.begin(), scatter.end());
        p2w.setFromTriplets(gather .begin(), gather .end());
    }





    void
    FullyImplicitTwophasePolymerSolver::
    step(const double   dt,
         PolymerState&  x,
         WellState&     xw,
         const std::vector<double>& polymer_inflow,
		 std::vector<double>& src)
    {
        
        V pvol(grid_.number_of_cells);
        // Pore volume
        const V::Index nc = grid_.number_of_cells;
        V rho = V::Constant(pvol.size(), 1, *fluid_.porosity());
        std::transform(grid_.cell_volumes, grid_.cell_volumes + nc,
                       rho.data(), pvol.data(),
                       std::multiplies<double>());

        const V pvdt = pvol / dt;
        
        std::cout << "Starting constantState:\n"; 
        const SolutionState old_state = constantState(x, xw);
        std::cout << "Starting cmax:\n";
		computeCmax(x, old_state.concentration);
        std::cout << "Starting Accum:\n";
		computeAccum(old_state, 0);
        const double atol  = 1.0e-12;
        const double rtol  = 5.0e-8;
        const int    maxit = 40;

        assemble(pvdt, x, xw, polymer_inflow, src);

        const double r0  = residualNorm();
        int          it  = 0;
        std::cout << "\nIteration         Residual\n"
                  << std::setw(9) << it << std::setprecision(9)
                  << std::setw(18) << r0 << std::endl;
        bool resTooLarge = r0 > atol;
        while (resTooLarge && (it < maxit)) {
            const V dx = solveJacobianSystem();
            updateState(dx, x, xw);

            assemble(pvdt, x, xw, polymer_inflow, src);

            const double r = residualNorm();

            resTooLarge = (r > atol) && (r > rtol*r0);

            it += 1;
            std::cout << std::setw(9) << it << std::setprecision(9)
                      << std::setw(18) << r << std::endl;
        }

        if (resTooLarge) {
            std::cerr << "Failed to compute converged solution in " << it << " iterations. Ignoring!\n";
            // OPM_THROW(std::runtime_error, "Failed to compute converged solution in " << it << " iterations.");
        }
    }





    FullyImplicitTwophasePolymerSolver::ReservoirResidualQuant::ReservoirResidualQuant()
        : accum(2, ADB::null())
        , mflux(   ADB::null())
        , b    (   ADB::null())
        , head (   ADB::null())
        , mob  (   ADB::null())
    {
    }





    FullyImplicitTwophasePolymerSolver::SolutionState::SolutionState(const int np)
        : pressure   (    ADB::null())
        , saturation (np, ADB::null())
        , concentration ( ADB::null())
        , shear_mult ( ADB::null())
        , qs            ( ADB::null())
        , bhp           ( ADB::null())
    {
    }





    FullyImplicitTwophasePolymerSolver::SolutionState
    FullyImplicitTwophasePolymerSolver::constantState(const PolymerState& x,
                                                      const WellState&    xw)
    {
        const int nc = grid_.number_of_cells;
        const int np = x.numPhases();
        // The block pattern assumes the following primary variables:
        //    pressure
        //    water saturation 
        //    polymer concentration
        //    shear_multiplier
        //    well surface rates
        //    well bottom-hole pressure
        // Note that oil is assumed to always be present, but is never
        // a primary variable.
        std::vector<int> bpat(np + 1, nc);
        bpat.push_back(xw.bhp().size() * np);
        bpat.push_back(xw.bhp().size());
        
        SolutionState state(np);

        // Pressure.
        assert (not x.pressure().empty());
        const V p = Eigen::Map<const V>(& x.pressure()[0], nc);
        state.pressure = ADB::constant(p);

        // Saturation.
        assert (not x.saturation().empty());
        const DataBlock s_all = Eigen::Map<const DataBlock>(& x.saturation()[0], nc, np);
        for (int phase = 0; phase < np; ++phase) {
            state.saturation[phase] = ADB::constant(s_all.col(phase));
        }

        // Concentration.
        assert(not x.concentration().empty());
        const V c = Eigen::Map<const V>(&x.concentration()[0], nc);
        state.concentration = ADB::constant(c);

        // Shear multiplier.
//        assert(not x.shear_mult().empty());
//        const V sm = Eigen::Map<const V>(&x.shear_mult()[0], nc);
//        state.shear_mult = ADB::constant(sm);

        // Well rates.
        assert (not xw.wellRates().empty());
        // Need to reshuffle well rates, from ordered by wells, then phase,
        // to ordered by phase, then wells.
        const int nw = wells_.number_of_wells;
        // The transpose() below switches the ordering.
        const DataBlock wrates = Eigen::Map<const DataBlock>(& xw.wellRates()[0], nw, np).transpose();
        const V qs = Eigen::Map<const V>(wrates.data(), nw * np);
        state.qs = ADB::constant(qs, bpat);

        // Bottom hole pressure.
        assert (not xw.bhp().empty());
        const V bhp = Eigen::Map<const V>(& xw.bhp()[0], xw.bhp().size());
        state.bhp = ADB::constant(bhp, bpat);
        return state;
    }





    FullyImplicitTwophasePolymerSolver::SolutionState
    FullyImplicitTwophasePolymerSolver::variableState(const PolymerState& x,
                                                      const WellState&    xw)
    {
        const int nc = grid_.number_of_cells;
        const int np = x.numPhases();

        std::vector<V> vars0;
        vars0.reserve(np + 3); 

        // Initial pressure.
        assert (not x.pressure().empty());
        const V p = Eigen::Map<const V>(& x.pressure()[0], nc);
        vars0.push_back(p);

        // Initial saturation.
        assert (not x.saturation().empty());
        const DataBlock s_all = Eigen::Map<const DataBlock>(& x.saturation()[0], nc, np);
        const V sw = s_all.col(0);
        vars0.push_back(sw);

        // Initial concentration.
        assert (not x.concentration().empty());
        const V c = Eigen::Map<const V>(& x.concentration()[0], nc);
        vars0.push_back(c);

//        assert (not x.shear_mult().empty());
//        const V sm = Eigen::Map<const V>(& x.shear_mult()[0], nc);
//        vars0.push_back(sm);       

        // Initial well rates.
        assert (not xw.wellRates().empty());
        // Need to reshuffle well rates, from ordered by wells, then phase,
        // to ordered by phase, then wells.
        const int nw = wells_.number_of_wells;
        // The transpose() below switches the ordering.
        const DataBlock wrates = Eigen::Map<const DataBlock>(& xw.wellRates()[0], nw, np).transpose();
        const V qs = Eigen::Map<const V>(wrates.data(), nw * np);
        vars0.push_back(qs);

        // Initial well bottom hole pressure.
        assert (not xw.bhp().empty());
        const V bhp = Eigen::Map<const V>(& xw.bhp()[0], xw.bhp().size());
        vars0.push_back(bhp);

        std::vector<ADB> vars = ADB::variables(vars0);

        SolutionState state(np);

        // Pressure.
        int nextvar = 0;
        state.pressure = vars[ nextvar++ ];

        // Saturation.
        const std::vector<int>& bpat = vars[0].blockPattern();
        {
            ADB so  = ADB::constant(V::Ones(nc, 1), bpat);
            ADB sw = vars[ nextvar++ ];
            state.saturation[0] = sw;
            so = so - sw;
            state.saturation[1] = so;
        }
        
        // Concentration.
        state.concentration = vars[nextvar++];

        // Shear multiplier.
//        state.shear_mult = vars[nextvar++];

        // Qs.
        state.qs = vars[ nextvar++ ];

        // BHP.
        state.bhp = vars[ nextvar++ ];
        assert(nextvar == int(vars.size()));

        return state;
    }
  




    void
    FullyImplicitTwophasePolymerSolver::
    computeCmax(PolymerState& state,
				const ADB& c)
    {
        const int nc = grid_.number_of_cells;
        for (int i = 0; i < nc; ++i) {
		    cmax_(i) = std::max(cmax_(i), c.value()(i));
        }
		std::copy(&cmax_[0], &cmax_[0] + nc, state.maxconcentration().begin());
    }





    void
    FullyImplicitTwophasePolymerSolver::
	computeAccum(const SolutionState& state,
                 const int            aix  )
    {

        const std::vector<ADB>& sat   = state.saturation;
        const ADB&              c     = state.concentration;

        rq_[0].accum[aix] = sat[0];
        rq_[1].accum[aix] = sat[1];

		const ADB cmax = ADB::constant(cmax_, state.concentration.blockPattern());
        const ADB ads = polymer_props_ad_.adsorption(state.concentration, cmax);
        const double rho_rock = polymer_props_ad_.rockDensity();
        const V phi = Eigen::Map<const V>(&fluid_.porosity()[0], grid_.number_of_cells, 1);
        const double dead_pore_vol = polymer_props_ad_.deadPoreVol();

        rq_[2].accum[aix] = sat[0] * c * (1. - dead_pore_vol) + rho_rock * (1. - phi) / phi * ads;
    }





    void
    FullyImplicitTwophasePolymerSolver::
    assemble(const V&             pvdt,
             PolymerState&  x,
             const WellState&     xw,
             const std::vector<double>& polymer_inflow,
			 std::vector<double>& src)
    {
        // Create the primary variables.
        const SolutionState state = variableState(x, xw);

        const V trans = subset(transmissibility(), ops_.internal_faces);
        // Compute Shear multiplier.
//        std::cout << "Starting u_sq:\n";
//        const ADB u_sq = computeSquareGradP(trans, state);
//        std::cout << "Starting ShearMult:\n";
//        computeShearMult(x, state, u_sq.value());
//        std::cout << "Starting variableState:\n";
//        const SolutionState state = variableState(x, xw);
//        std::cout << "shear mult value: \n" << state.shear_mult.value() << std::endl;
		computeAccum(state, 1);
        // -------- Mass balance equations for water and oil --------
        const std::vector<ADB> kr = computeRelPerm(state);
		const ADB cmax = ADB::constant(cmax_, state.concentration.blockPattern());
        const ADB krw_eff = polymer_props_ad_.effectiveRelPerm(state.concentration, cmax, kr[0], state.saturation[0]);
        const ADB mc = computeMc(state);
        const V sh = Eigen::Map<const V>(&x.shear_mult()[0], grid_.number_of_cells); 
        computeMassFlux(trans, mc, kr[1], krw_eff, state);
//        const V face_areas = Eigen::Map<const V>(&grid_.face_areas[0], grid_.number_of_faces);
//        const V internal_face_areas = subset(face_areas, ops_.internal_faces);
    //    std::cout << "Water Equation Face Flux:\n";
    //    std::cout << rq_[0].mflux.value() << std::endl;
   //     std::cout << "Water Equation Face Velocity:\n";
     //   std::cout << rq_[0].mflux.value() / internal_face_areas << std::endl;
   //     exit(1);
        residual_.mass_balance[0] = pvdt*(rq_[0].accum[1] - rq_[0].accum[0])
                                    + ops_.div*rq_[0].mflux;
        residual_.mass_balance[1] = pvdt*(rq_[1].accum[1] - rq_[1].accum[0])
                                    + ops_.div*rq_[1].mflux;
        residual_.mass_balance[2] = pvdt*(rq_[2].accum[1] - rq_[2].accum[0])
                                    + ops_.div*rq_[2].mflux;
        //Shear multiplier equation.
 //       residual_.mass_balance[3] = shearMultFunc(state, cmax, u_sq); 
//        std::cout << residual_.mass_balance[3] << std::endl;
        // -------- Well equation, and well contributions to the mass balance equations --------

        // Contribution to mass balance will have to wait.

        const int nc = grid_.number_of_cells;
        const int np = wells_.number_of_phases;
        const int nw = wells_.number_of_wells;
        const int nperf = wells_.well_connpos[nw];

        const std::vector<int> well_cells(wells_.well_cells, wells_.well_cells + nperf);
        const V transw = Eigen::Map<const V>(wells_.WI, nperf);

        const ADB& bhp = state.bhp;

        const DataBlock well_s = wops_.w2p * Eigen::Map<const DataBlock>(wells_.comp_frac, nw, np).matrix();

        // Extract variables for perforation cell pressures
        // and corresponding perforation well pressures.
        const ADB p_perfcell = subset(state.pressure, well_cells);
        // Finally construct well perforation pressures and well flows.

        // Compute well pressure differentials.
        // Construct pressure difference vector for wells.
        const int dim = grid_.dimensions;
        if (gravity_) {
            for (int dd = 0; dd < dim -1; ++dd) {
                assert(g[dd] == 0.0);
            }
        }
        ADB cell_rho_total = ADB::constant(V::Zero(nc), state.pressure.blockPattern());
        for (int phase = 0; phase < 2; ++phase) {
            // For incompressible flow cell rho is the same.
            const ADB cell_rho = fluidDensity(phase, state.pressure);
            cell_rho_total += state.saturation[phase] * cell_rho;
        }
        ADB inj_rho_total = ADB::constant(V::Zero(nperf), state.pressure.blockPattern());
        assert(np == wells_.number_of_phases);
        const DataBlock compi = Eigen::Map<const DataBlock>(wells_.comp_frac, nw, np);
        for (int phase = 0; phase < 2; ++phase) {
            const ADB cell_rho = fluidDensity(phase, state.pressure);
            const V fraction = compi.col(phase);
            inj_rho_total += (wops_.w2p * fraction.matrix()).array() * subset(cell_rho, well_cells);
        }
        const V rho_perf_cell = subset(cell_rho_total, well_cells).value();
        const V rho_perf_well = inj_rho_total.value();
        V prodperfs = V::Constant(nperf, -1.0);
        for (int w = 0; w < nw; ++w) {
            if (wells_.type[w] == PRODUCER) {
                std::fill(prodperfs.data() + wells_.well_connpos[w],
                          prodperfs.data() + wells_.well_connpos[w+1], 1.0);
            }
        }
        const Selector<double> producer(prodperfs);
        const V rho_perf = producer.select(rho_perf_cell, rho_perf_well);
        const V well_perf_dp = computePerfPress(grid_, wells_, rho_perf, gravity_ ? gravity_[dim - 1] : 0.0);

        const ADB p_perfwell = wops_.w2p * bhp + well_perf_dp;
        const ADB nkgradp_well = transw * (p_perfcell - p_perfwell);
        // DUMP(nkgradp_well);
        const Selector<double> cell_to_well_selector(nkgradp_well.value());
        ADB well_rates_all = ADB::constant(V::Zero(nw*np), state.bhp.blockPattern());

        ADB perf_total_mob = subset(rq_[0].mob, well_cells) + subset(rq_[1].mob, well_cells);

        std::vector<ADB> well_contribs(np, ADB::null());
        std::vector<ADB> well_perf_rates(np, ADB::null());
        for (int phase = 0; phase < np; ++phase) {
//            const ADB& cell_b = rq_[phase].b;
  //          const ADB perf_b = subset(cell_b, well_cells);
            const ADB& cell_mob = rq_[phase].mob;
            const V well_fraction = compi.col(phase);
            // Using total mobilities for all phases for injection.
            const ADB perf_mob_injector = (wops_.w2p * well_fraction.matrix()).array() * perf_total_mob;
            const ADB perf_mob = producer.select(subset(cell_mob, well_cells),
                                                 perf_mob_injector);
            const ADB perf_flux = perf_mob * (nkgradp_well); // No gravity term for perforations.
            well_perf_rates[phase] = perf_flux;
            const ADB well_rates = wops_.p2w * well_perf_rates[phase];
            well_rates_all += superset(well_rates, Span(nw, 1, phase*nw), nw*np);

            // const ADB well_contrib = superset(perf_flux*perf_b, well_cells, nc);
            well_contribs[phase] = superset(perf_flux, well_cells, nc);
            // DUMP(well_contribs[phase]);
            residual_.mass_balance[phase] += well_contribs[phase];
			for (int i = 0; i < nc; ++i) {
				src[i] += well_contribs[phase].value()[i];
			}
        }

        // well rates contribs to polymer mass balance eqn.
        // for injection wells.
        const V polyin = Eigen::Map<const V>(& polymer_inflow[0], nc);
        const V poly_in_perf = subset(polyin, well_cells);
        const V poly_mc_cell = subset(mc, well_cells).value();
        const V poly_c = producer.select(poly_mc_cell, poly_in_perf);
        residual_.mass_balance[2] += superset(well_perf_rates[0] * poly_c, well_cells, nc);

        // Set the well flux equation
        residual_.well_flux_eq = state.qs + well_rates_all;
        // DUMP(residual_.well_flux_eq);

        // Handling BHP and SURFACE_RATE wells.
        V bhp_targets(nw);
        V rate_targets(nw);
        M rate_distr(nw, np*nw);
        for (int w = 0; w < nw; ++w) {
            const WellControls* wc = wells_.ctrls[w];
			if (well_controls_get_current_type(wc) == BHP) {
				bhp_targets[w] = well_controls_get_current_target(wc);
                rate_targets[w] = -1e100;
            } else if (well_controls_get_current_type(wc) == SURFACE_RATE) {
                bhp_targets[w] = -1e100;
				rate_targets[w] = well_controls_get_current_target(wc);
				{	
					const double* distr = well_controls_get_current_distr(wc);
            	    for (int phase = 0; phase < np; ++phase) {
                	    rate_distr.insert(w, phase*nw + w) = distr[phase];
                	}
				}
            } else {
                OPM_THROW(std::runtime_error, "Can only handle BHP type controls.");
            }
        }
        const ADB bhp_residual = bhp - bhp_targets;
        const ADB rate_residual = rate_distr * state.qs - rate_targets;
        // Choose bhp residual for positive bhp targets.
        Selector<double> bhp_selector(bhp_targets);
        residual_.well_eq = bhp_selector.select(bhp_residual, rate_residual);
 //       residual_.well_eq = bhp_residual;

    }




   
    std::vector<ADB>
    FullyImplicitTwophasePolymerSolver::
	computePressures(const SolutionState& state) const
    {
        const ADB sw = state.saturation[0];
        const ADB so = state.saturation[1];

        // convert the pressure offsets to the capillary pressures
        std::vector<ADB> pressure = fluid_.capPress(sw, so, cells_);
		pressure[0] = pressure[0] - pressure[1];

        // add the total pressure to the capillary pressures
        for (int phaseIdx = 0; phaseIdx < 2; ++phaseIdx) {
            pressure[phaseIdx] += state.pressure;
        }

        return pressure;
    }





    void
    FullyImplicitTwophasePolymerSolver::computeMassFlux(const V&                trans,
                                               //         const V&                sh,
                                                        const ADB&              mc,
                                                        const ADB&              kro,
                                                        const ADB&              krw_eff,
                                                        const SolutionState&    state )
    {
        const double* mus = fluid_.viscosity();
 //       ADB inv_wat_eff_vis = polymer_props_ad_.effectiveInvWaterViscWithShear(state.concentration, sh, mus);
        ADB inv_wat_eff_vis = polymer_props_ad_.effectiveInvWaterVisc(state.concentration, mus);
        V   shear_mult = eclShearMult(trans, state);
//        std::cout << "shear_mult" << shear_mult << std::endl;
        rq_[0].mob = krw_eff * inv_wat_eff_vis / shear_mult ;
        rq_[1].mob = kro / V::Constant(kro.size(), 1, mus[1]);
        rq_[2].mob = mc * krw_eff * inv_wat_eff_vis / shear_mult;

        const int nc = grid_.number_of_cells; 
        V z(nc);
        // Compute z coordinates
        for (int c = 0; c < nc; ++c){
            z[c] = grid_.cell_centroids[c * 3 + 2];
        }
		std::vector<ADB> press = computePressures(state);
        for (int phase = 0; phase < 2; ++phase) {
            const ADB rho = fluidDensity(phase, state.pressure);
			ADB& head = rq_[phase].head;
            const ADB rhoavg = ops_.caver * rho;
            const ADB dp = ops_.ngrad * press[phase]
                           - gravity_[2] * (rhoavg * (ops_.ngrad * z.matrix()));
            head = trans * dp;
            UpwindSelector<double> upwind(grid_, ops_, head.value());
			const ADB& mob = rq_[phase].mob;
            rq_[phase].mflux = upwind.select(mob) * head;
        }
		rq_[2].head = rq_[0].head;
		UpwindSelector<double> upwind(grid_, ops_, rq_[2].head.value());
		rq_[2].mflux = upwind.select(rq_[2].mob) * rq_[2].head;
    }





    std::vector<ADB>
    FullyImplicitTwophasePolymerSolver::accumSource(const ADB&                 kro,
                                                    const ADB&                 krw_eff,
                                                    const ADB&                 c,
                                                    const std::vector<double>& src,
                                                    const std::vector<double>& polymer_inflow_c) const
    {
        //extract the source to out and in source.
        std::vector<double> outsrc;
        std::vector<double> insrc;
        std::vector<double>::const_iterator it;
        for (it = src.begin(); it != src.end(); ++it) {
            if (*it < 0) {
                outsrc.push_back(*it);
                insrc.push_back(0.0);
            } else if (*it > 0) {
                insrc.push_back(*it);
                outsrc.push_back(0.0);
            } else {
                outsrc.push_back(0);
                insrc.push_back(0);
            }
        }
        const V outSrc = Eigen::Map<const V>(& outsrc[0], grid_.number_of_cells);
        const V inSrc = Eigen::Map<const V>(& insrc[0], grid_.number_of_cells);
        const V polyin = Eigen::Map<const V>(& polymer_inflow_c[0], grid_.number_of_cells);
        // compute the out-fracflow.
        const std::vector<ADB> f = computeFracFlow();
        // compute the in-fracflow.
        V   zero = V::Zero(grid_.number_of_cells);
        V   one  = V::Ones(grid_.number_of_cells);

        std::vector<ADB> source;
        //water source
        source.push_back(f[0] * outSrc + one * inSrc);
        //oil source
        source.push_back(f[1] * outSrc + zero * inSrc);
        //polymer source
        source.push_back(f[0] * outSrc * c + one * inSrc * polyin);

        return source;
     }





    std::vector<ADB>
    FullyImplicitTwophasePolymerSolver::computeFracFlow() const
    {
        ADB total_mob = rq_[0].mob + rq_[1].mob;

        std::vector<ADB> fracflow;

        fracflow.push_back(rq_[0].mob / total_mob);
        fracflow.push_back(rq_[1].mob / total_mob);

        return fracflow;
    }





    V 
    FullyImplicitTwophasePolymerSolver::solveJacobianSystem() const
    {
        const int np = fluid_.numPhases();
    	if (np != 2) {
	        OPM_THROW(std::logic_error, "Only two-phase ok in FullyImplicitTwophasePolymerSolver.");
	    }
	    ADB mass_res = vertcat(residual_.mass_balance[0], residual_.mass_balance[1]);
        mass_res = vertcat(mass_res, residual_.mass_balance[2]);
        ADB well_res = vertcat(residual_.well_flux_eq, residual_.well_eq);
	    ADB total_res = collapseJacs(vertcat(mass_res, well_res));

        const Eigen::SparseMatrix<double, Eigen::RowMajor> matr = total_res.derivative()[0];
        V dx(V::Zero(total_res.size()));
        Opm::LinearSolverInterface::LinearSolverReport rep
            = linsolver_.solve(matr.rows(), matr.nonZeros(),
                               matr.outerIndexPtr(), matr.innerIndexPtr(), matr.valuePtr(),
                               total_res.value().data(), dx.data());
        if (!rep.converged) {
            OPM_THROW(std::runtime_error,
                      "FullyImplicitBlackoilSolver::solveJacobianSystem(): "
                      "Linear solver convergence failure.");
        }
        return dx;
    }





    void FullyImplicitTwophasePolymerSolver::updateState(const V& dx,
                                                  PolymerState& state,
                                                  WellState&    well_state) const
    {
        const int np = fluid_.numPhases();
        const int nc = grid_.number_of_cells;
        const int nw = wells_.number_of_wells;
        const V one = V::Constant(nc, 1.0);
        const V zero = V::Zero(nc);

        // Extract parts of dx corresponding to each part.
        const V dp = subset(dx, Span(nc));
        int varstart = nc;
        const V dsw = subset(dx, Span(nc, 1, varstart));
        varstart += dsw.size();
        const V dc = subset(dx, Span(nc, 1, varstart));
        varstart += dc.size();
//        const V dsm = subset(dx, Span(nc, 1, varstart));
//        varstart += dsm.size();
        const V dqs = subset(dx, Span(np*nw, 1, varstart));
        varstart += dqs.size();
        const V dbhp = subset(dx, Span(nw, 1, varstart));
        varstart += dbhp.size();

        assert(varstart == dx.size());

        // Pressure update.
        const V p_old = Eigen::Map<const V>(&state.pressure()[0], nc);
        const V p = p_old - dp;
        std::copy(&p[0], &p[0] + nc, state.pressure().begin());

        // Saturation updates.
        const double dsmax = 0.3;
        const DataBlock s_old = Eigen::Map<const DataBlock>(& state.saturation()[0], nc, np);
        V so = one;
        const V sw_old = s_old.col(0);
        const V dsw_limited = sign(dsw) * dsw.abs().min(dsmax);
        const V sw = (sw_old - dsw_limited).unaryExpr(Chop01());
        so -= sw;
        for (int c = 0; c < nc; ++c) {
            state.saturation()[c*np] = sw[c];
            state.saturation()[c*np + 1] = so[c];
        }
        
        // Concentration updates.
        const V c_old = Eigen::Map<const V>(&state.concentration()[0], nc);
        const V c = (c_old - dc).max(zero);
        std::copy(&c[0], &c[0] + nc, state.concentration().begin());

//        // Shear multiplier updates.
//        const V sm_old = Eigen::Map<const V>(&state.shear_mult()[0], nc);
//        const V sm = (sm_old - dsm).max(0.1); 
//        std::copy(&sm[0], &sm[0] + nc, state.shear_mult().begin());
        // Qs update.
        // Since we need to update the wellrates, that are ordered by wells,
        // from dqs which are ordered by phase, the simplest is to compute
        // dwr, which is the data from dqs but ordered by wells.
        const DataBlock wwr = Eigen::Map<const DataBlock>(dqs.data(), np, nw).transpose();
        const V dwr = Eigen::Map<const V>(wwr.data(), nw*np);
        const V wr_old = Eigen::Map<const V>(&well_state.wellRates()[0], nw*np);
        const V wr = wr_old - dwr;
        std::copy(&wr[0], &wr[0] + wr.size(), well_state.wellRates().begin());

        // Bhp update.
        const V bhp_old = Eigen::Map<const V>(&well_state.bhp()[0], nw, 1);
        const V bhp = bhp_old - dbhp;
        std::copy(&bhp[0], &bhp[0] + bhp.size(), well_state.bhp().begin());
    }
   




    std::vector<ADB>
    FullyImplicitTwophasePolymerSolver::computeRelPerm(const SolutionState& state) const
    {

        const ADB sw = state.saturation[0];
        const ADB so = state.saturation[1];

        return fluid_.relperm(sw, so, cells_);
    }
    
    
    
    
    
    double
    FullyImplicitTwophasePolymerSolver::residualNorm() const
    {
        double r = 0;
        for (std::vector<ADB>::const_iterator
                 b = residual_.mass_balance.begin(),
                 e = residual_.mass_balance.end();
             b != e; ++b)
        {
            r = std::max(r, (*b).value().matrix().lpNorm<Eigen::Infinity>());
        }

        r = std::max(r, residual_.well_flux_eq.value().matrix().lpNorm<Eigen::Infinity>());
        r = std::max(r, residual_.well_eq.value().matrix().lpNorm<Eigen::Infinity>());

        return r;
    }
   
   
   


    ADB
    FullyImplicitTwophasePolymerSolver::fluidDensity(const int phase,
		                                             const ADB p) const
    {
        const double* rhos = fluid_.surfaceDensity();
        ADB rho = ADB::constant(V::Constant(grid_.number_of_cells, 1, rhos[phase]),
                               p.blockPattern());
        
        return rho;
    }
  


 
   
    V
    FullyImplicitTwophasePolymerSolver::transmissibility() const
    {
        const V::Index nc = grid_.number_of_cells;
        V htrans(grid_.cell_facepos[nc]);
        V trans(grid_.cell_facepos[nc]);
        UnstructuredGrid* ug = const_cast<UnstructuredGrid*>(& grid_);
        tpfa_htrans_compute(ug, fluid_.permeability(), htrans.data());
        tpfa_trans_compute (ug, htrans.data(), trans.data());
        
        return trans;
    }

    // here mc means m(c) * c. 
    ADB
    FullyImplicitTwophasePolymerSolver::computeMc(const SolutionState& state) const
    {
        ADB c = state.concentration;
        return polymer_props_ad_.polymerWaterVelocityRatio(c);
    }

    V
    FullyImplicitTwophasePolymerSolver::eclShearMult(const V& trans, const SolutionState& state) const
    {
        const int nf = grid_.number_of_faces;
        const int nc = grid_.number_of_cells;
        V z(nc);
        ADB u_sq(ADB::null());
        ADB u(ADB::null());
        // Compute z coordinates
        for (int c = 0; c < nc; ++c){
            z[c] = grid_.cell_centroids[c * 3 + 2];
        }
		std::vector<ADB> press = computePressures(state);
        const ADB rhow = fluidDensity(0, state.pressure);
        const ADB rhoavg = ops_.caver * rhow;
        const ADB dp = ops_.ngrad * press[0]
                       - gravity_[2] * (rhoavg * (ops_.ngrad * z.matrix()));
        
        const V face_areas = Eigen::Map<const V>(&grid_.face_areas[0], nf);
        const V internal_face_areas = subset(face_areas, ops_.internal_faces);
        u = trans * dp / internal_face_areas;
        u_sq = ops_.div * u;

        V poro = V::Constant(grid_.number_of_cells, 1, *fluid_.porosity());
//        std::cout << "velocity:\n" << u_sq.value();
        V visc_mult = polymer_props_ad_.viscMult(state.concentration.value());
        V sh = polymer_props_ad_.shearMult(u_sq.value().abs()/poro);
//        std::cout << "sh:\n" << sh <<std::endl;
        const V one = V::Ones(grid_.number_of_cells);
        return (one + (visc_mult-one) * sh) / visc_mult;
    }
    
    ADB
    FullyImplicitTwophasePolymerSolver::computeSquareGradP(const V& trans, 
                                                           const SolutionState& state) const
    {
        const int nc = grid_.number_of_cells;
        const int nf = grid_.number_of_faces;
        V z(nc);
        ADB u_sq(ADB::null());
        ADB u(ADB::null());
        // Compute z coordinates
        for (int c = 0; c < nc; ++c){
            z[c] = grid_.cell_centroids[c * 3 + 2];
        }
		std::vector<ADB> press = computePressures(state);
        const ADB rhow = fluidDensity(0, state.pressure);
        const ADB rhoavg = ops_.caver * rhow;
        const ADB dp = ops_.ngrad * press[0]
                       - gravity_[2] * (rhoavg * (ops_.ngrad * z.matrix()));
        
        const V face_areas = Eigen::Map<const V>(&grid_.face_areas[0], nf);
        const V internal_face_areas = subset(face_areas, ops_.internal_faces);
        u = trans * dp / internal_face_areas;
        u_sq = ops_.div * (u * u);
        return u_sq;
    }

    ADB
    FullyImplicitTwophasePolymerSolver::
    shearMultFunc(const SolutionState& state, const ADB& cmax, const ADB& u_sq) const
    {
        const int nc = grid_.number_of_cells;
        // update parameters.
        ADB sw = state.saturation[0];
      
        ADB rk = polymer_props_ad_.permReduction(state.concentration, cmax);
        V poro = V::Constant(nc, 1, *fluid_.porosity());
        V perm = V::Constant(nc, 1, *fluid_.permeability());
        const std::vector<ADB> kr = computeRelPerm(state);
        const double* mus = fluid_.viscosity();
        const ADB inv_wat_eff_vis = polymer_props_ad_.effectiveInvWaterVisc(state.concentration, mus);
        ADB eps = ADB::constant(V::Constant(nc, 1e-9));
        ADB square_term = kr[0]*u_sq / ((sw)*poro*perm).pow(0.5);
        ADB sh = 6. * inv_wat_eff_vis * square_term / (state.shear_mult * rk.value());
//        std::cout << "sh:\n" << sh << std::endl;
        std::cout << "rk:\n" << rk << std::endl;
        V zeros = V::Zero(nc);
        return state.shear_mult - polymer_props_ad_.shearMult(sh) / polymer_props_ad_.shearMult(zeros);
    }

    ADB
    FullyImplicitTwophasePolymerSolver::
    shearMultFunc(const ADB& shear_mult, const SolutionState& state, const V& cmax_cells, const V& u_sq) const
    {
        const int nc = grid_.number_of_cells;
        // update parameters.
        V sw = state.saturation[0].value();
        V rk = polymer_props_ad_.permReduction(state.concentration.value(), cmax_cells);
        V poro = V::Constant(grid_.number_of_cells, 1, *fluid_.porosity());
        V perm = V::Constant(nc, 1, *fluid_.permeability());
        const std::vector<ADB> kr = computeRelPerm(state);
        const double* mus = fluid_.viscosity();
//        std::cout << "muw:  " << mus[0] << std::endl;
        const V inv_wat_eff_vis = polymer_props_ad_.effectiveInvWaterVisc(state.concentration.value(), mus);
        const V vw = polymer_props_ad_.shearWaterVelocity();
        const V vrf = polymer_props_ad_.shearVrf();
//        std::cout << "vw from table:\n";
//        std::cout << vw << std::endl;
//        std::cout << "vrf:\n";
//        std::cout << vrf << std::endl;
        V square_term(nc);
        square_term = (kr[0].value() * u_sq.abs() / ((sw ) * poro * perm)).pow(double(0.5));
//        std::cout << "square_term:\n" << square_term << std::endl; 
//        std::cout << "perm\n" << perm << std::endl;
//        std::cout << "inv_wat_eff_vis * 6:\n" <<  6. * inv_wat_eff_vis << std::endl;
        ADB sh = 6. * inv_wat_eff_vis * square_term / (shear_mult * rk);
        V zeros = V::Zero(nc);
        V new_shear_mult = polymer_props_ad_.shearMult(sh.value()) / polymer_props_ad_.shearMult(zeros); 
        
//        std::cout << "water velocity for shear:\n";
//        std::cout << fake_velocity.value() << std::endl;
//        std::cout << "shear mult after iterploation:\n";
//        std::cout << new_shear_mult << std::endl;

        return shear_mult - ADB::function(new_shear_mult, sh.derivative());
    }

    void
    FullyImplicitTwophasePolymerSolver::
    computeShearMult(PolymerState& x, const SolutionState& state, const V& u_sq)
    { 
        const int nc = grid_.number_of_cells;
        std::vector<int> blocksize;
        blocksize.push_back(nc);
//        std::cout << "starting creating shear mult variable:\n";
//        std::cout << state.shear_mult << std::endl;
        const V sh = Eigen::Map<const V>(&x.shear_mult()[0], nc);
        ADB shear_mult = ADB::variable(0, sh, blocksize);
//        std::cout << "ADB shear_mult:\n" << shear_mult << std::endl;
//        const V cmax = Eigen::Map<const V>(& x.maxconcentration()[0], nc);
//        std::cout << "compute shear mult func:\n"; 
        ADB eqs = shearMultFunc(shear_mult, state, cmax_, u_sq);
        double residual0 = eqs.value().matrix().lpNorm<Eigen::Infinity>();
        int iter = 0;
        double residual = residual0;
        const double dshear_max=1000;
        V shear_mult_value = shear_mult.value();
//        std::cout << "initial shear_mult value" << std::endl;
//        std::cout << shear_mult_value << std::endl;
        std::cout << "initial shear residual: " << residual0 << std::endl; 
        std::cout << "Starting computing shear multiplier.\n";
        std::cout << "Iter    Residual:\n";
        while (residual > 1e-4 && iter < 40) {
            iter = iter + 1;
            const Eigen::SparseMatrix<double, Eigen::RowMajor> matr = eqs.derivative()[0];
            V dshear(V::Zero(eqs.size()));
            Opm::LinearSolverInterface::LinearSolverReport rep 
                = linsolver_.solve(matr.rows(), matr.nonZeros(), matr.outerIndexPtr(), matr.innerIndexPtr(),
                                   matr.valuePtr(), eqs.value().data(), dshear.data());
            if (dshear.abs().maxCoeff() > dshear_max) {
                dshear = dshear_max /dshear.abs().maxCoeff() * dshear;
            }
            shear_mult_value = shear_mult.value() - dshear;
            shear_mult = ADB::variable(0, shear_mult_value, blocksize);
            eqs = shearMultFunc(shear_mult, state, cmax_, u_sq);
            residual = eqs.value().matrix().lpNorm<Eigen::Infinity>();
            std::cout << iter << "      "  << residual << std::endl;
        }
        std::cout << "Finish ShearMult Computing, total iteration is: " << iter << std::endl;
        if (iter >= 40) {
            std::cerr << "Failed to compute converged solution in" << iter << "iterations. Ignoring.\n";
        }
        if (residual < residual0) {
            std::copy(&shear_mult.value()[0], &shear_mult.value()[0] + nc, x.shear_mult().begin());
//            return shear_mult.value();
        } else {
            std::copy(&shear_mult_value[0], &shear_mult_value[0] + nc, x.shear_mult().begin());
//            return shear_mult_value;
        }
    }
} //namespace Opm 
