/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.
  Copyright 2014 STATOIL ASA.
  Copyright 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015 NTNU

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

#include <opm/polymer/fullyimplicit/FullyImplicitBlackoilPolymerSolver.hpp>
#include <opm/polymer/PolymerBlackoilState.hpp>

#include <opm/autodiff/AutoDiffBlock.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>
#include <opm/autodiff/GridHelpers.hpp>
#include <opm/autodiff/BlackoilPropsAdInterface.hpp>
#include <opm/autodiff/GeoProps.hpp>
#include <opm/autodiff/WellDensitySegmented.hpp>
#include <opm/autodiff/WellStateFullyImplicitBlackoil.hpp>

#include <opm/core/grid.h>
#include <opm/core/linalg/LinearSolverInterface.hpp>
#include <opm/core/linalg/ParallelIstlInformation.hpp>
#include <opm/core/props/rock/RockCompressibility.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/Exceptions.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/well_controls.h>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
//#include <fstream>

// A debugging utility.
#define DUMP(foo)                                                       \
    do {                                                                \
        std::cout << "==========================================\n"     \
                  << #foo ":\n"                                         \
                  << collapseJacs(foo) << std::endl;                    \
    } while (0)

#define DUMPVAL(foo)                                                    \
    do {                                                                \
        std::cout << "==========================================\n"     \
                  << #foo ":\n"                                         \
                  << foo.value() << std::endl;                          \
    } while (0)

#define DISKVAL(foo)                                                    \
    do {                                                                \
        std::ofstream os(#foo);                                         \
        os.precision(16);                                               \
        os << foo.value() << std::endl;                                 \
    } while (0)


namespace Opm {

typedef AutoDiffBlock<double> ADB;
typedef ADB::V V;
typedef ADB::M M;
typedef Eigen::Array<double,
                     Eigen::Dynamic,
                     Eigen::Dynamic,
                     Eigen::RowMajor> DataBlock;


namespace detail {


    std::vector<int>
    buildAllCells(const int nc)
    {
        std::vector<int> all_cells(nc);

        for (int c = 0; c < nc; ++c) { all_cells[c] = c; }

        return all_cells;
    }



    template <class PU>
    std::vector<bool>
    activePhases(const PU& pu)
    {
        const int maxnp = Opm::BlackoilPhases::MaxNumPhases;
        std::vector<bool> active(maxnp, false);

        for (int p = 0; p < pu.MaxNumPhases; ++p) {
            active[ p ] = pu.phase_used[ p ] != 0;
        }

        return active;
    }



    template <class PU>
    std::vector<int>
    active2Canonical(const PU& pu)
    {
        const int maxnp = Opm::BlackoilPhases::MaxNumPhases;
        std::vector<int> act2can(maxnp, -1);

        for (int phase = 0; phase < maxnp; ++phase) {
            if (pu.phase_used[ phase ]) {
                act2can[ pu.phase_pos[ phase ] ] = phase;
            }
        }

        return act2can;
    }



    template <class PU>
    int polymerPos(const PU& pu)
    {
        const int maxnp = Opm::BlackoilPhases::MaxNumPhases;
        int pos = 0;
        for (int phase = 0; phase < maxnp; ++phase) {
            if (pu.phase_used[phase]) {
                pos++;
            }
        }

        return pos;
    }

} // namespace detail

    template<class T>
    void FullyImplicitBlackoilPolymerSolver<T>::SolverParameter::
    reset()
    {
        // default values for the solver parameters
        dp_max_rel_      = 1.0e9;
        ds_max_          = 0.2;
        dr_max_rel_      = 1.0e9;
        relax_type_      = DAMPEN;
        relax_max_       = 0.5;
        relax_increment_ = 0.1;
        relax_rel_tol_   = 0.2;
        max_iter_        = 15;
        max_residual_allowed_  = std::numeric_limits< double >::max();
        tolerance_mb_    = 1.0e-7;
        tolerance_cnv_   = 1.0e-3;
        tolerance_wells_ = 1./Opm::unit::day;
    }

    template<class T>
    FullyImplicitBlackoilPolymerSolver<T>::SolverParameter::
    SolverParameter()
    {
        // set default values
        reset();
    }

    template<class T>
    FullyImplicitBlackoilPolymerSolver<T>::SolverParameter::
    SolverParameter( const parameter::ParameterGroup& param )
    {
        // set default values
        reset();

        // overload with given parameters
        dp_max_rel_  = param.getDefault("dp_max_rel", dp_max_rel_);
        ds_max_      = param.getDefault("ds_max", ds_max_);
        dr_max_rel_  = param.getDefault("dr_max_rel", dr_max_rel_);
        relax_max_   = param.getDefault("relax_max", relax_max_);
        max_iter_    = param.getDefault("max_iter", max_iter_);
        max_residual_allowed_ = param.getDefault("max_residual_allowed", max_residual_allowed_);

        tolerance_mb_    = param.getDefault("tolerance_mb", tolerance_mb_);
        tolerance_cnv_   = param.getDefault("tolerance_cnv", tolerance_cnv_);
        tolerance_wells_ = param.getDefault("tolerance_wells", tolerance_wells_ );

        std::string relaxation_type = param.getDefault("relax_type", std::string("dampen"));
        if (relaxation_type == "dampen") {
            relax_type_ = DAMPEN;
        } else if (relaxation_type == "sor") {
            relax_type_ = SOR;
        } else {
            OPM_THROW(std::runtime_error, "Unknown Relaxtion Type " << relaxation_type);
        }
    }


    template<class T>
    FullyImplicitBlackoilPolymerSolver<T>::
    FullyImplicitBlackoilPolymerSolver(const SolverParameter&          param,
                                       const Grid&                     grid ,
                                       const BlackoilPropsAdInterface& fluid,
                                       const DerivedGeology&           geo  ,
                                       const RockCompressibility*      rock_comp_props,
                                       const PolymerPropsAd&           polymer_props_ad,
                                       const Wells*                    wells,
                                       const NewtonIterationBlackoilInterface&    linsolver,
                                       const bool has_disgas,
                                       const bool has_vapoil,
                                       const bool has_polymer,
                                       const bool terminal_output)
        : grid_  (grid)
        , fluid_ (fluid)
        , geo_   (geo)
        , rock_comp_props_(rock_comp_props)
        , polymer_props_ad_(polymer_props_ad)
        , wells_ (wells)
        , linsolver_ (linsolver)
        , active_(detail::activePhases(fluid.phaseUsage()))
        , canph_ (detail::active2Canonical(fluid.phaseUsage()))
        , cells_ (detail::buildAllCells(Opm::AutoDiffGrid::numCells(grid)))
        , ops_   (grid)
        , wops_  (wells_)
        , cmax_(V::Zero(Opm::AutoDiffGrid::numCells(grid)))
        , has_disgas_(has_disgas)
        , has_vapoil_(has_vapoil)
        , has_polymer_(has_polymer)
        , poly_pos_(detail::polymerPos(fluid.phaseUsage()))
        , param_( param )
        , use_threshold_pressure_(false)
        , rq_    (fluid.numPhases())
        , phaseCondition_(AutoDiffGrid::numCells(grid))
        , residual_ ( { std::vector<ADB>(fluid.numPhases(), ADB::null()),
                        ADB::null(),
                        ADB::null() } )
        , terminal_output_ (terminal_output)
        , newtonIterations_( 0 )
        , linearIterations_( 0 )
    {
#if HAVE_MPI
        if ( terminal_output_ ) {
            if ( linsolver_.parallelInformation().type() == typeid(ParallelISTLInformation) )
            {
                const ParallelISTLInformation& info =
                    boost::any_cast<const ParallelISTLInformation&>(linsolver_.parallelInformation());
                // Only rank 0 does print to std::cout if terminal_output is enabled
                terminal_output_ = (info.communicator().rank()==0);
            }
        }
#endif
        if (has_polymer_) {
            if (!active_[Water]) {
                OPM_THROW(std::logic_error, "Polymer must solved in water!\n");
            }
            // If deck has polymer, residual_ should contain polymer equation.
            rq_.resize(fluid_.numPhases()+1);
            residual_.material_balance_eq.resize(fluid_.numPhases()+1, ADB::null());
            assert(poly_pos_ == fluid_.numPhases());
        }
    }



    template<class T>
    void
    FullyImplicitBlackoilPolymerSolver<T>::
    setThresholdPressures(const std::vector<double>& threshold_pressures_by_face)
    {
        const int num_faces = AutoDiffGrid::numFaces(grid_);
        if (int(threshold_pressures_by_face.size()) != num_faces) {
            OPM_THROW(std::runtime_error, "Illegal size of threshold_pressures_by_face input, must be equal to number of faces.");
        }
        use_threshold_pressure_ = true;
        // Map to interior faces.
        const int num_ifaces = ops_.internal_faces.size();
        threshold_pressures_by_interior_face_.resize(num_ifaces);
        for (int ii = 0; ii < num_ifaces; ++ii) {
            threshold_pressures_by_interior_face_[ii] = threshold_pressures_by_face[ops_.internal_faces[ii]];
        }
    }




    template<class T>
    int
    FullyImplicitBlackoilPolymerSolver<T>::
    step(const double   dt,
         PolymerBlackoilState& x ,
         WellStateFullyImplicitBlackoil& xw,
         const std::vector<double>& polymer_inflow)
    {
        const V pvdt = geo_.poreVolume() / dt;

        // Initial max concentration of this time step from PolymerBlackoilState.
        cmax_ = Eigen::Map<V>(&x.maxconcentration()[0], Opm::AutoDiffGrid::numCells(grid_));
        if (active_[Gas]) { updatePrimalVariableFromState(x); }

        // For each iteration we store in a vector the norms of the residual of
        // the mass balance for each active phase, the well flux and the well equations
        std::vector<std::vector<double>> residual_norms_history;

        assemble(pvdt, x, true, xw, polymer_inflow);


        bool converged = false;
        double omega = 1.;

        residual_norms_history.push_back(computeResidualNorms());

        int          it  = 0;
        converged = getConvergence(dt,it);
        const int sizeNonLinear = residual_.sizeNonLinear();

        V dxOld = V::Zero(sizeNonLinear);

        bool isOscillate = false;
        bool isStagnate = false;
        const enum RelaxType relaxtype = relaxType();
        int linearIterations = 0;

        while ((!converged) && (it < maxIter())) {
            V dx = solveJacobianSystem();

            // store number of linear iterations used
            linearIterations += linsolver_.iterations();

            detectNewtonOscillations(residual_norms_history, it, relaxRelTol(), isOscillate, isStagnate);

            if (isOscillate) {
                omega -= relaxIncrement();
                omega = std::max(omega, relaxMax());
                if (terminal_output_)
                {
                    std::cout << " Oscillating behavior detected: Relaxation set to " << omega << std::endl;
                }
            }

            stablizeNewton(dx, dxOld, omega, relaxtype);

            updateState(dx, x, xw);

            assemble(pvdt, x, false, xw, polymer_inflow);

            residual_norms_history.push_back(computeResidualNorms());

            // increase iteration counter
            ++it;

            converged = getConvergence(dt,it);
        }

        if (!converged) {
            // the runtime_error is caught by the AdaptiveTimeStepping
            OPM_THROW(std::runtime_error, "Failed to compute converged solution in " << it << " iterations.");
            return -1;
        }

        linearIterations_ += linearIterations;
        newtonIterations_ += it;
        // Update max concentration.
        computeCmax(x);
        return linearIterations;
    }





    template<class T>
    FullyImplicitBlackoilPolymerSolver<T>::ReservoirResidualQuant::ReservoirResidualQuant()
        : accum(2, ADB::null())
        , mflux(   ADB::null())
        , b    (   ADB::null())
        , head (   ADB::null())
        , mob  (   ADB::null())
    {
    }





    template<class T>
    FullyImplicitBlackoilPolymerSolver<T>::SolutionState::SolutionState(const int np)
        : pressure  (    ADB::null())
        , temperature(   ADB::null())
        , saturation(np, ADB::null())
        , rs        (    ADB::null())
        , rv        (    ADB::null())
        , concentration( ADB::null())
        , qs        (    ADB::null())
        , bhp       (    ADB::null())
        , canonical_phase_pressures(3, ADB::null())
    {
    }





    template<class T>
    FullyImplicitBlackoilPolymerSolver<T>::
    WellOps::WellOps(const Wells* wells)
      : w2p(),
        p2w()
    {
        if( wells )
        {
            w2p = M(wells->well_connpos[ wells->number_of_wells ], wells->number_of_wells);
            p2w = M(wells->number_of_wells, wells->well_connpos[ wells->number_of_wells ]);

            const int        nw   = wells->number_of_wells;
            const int* const wpos = wells->well_connpos;

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
    }





    template<class T>
    typename FullyImplicitBlackoilPolymerSolver<T>::SolutionState
    FullyImplicitBlackoilPolymerSolver<T>::constantState(const PolymerBlackoilState& x,
                                                         const WellStateFullyImplicitBlackoil&     xw) const
    {
        auto state = variableState(x, xw);
        makeConstantState(state);
        return state;
    }





    template<class T>
    void
    FullyImplicitBlackoilPolymerSolver<T>::makeConstantState(SolutionState& state) const
    {
        // HACK: throw away the derivatives. this may not be the most
        // performant way to do things, but it will make the state
        // automatically consistent with variableState() (and doing
        // things automatically is all the rage in this module ;)
        state.pressure = ADB::constant(state.pressure.value());
        state.temperature = ADB::constant(state.temperature.value());
        state.rs = ADB::constant(state.rs.value());
        state.rv = ADB::constant(state.rv.value());
        state.concentration = ADB::constant(state.concentration.value());
        const int num_phases = state.saturation.size();
        for (int phaseIdx = 0; phaseIdx < num_phases; ++ phaseIdx) {
            state.saturation[phaseIdx] = ADB::constant(state.saturation[phaseIdx].value());
        }
        state.qs = ADB::constant(state.qs.value());
        state.bhp = ADB::constant(state.bhp.value());
        assert(state.canonical_phase_pressures.size() == Opm::BlackoilPhases::MaxNumPhases);
        for (int canphase = 0; canphase < Opm::BlackoilPhases::MaxNumPhases; ++canphase) {
            ADB& pp = state.canonical_phase_pressures[canphase];
            pp = ADB::constant(pp.value());
        }
    }





    template<class T>
    typename FullyImplicitBlackoilPolymerSolver<T>::SolutionState
    FullyImplicitBlackoilPolymerSolver<T>::variableState(const PolymerBlackoilState& x,
                                                         const WellStateFullyImplicitBlackoil&     xw) const
    {
        using namespace Opm::AutoDiffGrid;
        const int nc = numCells(grid_);
        const int np = x.numPhases();

        std::vector<V> vars0;
        // p, Sw and Rs, Rv or Sg, concentration is used as primary depending on solution conditions
        vars0.reserve(np + 2);
        // Initial pressure.
        assert (not x.pressure().empty());
        const V p = Eigen::Map<const V>(& x.pressure()[0], nc, 1);
        vars0.push_back(p);

        // Initial saturation.
        assert (not x.saturation().empty());
        const DataBlock s = Eigen::Map<const DataBlock>(& x.saturation()[0], nc, np);
        const Opm::PhaseUsage pu = fluid_.phaseUsage();
        // We do not handle a Water/Gas situation correctly, guard against it.
        assert (active_[ Oil]);
        if (active_[ Water ]) {
            const V sw = s.col(pu.phase_pos[ Water ]);
            vars0.push_back(sw);
        }

        // store cell status in vectors
        V isRs = V::Zero(nc,1);
        V isRv = V::Zero(nc,1);
        V isSg = V::Zero(nc,1);

        if (active_[ Gas ]){
            for (int c = 0; c < nc ; c++ ) {
                switch (primalVariable_[c]) {
                case PrimalVariables::RS:
                    isRs[c] = 1;
                    break;

                case PrimalVariables::RV:
                    isRv[c] = 1;
                    break;

                default:
                    isSg[c] = 1;
                    break;
                }
            }


            // define new primary variable xvar depending on solution condition
            V xvar(nc);
            const V sg = s.col(pu.phase_pos[ Gas ]);
            const V rs = Eigen::Map<const V>(& x.gasoilratio()[0], x.gasoilratio().size());
            const V rv = Eigen::Map<const V>(& x.rv()[0], x.rv().size());
            xvar = isRs*rs + isRv*rv + isSg*sg;
            vars0.push_back(xvar);
        }
        
        // Initial polymer concentration.
        if (has_polymer_) {
            assert (not x.concentration().empty());
            const V c = Eigen::Map<const V>(& x.concentration()[0], nc, 1);
            vars0.push_back(c);
        }

        // Initial well rates.
        if( wellsActive() ) {
            // Need to reshuffle well rates, from phase running fastest
            // to wells running fastest.
            const int nw = wells().number_of_wells;

            // The transpose() below switches the ordering.
            const DataBlock wrates = Eigen::Map<const DataBlock>(& xw.wellRates()[0], nw, np).transpose();
            const V qs = Eigen::Map<const V>(wrates.data(), nw*np);
            vars0.push_back(qs);

            // Initial well bottom-hole pressure.
            assert (not xw.bhp().empty());
            const V bhp = Eigen::Map<const V>(& xw.bhp()[0], xw.bhp().size());
            vars0.push_back(bhp);
        } else {
            // push null states for qs and bhp
            vars0.push_back(V());
            vars0.push_back(V());
        }

        std::vector<ADB> vars = ADB::variables(vars0);

        SolutionState state(np);

        // Pressure.
        int nextvar = 0;
        state.pressure = vars[ nextvar++ ];

        // temperature
        const V temp = Eigen::Map<const V>(& x.temperature()[0], x.temperature().size());
        state.temperature = ADB::constant(temp);

        // Saturations
        const std::vector<int>& bpat = vars[0].blockPattern();
        {

            ADB so = ADB::constant(V::Ones(nc, 1), bpat);

            if (active_[ Water ]) {
                ADB& sw = vars[ nextvar++ ];
                state.saturation[pu.phase_pos[ Water ]] = sw;
                so -= sw;
            }

            if (active_[ Gas ]) {
                // Define Sg Rs and Rv in terms of xvar.
                // Xvar is only defined if gas phase is active
                const ADB& xvar = vars[ nextvar++ ];
                ADB& sg = state.saturation[ pu.phase_pos[ Gas ] ];
                sg = isSg*xvar + isRv* so;
                so -= sg;

                if (active_[ Oil ]) {
                    // RS and RV is only defined if both oil and gas phase are active.
                    const ADB& sw = (active_[ Water ]
                                             ? state.saturation[ pu.phase_pos[ Water ] ]
                                             : ADB::constant(V::Zero(nc, 1), bpat));
                    state.canonical_phase_pressures = computePressures(state.pressure, sw, so, sg);
                    const ADB rsSat = fluidRsSat(state.canonical_phase_pressures[ Oil ], so , cells_);
                    if (has_disgas_) {
                        state.rs = (1-isRs) * rsSat + isRs*xvar;
                    } else {
                        state.rs = rsSat;
                    }
                    const ADB rvSat = fluidRvSat(state.canonical_phase_pressures[ Gas ], so , cells_);
                    if (has_vapoil_) {
                        state.rv = (1-isRv) * rvSat + isRv*xvar;
                    } else {
                        state.rv = rvSat;
                    }
                }
            }

            if (active_[ Oil ]) {
                // Note that so is never a primary variable.
                state.saturation[pu.phase_pos[ Oil ]] = so;
            }
        }

        // Concentration.
        if (has_polymer_) {
            state.concentration = vars[nextvar++];
        }

        // Qs.
        state.qs = vars[ nextvar++ ];

        // Bhp.
        state.bhp = vars[ nextvar++ ];

        assert(nextvar == int(vars.size()));

        return state;
    }





    template<class T>
    void
    FullyImplicitBlackoilPolymerSolver<T>::computeAccum(const SolutionState& state,
                                                        const int            aix  )
    {
        const Opm::PhaseUsage& pu = fluid_.phaseUsage();

        const ADB&              press = state.pressure;
        const ADB&              temp  = state.temperature;
        const std::vector<ADB>& sat   = state.saturation;
        const ADB&              rs    = state.rs;
        const ADB&              rv    = state.rv;
        const ADB&              c     = state.concentration;

        const std::vector<PhasePresence> cond = phaseCondition();

        const ADB pv_mult = poroMult(press);

        const int maxnp = Opm::BlackoilPhases::MaxNumPhases;
        for (int phase = 0; phase < maxnp; ++phase) {
            if (active_[ phase ]) {
                const int pos = pu.phase_pos[ phase ];
                rq_[pos].b = fluidReciprocFVF(phase, state.canonical_phase_pressures[phase], temp, rs, rv, cond, cells_);
                rq_[pos].accum[aix] = pv_mult * rq_[pos].b * sat[pos];
                // DUMP(rq_[pos].b);
                // DUMP(rq_[pos].accum[aix]);
            }
        }

        if (active_[ Oil ] && active_[ Gas ]) {
            // Account for gas dissolved in oil and vaporized oil
            const int po = pu.phase_pos[ Oil ];
            const int pg = pu.phase_pos[ Gas ];

            // Temporary copy to avoid contribution of dissolved gas in the vaporized oil
            // when both dissolved gas and vaporized oil are present.
            const ADB accum_gas_copy =rq_[pg].accum[aix];

            rq_[pg].accum[aix] += state.rs * rq_[po].accum[aix];
            rq_[po].accum[aix] += state.rv * accum_gas_copy;
            //DUMP(rq_[pg].accum[aix]);
        }
                
        if (has_polymer_) {
            // compute polymer properties.
            const ADB cmax = ADB::constant(cmax_, state.concentration.blockPattern());
            const ADB ads  = polymer_props_ad_.adsorption(state.concentration, cmax);
            const double rho_rock = polymer_props_ad_.rockDensity();
            const V phi = Eigen::Map<const V>(& fluid_.porosity()[0], AutoDiffGrid::numCells(grid_), 1);
            const double dead_pore_vol = polymer_props_ad_.deadPoreVol();
            // compute total phases and determin polymer position.
            rq_[poly_pos_].accum[aix] = pv_mult * rq_[pu.phase_pos[Water]].b * sat[pu.phase_pos[Water]] * c * (1. - dead_pore_vol) 
                                        + pv_mult * rho_rock * (1. - phi) / phi * ads;
        }
 
    }






    template<class T>
    void FullyImplicitBlackoilPolymerSolver<T>::computeCmax(PolymerBlackoilState& state)
    {
        const int nc = AutoDiffGrid::numCells(grid_);
        V tmp = V::Zero(nc);
        for (int i = 0; i < nc; ++i) {
            tmp[i] = std::max(state.maxconcentration()[i], state.concentration()[i]);
        }
        std::copy(&tmp[0], &tmp[0] + nc, state.maxconcentration().begin());
    }





    template<class T>
    void FullyImplicitBlackoilPolymerSolver<T>::computeWellConnectionPressures(const SolutionState& state,
                                                                               const WellStateFullyImplicitBlackoil& xw)
    {
        if( ! wellsActive() ) return ;

        using namespace Opm::AutoDiffGrid;
        // 1. Compute properties required by computeConnectionPressureDelta().
        //    Note that some of the complexity of this part is due to the function
        //    taking std::vector<double> arguments, and not Eigen objects.
        const int nperf = wells().well_connpos[wells().number_of_wells];
        const int nw = wells().number_of_wells;
        const std::vector<int> well_cells(wells().well_cells, wells().well_cells + nperf);

        // Compute the average pressure in each well block
        const V perf_press = Eigen::Map<const V>(xw.perfPress().data(), nperf);
        V avg_press = perf_press*0;
        for (int w = 0; w < nw; ++w) {
            for (int perf = wells().well_connpos[w]; perf < wells().well_connpos[w+1]; ++perf) {
                const double p_above = perf == wells().well_connpos[w] ? state.bhp.value()[w] : perf_press[perf - 1];
                const double p_avg = (perf_press[perf] + p_above)/2;
                avg_press[perf] = p_avg;
            }
        }

        // Use cell values for the temperature as the wells don't knows its temperature yet.
        const ADB perf_temp = subset(state.temperature, well_cells);

        // Compute b, rsmax, rvmax values for perforations.
        // Evaluate the properties using average well block pressures
        // and cell values for rs, rv, phase condition and temperature.
        const ADB avg_press_ad = ADB::constant(avg_press);
        std::vector<PhasePresence> perf_cond(nperf);
        const std::vector<PhasePresence>& pc = phaseCondition();
        for (int perf = 0; perf < nperf; ++perf) {
            perf_cond[perf] = pc[well_cells[perf]];
        }
        const PhaseUsage& pu = fluid_.phaseUsage();
        DataBlock b(nperf, pu.num_phases);
        std::vector<double> rsmax_perf(nperf, 0.0);
        std::vector<double> rvmax_perf(nperf, 0.0);
        if (pu.phase_used[BlackoilPhases::Aqua]) {
            const V bw = fluid_.bWat(avg_press_ad, perf_temp, well_cells).value();
            b.col(pu.phase_pos[BlackoilPhases::Aqua]) = bw;
        }
        assert(active_[Oil]);
        const V perf_so =  subset(state.saturation[pu.phase_pos[Oil]].value(), well_cells);
        if (pu.phase_used[BlackoilPhases::Liquid]) {
            const ADB perf_rs = subset(state.rs, well_cells);
            const V bo = fluid_.bOil(avg_press_ad, perf_temp, perf_rs, perf_cond, well_cells).value();
            b.col(pu.phase_pos[BlackoilPhases::Liquid]) = bo;
            const V rssat = fluidRsSat(avg_press, perf_so, well_cells);
            rsmax_perf.assign(rssat.data(), rssat.data() + nperf);
        }
        if (pu.phase_used[BlackoilPhases::Vapour]) {
            const ADB perf_rv = subset(state.rv, well_cells);
            const V bg = fluid_.bGas(avg_press_ad, perf_temp, perf_rv, perf_cond, well_cells).value();
            b.col(pu.phase_pos[BlackoilPhases::Vapour]) = bg;
            const V rvsat = fluidRvSat(avg_press, perf_so, well_cells);
            rvmax_perf.assign(rvsat.data(), rvsat.data() + nperf);
        }
        // b is row major, so can just copy data.
        std::vector<double> b_perf(b.data(), b.data() + nperf * pu.num_phases);
        // Extract well connection depths.
        const V depth = cellCentroidsZToEigen(grid_);
        const V pdepth = subset(depth, well_cells);
        std::vector<double> perf_depth(pdepth.data(), pdepth.data() + nperf);
        // Surface density.
        std::vector<double> surf_dens(fluid_.surfaceDensity(), fluid_.surfaceDensity() + pu.num_phases);
        // Gravity
        double grav = 0.0;
        const double* g = geo_.gravity();
        const int dim = dimensions(grid_);
        if (g) {
            // Guard against gravity in anything but last dimension.
            for (int dd = 0; dd < dim - 1; ++dd) {
                assert(g[dd] == 0.0);
            }
            grav = g[dim - 1];
        }

        // 2. Compute pressure deltas, and store the results.
        std::vector<double> cdp = WellDensitySegmented
            ::computeConnectionPressureDelta(wells(), xw, fluid_.phaseUsage(),
                                             b_perf, rsmax_perf, rvmax_perf, perf_depth,
                                             surf_dens, grav);
        well_perforation_pressure_diffs_ = Eigen::Map<const V>(cdp.data(), nperf);
    }





    template<class T>
    void
    FullyImplicitBlackoilPolymerSolver<T>::
    assemble(const V&             pvdt,
             const PolymerBlackoilState& x   ,
             const bool initial_assembly,
             WellStateFullyImplicitBlackoil& xw,
             const std::vector<double>& polymer_inflow)
    {
        using namespace Opm::AutoDiffGrid;
        // Create the primary variables.
        SolutionState state = variableState(x, xw);

        if (initial_assembly) {
            // Create the (constant, derivativeless) initial state.
            SolutionState state0 = state;
            makeConstantState(state0);
            // Compute initial accumulation contributions
            // and well connection pressures.
            computeAccum(state0, 0);
            computeWellConnectionPressures(state0, xw);
        }

        // DISKVAL(state.pressure);
        // DISKVAL(state.saturation[0]);
        // DISKVAL(state.saturation[1]);
        // DISKVAL(state.saturation[2]);
        // DISKVAL(state.rs);
        // DISKVAL(state.rv);
        // DISKVAL(state.qs);
        // DISKVAL(state.bhp);

        // -------- Mass balance equations --------

        // Compute b_p and the accumulation term b_p*s_p for each phase,
        // except gas. For gas, we compute b_g*s_g + Rs*b_o*s_o.
        // These quantities are stored in rq_[phase].accum[1].
        // The corresponding accumulation terms from the start of
        // the timestep (b^0_p*s^0_p etc.) were already computed
        // on the initial call to assemble() and stored in rq_[phase].accum[0].
        computeAccum(state, 1);

        // Set up the common parts of the mass balance equations
        // for each active phase.
        const V transi = subset(geo_.transmissibility(), ops_.internal_faces);
        const std::vector<ADB> kr = computeRelPerm(state);
        computeMassFlux(transi, kr, state.canonical_phase_pressures, state);
        for (int phaseIdx = 0; phaseIdx < fluid_.numPhases(); ++phaseIdx) {
            //            computeMassFlux(phaseIdx, transi, kr[canph_[phaseIdx]], state.canonical_phase_pressures[canph_[phaseIdx]], state);
            residual_.material_balance_eq[ phaseIdx ] =
                pvdt*(rq_[phaseIdx].accum[1] - rq_[phaseIdx].accum[0])
                + ops_.div*rq_[phaseIdx].mflux;
        }

        // -------- Extra (optional) rs and rv contributions to the mass balance equations --------

        // Add the extra (flux) terms to the mass balance equations
        // From gas dissolved in the oil phase (rs) and oil vaporized in the gas phase (rv)
        // The extra terms in the accumulation part of the equation are already handled.
        if (active_[ Oil ] && active_[ Gas ]) {
            const int po = fluid_.phaseUsage().phase_pos[ Oil ];
            const int pg = fluid_.phaseUsage().phase_pos[ Gas ];

            const UpwindSelector<double> upwindOil(grid_, ops_,
                                                rq_[po].head.value());
            const ADB rs_face = upwindOil.select(state.rs);

            const UpwindSelector<double> upwindGas(grid_, ops_,
                                                rq_[pg].head.value());
            const ADB rv_face = upwindGas.select(state.rv);

            residual_.material_balance_eq[ pg ] += ops_.div * (rs_face * rq_[po].mflux);
            residual_.material_balance_eq[ po ] += ops_.div * (rv_face * rq_[pg].mflux);

            // DUMP(residual_.material_balance_eq[ Gas ]);

        }

        // Add polymer equation.
        if (has_polymer_) {
            residual_.material_balance_eq[poly_pos_] = pvdt*(rq_[poly_pos_].accum[1] - rq_[poly_pos_].accum[0])
                                               + ops_.div*rq_[poly_pos_].mflux;
        }
        // Note: updateWellControls() can change all its arguments if
        // a well control is switched.
        updateWellControls(state.bhp, state.qs, xw);
        V aliveWells;
        addWellEq(state, xw, aliveWells, polymer_inflow);
        addWellControlEq(state, xw, aliveWells);
    }





    template <class T>
    void FullyImplicitBlackoilPolymerSolver<T>::addWellEq(const SolutionState& state,
                                                          WellStateFullyImplicitBlackoil& xw,
                                                          V& aliveWells,
                                                          const std::vector<double>& polymer_inflow)
    {
        if( ! wellsActive() ) return ;

        const int nc = Opm::AutoDiffGrid::numCells(grid_);
        const int np = wells().number_of_phases;
        const int nw = wells().number_of_wells;
        const int nperf = wells().well_connpos[nw];
        const Opm::PhaseUsage& pu = fluid_.phaseUsage();
        V Tw = Eigen::Map<const V>(wells().WI, nperf);
        const std::vector<int> well_cells(wells().well_cells, wells().well_cells + nperf);

        // pressure diffs computed already (once per step, not changing per iteration)
        const V& cdp = well_perforation_pressure_diffs_;

        // Extract variables for perforation cell pressures
        // and corresponding perforation well pressures.
        const ADB p_perfcell = subset(state.pressure, well_cells);

        // DUMPVAL(p_perfcell);
        // DUMPVAL(state.bhp);
        // DUMPVAL(ADB::constant(cdp));

        // Perforation pressure
        const ADB perfpressure = (wops_.w2p * state.bhp) + cdp;
        std::vector<double> perfpressure_d(perfpressure.value().data(), perfpressure.value().data() + nperf);
        xw.perfPress() = perfpressure_d;

        // Pressure drawdown (also used to determine direction of flow)
        const ADB drawdown =  p_perfcell - perfpressure;

        // current injecting connections
        auto connInjInx = drawdown.value() < 0;

        // injector == 1, producer == 0
        V isInj = V::Zero(nw);
        for (int w = 0; w < nw; ++w) {
            if (wells().type[w] == INJECTOR) {
                isInj[w] = 1;
            }
        }

//        // A cross-flow connection is defined as a connection which has opposite
//        // flow-direction to the well total flow
//         V isInjPerf = (wops_.w2p * isInj);
//        auto crossFlowConns = (connInjInx != isInjPerf);

//        bool allowCrossFlow = true;

//        if (not allowCrossFlow) {
//            auto closedConns = crossFlowConns;
//            for (int c = 0; c < nperf; ++c) {
//                if (closedConns[c]) {
//                    Tw[c] = 0;
//                }
//            }
//            connInjInx = !closedConns;
//        }
//        TODO: not allow for crossflow


        V isInjInx = V::Zero(nperf);
        V isNotInjInx = V::Zero(nperf);
        for (int c = 0; c < nperf; ++c){
            if (connInjInx[c])
                isInjInx[c] = 1;
            else
                isNotInjInx[c] = 1;
        }


        // HANDLE FLOW INTO WELLBORE

        // compute phase volumerates standard conditions
        std::vector<ADB> cq_ps(np, ADB::null());
        for (int phase = 0; phase < np; ++phase) {
            const ADB& wellcell_mob = subset ( rq_[phase].mob, well_cells);
            const ADB cq_p = -(isNotInjInx * Tw) * (wellcell_mob * drawdown);
            cq_ps[phase] = subset(rq_[phase].b,well_cells) * cq_p;
        }
        if (active_[Oil] && active_[Gas]) {
            const int oilpos = pu.phase_pos[Oil];
            const int gaspos = pu.phase_pos[Gas];
            ADB cq_psOil = cq_ps[oilpos];
            ADB cq_psGas = cq_ps[gaspos];
            cq_ps[gaspos] += subset(state.rs,well_cells) * cq_psOil;
            cq_ps[oilpos] += subset(state.rv,well_cells) * cq_psGas;
        }

        // phase rates at std. condtions
        std::vector<ADB> q_ps(np, ADB::null());
        for (int phase = 0; phase < np; ++phase) {
            q_ps[phase] = wops_.p2w * cq_ps[phase];
        }

        // total rates at std
        ADB qt_s = ADB::constant(V::Zero(nw), state.bhp.blockPattern());
        for (int phase = 0; phase < np; ++phase) {
            qt_s += subset(state.qs, Span(nw, 1, phase*nw));
        }

        // compute avg. and total wellbore phase volumetric rates at std. conds
        const DataBlock compi = Eigen::Map<const DataBlock>(wells().comp_frac, nw, np);
        std::vector<ADB> wbq(np, ADB::null());
        ADB wbqt = ADB::constant(V::Zero(nw), state.pressure.blockPattern());
        for (int phase = 0; phase < np; ++phase) {
            const int pos = pu.phase_pos[phase];
            wbq[phase] = (isInj * compi.col(pos)) * qt_s - q_ps[phase];
            wbqt += wbq[phase];
        }
        // DUMPVAL(wbqt);

        // check for dead wells
        aliveWells = V::Constant(nw, 1.0);
        for (int w = 0; w < nw; ++w) {
            if (wbqt.value()[w] == 0) {
                aliveWells[w] = 0.0;
            }
        }
        // compute wellbore mixture at std conds
        Selector<double> notDeadWells_selector(wbqt.value(), Selector<double>::Zero);
        std::vector<ADB> mix_s(np, ADB::null());
        for (int phase = 0; phase < np; ++phase) {
            const int pos = pu.phase_pos[phase];
            mix_s[phase] = notDeadWells_selector.select(ADB::constant(compi.col(pos)), wbq[phase]/wbqt);
        }


        // HANDLE FLOW OUT FROM WELLBORE

        // Total mobilities
        ADB mt = subset(rq_[0].mob,well_cells);
        for (int phase = 1; phase < np; ++phase) {
            mt += subset(rq_[phase].mob,well_cells);
        }

        // DUMPVAL(ADB::constant(isInjInx));
        // DUMPVAL(ADB::constant(Tw));
        // DUMPVAL(mt);
        // DUMPVAL(drawdown);

        // injection connections total volumerates
        ADB cqt_i = -(isInjInx * Tw) * (mt * drawdown);

        // compute volume ratio between connection at standard conditions
        ADB volRat = ADB::constant(V::Zero(nperf), state.pressure.blockPattern());
        std::vector<ADB> cmix_s(np, ADB::null());
        for (int phase = 0; phase < np; ++phase) {
            cmix_s[phase] = wops_.w2p * mix_s[phase];
        }

        ADB well_rv = subset(state.rv,well_cells);
        ADB well_rs = subset(state.rs,well_cells);
        ADB d = V::Constant(nperf,1.0) -  well_rv * well_rs;

        for (int phase = 0; phase < np; ++phase) {
            ADB tmp = cmix_s[phase];

            if (phase == Oil && active_[Gas]) {
                const int gaspos = pu.phase_pos[Gas];
                tmp = tmp - subset(state.rv,well_cells) * cmix_s[gaspos] / d;
            }
            if (phase == Gas && active_[Oil]) {
                const int oilpos = pu.phase_pos[Oil];
                tmp = tmp - subset(state.rs,well_cells) * cmix_s[oilpos] / d;
            }
            volRat += tmp / subset(rq_[phase].b,well_cells);
        }

        // DUMPVAL(cqt_i);
        // DUMPVAL(volRat);

        // injecting connections total volumerates at std cond
        ADB cqt_is = cqt_i/volRat;

        // connection phase volumerates at std cond
        std::vector<ADB> cq_s(np, ADB::null());
        for (int phase = 0; phase < np; ++phase) {
            cq_s[phase] = cq_ps[phase] + (wops_.w2p * mix_s[phase])*cqt_is;
        }

        // DUMPVAL(mix_s[2]);
        // DUMPVAL(cq_ps[2]);

        // Add well contributions to mass balance equations
        for (int phase = 0; phase < np; ++phase) {
            residual_.material_balance_eq[phase] -= superset(cq_s[phase],well_cells,nc);
        }

        
        // Add well contributions to polymer mass balance equation
        if (has_polymer_) {
            const ADB mc = computeMc(state);
            const V polyin = Eigen::Map<const V>(&polymer_inflow[0], nc);
            const V poly_in_perf = subset(polyin, well_cells);
            const V poly_mc_perf = subset(mc, well_cells).value();
            residual_.material_balance_eq[poly_pos_] -= superset(cq_ps[pu.phase_pos[Water]]*poly_mc_perf
                                                     + (wops_.w2p * mix_s[pu.phase_pos[Water]])*cqt_is*poly_in_perf,
                                                     well_cells, nc);
        }


        // Add WELL EQUATIONS
        ADB qs = state.qs;
        for (int phase = 0; phase < np; ++phase) {
            qs -= superset(wops_.p2w * cq_s[phase], Span(nw, 1, phase*nw), nw*np);

        }


        V cq = superset(cq_s[0].value(), Span(nperf, np, 0), nperf*np);
        for (int phase = 1; phase < np; ++phase) {
            cq += superset(cq_s[phase].value(), Span(nperf, np, phase), nperf*np);
        }

        std::vector<double> cq_d(cq.data(), cq.data() + nperf*np);
        xw.perfPhaseRates() = cq_d;

        residual_.well_flux_eq = qs;
    }





    namespace detail
    {
        double rateToCompare(const ADB& well_phase_flow_rate,
                             const int well,
                             const int num_phases,
                             const double* distr)
        {
            const int num_wells = well_phase_flow_rate.size() / num_phases;
            double rate = 0.0;
            for (int phase = 0; phase < num_phases; ++phase) {
                // Important: well_phase_flow_rate is ordered with all rates for first
                // phase coming first, then all for second phase etc.
                rate += well_phase_flow_rate.value()[well + phase*num_wells] * distr[phase];
            }
            return rate;
        }

        bool constraintBroken(const ADB& bhp,
                              const ADB& well_phase_flow_rate,
                              const int well,
                              const int num_phases,
                              const WellType& well_type,
                              const WellControls* wc,
                              const int ctrl_index)
        {
            const WellControlType ctrl_type = well_controls_iget_type(wc, ctrl_index);
            const double target = well_controls_iget_target(wc, ctrl_index);
            const double* distr = well_controls_iget_distr(wc, ctrl_index);

            bool broken = false;

            switch (well_type) {
            case INJECTOR:
            {
                switch (ctrl_type) {
                case BHP:
                    broken = bhp.value()[well] > target;
                    break;

                case RESERVOIR_RATE: // Intentional fall-through
                case SURFACE_RATE:
                    broken = rateToCompare(well_phase_flow_rate,
                                           well, num_phases, distr) > target;
                    break;
                }
            }
            break;

            case PRODUCER:
            {
                switch (ctrl_type) {
                case BHP:
                    broken = bhp.value()[well] < target;
                    break;

                case RESERVOIR_RATE: // Intentional fall-through
                case SURFACE_RATE:
                    // Note that the rates compared below are negative,
                    // so breaking the constraints means: too high flow rate
                    // (as for injection).
                    broken = rateToCompare(well_phase_flow_rate,
                                           well, num_phases, distr) < target;
                    break;
                }
            }
            break;

            default:
                OPM_THROW(std::logic_error, "Can only handle INJECTOR and PRODUCER wells.");
            }

            return broken;
        }
    } // namespace detail




    template<class T>
    void FullyImplicitBlackoilPolymerSolver<T>::updateWellControls(ADB& bhp,
                                                                   ADB& well_phase_flow_rate,
                                                                   WellStateFullyImplicitBlackoil& xw) const
    {
        if( ! wellsActive() ) return ;

        std::string modestring[3] = { "BHP", "RESERVOIR_RATE", "SURFACE_RATE" };
        // Find, for each well, if any constraints are broken. If so,
        // switch control to first broken constraint.
        const int np = wells().number_of_phases;
        const int nw = wells().number_of_wells;
        bool bhp_changed = false;
        bool rates_changed = false;
        for (int w = 0; w < nw; ++w) {
            const WellControls* wc = wells().ctrls[w];
            // The current control in the well state overrides
            // the current control set in the Wells struct, which
            // is instead treated as a default.
            const int current = xw.currentControls()[w];
            // Loop over all controls except the current one, and also
            // skip any RESERVOIR_RATE controls, since we cannot
            // handle those.
            const int nwc = well_controls_get_num(wc);
            int ctrl_index = 0;
            for (; ctrl_index < nwc; ++ctrl_index) {
                if (ctrl_index == current) {
                    // This is the currently used control, so it is
                    // used as an equation. So this is not used as an
                    // inequality constraint, and therefore skipped.
                    continue;
                }
                if (detail::constraintBroken(bhp, well_phase_flow_rate, w, np, wells().type[w], wc, ctrl_index)) {
                    // ctrl_index will be the index of the broken constraint after the loop.
                    break;
                }
            }
            if (ctrl_index != nwc) {
                // Constraint number ctrl_index was broken, switch to it.
                if (terminal_output_)
                {
                    std::cout << "Switching control mode for well " << wells().name[w]
                              << " from " << modestring[well_controls_iget_type(wc, current)]
                              << " to " << modestring[well_controls_iget_type(wc, ctrl_index)] << std::endl;
                }
                xw.currentControls()[w] = ctrl_index;
                // Also updating well state and primary variables.
                // We can only be switching to BHP and SURFACE_RATE
                // controls since we do not support RESERVOIR_RATE.
                const double target = well_controls_iget_target(wc, ctrl_index);
                const double* distr = well_controls_iget_distr(wc, ctrl_index);
                switch (well_controls_iget_type(wc, ctrl_index)) {
                case BHP:
                    xw.bhp()[w] = target;
                    bhp_changed = true;
                    break;

                case RESERVOIR_RATE:
                    // No direct change to any observable quantity at
                    // surface condition.  In this case, use existing
                    // flow rates as initial conditions as reservoir
                    // rate acts only in aggregate.
                    //
                    // Just record the fact that we need to recompute
                    // the 'well_phase_flow_rate'.
                    rates_changed = true;
                    break;

                case SURFACE_RATE:
                    for (int phase = 0; phase < np; ++phase) {
                        if (distr[phase] > 0.0) {
                            xw.wellRates()[np*w + phase] = target * distr[phase];
                        }
                    }
                    rates_changed = true;
                    break;
                }
            }
        }

        // Update primary variables, if necessary.
        if (bhp_changed) {
            ADB::V new_bhp = Eigen::Map<ADB::V>(xw.bhp().data(), nw);
            // Avoiding the copy below would require a value setter method
            // in AutoDiffBlock.
            std::vector<ADB::M> old_derivs = bhp.derivative();
            bhp = ADB::function(std::move(new_bhp), std::move(old_derivs));
        }
        if (rates_changed) {
            // Need to reshuffle well rates, from phase running fastest
            // to wells running fastest.
            // The transpose() below switches the ordering.
            const DataBlock wrates = Eigen::Map<const DataBlock>(xw.wellRates().data(), nw, np).transpose();
            ADB::V new_qs = Eigen::Map<const V>(wrates.data(), nw*np);
            std::vector<ADB::M> old_derivs = well_phase_flow_rate.derivative();
            well_phase_flow_rate = ADB::function(std::move(new_qs), std::move(old_derivs));
        }
    }




    template<class T>
    void FullyImplicitBlackoilPolymerSolver<T>::addWellControlEq(const SolutionState& state,
                                                                 const WellStateFullyImplicitBlackoil& xw,
                                                                 const V& aliveWells)
    {
        if( ! wellsActive() ) return;

        const int np = wells().number_of_phases;
        const int nw = wells().number_of_wells;

        V bhp_targets  = V::Zero(nw);
        V rate_targets = V::Zero(nw);
        M rate_distr(nw, np*nw);
        for (int w = 0; w < nw; ++w) {
            const WellControls* wc = wells().ctrls[w];
            // The current control in the well state overrides
            // the current control set in the Wells struct, which
            // is instead treated as a default.
            const int current = xw.currentControls()[w];

            switch (well_controls_iget_type(wc, current)) {
            case BHP:
            {
                bhp_targets (w) = well_controls_iget_target(wc, current);
                rate_targets(w) = -1e100;
            }
            break;

            case RESERVOIR_RATE: // Intentional fall-through
            case SURFACE_RATE:
            {
                // RESERVOIR and SURFACE rates look the same, from a
                // high-level point of view, in the system of
                // simultaneous linear equations.

                const double* const distr =
                    well_controls_iget_distr(wc, current);

                for (int p = 0; p < np; ++p) {
                    rate_distr.insert(w, p*nw + w) = distr[p];
                }

                bhp_targets (w) = -1.0e100;
                rate_targets(w) = well_controls_iget_target(wc, current);
            }
            break;
            }
        }
        const ADB bhp_residual = state.bhp - bhp_targets;
        const ADB rate_residual = rate_distr * state.qs - rate_targets;
        // Choose bhp residual for positive bhp targets.
        Selector<double> bhp_selector(bhp_targets);
        residual_.well_eq = bhp_selector.select(bhp_residual, rate_residual);
        // For wells that are dead (not flowing), and therefore not communicating
        // with the reservoir, we set the equation to be equal to the well's total
        // flow. This will be a solution only if the target rate is also zero.
        M rate_summer(nw, np*nw);
        for (int w = 0; w < nw; ++w) {
            for (int phase = 0; phase < np; ++phase) {
                rate_summer.insert(w, phase*nw + w) = 1.0;
            }
        }
        Selector<double> alive_selector(aliveWells, Selector<double>::NotEqualZero);
        residual_.well_eq = alive_selector.select(residual_.well_eq, rate_summer * state.qs);
        // DUMP(residual_.well_eq);
    }





    template<class T>
    V FullyImplicitBlackoilPolymerSolver<T>::solveJacobianSystem() const
    {
        return linsolver_.computeNewtonIncrement(residual_);
    }





    namespace detail
    {

        double infinityNorm( const ADB& a )
        {
            if( a.value().size() > 0 ) {
                return a.value().matrix().lpNorm<Eigen::Infinity> ();
            }
            else { // this situation can occur when no wells are present
                return 0.0;
            }
        }

    } // namespace detail





    template<class T>
    void FullyImplicitBlackoilPolymerSolver<T>::updateState(const V& dx,
                                                            PolymerBlackoilState& state,
                                                            WellStateFullyImplicitBlackoil& well_state)
    {
        using namespace Opm::AutoDiffGrid;
        const int np = fluid_.numPhases();
        const int nc = numCells(grid_);
        const int nw = wellsActive() ? wells().number_of_wells : 0;
        const V null;
        assert(null.size() == 0);
        const V zero = V::Zero(nc);

        // store cell status in vectors
        V isRs = V::Zero(nc,1);
        V isRv = V::Zero(nc,1);
        V isSg = V::Zero(nc,1);
        if (active_[Gas]) {
            for (int c = 0; c < nc; ++c) {
                switch (primalVariable_[c]) {
                case PrimalVariables::RS:
                    isRs[c] = 1;
                    break;

                case PrimalVariables::RV:
                    isRv[c] = 1;
                    break;

                default:
                    isSg[c] = 1;
                    break;
                }
            }
        }

        // Extract parts of dx corresponding to each part.
        const V dp = subset(dx, Span(nc));
        int varstart = nc;
        const V dsw = active_[Water] ? subset(dx, Span(nc, 1, varstart)) : null;
        varstart += dsw.size();

        const V dxvar = active_[Gas] ? subset(dx, Span(nc, 1, varstart)): null;
        varstart += dxvar.size();
        
        const V dc = (has_polymer_) ? subset(dx, Span(nc, 1, varstart)) : null;
        varstart += dc.size();
        
        const V dqs = subset(dx, Span(np*nw, 1, varstart));
        varstart += dqs.size();
        const V dbhp = subset(dx, Span(nw, 1, varstart));
        varstart += dbhp.size();
        assert(varstart == dx.size());

        // Pressure update.
        const double dpmaxrel = dpMaxRel();
        const V p_old = Eigen::Map<const V>(&state.pressure()[0], nc, 1);
        const V absdpmax = dpmaxrel*p_old.abs();
        const V dp_limited = sign(dp) * dp.abs().min(absdpmax);
        const V p = (p_old - dp_limited).max(zero);
        std::copy(&p[0], &p[0] + nc, state.pressure().begin());


        // Saturation updates.
        const Opm::PhaseUsage& pu = fluid_.phaseUsage();
        const DataBlock s_old = Eigen::Map<const DataBlock>(& state.saturation()[0], nc, np);
        const double dsmax = dsMax();

        V so;
        V sw;
        V sg;

        {
            V maxVal = zero;
            V dso = zero;
            if (active_[Water]){
                maxVal = dsw.abs().max(maxVal);
                dso = dso - dsw;
            }

            V dsg;
            if (active_[Gas]){
                dsg = isSg * dxvar - isRv * dsw;
                maxVal = dsg.abs().max(maxVal);
                dso = dso - dsg;
            }

            maxVal = dso.abs().max(maxVal);

            V step = dsmax/maxVal;
            step = step.min(1.);

            if (active_[Water]) {
                const int pos = pu.phase_pos[ Water ];
                const V sw_old = s_old.col(pos);
                sw = sw_old - step * dsw;
            }

            if (active_[Gas]) {
                const int pos = pu.phase_pos[ Gas ];
                const V sg_old = s_old.col(pos);
                sg = sg_old - step * dsg;
            }

            const int pos = pu.phase_pos[ Oil ];
            const V so_old = s_old.col(pos);
            so = so_old - step * dso;
        }

        // Appleyard chop process.
        auto ixg = sg < 0;
        for (int c = 0; c < nc; ++c) {
            if (ixg[c]) {
                sw[c] = sw[c] / (1-sg[c]);
                so[c] = so[c] / (1-sg[c]);
                sg[c] = 0;
            }
        }


        auto ixo = so < 0;
        for (int c = 0; c < nc; ++c) {
            if (ixo[c]) {
                sw[c] = sw[c] / (1-so[c]);
                sg[c] = sg[c] / (1-so[c]);
                so[c] = 0;
            }
        }

        auto ixw = sw < 0;
        for (int c = 0; c < nc; ++c) {
            if (ixw[c]) {
                so[c] = so[c] / (1-sw[c]);
                sg[c] = sg[c] / (1-sw[c]);
                sw[c] = 0;
            }
        }

        const V sumSat = sw + so + sg;
        sw = sw / sumSat;
        so = so / sumSat;
        sg = sg / sumSat;

        // Update the state
        for (int c = 0; c < nc; ++c) {
            state.saturation()[c*np + pu.phase_pos[ Water ]] = sw[c];
        }

        for (int c = 0; c < nc; ++c) {
            state.saturation()[c*np + pu.phase_pos[ Gas ]] = sg[c];
        }

        if (active_[ Oil ]) {
            const int pos = pu.phase_pos[ Oil ];
            for (int c = 0; c < nc; ++c) {
                state.saturation()[c*np + pos] = so[c];
            }
        }

        // Update rs and rv
        const double drmaxrel = drMaxRel();
        V rs;
        if (has_disgas_) {
            const V rs_old = Eigen::Map<const V>(&state.gasoilratio()[0], nc);
            const V drs = isRs * dxvar;
            const V drs_limited = sign(drs) * drs.abs().min(rs_old.abs()*drmaxrel);
            rs = rs_old - drs_limited;
        }
        V rv;
        if (has_vapoil_) {
            const V rv_old = Eigen::Map<const V>(&state.rv()[0], nc);
            const V drv = isRv * dxvar;
            const V drv_limited = sign(drv) * drv.abs().min(rv_old.abs()*drmaxrel);
            rv = rv_old - drv_limited;
        }


        // Sg is used as primal variable for water only cells.
        const double epsilon = std::sqrt(std::numeric_limits<double>::epsilon());
        auto watOnly = sw >  (1 - epsilon);

        // phase translation sg <-> rs
        std::fill(primalVariable_.begin(), primalVariable_.end(), PrimalVariables::Sg);

        if (has_disgas_) {
            const V rsSat0 = fluidRsSat(p_old, s_old.col(pu.phase_pos[Oil]), cells_);
            const V rsSat = fluidRsSat(p, so, cells_);
            // The obvious case
            auto hasGas = (sg > 0 && isRs == 0);

            // Set oil saturated if previous rs is sufficiently large
            const V rs_old = Eigen::Map<const V>(&state.gasoilratio()[0], nc);
            auto gasVaporized =  ( (rs > rsSat * (1+epsilon) && isRs == 1 ) && (rs_old > rsSat0 * (1-epsilon)) );
            auto useSg = watOnly || hasGas || gasVaporized;
            for (int c = 0; c < nc; ++c) {
                if (useSg[c]) {
                    rs[c] = rsSat[c];
                } else {
                    primalVariable_[c] = PrimalVariables::RS;
                }
            }

        }

        // phase transitions so <-> rv
        if (has_vapoil_) {

            // The gas pressure is needed for the rvSat calculations
            const V gaspress_old = computeGasPressure(p_old, s_old.col(Water), s_old.col(Oil), s_old.col(Gas));
            const V gaspress = computeGasPressure(p, sw, so, sg);
            const V rvSat0 = fluidRvSat(gaspress_old, s_old.col(pu.phase_pos[Oil]), cells_);
            const V rvSat = fluidRvSat(gaspress, so, cells_);

            // The obvious case
            auto hasOil = (so > 0 && isRv == 0);

            // Set oil saturated if previous rv is sufficiently large
            const V rv_old = Eigen::Map<const V>(&state.rv()[0], nc);
            auto oilCondensed = ( (rv > rvSat * (1+epsilon) && isRv == 1) && (rv_old > rvSat0 * (1-epsilon)) );
            auto useSg = watOnly || hasOil || oilCondensed;
            for (int c = 0; c < nc; ++c) {
                if (useSg[c]) {
                    rv[c] = rvSat[c];
                } else {
                    primalVariable_[c] = PrimalVariables::RV;
                }
            }

        }

        // Update the state
        if (has_disgas_) {
            std::copy(&rs[0], &rs[0] + nc, state.gasoilratio().begin());
        }

        if (has_vapoil_) {
            std::copy(&rv[0], &rv[0] + nc, state.rv().begin());
        }

        //Polymer concentration updates.
        if (has_polymer_) {
            const V c_old = Eigen::Map<const V>(&state.concentration()[0], nc, 1);
            const V c = (c_old - dc).max(zero);
            std::copy(&c[0], &c[0] + nc, state.concentration().begin());
        }
        
        if( wellsActive() )
        {
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
            const V dbhp_limited = sign(dbhp) * dbhp.abs().min(bhp_old.abs()*dpmaxrel);
            const V bhp = bhp_old - dbhp_limited;
            std::copy(&bhp[0], &bhp[0] + bhp.size(), well_state.bhp().begin());
        }

        // Update phase conditions used for property calculations.
        updatePhaseCondFromPrimalVariable();
    }





    template<class T>
    std::vector<ADB>
    FullyImplicitBlackoilPolymerSolver<T>::computeRelPerm(const SolutionState& state) const
    {
        using namespace Opm::AutoDiffGrid;
        const int               nc   = numCells(grid_);
        const std::vector<int>& bpat = state.pressure.blockPattern();

        const ADB null = ADB::constant(V::Zero(nc, 1), bpat);

        const Opm::PhaseUsage& pu = fluid_.phaseUsage();
        const ADB sw = (active_[ Water ]
                        ? state.saturation[ pu.phase_pos[ Water ] ]
                        : null);

        const ADB so = (active_[ Oil ]
                        ? state.saturation[ pu.phase_pos[ Oil ] ]
                        : null);

        const ADB sg = (active_[ Gas ]
                        ? state.saturation[ pu.phase_pos[ Gas ] ]
                        : null);

        return fluid_.relperm(sw, so, sg, cells_);
    }





    template<class T>
    std::vector<ADB>
    FullyImplicitBlackoilPolymerSolver<T>::computePressures(const SolutionState& state) const
    {
        using namespace Opm::AutoDiffGrid;
        const int               nc   = numCells(grid_);
        const std::vector<int>& bpat = state.pressure.blockPattern();

        const ADB null = ADB::constant(V::Zero(nc, 1), bpat);

        const Opm::PhaseUsage& pu = fluid_.phaseUsage();
        const ADB& sw = (active_[ Water ]
                        ? state.saturation[ pu.phase_pos[ Water ] ]
                        : null);

        const ADB& so = (active_[ Oil ]
                        ? state.saturation[ pu.phase_pos[ Oil ] ]
                        : null);

        const ADB& sg = (active_[ Gas ]
                        ? state.saturation[ pu.phase_pos[ Gas ] ]
                        : null);
        return computePressures(state.pressure, sw, so, sg);


    }





    template <class T>
    std::vector<ADB>
    FullyImplicitBlackoilPolymerSolver<T>::
    computePressures(const ADB& po,
                     const ADB& sw,
                     const ADB& so,
                     const ADB& sg) const
    {
        // convert the pressure offsets to the capillary pressures
        std::vector<ADB> pressure = fluid_.capPress(sw, so, sg, cells_);
        for (int phaseIdx = 0; phaseIdx < BlackoilPhases::MaxNumPhases; ++phaseIdx) {
            // The reference pressure is always the liquid phase (oil) pressure.
            if (phaseIdx == BlackoilPhases::Liquid)
                continue;
            pressure[phaseIdx] = pressure[phaseIdx] - pressure[BlackoilPhases::Liquid];
        }

        // Since pcow = po - pw, but pcog = pg - po,
        // we have
        //   pw = po - pcow
        //   pg = po + pcgo
        // This is an unfortunate inconsistency, but a convention we must handle.
        for (int phaseIdx = 0; phaseIdx < BlackoilPhases::MaxNumPhases; ++phaseIdx) {
            if (phaseIdx == BlackoilPhases::Aqua) {
                pressure[phaseIdx] = po - pressure[phaseIdx];
            } else {
                pressure[phaseIdx] += po;
            }
        }

        return pressure;
    }





    template<class T>
    V
    FullyImplicitBlackoilPolymerSolver<T>::computeGasPressure(const V& po,
                                                       const V& sw,
                                                       const V& so,
                                                       const V& sg) const
    {
        assert (active_[Gas]);
        std::vector<ADB> cp = fluid_.capPress(ADB::constant(sw),
                                              ADB::constant(so),
                                              ADB::constant(sg),
                                              cells_);
        return cp[Gas].value() + po;
    }





    template<class T>
    std::vector<ADB>
    FullyImplicitBlackoilPolymerSolver<T>::computeRelPermWells(const SolutionState& state,
                                                               const DataBlock& well_s,
                                                               const std::vector<int>& well_cells) const
    {
        const int nw = wells().number_of_wells;
        const int nperf = wells().well_connpos[nw];
        const std::vector<int>& bpat = state.pressure.blockPattern();

        const ADB null = ADB::constant(V::Zero(nperf), bpat);

        const Opm::PhaseUsage& pu = fluid_.phaseUsage();
        const ADB sw = (active_[ Water ]
                        ? ADB::constant(well_s.col(pu.phase_pos[ Water ]), bpat)
                        : null);

        const ADB so = (active_[ Oil ]
                        ? ADB::constant(well_s.col(pu.phase_pos[ Oil ]), bpat)
                        : null);

        const ADB sg = (active_[ Gas ]
                        ? ADB::constant(well_s.col(pu.phase_pos[ Gas ]), bpat)
                        : null);

        return fluid_.relperm(sw, so, sg, well_cells);
    }





    template<class T>
    void
    FullyImplicitBlackoilPolymerSolver<T>::computeMassFlux(const V& transi,
                                                           const std::vector<ADB>& kr    ,
                                                           const std::vector<ADB>& phasePressure,
                                                           const SolutionState&    state)
    {
        for (int phase = 0; phase < fluid_.numPhases(); ++phase) {
            const int canonicalPhaseIdx = canph_[phase];
            const std::vector<PhasePresence> cond = phaseCondition();

            const ADB tr_mult = transMult(state.pressure);
            const ADB mu    = fluidViscosity(canonicalPhaseIdx, phasePressure[canonicalPhaseIdx], state.temperature, state.rs, state.rv,cond, cells_);

            rq_[phase].mob = tr_mult * kr[canonicalPhaseIdx] / mu;

            const ADB rho   = fluidDensity(canonicalPhaseIdx, phasePressure[canonicalPhaseIdx], state.temperature, state.rs, state.rv,cond, cells_);

            ADB& head = rq_[phase].head;

            // compute gravity potensial using the face average as in eclipse and MRST
            const ADB rhoavg = ops_.caver * rho;

            ADB dp = ops_.ngrad * phasePressure[canonicalPhaseIdx] - geo_.gravity()[2] * (rhoavg * (ops_.ngrad * geo_.z().matrix()));

            if (use_threshold_pressure_) {
                applyThresholdPressures(dp);
            }

            head = transi*dp;
            if (canonicalPhaseIdx == Water) {
                if(has_polymer_) {
                    const ADB cmax = ADB::constant(cmax_, state.concentration.blockPattern());
                    const ADB mc = computeMc(state);
                    ADB krw_eff = polymer_props_ad_.effectiveRelPerm(state.concentration,
                                                                     cmax,
                                                                     kr[canonicalPhaseIdx],
                                                                     state.saturation[canonicalPhaseIdx]);
                    ADB inv_wat_eff_visc = polymer_props_ad_.effectiveInvWaterVisc(state.concentration, mu.value().data());
                    rq_[phase].mob = tr_mult * krw_eff * inv_wat_eff_visc;
                    rq_[poly_pos_].mob = tr_mult * mc * krw_eff * inv_wat_eff_visc; 
                    rq_[poly_pos_].b = rq_[phase].b;
                    rq_[poly_pos_].head = rq_[phase].head;
                    UpwindSelector<double> upwind(grid_, ops_, rq_[poly_pos_].head.value());
                    rq_[poly_pos_].mflux = upwind.select(rq_[poly_pos_].b * rq_[poly_pos_].mob) * rq_[poly_pos_].head;
                }
            }
            //head      = transi*(ops_.ngrad * phasePressure) + gflux;

            UpwindSelector<double> upwind(grid_, ops_, head.value());

            const ADB& b       = rq_[phase].b;
            const ADB& mob     = rq_[phase].mob;
            rq_[phase].mflux = upwind.select(b * mob) * head;
        }       
    }





    template<class T>
    void
    FullyImplicitBlackoilPolymerSolver<T>::applyThresholdPressures(ADB& dp)
    {
        // We support reversible threshold pressures only.
        // Method: if the potential difference is lower (in absolute
        // value) than the threshold for any face, then the potential
        // (and derivatives) is set to zero. If it is above the
        // threshold, the threshold pressure is subtracted from the
        // absolute potential (the potential is moved towards zero).

        // Identify the set of faces where the potential is under the
        // threshold, that shall have zero flow. Storing the bool
        // Array as a V (a double Array) with 1 and 0 elements, a
        // 1 where flow is allowed, a 0 where it is not.
        const V high_potential = (dp.value().abs() >= threshold_pressures_by_interior_face_).template cast<double>();

        // Create a sparse vector that nullifies the low potential elements.
        const M keep_high_potential = spdiag(high_potential);

        // Find the current sign for the threshold modification
        const V sign_dp = sign(dp.value());
        const V threshold_modification = sign_dp * threshold_pressures_by_interior_face_;

        // Modify potential and nullify where appropriate.
        dp = keep_high_potential * (dp - threshold_modification);
    }





    template<class T>
    double
    FullyImplicitBlackoilPolymerSolver<T>::residualNorm() const
    {
        double globalNorm = 0;
        std::vector<ADB>::const_iterator quantityIt = residual_.material_balance_eq.begin();
        const std::vector<ADB>::const_iterator endQuantityIt = residual_.material_balance_eq.end();
        for (; quantityIt != endQuantityIt; ++quantityIt) {
            const double quantityResid = (*quantityIt).value().matrix().norm();
            if (!std::isfinite(quantityResid)) {
                const int trouble_phase = quantityIt - residual_.material_balance_eq.begin();
                OPM_THROW(Opm::NumericalProblem,
                          "Encountered a non-finite residual in material balance equation "
                          << trouble_phase);
            }
            globalNorm = std::max(globalNorm, quantityResid);
        }
        globalNorm = std::max(globalNorm, residual_.well_flux_eq.value().matrix().norm());
        globalNorm = std::max(globalNorm, residual_.well_eq.value().matrix().norm());

        return globalNorm;
    }


    template<class T>
    std::vector<double>
    FullyImplicitBlackoilPolymerSolver<T>::computeResidualNorms() const
    {
        std::vector<double> residualNorms;

        std::vector<ADB>::const_iterator massBalanceIt = residual_.material_balance_eq.begin();
        const std::vector<ADB>::const_iterator endMassBalanceIt = residual_.material_balance_eq.end();

        for (; massBalanceIt != endMassBalanceIt; ++massBalanceIt) {
            const double massBalanceResid = detail::infinityNorm( (*massBalanceIt) );
            if (!std::isfinite(massBalanceResid)) {
                OPM_THROW(Opm::NumericalProblem,
                          "Encountered a non-finite residual");
            }
            residualNorms.push_back(massBalanceResid);
        }

        // the following residuals are not used in the oscillation detection now
        const double wellFluxResid = detail::infinityNorm( residual_.well_flux_eq );
        if (!std::isfinite(wellFluxResid)) {
            OPM_THROW(Opm::NumericalProblem,
               "Encountered a non-finite residual");
        }
        residualNorms.push_back(wellFluxResid);

        const double wellResid = detail::infinityNorm( residual_.well_eq );
        if (!std::isfinite(wellResid)) {
           OPM_THROW(Opm::NumericalProblem,
               "Encountered a non-finite residual");
        }
        residualNorms.push_back(wellResid);

        return residualNorms;
    }

    template<class T>
    void
    FullyImplicitBlackoilPolymerSolver<T>::detectNewtonOscillations(const std::vector<std::vector<double>>& residual_history,
                                                                    const int it, const double relaxRelTol,
                                                                    bool& oscillate, bool& stagnate) const
    {
        // The detection of oscillation in two primary variable results in the report of the detection
        // of oscillation for the solver.
        // Only the saturations are used for oscillation detection for the black oil model.
        // Stagnate is not used for any treatment here.

        if ( it < 2 ) {
            oscillate = false;
            stagnate = false;
            return;
        }

        stagnate = true;
        int oscillatePhase = 0;
        const std::vector<double>& F0 = residual_history[it];
        const std::vector<double>& F1 = residual_history[it - 1];
        const std::vector<double>& F2 = residual_history[it - 2];
        for (int p= 0; p < fluid_.numPhases(); ++p){
            const double d1 = std::abs((F0[p] - F2[p]) / F0[p]);
            const double d2 = std::abs((F0[p] - F1[p]) / F0[p]);

            oscillatePhase += (d1 < relaxRelTol) && (relaxRelTol < d2);

            // Process is 'stagnate' unless at least one phase
            // exhibits significant residual change.
            stagnate = (stagnate && !(std::abs((F1[p] - F2[p]) / F2[p]) > 1.0e-3));
        }

        oscillate = (oscillatePhase > 1);
    }


    template<class T>
    void
    FullyImplicitBlackoilPolymerSolver<T>::stablizeNewton(V& dx, V& dxOld, const double omega,
                                                          const RelaxType relax_type) const
    {
        // The dxOld is updated with dx.
        // If omega is equal to 1., no relaxtion will be appiled.

        const V tempDxOld = dxOld;
        dxOld = dx;

        switch (relax_type) {
            case DAMPEN:
                if (omega == 1.) {
                    return;
                }
                dx = dx*omega;
                return;
            case SOR:
                if (omega == 1.) {
                    return;
                }
                dx = dx*omega + (1.-omega)*tempDxOld;
                return;
            default:
                OPM_THROW(std::runtime_error, "Can only handle DAMPEN and SOR relaxation type.");
        }

        return;
    }

    template<class T>
    double
    FullyImplicitBlackoilPolymerSolver<T>::convergenceReduction(const Eigen::Array<double, Eigen::Dynamic, MaxNumPhases+1>& B,
                                                                const Eigen::Array<double, Eigen::Dynamic, MaxNumPhases+1>& tempV,
                                                                const Eigen::Array<double, Eigen::Dynamic, MaxNumPhases+1>& R,
                                                                std::array<double,MaxNumPhases+1>& R_sum,
                                                                std::array<double,MaxNumPhases+1>& maxCoeff,
                                                                std::array<double,MaxNumPhases+1>& B_avg,
                                                                int nc) const
    {
        // Do the global reductions
#if HAVE_MPI
        if ( linsolver_.parallelInformation().type() == typeid(ParallelISTLInformation) )
        {
            const ParallelISTLInformation& info =
                boost::any_cast<const ParallelISTLInformation&>(linsolver_.parallelInformation());

            // Compute the global number of cells and porevolume
            std::vector<int> v(nc, 1);
            auto nc_and_pv = std::tuple<int, double>(0, 0.0);
            auto nc_and_pv_operators = std::make_tuple(Opm::Reduction::makeGlobalSumFunctor<int>(),
                                                        Opm::Reduction::makeGlobalSumFunctor<double>());
            auto nc_and_pv_containers  = std::make_tuple(v, geo_.poreVolume());
            info.computeReduction(nc_and_pv_containers, nc_and_pv_operators, nc_and_pv);

            for ( int idx=0; idx<MaxNumPhases; ++idx )
            {
                if (active_[idx]) {
                    auto values     = std::tuple<double,double,double>(0.0 ,0.0 ,0.0);
                    auto containers = std::make_tuple(B.col(idx),
                                                      tempV.col(idx),
                                                      R.col(idx));
                    auto operators  = std::make_tuple(Opm::Reduction::makeGlobalSumFunctor<double>(),
                                                      Opm::Reduction::makeGlobalMaxFunctor<double>(),
                                                      Opm::Reduction::makeGlobalSumFunctor<double>());
                    info.computeReduction(containers, operators, values);
                    B_avg[idx]    = std::get<0>(values)/std::get<0>(nc_and_pv);
                    maxCoeff[idx] = std::get<1>(values);
                    R_sum[idx]    = std::get<2>(values);
                }
                else
                {
                    R_sum[idx] = B_avg[idx] = maxCoeff[idx] = 0.0;
                }
            }
            // Compute pore volume
            return std::get<1>(nc_and_pv);
        }
        else
#endif
        {
            for ( int idx=0; idx<MaxNumPhases; ++idx )
            {
                if (active_[idx]) {
                    B_avg[idx] = B.col(idx).sum()/nc;
                    maxCoeff[idx]=tempV.col(idx).maxCoeff();
                    R_sum[idx] = R.col(idx).sum();
                }
                else
                {
                    R_sum[idx] = B_avg[idx] = maxCoeff[idx] =0.0;
                }
            }
            if (has_polymer_) {
                B_avg[MaxNumPhases] = B.col(MaxNumPhases).sum()/nc;
                maxCoeff[MaxNumPhases]=tempV.col(MaxNumPhases).maxCoeff();
                R_sum[MaxNumPhases] = R.col(MaxNumPhases).sum();
            }
            // Compute total pore volume
            return geo_.poreVolume().sum();
        }
    }

    template<class T>
    bool
    FullyImplicitBlackoilPolymerSolver<T>::getConvergence(const double dt, const int iteration)
    {
        const double tol_mb    = param_.tolerance_mb_;
        const double tol_cnv   = param_.tolerance_cnv_;
        const double tol_wells = param_.tolerance_wells_;

        const int nc = Opm::AutoDiffGrid::numCells(grid_);
        const Opm::PhaseUsage& pu = fluid_.phaseUsage();

        const V pv = geo_.poreVolume();

        const std::vector<PhasePresence> cond = phaseCondition();

        std::array<double,MaxNumPhases+1> CNV                   = {{0., 0., 0., 0.}};
        std::array<double,MaxNumPhases+1> R_sum                 = {{0., 0., 0., 0.}};
        std::array<double,MaxNumPhases+1> B_avg                 = {{0., 0., 0., 0.}};
        std::array<double,MaxNumPhases+1> maxCoeff              = {{0., 0., 0., 0.}};
        std::array<double,MaxNumPhases+1> mass_balance_residual = {{0., 0., 0., 0.}};
        std::size_t cols = MaxNumPhases+1; // needed to pass the correct type to Eigen
        Eigen::Array<V::Scalar, Eigen::Dynamic, MaxNumPhases+1> B(nc, cols);
        Eigen::Array<V::Scalar, Eigen::Dynamic, MaxNumPhases+1> R(nc, cols);
        Eigen::Array<V::Scalar, Eigen::Dynamic, MaxNumPhases+1> tempV(nc, cols);

        for ( int idx=0; idx<MaxNumPhases; ++idx )
        {
            if (active_[idx]) {
                const int pos    = pu.phase_pos[idx];
                const ADB& tempB = rq_[pos].b;
                B.col(idx)       = 1./tempB.value();
                R.col(idx)       = residual_.material_balance_eq[idx].value();
                tempV.col(idx)   = R.col(idx).abs()/pv;
            }
        }
        if (has_polymer_) {
            const ADB& tempB = rq_[poly_pos_].b;
            B.col(MaxNumPhases) = 1. / tempB.value();
            R.col(MaxNumPhases) = residual_.material_balance_eq[poly_pos_].value();
            tempV.col(MaxNumPhases) = R.col(MaxNumPhases).abs()/pv;
        }

        const double pvSum = convergenceReduction(B, tempV, R, R_sum, maxCoeff, B_avg, nc);

        bool converged_MB = true;
        bool converged_CNV = true;
        // Finish computation
        for ( int idx=0; idx<MaxNumPhases+1; ++idx )
        {
            CNV[idx]                   = B_avg[idx] * dt * maxCoeff[idx];
            mass_balance_residual[idx] = std::abs(B_avg[idx]*R_sum[idx]) * dt / pvSum;
            converged_MB               = converged_MB && (mass_balance_residual[idx] < tol_mb);
            converged_CNV              = converged_CNV && (CNV[idx] < tol_cnv);
        }

        const double residualWellFlux = detail::infinityNorm(residual_.well_flux_eq);
        const double residualWell     = detail::infinityNorm(residual_.well_eq);
        const bool   converged_Well   = (residualWellFlux < tol_wells) && (residualWell < Opm::unit::barsa);
        const bool   converged        = converged_MB && converged_CNV && converged_Well;

        // if one of the residuals is NaN, throw exception, so that the solver can be restarted
        if ( std::isnan(mass_balance_residual[Water]) || mass_balance_residual[Water] > maxResidualAllowed() ||
            std::isnan(mass_balance_residual[Oil])   || mass_balance_residual[Oil]   > maxResidualAllowed() ||
            std::isnan(mass_balance_residual[Gas])   || mass_balance_residual[Gas]   > maxResidualAllowed() ||
            std::isnan(mass_balance_residual[MaxNumPhases])   || mass_balance_residual[MaxNumPhases]   > maxResidualAllowed() ||             
            std::isnan(CNV[Water]) || CNV[Water] > maxResidualAllowed() ||
            std::isnan(CNV[Oil]) || CNV[Oil] > maxResidualAllowed() ||
            std::isnan(CNV[Gas]) || CNV[Gas] > maxResidualAllowed() ||
            std::isnan(CNV[MaxNumPhases]) || CNV[MaxNumPhases] > maxResidualAllowed() ||
            std::isnan(residualWellFlux) || residualWellFlux > maxResidualAllowed() ||
            std::isnan(residualWell)     || residualWell     > maxResidualAllowed() )
        {
            OPM_THROW(Opm::NumericalProblem,"One of the residuals is NaN or to large!");
        }

        if ( terminal_output_ )
        {
            // Only rank 0 does print to std::cout
            if (iteration == 0) {
                std::cout << "\nIter    MB(OIL)  MB(WATER)    MB(GAS)       CNVW       CNVO       CNVG    CNVP   WELL-FLOW WELL-CNTRL\n";
            }
            const std::streamsize oprec = std::cout.precision(3);
            const std::ios::fmtflags oflags = std::cout.setf(std::ios::scientific);
            std::cout << std::setw(4) << iteration
                      << std::setw(11) << mass_balance_residual[Water]
                      << std::setw(11) << mass_balance_residual[Oil]
                      << std::setw(11) << mass_balance_residual[Gas]
                      << std::setw(11) << mass_balance_residual[MaxNumPhases]
                      << std::setw(11) << CNV[Water]
                      << std::setw(11) << CNV[Oil]
                      << std::setw(11) << CNV[Gas]
                      << std::setw(11) << CNV[MaxNumPhases]
                      << std::setw(11) << residualWellFlux
                      << std::setw(11) << residualWell
                      << std::endl;
            std::cout.precision(oprec);
            std::cout.flags(oflags);
        }
        return converged;
    }


    template<class T>
    ADB
    FullyImplicitBlackoilPolymerSolver<T>::fluidViscosity(const int               phase,
                                                          const ADB&              p    ,
                                                          const ADB&              temp ,
                                                          const ADB&              rs   ,
                                                          const ADB&              rv   ,
                                                          const std::vector<PhasePresence>& cond,
                                                          const std::vector<int>& cells) const
    {
        switch (phase) {
        case Water:
            return fluid_.muWat(p, temp, cells);
        case Oil: {
            return fluid_.muOil(p, temp, rs, cond, cells);
        }
        case Gas:
            return fluid_.muGas(p, temp, rv, cond, cells);
        default:
            OPM_THROW(std::runtime_error, "Unknown phase index " << phase);
        }
    }





    template<class T>
    ADB
    FullyImplicitBlackoilPolymerSolver<T>::fluidReciprocFVF(const int               phase,
                                                            const ADB&              p    ,
                                                            const ADB&              temp ,
                                                            const ADB&              rs   ,
                                                            const ADB&              rv   ,
                                                            const std::vector<PhasePresence>& cond,
                                                            const std::vector<int>& cells) const
    {
        switch (phase) {
        case Water:
            return fluid_.bWat(p, temp, cells);
        case Oil: {
            return fluid_.bOil(p, temp, rs, cond, cells);
        }
        case Gas:
            return fluid_.bGas(p, temp, rv, cond, cells);
        default:
            OPM_THROW(std::runtime_error, "Unknown phase index " << phase);
        }
    }





    template<class T>
    ADB
    FullyImplicitBlackoilPolymerSolver<T>::fluidDensity(const int               phase,
                                                        const ADB&              p    ,
                                                        const ADB&              temp ,
                                                        const ADB&              rs   ,
                                                        const ADB&              rv   ,
                                                        const std::vector<PhasePresence>& cond,
                                                        const std::vector<int>& cells) const
    {
        const double* rhos = fluid_.surfaceDensity();
        ADB b = fluidReciprocFVF(phase, p, temp, rs, rv, cond, cells);
        ADB rho = V::Constant(p.size(), 1, rhos[phase]) * b;
        if (phase == Oil && active_[Gas]) {
            // It is correct to index into rhos with canonical phase indices.
            rho += V::Constant(p.size(), 1, rhos[Gas]) * rs * b;
        }
        if (phase == Gas && active_[Oil]) {
            // It is correct to index into rhos with canonical phase indices.
            rho += V::Constant(p.size(), 1, rhos[Oil]) * rv * b;
        }
        return rho;
    }





    template<class T>
    V
    FullyImplicitBlackoilPolymerSolver<T>::fluidRsSat(const V&                p,
                                                      const V&                satOil,
                                                      const std::vector<int>& cells) const
    {
        return fluid_.rsSat(ADB::constant(p), ADB::constant(satOil), cells).value();
    }





    template<class T>
    ADB
    FullyImplicitBlackoilPolymerSolver<T>::fluidRsSat(const ADB&              p,
                                                      const ADB&              satOil,
                                                      const std::vector<int>& cells) const
    {
        return fluid_.rsSat(p, satOil, cells);
    }


    template<class T>
    V
    FullyImplicitBlackoilPolymerSolver<T>::fluidRvSat(const V&                p,
                                                      const V&              satOil,
                                                      const std::vector<int>& cells) const
    {
        return fluid_.rvSat(ADB::constant(p), ADB::constant(satOil), cells).value();
    }





    template<class T>
    ADB
    FullyImplicitBlackoilPolymerSolver<T>::fluidRvSat(const ADB&              p,
                                                      const ADB&              satOil,
                                                      const std::vector<int>& cells) const
    {
        return fluid_.rvSat(p, satOil, cells);
    }

    template<class T>
    ADB
    FullyImplicitBlackoilPolymerSolver<T>::computeMc(const SolutionState& state) const
    {
        return polymer_props_ad_.polymerWaterVelocityRatio(state.concentration);
    }



    template<class T>
    ADB
    FullyImplicitBlackoilPolymerSolver<T>::poroMult(const ADB& p) const
    {
        const int n = p.size();
        if (rock_comp_props_ && rock_comp_props_->isActive()) {
            V pm(n);
            V dpm(n);
            for (int i = 0; i < n; ++i) {
                pm[i] = rock_comp_props_->poroMult(p.value()[i]);
                dpm[i] = rock_comp_props_->poroMultDeriv(p.value()[i]);
            }
            ADB::M dpm_diag = spdiag(dpm);
            const int num_blocks = p.numBlocks();
            std::vector<ADB::M> jacs(num_blocks);
            for (int block = 0; block < num_blocks; ++block) {
                fastSparseProduct(dpm_diag, p.derivative()[block], jacs[block]);
            }
            return ADB::function(std::move(pm), std::move(jacs));
         } else {
            return ADB::constant(V::Constant(n, 1.0), p.blockPattern());
        }
    }





    template<class T>
    ADB
    FullyImplicitBlackoilPolymerSolver<T>::transMult(const ADB& p) const
    {
        const int n = p.size();
        if (rock_comp_props_ && rock_comp_props_->isActive()) {
            V tm(n);
            V dtm(n);
            for (int i = 0; i < n; ++i) {
                tm[i] = rock_comp_props_->transMult(p.value()[i]);
                dtm[i] = rock_comp_props_->transMultDeriv(p.value()[i]);
            }
            ADB::M dtm_diag = spdiag(dtm);
            const int num_blocks = p.numBlocks();
            std::vector<ADB::M> jacs(num_blocks);
            for (int block = 0; block < num_blocks; ++block) {
                fastSparseProduct(dtm_diag, p.derivative()[block], jacs[block]);
            }
            return ADB::function(std::move(tm), std::move(jacs));
        } else {
            return ADB::constant(V::Constant(n, 1.0), p.blockPattern());
        }
    }


    /*
    template<class T>
    void
    FullyImplicitBlackoilSolver<T>::
    classifyCondition(const SolutionState&        state,
                      std::vector<PhasePresence>& cond ) const
    {
        const PhaseUsage& pu = fluid_.phaseUsage();

        if (active_[ Gas ]) {
            // Oil/Gas or Water/Oil/Gas system
            const int po = pu.phase_pos[ Oil ];
            const int pg = pu.phase_pos[ Gas ];

            const V&  so = state.saturation[ po ].value();
            const V&  sg = state.saturation[ pg ].value();

            cond.resize(sg.size());

            for (V::Index c = 0, e = sg.size(); c != e; ++c) {
                if (so[c] > 0)        { cond[c].setFreeOil  (); }
                if (sg[c] > 0)        { cond[c].setFreeGas  (); }
                if (active_[ Water ]) { cond[c].setFreeWater(); }
            }
        }
        else {
            // Water/Oil system
            assert (active_[ Water ]);

            const int po = pu.phase_pos[ Oil ];
            const V&  so = state.saturation[ po ].value();

            cond.resize(so.size());

            for (V::Index c = 0, e = so.size(); c != e; ++c) {
                cond[c].setFreeWater();

                if (so[c] > 0) { cond[c].setFreeOil(); }
            }
        }
    } */


    template<class T>
    void
    FullyImplicitBlackoilPolymerSolver<T>::classifyCondition(const PolymerBlackoilState& state)
    {
        using namespace Opm::AutoDiffGrid;
        const int nc = numCells(grid_);
        const int np = state.numPhases();

        const PhaseUsage& pu = fluid_.phaseUsage();
        const DataBlock s = Eigen::Map<const DataBlock>(& state.saturation()[0], nc, np);
        if (active_[ Gas ]) {
            // Oil/Gas or Water/Oil/Gas system
            const V so = s.col(pu.phase_pos[ Oil ]);
            const V sg = s.col(pu.phase_pos[ Gas ]);

            for (V::Index c = 0, e = sg.size(); c != e; ++c) {
                if (so[c] > 0)        { phaseCondition_[c].setFreeOil  (); }
                if (sg[c] > 0)        { phaseCondition_[c].setFreeGas  (); }
                if (active_[ Water ]) { phaseCondition_[c].setFreeWater(); }
            }
        }
        else {
            // Water/Oil system
            assert (active_[ Water ]);

            const V so = s.col(pu.phase_pos[ Oil ]);


            for (V::Index c = 0, e = so.size(); c != e; ++c) {
                phaseCondition_[c].setFreeWater();

                if (so[c] > 0) { phaseCondition_[c].setFreeOil(); }
            }
        }


    }

    template<class T>
    void
    FullyImplicitBlackoilPolymerSolver<T>::updatePrimalVariableFromState(const PolymerBlackoilState& state)
    {
        using namespace Opm::AutoDiffGrid;
        const int nc = numCells(grid_);
        const int np = state.numPhases();

        const PhaseUsage& pu = fluid_.phaseUsage();
        const DataBlock s = Eigen::Map<const DataBlock>(& state.saturation()[0], nc, np);

        // Water/Oil/Gas system
        assert (active_[ Gas ]);

        // reset the primary variables if RV and RS is not set Sg is used as primary variable.
        primalVariable_.resize(nc);
        std::fill(primalVariable_.begin(), primalVariable_.end(), PrimalVariables::Sg);

        const V sg = s.col(pu.phase_pos[ Gas ]);
        const V so = s.col(pu.phase_pos[ Oil ]);
        const V sw = s.col(pu.phase_pos[ Water ]);

        const double epsilon = std::sqrt(std::numeric_limits<double>::epsilon());
        auto watOnly = sw >  (1 - epsilon);
        auto hasOil = so > 0;
        auto hasGas = sg > 0;

        // For oil only cells Rs is used as primal variable. For cells almost full of water
        // the default primal variable (Sg) is used.
        if (has_disgas_) {
            for (V::Index c = 0, e = sg.size(); c != e; ++c) {
                if ( !watOnly[c] && hasOil[c] && !hasGas[c] ) {primalVariable_[c] = PrimalVariables::RS; }
            }
        }

        // For gas only cells Rv is used as primal variable. For cells almost full of water
        // the default primal variable (Sg) is used.
        if (has_vapoil_) {
            for (V::Index c = 0, e = so.size(); c != e; ++c) {
                if ( !watOnly[c] && hasGas[c] && !hasOil[c] ) {primalVariable_[c] = PrimalVariables::RV; }
            }
        }
        updatePhaseCondFromPrimalVariable();
    }





    /// Update the phaseCondition_ member based on the primalVariable_ member.
    template<class T>
    void
    FullyImplicitBlackoilPolymerSolver<T>::updatePhaseCondFromPrimalVariable()
    {
        if (! active_[Gas]) {
            OPM_THROW(std::logic_error, "updatePhaseCondFromPrimarVariable() logic requires active gas phase.");
        }
        const int nc = primalVariable_.size();
        for (int c = 0; c < nc; ++c) {
            phaseCondition_[c] = PhasePresence(); // No free phases.
            phaseCondition_[c].setFreeWater(); // Not necessary for property calculation usage.
            switch (primalVariable_[c]) {
            case PrimalVariables::Sg:
                phaseCondition_[c].setFreeOil();
                phaseCondition_[c].setFreeGas();
                break;
            case PrimalVariables::RS:
                phaseCondition_[c].setFreeOil();
                break;
            case PrimalVariables::RV:
                phaseCondition_[c].setFreeGas();
                break;
            default:
                OPM_THROW(std::logic_error, "Unknown primary variable enum value in cell " << c << ": " << primalVariable_[c]);
            }
        }
    }




} // namespace Opm
