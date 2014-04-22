#include <cmath>
#include <vector>
#include <opm/polymer/fullyimplicit/AutoDiffBlock.hpp>
#include <opm/polymer/fullyimplicit/AutoDiffHelpers.hpp>
#include <opm/polymer/fullyimplicit/PolymerPropsAd.hpp>

namespace Opm {





    typedef PolymerPropsAd::ADB ADB;
    typedef PolymerPropsAd::V V;





    double
    PolymerPropsAd::rockDensity() const
    {
        return polymer_props_.rockDensity();
    }





    double
    PolymerPropsAd::deadPoreVol() const
    {
        return polymer_props_.deadPoreVol();
    }





	double
	PolymerPropsAd::cMax() const
	{
		return polymer_props_.cMax();
	}



	

    PolymerPropsAd::PolymerPropsAd(const PolymerProperties& polymer_props)
        : polymer_props_ (polymer_props)
    {
    }





    PolymerPropsAd::~PolymerPropsAd()
    {
    }





    V PolymerPropsAd::effectiveInvWaterVisc(const V& c,
                                            const double* visc) const
    {
        const int nc = c.size();
        V inv_mu_w_eff(nc);
        for (int i = 0; i < nc; ++i) {
            double im = 0;
            polymer_props_.effectiveInvVisc(c(i), visc, im);
            inv_mu_w_eff(i) = im;
        }

        return inv_mu_w_eff;
    }




    ADB PolymerPropsAd::effectiveInvWaterVisc(const ADB& c,
	                    				      const double* visc) const
    {
	    const int nc = c.size();
    	V inv_mu_w_eff(nc);
    	V dinv_mu_w_eff(nc);
    	for (int i = 0; i < nc; ++i) {
    	    double im = 0, dim = 0;
    	    polymer_props_.effectiveInvViscWithDer(c.value()(i), visc, im, dim);
    	    inv_mu_w_eff(i) = im;
    	    dinv_mu_w_eff(i) = dim;
    	}
        ADB::M dim_diag = spdiag(dinv_mu_w_eff);
        const int num_blocks = c.numBlocks();
        std::vector<ADB::M> jacs(num_blocks);
        for (int block = 0; block < num_blocks; ++block) {
            jacs[block] = dim_diag * c.derivative()[block];
        }
        return ADB::function(inv_mu_w_eff, jacs);
    }


    V PolymerPropsAd::effectiveInvWaterViscWithShear(const V& c,
                                                     const V& shear_mult,
                                                     const double* visc) const
    {
       return effectiveInvWaterVisc(c, visc) / shear_mult;

    }


    ADB PolymerPropsAd::effectiveInvWaterViscWithShear(const ADB& c,
                                                       const ADB& shear_mult,
                                                       const double* visc) const
    {
       return effectiveInvWaterVisc(c, visc) / shear_mult;

    }




    V PolymerPropsAd::polymerWaterVelocityRatio(const V& c) const
    {
        const int nc = c.size();
        V mc(nc);

        for (int i = 0; i < nc; ++i) {
            double m = 0;
            polymer_props_.computeMc(c(i), m);
            mc(i) = m;
        }
       
       return mc;
    }





    ADB PolymerPropsAd::polymerWaterVelocityRatio(const ADB& c) const
    {
    
        const int nc = c.size();
        V mc(nc);
        V dmc(nc);
        
        for (int i = 0; i < nc; ++i) {
            double m = 0;
            double dm = 0;
            polymer_props_.computeMcWithDer(c.value()(i), m, dm);

            mc(i) = m;
            dmc(i) = dm;
        }

        ADB::M dmc_diag = spdiag(dmc);
        const int num_blocks = c.numBlocks();
        std::vector<ADB::M> jacs(num_blocks);
        for (int block = 0; block < num_blocks; ++block) {
            jacs[block] = dmc_diag * c.derivative()[block];
        }

        return ADB::function(mc, jacs);
    }





    V PolymerPropsAd::adsorption(const V& c, const V& cmax_cells) const
    {
        const int nc = c.size();
        V ads(nc);

        for (int i = 0; i < nc; ++i) {
            double c_ads = 0;
            polymer_props_.adsorption(c(i), cmax_cells(i), c_ads);

            ads(i) = c_ads;
        }

        return ads;
    }





    ADB PolymerPropsAd::adsorption(const ADB& c, const ADB& cmax_cells) const
    {
        const int nc = c.value().size();

        V ads(nc);
        V dads(nc);

        for (int i = 0; i < nc; ++i) {
            double c_ads = 0;
            double dc_ads = 0;
            polymer_props_.adsorptionWithDer(c.value()(i), cmax_cells.value()(i), c_ads, dc_ads);
            ads(i) = c_ads;
            dads(i) = dc_ads;
        }

        ADB::M dads_diag = spdiag(dads);
        int num_blocks = c.numBlocks();
        std::vector<ADB::M> jacs(num_blocks);
        for (int block = 0; block < num_blocks; ++block) {
            jacs[block] = dads_diag * c.derivative()[block];
        }

        return ADB::function(ads, jacs);
    }

    V
    PolymerPropsAd::permReduction(const V& c, 
                                  const V& cmax_cells) const
    {

        const int nc = c.size();

        V one  = V::Ones(nc);
        V ads = adsorption(c, cmax_cells);
        double max_ads = polymer_props_.cMaxAds();
        double res_factor = polymer_props_.resFactor();
        double factor = (res_factor -1.) / max_ads;
        V rk = one + factor * ads;
        
        return rk;
    }


    ADB
    PolymerPropsAd::permReduction(const ADB& c,
                                  const ADB& cmax_cells) const
    {
        const int nc = c.size();
        V one = V::Ones(nc);
        ADB ads = adsorption(c, cmax_cells);

        double max_ads = polymer_props_.cMaxAds();
        double res_factor = polymer_props_.resFactor();
        double factor = (res_factor - 1.) / max_ads;
        ADB rk = one + ads * factor; 

        return rk;
    }



    V
    PolymerPropsAd::effectiveRelPerm(const V& c, 
                                     const V& cmax_cells,
                                     const V& krw) const
    {
        const int nc = c.size();

        V one  = V::Ones(nc);
        V ads = adsorption(c, cmax_cells);
        double max_ads = polymer_props_.cMaxAds();
        double res_factor = polymer_props_.resFactor();
        double factor = (res_factor -1.) / max_ads;
        V rk = one + factor * ads;
        V krw_eff = krw / rk;

        return krw_eff;
    }





    ADB
    PolymerPropsAd::effectiveRelPerm(const ADB& c,
                                     const ADB& cmax_cells,
                                     const ADB& krw,
                                     const ADB& sw) const
    {
        const int nc = c.value().size();
        V one = V::Ones(nc);
        ADB ads = adsorption(c, cmax_cells);
        V krw_eff = effectiveRelPerm(c.value(), cmax_cells.value(), krw.value());

        double max_ads = polymer_props_.cMaxAds();
        double res_factor = polymer_props_.resFactor();
        double factor = (res_factor - 1.) / max_ads;
        ADB rk = one + ads * factor; 
		ADB dkrw_ds = krw / rk;
        ADB dkrw_dc = -factor * krw / (rk * rk) * ads ;

        const int num_blocks = c.numBlocks();
        std::vector<ADB::M> jacs(num_blocks);
        for (int block = 0; block < num_blocks; ++block) {
            jacs[block] = dkrw_ds.derivative()[block] * sw.derivative()[block] 
						+ dkrw_dc.derivative()[block] * c.derivative()[block];
        }

        return ADB::function(krw_eff, jacs);
    }
    
    const V
    PolymerPropsAd::shearWaterVelocity() const
    {
        const std::vector<double> vw = polymer_props_.shearWaterVelocity();
        const V watervelo = Eigen::Map<const V>(&vw[0], vw.size());
        return watervelo;
    }
    const V
    PolymerPropsAd::shearVrf() const
    {
        const std::vector<double>vrf = polymer_props_.shearVrf();
        const V shear_vrf = Eigen::Map<const V>(&vrf[0], vrf.size());
        return shear_vrf;
    }
    V
    PolymerPropsAd::shearMult(const V& velocity) const
    {
        const int nc = velocity.size();
        V shear_mult(nc);
        for (int i = 0; i < nc; ++i) {
            shear_mult(i) = polymer_props_.shearVrf(velocity(i));
        }
        return shear_mult;
    }

    ADB
    PolymerPropsAd::shearMult(const ADB& velocity) const
    {
        const int nc = velocity.size();
        V shear_mult(nc);
        V dshear_mult(nc);
        for (int i = 0; i < nc; ++i) {
            double sm = 0, dsm = 0;
            polymer_props_.shearVrfWithDer(sm, dsm);
            shear_mult(i) = sm;
            dshear_mult(i) = dsm;
        }
        ADB::M dsm_diag = spdiag(dshear_mult);
        const int num_blocks = velocity.numBlocks();
        std::vector<ADB::M> jacs(num_blocks);
        for ( int block = 0; block < num_blocks; ++block) {
            jacs[block] = dsm_diag * velocity.derivative()[block];
        }
        return ADB::function(shear_mult, jacs);
    }
}// namespace Opm
