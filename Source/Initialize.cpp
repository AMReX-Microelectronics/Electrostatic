#include "Initialize.H"


void SetupSphereParams(Real& R,
                       GpuArray<Real, AMREX_SPACEDIM>& center,
                       GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                       GpuArray<Real, AMREX_SPACEDIM> prob_hi)
{
    amrex::Vector<Real> domainLen(AMREX_SPACEDIM); //domain length
    for (auto i=0; i<AMREX_SPACEDIM; ++i) {
        center[i] = prob_lo[i] + 0.5*(prob_hi[i] - prob_lo[i]);
        domainLen[i] = std::abs(center[i] - prob_lo[i]);
    }
    //define radius of the charge sphere=(1/8)*minimum domain length 
    R = *std::min_element(domainLen.begin(), domainLen.end())/8.0;

    Print() << "center: " << center[0] << std::setw(5) << center[1] << std:: setw(5) << center[2] << "\n";
    Print() << "domain length: " << domainLen[0] << std::setw(5) << domainLen[1] << std:: setw(5) << domainLen[2] << "\n";
    Print() << "Radius of the charge: " << R << "\n";
}

//Initialize charge
void InitializeCharge(MultiFab& rho,
                   GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                   GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                   const       Geometry& geom
                   )
{
   amrex::Print() << "Initializing charge sphere in the center of the domain.\n";
   // loop over boxes
	   for (MFIter mfi(rho); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();
 
        Real R;
        amrex::GpuArray<Real, AMREX_SPACEDIM> center; // physical domain center
        SetupSphereParams(R, center, prob_lo, prob_hi);

        const Array4<Real>& chargeDensity_arr = rho.array(mfi);
 
        amrex::ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
             GpuArray<Real,AMREX_SPACEDIM> diff{AMREX_D_DECL(std::pow(prob_lo[0] + (i+0.5)*dx[0] - center[0],2),
                                                             std::pow(prob_lo[1] + (j+0.5)*dx[1] - center[1],2),
                                                             std::pow(prob_lo[2] + (k+0.5)*dx[2] - center[2],2) )};

             auto rad = std::pow(diff[0] + diff[1] + diff[2], 0.5);

             if(rad <= R){ 

                  chargeDensity_arr(i,j,k) = q;

             } else {

                  chargeDensity_arr(i,j,k) = 0.0;

             }
        });
    }
}

void InitializePermittivity(MultiFab& eps_mfab,
                   GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                   GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                   const       Geometry& geom,
                   Real eps_0,
                   Real eps_r
                   ) {
    amrex::Print() << "Initializing permittivity in the domain.\n";
	   for (MFIter mfi(eps_mfab); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();

        Real R;
        amrex::GpuArray<Real, AMREX_SPACEDIM> center; // physical domain center
        SetupSphereParams(R, center, prob_lo, prob_hi);

        const Array4<Real>& domain_eps = eps_mfab.array(mfi);
 
        amrex::ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
             GpuArray<Real,AMREX_SPACEDIM> diff{AMREX_D_DECL(std::pow(prob_lo[0] + (i+0.5)*dx[0] - center[0],2),
                                                             std::pow(prob_lo[1] + (j+0.5)*dx[1] - center[1],2),
                                                             std::pow(prob_lo[2] + (k+0.5)*dx[2] - center[2],2) )};

             auto rad = std::pow(diff[0] + diff[1] + diff[2], 0.5);


             if(rad <= R){ 

                  domain_eps(i,j,k) = eps_r*eps_0;

             } else {

                  domain_eps(i,j,k) = eps_0;

             }
        });
    }
    eps_mfab.FillBoundary(geom.periodicity());
}


