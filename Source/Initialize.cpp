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
                   const       Geometry& geom)
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

void InitializePermittivity(MultiFab& eps_cc,
                   GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                   GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                   const       Geometry& geom,
                   Real eps_0,
                   Real eps_r
                   ) 
{
    amrex::Print() << "Initializing permittivity at cell-center values in the domain.\n";
	   for (MFIter mfi(eps_cc); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();

        Real R;
        amrex::GpuArray<Real, AMREX_SPACEDIM> center; // physical domain center
        SetupSphereParams(R, center, prob_lo, prob_hi);

        const Array4<Real>& domain_eps = eps_cc.array(mfi);
 
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
    eps_cc.FillBoundary(geom.periodicity());
}


void AveragePermittivityToCellFaces(MultiFab& cc_arr,
                   std::array< MultiFab, AMREX_SPACEDIM >& face_arr) 
{
	   for (MFIter mfi(cc_arr, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Array4<Real const> & cc = cc_arr.array(mfi);
        AMREX_D_TERM(const Array4<Real> & facex = face_arr[0].array(mfi);,
                     const Array4<Real> & facey = face_arr[1].array(mfi);,
                     const Array4<Real> & facez = face_arr[2].array(mfi););

        AMREX_D_TERM(const Box & bx_x = mfi.nodaltilebox(0);,
                     const Box & bx_y = mfi.nodaltilebox(1);,
                     const Box & bx_z = mfi.nodaltilebox(2););
 
        amrex::ParallelFor(bx_x, bx_y, bx_z,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            facex(i,j,k) = 0.5*(cc(i,j,k)+cc(i-1,j,k));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            facey(i,j,k) = 0.5*(cc(i,j,k)+cc(i,j-1,k));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            facez(i,j,k) = 0.5*(cc(i,j,k)+cc(i,j,k-1));
        });
    }
}
