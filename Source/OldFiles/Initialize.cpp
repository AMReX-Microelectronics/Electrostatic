#include "Initialize.H"

void InitializeBoxArrayAndGeom(BoxArray& ba, 
                               Geometry& geom, 
                               GpuArray<Real, AMREX_SPACEDIM> const& prob_lo,
                               GpuArray<Real, AMREX_SPACEDIM> const& prob_hi,
                               GpuArray<int, AMREX_SPACEDIM> const& n_cell,
                               int max_grid_size)
{

    IntVect dom_lo(AMREX_D_DECL(       0,        0,        0)); // domain low 
    IntVect dom_hi(AMREX_D_DECL(n_cell[0]-1, n_cell[1]-1, n_cell[2]-1)); // domain high

    Box domain(dom_lo, dom_hi); // Make a single box that is the entire domain

    ba.define(domain); // initialize the boxarray 'ba' from the single box 'domain'

    ba.maxSize(max_grid_size); // break up ba into chunks no larger than 'max_grid_size' along a direction

    RealBox real_box({AMREX_D_DECL( prob_lo[0], prob_lo[1], prob_lo[2])},
                     {AMREX_D_DECL( prob_hi[0], prob_hi[1], prob_hi[2])});  //physical domain 

    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)}; // 0: not periodic, 1: periodic

    geom.define(domain, real_box, CoordSys::cartesian, is_periodic); //define the geom object

    GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray(); // obtain cell size array dx from the geom object
    Print() << "dx: " << dx[0] << " " << dx[1] << " " << dx[2] << "\n";

}

void SetupSphereParams(Real& R,
                       GpuArray<Real, AMREX_SPACEDIM>& center,
                       GpuArray<Real, AMREX_SPACEDIM> const& prob_lo,
                       GpuArray<Real, AMREX_SPACEDIM> const& prob_hi)
{
    amrex::Vector<Real> domainLen(AMREX_SPACEDIM); //domain length
    for (auto i=0; i<AMREX_SPACEDIM; ++i) {
        center[i] = prob_lo[i] + 0.5*(prob_hi[i] - prob_lo[i]);
        domainLen[i] = std::abs(center[i] - prob_lo[i]);
    }
    //define radius of the charge sphere=(1/4)*minimum domain length 
    R = *std::min_element(domainLen.begin(), domainLen.end())/4.0;

    //Print() << "center: " << center[0] << std::setw(5) << center[1] << std:: setw(5) << center[2] << "\n";
    //Print() << "domain length: " << domainLen[0] << std::setw(5) << domainLen[1] << std:: setw(5) << domainLen[2] << "\n";
    //Print() << "Radius of the charge: " << R << "\n";
}

//Initialize charge
void InitializeCharge(MultiFab& rho,
                   GpuArray<Real, AMREX_SPACEDIM> const& prob_lo,
                   GpuArray<Real, AMREX_SPACEDIM> const& prob_hi,
                   const Geometry& geom)
{
   amrex::Print() << "Initializing charge sphere in the center of the domain.\n";
   // loop over boxex
    for (MFIter mfi(rho); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();
 
        Real R, V;
        amrex::GpuArray<Real, AMREX_SPACEDIM> center; // physical domain center
        SetupSphereParams(R, center, prob_lo, prob_hi);
        Print() << "**************CHARGED_SPHERE_PROBLEM_SPECIFICATIONS************** \n";
        Print() << "Radius of the charged sphere, R: " << R << "\n";
        Print() << "Center of the charged sphere: " << center[0] << " " << center[1] << " " << center[2] << "\n";
        V=(4./3.)*pi*std::pow(R,3);
        Print() << "Volume of the sphere, V: " << V << "\n";
        Real Q = 100*q;
        Print() << "Total charge assumed inside the sphere, Q=100q: " << Q << "\n";
        Real rho_C = Q/V;
        Print() << "Uniform charge density, rho_C=Q/V: " << rho_C << "\n";
        Print() << "We expect: div(-epsilon*gradPhi)=rho, i.e., OVER_A_CLOSED_SURFACE_NET_SUM_OF_OUTWARD_FLUXES_OF(-epsilon*gradPhi)=Q " << "\n";
        Print() << "***************************************************************** \n";

        const Array4<Real>& chargeDensity_arr = rho.array(mfi);
 
        amrex::ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
             GpuArray<Real,AMREX_SPACEDIM> diff{AMREX_D_DECL(std::pow(prob_lo[0] + (i+0.5)*dx[0] - center[0],2),
                                                             std::pow(prob_lo[1] + (j+0.5)*dx[1] - center[1],2),
                                                             std::pow(prob_lo[2] + (k+0.5)*dx[2] - center[2],2) )};

             auto rad = std::pow(diff[0] + diff[1] + diff[2], 0.5);

             if(rad <= R){ 

                  chargeDensity_arr(i,j,k) = rho_C;

             } else {

                  chargeDensity_arr(i,j,k) = 0.0;

             }
        });
    }
}

void InitializePermittivity(MultiFab& eps_cc,
                   GpuArray<Real, AMREX_SPACEDIM> const& prob_lo,
                   GpuArray<Real, AMREX_SPACEDIM> const& prob_hi,
                   const Geometry& geom,
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
