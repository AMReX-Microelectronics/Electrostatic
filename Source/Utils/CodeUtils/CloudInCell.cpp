#include "CloudInCell.H"

AMREX_GPU_HOST_DEVICE
amrex::Real 
CloudInCell::Gather_Trilinear(const amrex::RealVect& pos, 
                              const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& plo, 
                              const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx, 
                              const amrex::Array4<amrex::Real>& phi) 
{
    amrex::Real lx = (pos[0] - plo[0] - dx[0]*0.5)/dx[0];
    amrex::Real ly = (pos[1] - plo[1] - dx[1]*0.5)/dx[1];
    amrex::Real lz = (pos[2] - plo[2] - dx[2]*0.5)/dx[2];

    int i = static_cast<int>(amrex::Math::floor(lx));
    int j = static_cast<int>(amrex::Math::floor(ly));
    int k = static_cast<int>(amrex::Math::floor(lz));

    amrex::Real wx_hi = lx - i; 
    amrex::Real wy_hi = ly - j; 
    amrex::Real wz_hi = lz - k; 

    amrex::Real wx_lo = amrex::Real(1.0) - wx_hi;
    amrex::Real wy_lo = amrex::Real(1.0) - wy_hi;
    amrex::Real wz_lo = amrex::Real(1.0) - wz_hi;

    return ( wx_lo*wy_lo*wz_lo*phi(i  , j  , k  , 0)

           + wx_hi*wy_lo*wz_lo*phi(i+1, j  , k  , 0)
           + wx_lo*wy_hi*wz_lo*phi(i  , j+1, k  , 0)
           + wx_lo*wy_lo*wz_hi*phi(i  , j  , k+1, 0)

           + wx_hi*wy_hi*wz_lo*phi(i+1, j+1, k  , 0)
           + wx_lo*wy_hi*wz_hi*phi(i  , j+1, k+1, 0)
           + wx_hi*wy_lo*wz_hi*phi(i+1, j  , k+1, 0)

           + wx_hi*wy_hi*wz_hi*phi(i+1, j+1, k+1, 0) );
}


AMREX_GPU_HOST_DEVICE
void
CloudInCell::Deposit_Trilinear(const amrex::RealVect& pos, 
                               const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& plo, 
                               const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx, 
                               const amrex::Array4<amrex::Real>& rho,
                               const amrex::Real qp) 
{
    amrex::Real lx = (pos[0] - plo[0] - dx[0]*0.5)/dx[0];
    amrex::Real ly = (pos[1] - plo[1] - dx[1]*0.5)/dx[1];
    amrex::Real lz = (pos[2] - plo[2] - dx[2]*0.5)/dx[2];

    int i = static_cast<int>(amrex::Math::floor(lx));
    int j = static_cast<int>(amrex::Math::floor(ly));
    int k = static_cast<int>(amrex::Math::floor(lz));

    amrex::Real wx_hi = lx - i;
    amrex::Real wy_hi = ly - j;
    amrex::Real wz_hi = lz - k;

    amrex::Real wx_lo = amrex::Real(1.0) - wx_hi;
    amrex::Real wy_lo = amrex::Real(1.0) - wy_hi;
    amrex::Real wz_lo = amrex::Real(1.0) - wz_hi;

    amrex::Gpu::Atomic::AddNoRet(&rho(i  , j  , k  , 0), wx_lo*wy_lo*wz_lo*qp);

    amrex::Gpu::Atomic::AddNoRet(&rho(i+1, j  , k  , 0), wx_hi*wy_lo*wz_lo*qp);
    amrex::Gpu::Atomic::AddNoRet(&rho(i  , j+1, k  , 0), wx_lo*wy_hi*wz_lo*qp);
    amrex::Gpu::Atomic::AddNoRet(&rho(i  , j  , k+1, 0), wx_lo*wy_lo*wz_hi*qp);

    amrex::Gpu::Atomic::AddNoRet(&rho(i+1, j+1, k  , 0), wx_hi*wy_hi*wz_lo*qp);
    amrex::Gpu::Atomic::AddNoRet(&rho(i  , j+1, k+1, 0), wx_lo*wy_hi*wz_hi*qp);
    amrex::Gpu::Atomic::AddNoRet(&rho(i+1, j  , k+1, 0), wx_hi*wy_lo*wz_hi*qp);

    amrex::Gpu::Atomic::AddNoRet(&rho(i+1, j+1, k+1, 0), wx_hi*wy_hi*wz_hi*qp);
}
