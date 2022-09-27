#include "MLMG.H"

#include "../Utils/WarpXUtil.H"
#include "Code.H"
#include "Input/GeometryProperties.H"
#include "Input/MacroscopicProperties.H"

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Parser.H>
#include <AMReX_RealBox.H>
#include <AMReX_MultiFab.H>
//#include <AMReX_MultiFabUtil.H>


using namespace amrex;

c_MLMGSolver::c_MLMGSolver ()
{
    ReadData();
} 


void
c_MLMGSolver::ReadData() 
{

    ParmParse pp_mlmg("mlmg");
    pp_mlmg.get("ascalar",ascalar);
    pp_mlmg.get("bscalar",bscalar);

    if (queryWithParser(pp_mlmg, "set_verbose" , set_verbose) ) {
    }
    else {
        set_verbose = 0;
        std::stringstream warnMsg;
        warnMsg << "MLMG parameter 'set_verbose'" << " is not specified in the input file. The default value of "
                <<  set_verbose
                << " is used.";
        c_Code::GetInstance().RecordWarning("MLMG properties", warnMsg.str());
    }

    if (queryWithParser(pp_mlmg, "max_order" , set_verbose) ) {
    }
    else {
        max_order = 2;
        std::stringstream warnMsg;
        warnMsg << "MLMG parameter 'max_order'" << " is not specified in the input file. The default value of "
                <<  set_verbose
                << " is used.";
        c_Code::GetInstance().RecordWarning("MLMG properties", warnMsg.str());
    }

    if (queryWithParser(pp_mlmg, "relative_tolerance" , relative_tolerance) ) {
    }
    else {
        relative_tolerance = 1.0e-10;
        std::stringstream warnMsg;
        warnMsg << "MLMG parameter 'relative_tolerance'" << " is not specified in the input file. The default value of "
                <<  relative_tolerance
                << " is used.";
        c_Code::GetInstance().RecordWarning("MLMG properties", warnMsg.str());
    }

    if (queryWithParser(pp_mlmg, "absolute_tolerance" , absolute_tolerance) ) {
    }
    else {
        absolute_tolerance = 0.;
        std::stringstream warnMsg;
        warnMsg << "MLMG parameter 'absolute_tolerance'" << " is not specified in the input file. The default value of "
                <<  absolute_tolerance
                << " is used.";
        c_Code::GetInstance().RecordWarning("MLMG properties", warnMsg.str());
    }

}


void
c_MLMGSolver:: Setup_MLABecLaplacian_ForPoissonEqn()
{

    auto& rCode = c_Code::GetInstance();
    auto& rGprop = rCode.get_GeometryProperties();
    auto& ba = rGprop.ba;
    auto& dm = rGprop.dm;
    auto& geom = rGprop.geom;
    auto& rMprop = rCode.get_MacroscopicProperties();
    auto& eps_cc = rMprop.get_mf("epsilon");
    auto& phi_cc = rMprop.get_mf("phi");

    mlabec.define({geom}, {ba}, {dm}, info); // Implicit solve using MLABecLaplacian class

    // Force singular system to be solvable
    mlabec.setEnforceSingularSolvable(false);

    // set order of stencil
    mlabec.setMaxOrder(max_order);

    // build array of boundary conditions needed by MLABecLaplacian
    // see Src/Boundary/AMReX_LO_BCTYPES.H for supported types
    std::array<LinOpBCType, AMREX_SPACEDIM> bc_lo;
    std::array<LinOpBCType, AMREX_SPACEDIM> bc_hi;

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
          bc_lo[idim] = bc_hi[idim] = LinOpBCType::Dirichlet;
    }

    // tell the solver what the domain boundary conditions are
    mlabec.setDomainBC(bc_lo, bc_hi);

    // initialize phi to zero including the ghost cells
    // PoissonPhi.setVal(0.);

    // set Dirichlet BC by reading in the ghost cell values
    int amrlev = 0; //refers to the setcoarsest level of the solve
    mlabec.setLevelBC(amrlev, &phi_cc);

    // set scaling factors 
    mlabec.setScalars(ascalar, bscalar);

    int Ncomp1=1;
    int Nghost0=0;
    // set alpha_cc, and beta_fc coefficients
    // alpha_cc is a cell-centered multifab
    alpha_cc =  std::make_unique<amrex::MultiFab>(ba, dm, Ncomp1, Nghost0); 
    alpha_cc->setVal(0.); // fill in alpha_cc multifab to the value of 0.0

    mlabec.setACoeffs(amrlev, *alpha_cc);

 //   beta_fc =  std::make_unique<amrex::FArrayBox,AMREX_SPACEDIM>; 
    // beta_fc is a face-centered multifab
    AMREX_D_TERM(beta_fc[0].define(convert(ba,IntVect(AMREX_D_DECL(1,0,0))), dm, 1, 0);,
                 beta_fc[1].define(convert(ba,IntVect(AMREX_D_DECL(0,1,0))), dm, 1, 0);,
                 beta_fc[2].define(convert(ba,IntVect(AMREX_D_DECL(0,0,1))), dm, 1, 0););

    AverageCellCenteredMultiFabToCellFaces(eps_cc, beta_fc); //converts from cell-centered permittivity to face-center and store in beta_fc

    mlabec.setBCoeffs(0, amrex::GetArrOfConstPtrs(beta_fc));

    pMLMG = std::make_unique<MLMG>(mlabec);

    pMLMG->setVerbose(set_verbose);

}

void
c_MLMGSolver:: Solve_PoissonEqn()
{

    auto& rCode = c_Code::GetInstance();
    auto& rMprop = rCode.get_MacroscopicProperties();
    auto& rho_cc = rMprop.get_mf("charge_density");
    auto& phi_cc = rMprop.get_mf("phi");

    pMLMG->solve({&phi_cc}, {&rho_cc},
                 relative_tolerance,
                 absolute_tolerance);

}


void
c_MLMGSolver:: Compute_vecE(std::array<amrex::MultiFab, AMREX_SPACEDIM>& vecE)
{
    /** First the gradient of phi is computed then it is multiplied by one. */
    amrex::Print() << "Computing vecE using MLMG. \n";
    auto& rCode = c_Code::GetInstance();
    auto& rGprop = rCode.get_GeometryProperties();
    auto& ba = rGprop.ba;
    auto& dm = rGprop.dm;

    amrex::Vector<Array<MultiFab, AMREX_SPACEDIM>> gradPhi;
    gradPhi.resize(1);
    //TD<decltype(gradPhi)> gradPhi_type; //amrex::Vector<std::array<amrex::MultiFab, 3> >
    //TD<decltype(gradPhi[0])> gradPhi_0_type; //std::array<amrex::MultiFab, 3>&
    //TD<decltype(gradPhi[0][0])> gradPhi_00_type; //amrex::MultiFab&
    //TD<decltype(gradPhi[0][0][0])> gradPhi_000_type; //amrex::FArrayBox&
    //TD<decltype(gradPhi[0][0][0].array())> gradPhi_000array_type; //amrex::Array4<double>
    //TD<decltype(vecE)> vecE_type; //std::array<amrex::MultiFab, 3>&
    //TD<decltype(vecE[0])> vecE_0_type; //amrex::MultiFab&
    //TD<decltype(vecE[0][0])> vecE_00_type;//amrex::FArrayBox&

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        
        gradPhi[0][idim].define(amrex::convert(ba,amrex::IntVect::TheDimensionVector(idim)),dm, 1, 0);

        vecE[idim].define(amrex::convert(ba,IntVect::TheDimensionVector(idim)),dm, 1, 0);

    }
    pMLMG->getGradSolution(amrex::GetVecOfArrOfPtrs(gradPhi));

    /*Saxpy does E[x] += -1*gradPhi[0][x], where initial value of E[x]=0*/
    AMREX_D_TERM ( MultiFab::Saxpy (vecE[0], -1.0, gradPhi[0][0],0,0,1,0); ,
                   MultiFab::Saxpy (vecE[1], -1.0, gradPhi[0][1],0,0,1,0); ,
                   MultiFab::Saxpy (vecE[2], -1.0, gradPhi[0][2],0,0,1,0); );



}


void
c_MLMGSolver:: Compute_vecFlux(std::array<amrex::MultiFab, AMREX_SPACEDIM>& vecFlux)
{

    /** Flux (-beta*grad phi)  is computed. */
    amrex::Print() << "Computing vecFlux using MLMG. \n";

    auto& rCode = c_Code::GetInstance();
    auto& rGprop = rCode.get_GeometryProperties();
    auto& ba = rGprop.ba;
    auto& dm = rGprop.dm;
    amrex::Vector<Array<MultiFab, AMREX_SPACEDIM>> minus_epsGradPhi;
    minus_epsGradPhi.resize(1);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        minus_epsGradPhi[0][idim].define(amrex::convert(ba,amrex::IntVect::TheDimensionVector(idim)),dm, 1, 0);
        vecFlux[idim].define(amrex::convert(ba,amrex::IntVect::TheDimensionVector(idim)),dm, 1, 0);
    }
     pMLMG->getFluxes(amrex::GetVecOfArrOfPtrs(minus_epsGradPhi));
     
     /*Copy copies minus_epsGradPhi[0][x] multifab to vecFlux[x] multifab */
     AMREX_D_TERM ( MultiFab::Copy (vecFlux[0], minus_epsGradPhi[0][0],0,0,1,0); ,
                    MultiFab::Copy (vecFlux[1], minus_epsGradPhi[0][1],0,0,1,0); ,
                    MultiFab::Copy (vecFlux[2], minus_epsGradPhi[0][2],0,0,1,0); );

}


void 
c_MLMGSolver:: AverageCellCenteredMultiFabToCellFaces(const amrex::MultiFab& cc_arr,
                   std::array< amrex::MultiFab, AMREX_SPACEDIM >& face_arr)
{
    for (MFIter mfi(cc_arr, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Array4<Real const> & cc = cc_arr.array(mfi);
        AMREX_D_TERM(const Array4<Real> & facex = face_arr[0].array(mfi);,
                     const Array4<Real> & facey = face_arr[1].array(mfi);,
                     const Array4<Real> & facez = face_arr[2].array(mfi););

        AMREX_D_TERM(const Box & nodal_x = mfi.nodaltilebox(0);,
                     const Box & nodal_y = mfi.nodaltilebox(1);,
                     const Box & nodal_z = mfi.nodaltilebox(2););

        amrex::ParallelFor(nodal_x, nodal_y, nodal_z,
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

