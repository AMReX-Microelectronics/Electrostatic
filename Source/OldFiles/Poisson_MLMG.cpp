#include "Poisson_MLMG.H"

void setupMLABecForPoissonProblem(MLABecLaplacian& mlabec,
                                  Real mlmg_ascalar,
                                  Real mlmg_bscalar, 
                                  int mlmg_max_order,
                                  const BoxArray& ba,
                                  const DistributionMapping& dm,
                                  MultiFab& PoissonPhi, 
                                  MultiFab& eps_cc, 
                                  MultiFab& rho) 
{
    // Multi-Level Multi-Grid solver is used which solves for phi in (A*alpha_cc - B * div (beta_face grad)) phi = rho
    // For Poisson solver, A and alpha_cc are set to zero.
    // see AMReX_MLLinOp.H for the defaults, accessors, and mutators

    // Force singular system to be solvable
    mlabec.setEnforceSingularSolvable(false);

    // set order of stencil
    mlabec.setMaxOrder(mlmg_max_order);

    // build array of boundary conditions needed by MLABecLaplacian
    // see Src/Boundary/AMReX_LO_BCTYPES.H for supported types
    std::array<LinOpBCType, AMREX_SPACEDIM> mlmg_bc_lo;
    std::array<LinOpBCType, AMREX_SPACEDIM> mlmg_bc_hi;

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
          mlmg_bc_lo[idim] = mlmg_bc_hi[idim] = LinOpBCType::Dirichlet;
    }

    // tell the solver what the domain boundary conditions are
    mlabec.setDomainBC(mlmg_bc_lo,mlmg_bc_hi);

    // initialize phi to zero including the ghost cells
    PoissonPhi.setVal(0.);
    // set Dirichlet BC by reading in the ghost cell values
    int amrlev = 0; //refers to the setcoarsest level of the solve
    mlabec.setLevelBC(amrlev, &PoissonPhi);

    // set scaling factors 
    mlabec.setScalars(mlmg_ascalar, mlmg_bscalar);

    // set alpha_cc, and beta_face coefficients
    MultiFab alpha_cc(ba, dm, 1, 0); // cell-centered 
    alpha_cc.setVal(0.); // fill in alpha_cc multifab to the value of 0.0
    mlabec.setACoeffs(amrlev, alpha_cc);

    std::array< MultiFab, AMREX_SPACEDIM > beta_face; // beta_face is the permittivity in Poisson's equation
    // beta_face lives on faces so we make an array of face-centered MultiFabs
    AMREX_D_TERM(beta_face[0].define(convert(ba,IntVect(AMREX_D_DECL(1,0,0))), dm, 1, 0);,
                 beta_face[1].define(convert(ba,IntVect(AMREX_D_DECL(0,1,0))), dm, 1, 0);,
                 beta_face[2].define(convert(ba,IntVect(AMREX_D_DECL(0,0,1))), dm, 1, 0););

    AverageCellCenteredMultiFabToCellFaces(eps_cc, beta_face); //converts from cell-centered permittivity to face-center and store in beta_face

    mlabec.setBCoeffs(0, amrex::GetArrOfConstPtrs(beta_face));
}

template<typename T>
class TD;

void GetEFieldUsingMLMG(MLMG& mlmg,
                        Array< MultiFab, AMREX_SPACEDIM >& E,
                        const BoxArray& ba,
                        const DistributionMapping& dm)
{

    amrex::Vector<Array<MultiFab, AMREX_SPACEDIM>> gradPhi;
    gradPhi.resize(1);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        gradPhi[0][idim].define(amrex::convert(ba,IntVect::TheDimensionVector(idim)),dm, 1, 0);
        E[idim].define(amrex::convert(ba,IntVect::TheDimensionVector(idim)),dm, 1, 0);
    }
    mlmg.getGradSolution(GetVecOfArrOfPtrs(gradPhi));
    
    //TD<decltype(GetVecOfArrOfPtrs(gradPhi))> type;
    //TD<decltype(gradPhi[0][0])> gradPhi_00_type;
    //TD<decltype(gradPhi[0][1][0].array())> gradPhi_0_type;
    //TD<decltype(E)> E_type;
    //TD<decltype(E[1])> E0_type;
    //TD<decltype(E[0][0])> E00_type;
    
    MultiFab::Saxpy (E[0], -1.0, gradPhi[0][0],0,0,1,0); /*E[0] += -1*gradPhi[0][0], where initial value of E[0]=0*/
    MultiFab::Saxpy (E[1], -1.0, gradPhi[0][1],0,0,1,0);
    MultiFab::Saxpy (E[2], -1.0, gradPhi[0][2],0,0,1,0);

    //Print() << "gradPhi: " << gradPhi[0][0][0].array()(0,0,0);
    //Print() << "E: " << E[0][0].array()(0,0,0);
}

void GetFluxesUsingMLMG(MLMG& mlmg,
                        amrex::Vector<Array< MultiFab, AMREX_SPACEDIM >>& minus_epsGradPhi,
                        const BoxArray& ba,
                        const DistributionMapping& dm)
{
    /*Flux (-beta*grad phi)  is computed*/
    minus_epsGradPhi.resize(1);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        minus_epsGradPhi[0][idim].define(amrex::convert(ba,IntVect::TheDimensionVector(idim)),dm, 1, 0);
    }
    mlmg.getFluxes(GetVecOfArrOfPtrs(minus_epsGradPhi));
}
