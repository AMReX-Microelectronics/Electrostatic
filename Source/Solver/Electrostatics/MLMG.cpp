#include "MLMG.H"

#include "../../Utils/SelectWarpXUtils/WarpXUtil.H"
#include "../../Utils/SelectWarpXUtils/TextMsg.H"
#include "../../Utils/CodeUtils/CodeUtil.H"
#include "Code.H"
#include "Input/GeometryProperties/GeometryProperties.H"
#include "Input/MacroscopicProperties/MacroscopicProperties.H"
#include "Input/BoundaryConditions/BoundaryConditions.H"

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Parser.H>
#include <AMReX_RealBox.H>
#include <AMReX_MultiFab.H>
//#include <AMReX_MultiFabUtil.H>


using namespace amrex;

c_MLMGSolver::c_MLMGSolver ()
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t{************************c_MLMGSolver Constructor()************************\n";
    amrex::Print() << "\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    ReadData();

    all_homogeneous_boundaries=true;
    some_functionbased_inhomogeneous_boundaries=false;
    some_constant_inhomogeneous_boundaries=false;

    AssignLinOpBCTypeToBoundaries();

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t}************************c_MLMGSolver Constructor()************************\n";
#endif
} 


c_MLMGSolver::~c_MLMGSolver ()
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t{************************c_MLMGSolver Destructor()************************\n";
    amrex::Print() << "\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    alpha_cc=nullptr;
    beta_cc=nullptr;
    soln_cc=nullptr;
    rhs_cc=nullptr;
    robin_a=nullptr;
    robin_b=nullptr;
    robin_f=nullptr;

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t}************************c_MLMGSolver Destructor()************************\n";
#endif
} 


void
c_MLMGSolver::ReadData() 
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t\t{************************c_MLMGSolver::ReadData()************************\n";
    amrex::Print() << "\t\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif


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

    if (queryWithParser(pp_mlmg, "max_order" , max_order) ) {
    }
    else {
        max_order = 2;
        std::stringstream warnMsg;
        warnMsg << "MLMG parameter 'max_order'" << " is not specified in the input file. The default value of "
                <<  max_order
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
     


    auto& rCode = c_Code::GetInstance();
    auto& rMprop = rCode.get_MacroscopicProperties();
    std::map<std::string,int>::iterator it_Mprop;

    pp_mlmg.get("alpha",alpha_cc_str);

    //Assert that the multifab name is specified through 'macroscopic.fields_to_define' in the input file.
    it_Mprop = rMprop.map_param_all.find(alpha_cc_str);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(it_Mprop != rMprop.map_param_all.end(), 
                                     "'" + alpha_cc_str + "' is not specified through 'macroscopic.fields_to_define' in the input file.");

    pp_mlmg.get("beta" ,beta_cc_str);

    it_Mprop = rMprop.map_param_all.find(beta_cc_str);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(it_Mprop != rMprop.map_param_all.end(), 
                                     beta_cc_str + " is not specified through 'macroscopic.fields_to_define' in the input file.");

    pp_mlmg.get("soln" ,soln_cc_str);

    it_Mprop = rMprop.map_param_all.find(soln_cc_str);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(it_Mprop != rMprop.map_param_all.end(), 
                                     soln_cc_str + " is not specified through 'macroscopic.fields_to_define' in the input file.");

    pp_mlmg.get("rhs"  ,rhs_cc_str);

    it_Mprop = rMprop.map_param_all.find(rhs_cc_str);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(it_Mprop != rMprop.map_param_all.end(), 
                                     rhs_cc_str + " is not specified through 'macroscopic.fields_to_define' in the input file.");

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t\t}************************c_MLMGSolver::ReadData()************************\n";
#endif
}


void
c_MLMGSolver::AssignLinOpBCTypeToBoundaries ()
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t\t{************************c_MLMGSolver::AssignLinOpBCTypeToBoundaries()************************\n";
    amrex::Print() << "\t\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
    std::string prt =  "\t\t\t\t";
#endif
/*
    LinOpBCType::interior:
    LinOpBCType::Dirichlet:
    LinOpBCType::Neumann:
    LinOpBCType::reflect_odd:
    LinOpBCType::Marshak:
    LinOpBCType::SanchezPomraning:
    LinOpBCType::inflow:
    LinOpBCType::inhomogNeumann:
    LinOpBCType::Robin:
    LinOpBCType::Periodic:
*/
    
    auto& rCode = c_Code::GetInstance();
    auto& rBC = rCode.get_BoundaryConditions();
    auto& map_boundary_type = rBC.map_boundary_type;
    auto& bcType_2d = rBC.bcType_2d;
    auto& map_bcAny_2d = rBC.map_bcAny_2d;
   

    for (std::size_t i = 0; i < 2; ++i)
    {
        for (std::size_t j = 0; j < AMREX_SPACEDIM; ++j)
        {

            switch(map_boundary_type[ bcType_2d[i][j] ])
            {
                case s_BoundaryConditions::dir :
                {
                     LinOpBCType_2d[i][j] = LinOpBCType::Dirichlet;

                     if(map_bcAny_2d[i][j] == "inhomogeneous_constant") 
                     {
                         all_homogeneous_boundaries = false;
                         some_constant_inhomogeneous_boundaries = true;
                     }
                     if(map_bcAny_2d[i][j] == "inhomogeneous_function") 
                     {
                         all_homogeneous_boundaries = false;
                         some_functionbased_inhomogeneous_boundaries = true;
                     }
                     break;
                }
                case s_BoundaryConditions::neu :
                {
                     if(map_bcAny_2d[i][j] == "homogeneous")
                     {
                         LinOpBCType_2d[i][j] = LinOpBCType::Neumann;
                     }
                     else if(map_bcAny_2d[i][j] == "inhomogeneous_constant") 
                     {
                         LinOpBCType_2d[i][j] = LinOpBCType::inhomogNeumann;
                         all_homogeneous_boundaries = false;
                         some_constant_inhomogeneous_boundaries = true;
                     }
                     else if(map_bcAny_2d[i][j] == "inhomogeneous_function") 
                     {
                         LinOpBCType_2d[i][j] = LinOpBCType::inhomogNeumann;
                         all_homogeneous_boundaries = false;
                         some_functionbased_inhomogeneous_boundaries = true;
                     }
                     break;
                }
                case s_BoundaryConditions::per :
                {
                    LinOpBCType_2d[i][j] = LinOpBCType::Periodic;
                    break;
                }
                case s_BoundaryConditions::rob :
                {
                    LinOpBCType_2d[i][j] = LinOpBCType::Robin;
                    all_homogeneous_boundaries = false;
                    some_functionbased_inhomogeneous_boundaries = true;
                    break;
                }
                case s_BoundaryConditions::ref :
                {
                    LinOpBCType_2d[i][j] = LinOpBCType::reflect_odd;
                    break;
                }
            }

        }
    }

#ifdef PRINT_LOW
    for (std::size_t i = 0; i < 2; ++i)
    {
        for (std::size_t j = 0; j < AMREX_SPACEDIM; ++j)
        {
            amrex::Print() << prt << "side: " << i << " direction: " << j << " LinOpBCType: " << LinOpBCType_2d[i][j] << " type: " << map_bcAny_2d[i][j]  << "\n";
        }
    }
#endif 
#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t\t}************************c_MLMGSolver::AssignLinOpBCTypeToBoundaries()************************\n";
#endif
} 


void
c_MLMGSolver:: InitData()
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t{************************c_MLMGSolver::InitData()************************\n";
    amrex::Print() << "\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    Set_MLMG_CellCentered_Multifabs();
    Setup_MLABecLaplacian_ForPoissonEqn();

#ifdef PRINT_NAME
    amrex::Print() << "\t\t}************************c_MLMGSolver::InitData()************************\n";
#endif
}


void
c_MLMGSolver:: Set_MLMG_CellCentered_Multifabs()
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t{************************c_MLMGSolver::Set_MLMG_CellCentered_Multifabs()************************\n";
    amrex::Print() << "\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    auto& rCode = c_Code::GetInstance();
    auto& rMprop = rCode.get_MacroscopicProperties();

    alpha_cc = rMprop.get_p_mf(alpha_cc_str);
    beta_cc = rMprop.get_p_mf(beta_cc_str);
    soln_cc = rMprop.get_p_mf(soln_cc_str);
    rhs_cc  = rMprop.get_p_mf(rhs_cc_str);

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t}************************c_MLMGSolver::Set_MLMG_CellCentered_Multifabs()************************\n";
#endif
}


void
c_MLMGSolver:: Setup_MLABecLaplacian_ForPoissonEqn()
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t{************************c_MLMGSolver::Setup_MLABecLaplacian_ForPoissonEqn()************************\n";
    amrex::Print() << "\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif


    auto& rCode = c_Code::GetInstance();
    auto& rGprop = rCode.get_GeometryProperties();
    auto& ba = rGprop.ba;
    auto& dm = rGprop.dm;
    auto& geom = rGprop.geom;

    mlabec.define({geom}, {ba}, {dm}, info); // Implicit solve using MLABecLaplacian class

    // Force singular system to be solvable
    mlabec.setEnforceSingularSolvable(false);

    // set order of stencil
    mlabec.setMaxOrder(max_order);

    // assign domain boundary conditions to the solver
    // see Src/Boundary/AMReX_LO_BCTYPES.H for supported types
    mlabec.setDomainBC(LinOpBCType_2d[0], LinOpBCType_2d[1]);

    auto& rBC = rCode.get_BoundaryConditions();
    auto& map_boundary_type = rBC.map_boundary_type;
    auto& bcType_2d = rBC.bcType_2d;
    auto& map_bcAny_2d = rBC.map_bcAny_2d;

    // Fill the ghost cells of each grid from the other grids
    // includes periodic domain boundaries

    int amrlev = 0; //refers to the setcoarsest level of the solve

    if(some_constant_inhomogeneous_boundaries)
    {
        Fill_Constant_Inhomogeneous_Boundaries();
    }
    if(some_functionbased_inhomogeneous_boundaries) 
    {
        Fill_FunctionBased_Inhomogeneous_Boundaries();
    }
    soln_cc->FillBoundary(geom.periodicity());

    mlabec.setLevelBC(amrlev, soln_cc, robin_a, robin_b, robin_f);


    // set scaling factors 
    mlabec.setScalars(ascalar, bscalar);

    // set alpha_cc, and beta_fc coefficients
    mlabec.setACoeffs(amrlev, *alpha_cc);

    // beta_fc =  std::make_unique<amrex::FArrayBox,AMREX_SPACEDIM>; 
    // beta_fc is a face-centered multifab
    AMREX_D_TERM(beta_fc[0].define(convert(ba,IntVect(AMREX_D_DECL(1,0,0))), dm, 1, 0);,
                 beta_fc[1].define(convert(ba,IntVect(AMREX_D_DECL(0,1,0))), dm, 1, 0);,
                 beta_fc[2].define(convert(ba,IntVect(AMREX_D_DECL(0,0,1))), dm, 1, 0););

    Multifab_Manipulation::AverageCellCenteredMultiFabToCellFaces(*beta_cc, beta_fc); //converts from cell-centered permittivity to face-center and store in beta_fc

    mlabec.setBCoeffs(0, amrex::GetArrOfConstPtrs(beta_fc));

    pMLMG = std::make_unique<MLMG>(mlabec);

    pMLMG->setVerbose(set_verbose);

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t}************************c_MLMGSolver::Setup_MLABecLaplacian_ForPoissonEqn()************************\n";
#endif
}


void
c_MLMGSolver:: Fill_Constant_Inhomogeneous_Boundaries()
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t\t{************************c_MLMGSolver::Fill_Constant_Inhomogeneous_Boundaries()************************\n";
    amrex::Print() << "\t\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
    std::string prt = "\t\t\t\t";
#endif

    auto& rCode = c_Code::GetInstance();
    auto& rGprop = rCode.get_GeometryProperties();
    auto& geom = rGprop.geom;
    auto& real_box = geom.ProbDomain();
    auto dx = geom.CellSizeArray();
    Box const& domain = geom.Domain();

    auto& rBC = rCode.get_BoundaryConditions();
    auto& map_boundary_type = rBC.map_boundary_type;
    auto& bcType_2d = rBC.bcType_2d;
    auto& bcAny_2d = rBC.bcAny_2d;
    auto& map_bcAny_2d = rBC.map_bcAny_2d;

    std::vector<int> dir_inhomo_const_lo;
    std::string value = "inhomogeneous_constant";
    bool result = findByValue(dir_inhomo_const_lo, map_bcAny_2d[0], value);

#ifdef PRINT_LOW
    if(result)
    {
        amrex::Print() << "\n"<< prt <<"Low directions with value `"<< value << "' are:\n";
        for(auto dir : dir_inhomo_const_lo)  amrex::Print() << prt << "direction: " << dir << " boundary value: " << std::any_cast<amrex::Real>(bcAny_2d[0][dir]) << "\n";
    }
#endif

    std::vector<int> dir_inhomo_const_hi;
    result = findByValue(dir_inhomo_const_hi, map_bcAny_2d[1], value);

#ifdef PRINT_LOW
    if(result)
    {
        amrex::Print() << "\n"<< prt <<"High directions with value `"<< value << "' are:\n";
        for(auto dir : dir_inhomo_const_hi)  amrex::Print() << prt << "direction: " << dir  << " boundary value: " << std::any_cast<amrex::Real>(bcAny_2d[1][dir]) << "\n";
    }
#endif

    int len = 1;
    for (MFIter mfi(*soln_cc, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const auto& soln_arr = soln_cc->array(mfi);

        const auto& bx = mfi.tilebox();
        //const auto& gbx = mfi.growntilebox(1);


        for (auto dir : dir_inhomo_const_lo) {
            if (bx.smallEnd(dir) == domain.smallEnd(dir)) {
                Box const& bxlo = amrex::adjCellLo(bx, dir,len);
                amrex::ParallelFor(bxlo,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    soln_arr(i,j,k) = std::any_cast<amrex::Real>(bcAny_2d[0][dir]);
                });
            }
        }
        for (auto dir : dir_inhomo_const_hi) {
            if (bx.bigEnd(dir) == domain.bigEnd(dir)) {
                Box const& bxhi = amrex::adjCellHi(bx, dir,len);
                amrex::ParallelFor(bxhi,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    soln_arr(i,j,k) = std::any_cast<amrex::Real>(bcAny_2d[1][dir]);
                });
            }
        }
    }

#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t\t}************************c_MLMGSolver::Fill_Constant_Inhomogeneous_Boundaries()************************\n";
#endif
}


void
c_MLMGSolver:: Fill_FunctionBased_Inhomogeneous_Boundaries()
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t\t{************************c_MLMGSolver::Fill_FunctionBased_Inhomogeneous_Boundaries()************************\n";
    amrex::Print() << "\t\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

//TO BE CODED

#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t\t}************************c_MLMGSolver::Fill_FunctionBased_Inhomogeneous_Boundaries()************************\n";
#endif
}


void
c_MLMGSolver:: Solve_PoissonEqn()
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t{************************c_MLMGSolver::Solve_PoissonEqn()************************\n";
    amrex::Print() << "\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    pMLMG->solve({soln_cc}, {rhs_cc},
                 relative_tolerance,
                 absolute_tolerance);

#ifdef PRINT_NAME
    amrex::Print() << "\t\t}************************c_MLMGSolver::Solve_PoissonEqn()************************\n";
#endif
}


void
c_MLMGSolver:: Compute_vecField(std::array<amrex::MultiFab, AMREX_SPACEDIM>& vecField)
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t{************************c_MLMGSolver::Compute_vecField(*)************************\n";
    amrex::Print() << "\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    /** First the gradient of phi is computed then it is multiplied by one. */
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
    //TD<decltype(vecField)> vecField_type; //std::array<amrex::MultiFab, 3>&
    //TD<decltype(vecField[0])> vecField_0_type; //amrex::MultiFab&
    //TD<decltype(vecField[0][0])> vecField_00_type;//amrex::FArrayBox&

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        
        gradPhi[0][idim].define(amrex::convert(ba,amrex::IntVect::TheDimensionVector(idim)),dm, 1, 0);

        vecField[idim].define(amrex::convert(ba,IntVect::TheDimensionVector(idim)),dm, 1, 0);

    }
    pMLMG->getGradSolution(amrex::GetVecOfArrOfPtrs(gradPhi));

    /*Saxpy does E[x] += -1*gradPhi[0][x], where initial value of E[x]=0*/
    AMREX_D_TERM ( MultiFab::Saxpy (vecField[0], -1.0, gradPhi[0][0],0,0,1,0); ,
                   MultiFab::Saxpy (vecField[1], -1.0, gradPhi[0][1],0,0,1,0); ,
                   MultiFab::Saxpy (vecField[2], -1.0, gradPhi[0][2],0,0,1,0); );


#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t}************************c_MLMGSolver::Compute_vecField(*)************************\n";
#endif
}


void
c_MLMGSolver:: Compute_vecFlux(std::array<amrex::MultiFab, AMREX_SPACEDIM>& vecFlux)
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t{************************c_MLMGSolver::Compute_vecFlux(*)************************\n";
    amrex::Print() << "\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    /** Flux (-beta*grad phi)  is computed. */

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

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t}************************c_MLMGSolver::Compute_vecFlux(*)************************\n";
#endif
}

