#include <AMReX_ParmParse.H>
#include <map>
#include "Test.H"


using namespace amrex;

enum class Test: int
{
    MM_Mul, DM_Inv, ConjT, DM_Eig
};

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    int num_proc = ParallelDescriptor::NProcs(); 
    int my_rank = ParallelDescriptor::MyProc();
    amrex::Print() << "total number of procs: " << num_proc << "\n";

    amrex::ParmParse pp;

    int A_rows=4, A_cols=4, B_cols = 4; //assume, B_rows = A_cols
    pp.query("A_rows", A_rows);
    pp.query("A_cols", A_cols);
    pp.query("B_cols", B_cols);
    
    std::vector<int> print_flags(2,0);
    pp.query("print_input", print_flags[0]);
    pp.query("print_output", print_flags[1]);

    std::unordered_map<Test, bool> flag_map{
        {Test::MM_Mul, true}, 
        {Test::DM_Inv, true},
        {Test::ConjT,  true},
        {Test::DM_Eig, true}
    };
    amrex::ParmParse pp_test("test");
    pp_test.query("mm_mul" , flag_map[ Test::MM_Mul ]);
    pp_test.query("dm_inv" , flag_map[ Test::DM_Inv ]);
    pp_test.query("conjT"  , flag_map[ Test::ConjT  ]);
    pp_test.query("dm_eig" , flag_map[ Test::DM_Eig ]);
 
    {

        //Matrix-Matrix Multiplication
        if(flag_map[ Test::MM_Mul ]) 
        {
            Test_MathLibrary<Test_MM_Mul> 
                mm_mul(Test_MM_Mul(A_rows, A_cols, B_cols), print_flags);   

            mm_mul.Perform_Test();
        }

        //Dense Matrix Inversion
        if(flag_map[ Test::DM_Inv ]) 
        {
            Test_MathLibrary<Test_DM_Inv> 
                dm_inv(Test_DM_Inv(A_rows, A_cols), print_flags);   

            dm_inv.Perform_Test();
        }

        ////Conjugate Transpose 
        //if(flag_map[ Test::ConjT ]) 
        //{
        //    Test_MathLibrary<Test_ConjT> 
        //        conjT(Test_DM_Inv(A_rows, A_cols), print_flags);   

        //    conjT.Perform_Test();
        //}

        ////Dense Matrix Eigendecomposition
        //if(flag_map[ Test::DM_Eig ]) 
        //{
        //    Test_MathLibrary<Test_DM_Eig> 
        //        dm_eig(Test_DM_Eig(A_rows, A_cols), print_flags);   

        //    dm_eig.Perform_Test();
        //}
    }

    amrex::Finalize();
}
