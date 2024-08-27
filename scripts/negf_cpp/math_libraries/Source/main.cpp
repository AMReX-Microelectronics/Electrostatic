#include <AMReX_ParmParse.H>

#include "Test.H"

using namespace amrex;

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
 
    {
        Test_MathLibrary<Test_MM_Mul> 
            mm_mul(Test_MM_Mul(A_rows, A_cols, B_cols, print_flags));   

        mm_mul.Perform_Test();

        mm_mul.Print_Output();
    }

    amrex::Finalize();
}
