#include "Code.H"

#include "Utils/MsgLogger/MsgLogger.H"
#include "Utils/WarnManager.H"
#include "Utils/WarpXUtil.H"
#include "Utils/WarpXProfilerWrapper.H"

#include "Input/GeometryProperties.H"
#include "Input/MacroscopicProperties.H"
#include "Solver/Electrostatics/MLMG.H"
#include "PostProcessor/PostProcessor.H"
#include "Output/Output.H"



c_Code* c_Code::m_instance = nullptr;
#ifdef AMREX_USE_GPU
bool c_Code::do_device_synchronize = true;
#else
bool c_Code::do_device_synchronize = false;
#endif

c_Code& c_Code::GetInstance() 
{

    if (!m_instance) {
        m_instance = new c_Code();
    }
    return *m_instance;

}


void
c_Code::ResetInstance ()
{
    delete m_instance;
    m_instance = nullptr;
}


c_Code::c_Code ()
{
    m_instance = this;
    m_p_warn_manager = std::make_unique<Utils::WarnManager>();

    ReadData();

}


c_Code::~c_Code ()
{
//
}


void
c_Code::RecordWarning(
        std::string topic,
        std::string text,
        WarnPriority priority)
{
    WARPX_PROFILE("WarpX::RecordWarning");

    auto msg_priority = Utils::MsgLogger::Priority::high;
    if(priority == WarnPriority::low)
        msg_priority = Utils::MsgLogger::Priority::low;
    else if(priority == WarnPriority::medium)
        msg_priority = Utils::MsgLogger::Priority::medium;

    if(m_always_warn_immediately){
        amrex::Warning(
            "!!!!!! WARNING: ["
            + std::string(Utils::MsgLogger::PriorityToString(msg_priority))
            + "][" + topic + "] " + text);
    }

#ifdef AMREX_USE_OMP
    #pragma omp critical
#endif
    {
        m_p_warn_manager->record_warning(topic, text, msg_priority);
    }
}


void
c_Code::PrintLocalWarnings(const std::string& when)
{

    WARPX_PROFILE("WarpX::PrintLocalWarnings");
    const std::string warn_string = m_p_warn_manager->print_local_warnings(when);
    amrex::AllPrint() << warn_string;

}


void
c_Code::PrintGlobalWarnings(const std::string& when)
{

    WARPX_PROFILE("WarpX::PrintGlobalWarnings");
    const std::string warn_string = m_p_warn_manager->print_global_warnings(when);
    amrex::Print() << warn_string;

}


void 
c_Code::ReadData ()
{

     m_pGeometryProperties = std::make_unique<c_GeometryProperties>();
     
     m_pMacroscopicProperties = std::make_unique<c_MacroscopicProperties>();
     
     m_pMLMGSolver = std::make_unique<c_MLMGSolver>();
     
     m_pPostProcessor = std::make_unique<c_PostProcessor>();

     m_pOutput = std::make_unique<c_Output>();

}


void 
c_Code::InitData ()
{
 
    m_pGeometryProperties->InitData();

    m_pMacroscopicProperties->InitData();

    m_pMLMGSolver->Setup_MLABecLaplacian_ForPoissonEqn();

    m_pPostProcessor->InitData();

    m_pOutput->InitData();

}


void 
c_Code::Solve()
{
//    auto& phi = m_pMacroscopicProperties->get_mf("phi");
//    const auto& phi_arr = phi[0].array();
//    amrex::Print() << "\nphi BEFORE Poisson solve: \n";
//    amrex::Print() << "phi_0,0,0 :  "   << phi_arr(0,0,0) << "\n";
//    amrex::Print() << "phi_15,49,49:  " << phi_arr(15,49,49) << "\n";
//    amrex::Print() << "phi_24,49,49:  " << phi_arr(24,49,49) << "\n";
//    amrex::Print() << "phi_25,49,49:  " << phi_arr(25,49,49) << "\n";
//    amrex::Print() << "phi_49,49,49:  " << phi_arr(49,49,49) << "\n";
//    amrex::Print() << "phi_74,49,49:  " << phi_arr(74,49,49) << "\n";
//    amrex::Print() << "phi_75,49,49:  " << phi_arr(75,49,49) << "\n";
//    amrex::Print() << "phi_85,49,49:  " << phi_arr(85,49,49) << "\n";

    m_pMLMGSolver->Solve_PoissonEqn();
//    amrex::Print() << "\nphi AFTER Poisson solve: \n";
//    amrex::Print() << "phi_0,0,0 :  "   << phi_arr(0,0,0) << "\n";
//    amrex::Print() << "phi_15,49,49:  " << phi_arr(15,49,49) << "\n";
//    amrex::Print() << "phi_24,49,49:  " << phi_arr(24,49,49) << "\n";
//    amrex::Print() << "phi_25,49,49:  " << phi_arr(25,49,49) << "\n";
//    amrex::Print() << "phi_49,49,49:  " << phi_arr(49,49,49) << "\n";
//    amrex::Print() << "phi_74,49,49:  " << phi_arr(74,49,49) << "\n";
//    amrex::Print() << "phi_75,49,49:  " << phi_arr(75,49,49) << "\n";
//    amrex::Print() << "phi_85,49,49:  " << phi_arr(85,49,49) << "\n";
}


void 
c_Code::PostProcess()
{
    
    m_pPostProcessor->Compute("vecE"); 
    m_pPostProcessor->Compute("vecFlux"); 

//    auto& Ex = m_pPostProcessor->get_array_mf_component("vecE", 0);
//    const auto& Ex_arr = Ex[0].array();
//    amrex::Print() << "\nEx field: \n";
//    amrex::Print() << "Ex_0,0,0 :  "   << Ex_arr(0,0,0) << "\n";
//    amrex::Print() << "Ex_15,49,49:  " << Ex_arr(15,49,49) << "\n";
//    amrex::Print() << "Ex_24,49,49:  " << Ex_arr(24,49,49) << "\n";
//    amrex::Print() << "Ex_49,49,49:  " << Ex_arr(49,49,49) << "\n";
//    amrex::Print() << "Ex_76,49,49:  " << Ex_arr(76,49,49) << "\n";
//    amrex::Print() << "Ex_85,49,49:  " << Ex_arr(85,49,49) << "\n";
//    amrex::Print() << "Ex_100,49,49:  " << Ex_arr(100,49,49) << "\n";
//    amrex::Print() << "Ex_101,49,49:  " << Ex_arr(101,49,49) << "\n";

}


void 
c_Code::Output()
{

    int step=0;
    amrex::Real time=0.0;
    m_pOutput->WriteSingleLevelPlotFile(step, time);

}
