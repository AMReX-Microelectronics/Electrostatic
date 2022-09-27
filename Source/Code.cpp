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
    m_pMLMGSolver->Solve_PoissonEqn();
}


void 
c_Code::PostProcess()
{
    
    m_pPostProcessor->Compute("vecE"); 
    m_pPostProcessor->Compute("vecFlux"); 

}


void 
c_Code::Output()
{

    int step=0;
    amrex::Real time=0.0;
    m_pOutput->WriteSingleLevelPlotFile(step, time);

}
