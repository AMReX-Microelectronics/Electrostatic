#include "Code.H"

#include "Utils/MsgLogger/MsgLogger.H"
#include "Utils/WarnManager.H"
#include "Utils/WarpXUtil.H"
#include "Utils/WarpXProfilerWrapper.H"

#include "Input/MaterialProperties.H"
#include "Input/GeometryProperties.H"



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

     m_pGeometryProperties.reset(); 
     m_pMaterialProperties.reset();

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
     m_pMaterialProperties = std::make_unique<c_MaterialProperties>();

}


void 
c_Code::InitData ()
{
 
    m_pGeometryProperties->InitData();
    m_pMaterialProperties->InitData();

}

