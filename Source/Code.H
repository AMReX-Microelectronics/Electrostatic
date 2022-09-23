#ifndef CODE_H_
#define CODE_H_

#include "Input/GeometryProperties_fwd.H"
#include "Input/MaterialProperties_fwd.H"
#include "Utils/WarnManager_fwd.H"


#include <AMReX.H>
//#include <AMReX_AmrCore.H>
//#include <AMReX_Array.H>
//#include <AMReX_Config.H>
#include <AMReX_GpuContainers.H>
//#include <AMReX_IntVect.H>
//#include <AMReX_LayoutData.H>
//#include <AMReX_Parser.H>
#include <AMReX_REAL.H>
//#include <AMReX_RealBox.H>
//#include <AMReX_RealVect.H>
//#include <AMReX_Vector.H>
//#include <AMReX_VisMF.H>

//#include <array>
//#include <iostream>
//#include <limits>
//#include <memory>
//#include <optional>
//#include <string>
//#include <vector>


//#include <AMReX_AmrCoreFwd.H>


enum class WarnPriority
{
    /** Low priority warning:
     * essentially an informative message
     */
    low,
    /** Medium priority warning:
     * a bug or a performance issue may affect the simulation
     */
    medium,
    /** High priority warning:
     * a very serious bug or performance issue
     * almost certainly affects the simulation
     */
    high
};

class c_Code
{

public:

    static c_Code& GetInstance ();
    static void ResetInstance ();

    c_Code ();
    ~c_Code ();


    /**
     * \brief This function is borrowed from WarpX. It records a warning message.
     * RecordWarning is thread safe: it can be used within OpenMP parallel loops.
     *
     * @param[in] topic a string to identify the topic of the warning (e.g., "parallelization", "pbc", "particles"...)
     * @param[in] text the text of the warning message
     * @param[in] priority priority of the warning message ("medium" by default)
     */
    void RecordWarning(
        std::string topic,
        std::string text,
        WarnPriority priority = WarnPriority::medium);

    /**
     * \brief This function prints all the warning messages collected on the present MPI rank
     * (i.e., this is not a collective call). This function is mainly intended for debug purposes.
     *
     * @param[in] when a string to mark when the warnings are printed out (it appears in the warning list)
     */
    void PrintLocalWarnings(const std::string& when);

    /**
     * \brief This function prints all the warning messages collected by all the MPI ranks
     * (i.e., this is a collective call). Only the I/O rank prints the message.
     *
     * @param[in] when a string to mark when the warnings are printed out (it appears in the warning list)
     */
    void PrintGlobalWarnings(const std::string& when);


    void ReadData ();
    void InitData ();

    static bool do_device_synchronize;

    c_GeometryProperties& get_GeometryProperties () { return *m_pGeometryProperties;}
    c_MaterialProperties& get_MaterialProperties () { return *m_pMaterialProperties;}


    // Warning manager: it allows recording and printing error messages
    std::unique_ptr<Utils::WarnManager> m_p_warn_manager;
    // Flag to control if WarpX has to emit a warning message as soon as a warning is recorded
    bool m_always_warn_immediately = false;

private:
    static c_Code* m_instance;

    std::unique_ptr<c_GeometryProperties> m_pGeometryProperties;
    std::unique_ptr<c_MaterialProperties> m_pMaterialProperties;

};
#endif
