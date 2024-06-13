#include "MacroscopicProperties.H"

#include "../../Utils/SelectWarpXUtils/WarpXUtil.H"
#include "../../Utils/CodeUtils/CodeUtil.H"

#include "Code.H"
#include "GeometryProperties.H"

#include <AMReX_ParmParse.H>
#include <AMReX_Parser.H>

#include <ctype.h>

using namespace amrex;

c_MacroscopicProperties::c_MacroscopicProperties ()
{

    DefineParameterNameMap();
    DefineDefaultValueMap();
    ReadData();

} 


void 
c_MacroscopicProperties::ReadData()
{ 
    num_params = ReadParameterMapAndNumberOfGhostCells();

    DefineMacroVariableVectorSizes();
    std::map<std::string,amrex::Real>::iterator it_default;

    for (auto it: map_param_all)
    {
        amrex::Real default_val;

        it_default = map_default_value.find(it.first);

        if (it_default == map_default_value.end()) {
            default_val = 0.0;
        }
        else {
            default_val = map_default_value[it.first];
        }
        ReadMacroparam(it.first, default_val);
    }

    Define_ExternalChargeDensitySources();
}


int 
c_MacroscopicProperties::ReadParameterMapAndNumberOfGhostCells()
{ 

    amrex::Vector< std::string > fields_to_define;
    amrex::Vector< std::string > ghostcells_for_fields;

    amrex::ParmParse pp_macro("macroscopic");

    bool varnames_specified = pp_macro.queryarr("fields_to_define", fields_to_define);
    bool ghostcells_specified = pp_macro.queryarr("ghostcells_for_fields", ghostcells_for_fields);

    std::map<std::string,int>::iterator it_map_param_all;

    int c=0;
    for (auto it: fields_to_define)
    {

        it_map_param_all = map_param_all.find(it);

        if (it_map_param_all == map_param_all.end()) {
                map_param_all[it] = c;
                ++c;
        }
    }
    fields_to_define.clear();

    for (auto it: map_param_all)
    {
        std::string str = it.first;
        bool comparison_true=false;
        std::string str_with_ghost;
        for (auto str_ghost: ghostcells_for_fields) {
             int compare = strncmp(str.c_str(), str_ghost.c_str(), str.length());
             if(compare == 0){
               comparison_true = true; 
               str_with_ghost = str_ghost;
               break;
             }   
        }
        if (comparison_true) {
            std::string num_ghost_str = str_with_ghost.substr(str.length()+1);

            if(std::isdigit(*num_ghost_str.c_str())) {
                amrex::Real num_ghost = std::stod(num_ghost_str);
                map_num_ghostcell[str] = num_ghost;
            } 
            else {
                map_num_ghostcell[str] = 0;
                std::stringstream warnMsg;
                warnMsg << "Macroscopic specification: '" << str_with_ghost  << "' does not contain number of ghost cells.\n"
                << "Default ghost cell value of 0 is chosen\n";
                c_Code::GetInstance().RecordWarning("Macroscopic properties", warnMsg.str());
            }
        }
        else {
            map_num_ghostcell[str] = 0;
            std::stringstream warnMsg;
            warnMsg << "For macroscopic property: '" << str  << "' ghost cell value is not specified.\n"
            << "Default ghost cell value of 0 is chosen\n";
            c_Code::GetInstance().RecordWarning("Macroscopic properties", warnMsg.str());
        }
        
    }
    ghostcells_for_fields.clear();

    amrex::Print() <<  "\n##### Macroscopic Properties #####\n\n";
    amrex::Print() <<  "##### fields_to_define: number, name,  ghost_cells\n";
    amrex::Print() << "##### " <<  "number" << std::setw(5) << "name" << std::setw(15) << "ghost_cells" << "\n";
    for (auto it: map_param_all)
    {
        amrex::Print() << "##### " <<  it.second << std::setw(20) << it.first << std::setw(5) << map_num_ghostcell[it.first] << "\n";
    }
    amrex::Print() << "##### Total number of parameters: " << map_param_all.size() << "\n\n";


    return map_param_all.size();

}


void 
c_MacroscopicProperties::DefineMacroVariableVectorSizes()
{ 
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t\t\t{************************c_MacroscopicProperties::DefineMacroVariableVectorSizes()************************\n";
    amrex::Print() << "\t\t\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif
    
    m_macro_type.resize(num_params, "constant"); //initialized to constant
    m_macro_value.resize(num_params);
    m_macro_str_function.resize(num_params);
    m_p_macro_parser.resize(num_params);
    m_p_mf.resize(num_params);

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t\t\t}************************c_MacroscopicProperties::DefineMacroVariableVectorSizes()************************\n";
#endif
}


void 
c_MacroscopicProperties::InitData()
{
    const int Ncomp1=1;

    auto& rCode = c_Code::GetInstance();
    auto& rGprop = rCode.get_GeometryProperties();
    auto& ba = rGprop.ba;
    auto& dm = rGprop.dm;
    auto& geom = rGprop.geom;

    amrex::Print()  << "\n##### MACROSCOPIC INITIALIZATION #####\n\n";
    for (auto it: map_param_all)
    {
        auto macro_str = it.first;
        auto macro_num = it.second;
        DefineAndInitializeMacroparam(macro_str, macro_num, ba, dm, geom, Ncomp1, map_num_ghostcell[macro_str]);
    }
    
    Deposit_ExternalChargeDensitySources();
}


template < class T >
void 
c_MacroscopicProperties::ReadMacroparam(std::string macro_str, 
                                        T default_value)
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t\t\t{************************c_MacroscopicProperties::ReadMacroparam()************************\n";
    amrex::Print() << "\t\t\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    auto macro_num = map_param_all[macro_str];
    m_macro_value[macro_num] = default_value;

    ParmParse pp_macroscopic("macroscopic");

    bool specified = false; 
    std::string macro_functionXYZ = macro_str+"_function";
    if (queryWithParser(pp_macroscopic, macro_str.c_str() , m_macro_value[macro_num]) ) {
        m_macro_type[macro_num] = "constant";
        specified = true;
    }
    if (pp_macroscopic.query( macro_functionXYZ.c_str(), m_macro_str_function[macro_num]) ) {
        m_macro_type[macro_num] = "parse_" + macro_str + "_function";
        specified = true;
    }
    if (!specified) {
        std::stringstream warnMsg;
        warnMsg << "Macroscopic parameter '" << macro_str << "' is not specified in the input file. The default value of " 
                <<  m_macro_value[macro_num]
                << " is used.";
        c_Code::GetInstance().RecordWarning("Macroscopic properties", warnMsg.str());
    }

    if (m_macro_type[macro_num] == "parse_" + macro_str + "_function") 
    { 
        Store_parserString(pp_macroscopic, macro_functionXYZ.c_str(),  m_macro_str_function[macro_num]);

        m_p_macro_parser[macro_num] = std::make_unique<amrex::Parser>(
                                           makeParser( m_macro_str_function[macro_num], m_parser_varname_vector));
    }

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t\t\t}************************c_MacroscopicProperties::ReadMacroparam()************************\n";
#endif
}


void 
c_MacroscopicProperties::DefineAndInitializeMacroparam(std::string macro_str, 
                                                       int num,
                                                       amrex::BoxArray& ba, 
                                                       amrex::DistributionMapping& dm, 
                                                       amrex::Geometry& geom, 
                                                       int Ncomp, 
                                                       int Nghost)
{

    auto macro_num = map_param_all[macro_str];
    //amrex::Print()  << " Initializing macro_str: " << macro_str << " macro_num: " << macro_num << " macro_type: " << m_macro_type[macro_num] << "\n";

    m_p_mf[macro_num] = std::make_unique<amrex::MultiFab>(ba, dm, Ncomp, Nghost); //cell-centered multifab

    if (m_macro_type[macro_num] == "constant") {
        amrex::Print()  << "##### " << macro_str << " initialized with a constant value: " << m_macro_value[macro_num] << "\n";
        m_p_mf[macro_num] -> setVal(m_macro_value[macro_num]);

    } else if (m_macro_type[macro_num] == "parse_" + macro_str + "_function") {
        amrex::Print() << "##### " << macro_str << " initialized with a parser function with name: " << m_macro_type[macro_num] << "\n";

        auto& rCode = c_Code::GetInstance();
        const amrex::Real time = rCode.get_time();

        #ifdef TIME_DEPENDENT
           Multifab_Manipulation::InitializeMacroMultiFabUsingParser_4vars(m_p_mf[macro_num].get(), m_p_macro_parser[macro_num]->compile<4>(), geom, time);
       #else
           Multifab_Manipulation::InitializeMacroMultiFabUsingParser_3vars(m_p_mf[macro_num].get(), m_p_macro_parser[macro_num]->compile<3>(), geom);
       #endif 
    }

}


void 
c_MacroscopicProperties::ReInitializeMacroparam(std::string macro_str)
{
    auto macro_num = map_param_all[macro_str];
    auto& rCode = c_Code::GetInstance();
    auto& rGprop = rCode.get_GeometryProperties();
    auto& geom = rGprop.geom;
    auto& ba = rGprop.ba;
    auto& dm = rGprop.dm;
    auto Nghost = map_num_ghostcell[macro_str];
    auto Ncomp = 1;

    if (m_macro_type[macro_num] == "constant") 
    {
        m_p_mf[macro_num] -> setVal(m_macro_value[macro_num]);
    } else 
    if (m_macro_type[macro_num] == "parse_" + macro_str + "_function") {

        auto& rCode = c_Code::GetInstance();
        const amrex::Real time = rCode.get_time();

        #ifdef TIME_DEPENDENT
           Multifab_Manipulation::InitializeMacroMultiFabUsingParser_4vars(m_p_mf[macro_num].get(), m_p_macro_parser[macro_num]->compile<4>(), geom, time);
       #else
           Multifab_Manipulation::InitializeMacroMultiFabUsingParser_3vars(m_p_mf[macro_num].get(), m_p_macro_parser[macro_num]->compile<3>(), geom);
       #endif 
    }

    m_p_mf[macro_num] -> FillBoundary(geom.periodicity());
}


void
c_MacroscopicProperties::Define_ExternalChargeDensitySources()
{

    amrex::Print()  << "\n##### EXTERNAL CHARGE DENSITY SOURCES #####\n\n";

    amrex::ParmParse pp_cds("charge_density_source");

    int num = 0;
    std::string type_str;
    auto num_isDefined  = pp_cds.query("num", num);
    auto type_isDefined = pp_cds.query("type", type_str);
    
    amrex::Print() << " charge_density_source.num " << num << "\n";
    amrex::Print() << " charge_density_source.type " << type_str << "\n";
    if(type_isDefined) 
    {
        if(type_str == "point_charge") 
        {
            if(p_pointChargeSource==nullptr) 
                p_pointChargeSource = std::make_unique<PointChargeSource>(num);

            for(size_t s=0; s<num; ++s) 
            {
                amrex::ParmParse pp_pc("pc_" + std::to_string(s+1));

                amrex::Vector<amrex::Real> pos(AMREX_SPACEDIM);

                getArrWithParser(pp_pc, "location",
                                 pos, 
                                 0, AMREX_SPACEDIM);
    
                amrex::Real sigma = 2.e-10;
                pp_pc.query("sigma", sigma);
                
                int charge_unit = 1;
                pp_pc.query("charge_unit", charge_unit);

                p_pointChargeSource->Define_PointCharge(s, 
                                       PointCharge(pos.data(), sigma, charge_unit));
                       
            }
            #ifdef AMREX_USE_GPU
            p_pointChargeSource->Copy_HostVectorToGPU();
            #endif
            p_pointChargeSource->Print_PointCharge();
        }
    }
}


void
c_MacroscopicProperties::Deposit_ExternalChargeDensitySources()
{
    amrex::Print() << "In Deposit_External\n";
    auto& rCode  = c_Code::GetInstance();
    auto& rGprop = rCode.get_GeometryProperties();
    auto& geom = rGprop.geom;
    auto dx = geom.CellSizeArray();
    auto& real_box = geom.ProbDomain();

    amrex::MultiFab* const p_charge_density_mf = get_p_mf("charge_density");
    
    if(p_charge_density_mf == nullptr) 
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(p_charge_density_mf != nullptr, 
                "charge_density mf string wrong!");


    auto iv = p_charge_density_mf->ixType().toIntVect();

    for ( amrex::MFIter mfi(*p_charge_density_mf, amrex::TilingIfNotGPU()); 
                        mfi.isValid(); 
                        ++mfi )
    {
        const auto& tb = mfi.tilebox( iv, p_charge_density_mf->nGrowVect() );
        auto const& mf_array = p_charge_density_mf->array(mfi);

        auto num_source = p_pointChargeSource->Get_NumSources();
        //auto* p_vec_source = p_pointChargeSource->Get_pVecSources();
        auto* p_vec_source = p_pointChargeSource->d_vec_source.dataPtr();
        auto get_charge = get_charge_sum();

        //auto* const p_charge_source = this->p_pointChargeSource.get();

        amrex::ParallelFor (tb,
            [=] 
            AMREX_GPU_DEVICE (int i, int j, int k)
        {
                amrex::Real fac_x = (1._rt - iv[0]) * dx[0] * 0.5_rt;
                amrex::Real x = i * dx[0] + real_box.lo(0) + fac_x;

                amrex::Real fac_y = (1._rt - iv[1]) * dx[1] * 0.5_rt;
                amrex::Real y = j * dx[1] + real_box.lo(1) + fac_y;

                amrex::Real fac_z = (1._rt - iv[2]) * dx[2] * 0.5_rt;
                amrex::Real z = k * dx[2] + real_box.lo(2) + fac_z;

                for(int s=0; s< num_source; ++s) {
                    mf_array(i,j,k) += get_charge(p_vec_source[s],x,y,z);
                }
                //mf_array(i,j,k) += p_charge_source->get_charge_sum(x,y,z);
        });
    }
//    amrex::Gpu::streamSynchronize();
//    p_charge_density_mf->FillBoundary(rGprop.geom.periodicity());
}
