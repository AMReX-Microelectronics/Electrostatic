#include "PostProcessor.H"

#include "Code.H"
#include "Solver/Electrostatics/MLMG.H"



using namespace amrex;


c_PostProcessor::c_PostProcessor ()
{
    DefineParameterNameMap();
    num_all_params = ReadData();
    num_arraymf_params = SortParameterType();
    num_mf_params = num_all_params - num_arraymf_params;

    //amrex::Print() << "num_all_params: " << num_all_params << "\n";
    //amrex::Print() << "num_arraymf_params: " << num_arraymf_params << "\n";
    //amrex::Print() << "num_mf_params: " << num_mf_params << "\n";

    m_p_mf.resize(num_mf_params);
    m_p_array_mf.resize(num_arraymf_params);

} 


c_PostProcessor::~c_PostProcessor ()
{
//    for (auto& elem : m_p_mf) {
//        elem.release();
//    }
//    for (auto& elem : m_p_array_mf) {
//        elem.release();
//    }
} 

int 
c_PostProcessor::ReadData()
{ 

   map_param_all["vecE"] = 0;
   map_param_all["vecFlux"] = 1;

   return map_param_all.size();

}


int 
c_PostProcessor::SortParameterType()
{ 

  int i=0;
  int j=0;
  for (auto it : map_param_all) {
       std::string str = it.first;
       std::string first_three_char_str = str.substr(0, 3);

       if(first_three_char_str == "vec"){
           map_param_arraymf[it.first] = i;
           //amrex::Print() << "key: " << it.first << "  value: " << it.second << "  dimension: " << AMREX_SPACEDIM << "  counter, i: "<< i<< "\n";
           ++i;
       }
       else {
           map_param_mf[it.first] = j;
           //amrex::Print() << "key: " << it.first << "  value: " << it.second << "  dimension: " << AMREX_SPACEDIM << "  counter, j: "<< j<< "\n";
           ++j;
       }
   }

   return i;

}


void 
c_PostProcessor::InitData()
{

    for (auto it : map_param_mf) {
        m_p_mf[it.second] = std::make_unique<amrex::MultiFab>(); //cell-centered multifab
    }
    for (auto it : map_param_arraymf) {
        m_p_array_mf[it.second] = std::make_unique<  std::array<amrex::MultiFab, AMREX_SPACEDIM>  >();
    }

}


void 
c_PostProcessor::Compute(std::string macro_str)
{
 
    auto macro_num = map_param_all[macro_str];
    auto& rCode = c_Code::GetInstance();

    std::map<std::string,s_PostProcessMacroName::macro_name>::iterator it_Post;

    it_Post = map_macro_name.find(macro_str);

    if(it_Post == map_macro_name.end())
    {
        amrex::Print() << "Computation of "<< macro_str << " is not implemented at present.\n";
    }
    else {

        switch(map_macro_name[macro_str])
        {
            case s_PostProcessMacroName::vecE : 
            {
                auto val = map_param_arraymf[macro_str];    
                auto& rMLMGsolver = rCode.get_MLMGSolver();
                rMLMGsolver.Compute_vecE( *m_p_array_mf[val] );
                break;
            }
            case s_PostProcessMacroName::vecFlux : 
            {
                auto val = map_param_arraymf[macro_str];    
                auto& rMLMGsolver = rCode.get_MLMGSolver();
                rMLMGsolver.Compute_vecFlux( *m_p_array_mf[val] );
                break;
            }
        }
    }
 
}
