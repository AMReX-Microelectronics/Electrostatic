#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Parser.H>


#include <ctype.h>
#include <typeindex>
#include <typeinfo>
#include <unordered_map>
#include <string>
#include <vector>
#include <any>


namespace Multifab_Manipulation
{
void InitializeMacroMultiFabUsingParser_3vars (
                       amrex::MultiFab *macro_mf,
                       amrex::ParserExecutor<3> const& macro_parser,
                       amrex::Geometry& geom);

void AverageCellCenteredMultiFabToCellFaces(const amrex::MultiFab& cc_arr,
                                            std::array< amrex::MultiFab, 
                                            AMREX_SPACEDIM >& face_arr);
}

void PrintRunDiagnostics(amrex::Real initial_time);


template<class T, class F>
inline std::pair<const std::type_index, std::function<void(std::any const&)>>
    to_any_visitor(F const &f)
{
    return {
        std::type_index(typeid(T)),
        [g = f](std::any const &a)
        {
            if constexpr (std::is_void_v<T>)
                g();
            else
                g(std::any_cast<T const&>(a));
        }
    };
}


static std::unordered_map<
    std::type_index, std::function<void(std::any const&)>>
    any_visitor {
        to_any_visitor<float>([](float x){ amrex::Print() << "contains a value: " << x; }),
        to_any_visitor<double>([](double x){ amrex::Print() << "contains a value: " << x; }),
        to_any_visitor<std::string>([](std::string s)
            { amrex::Print() << "contains a function parser name: "<< s; })
        // ... add more handlers for your types ...
    };

inline void process_std_any(const std::any& a)
{
    if (const auto it = any_visitor.find(std::type_index(a.type()));
        it != any_visitor.cend()) {
        it->second(a);
    }
    else {
//        std::cout << "unregistered type! "<< std::quoted(a.type().name());
        amrex::Print() << "contents do not matter! ";
    }
}

template<typename K, typename V> bool findByValue(std::vector<K> & vec, std::map<K, V> mapOfElemen, V value)
{
    bool bResult = false;
    auto it = mapOfElemen.begin();

    while(it != mapOfElemen.end())
    {
        if(it->second == value)
        {
            bResult = true;
            vec.push_back(it->first);
        }
        it++;
    }
    return bResult;
}

