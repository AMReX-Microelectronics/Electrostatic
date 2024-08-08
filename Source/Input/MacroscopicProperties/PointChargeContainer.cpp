#include "PointChargeContainer.H"
#include "Code.H"
#include "GeometryProperties.H"
#include "../../Utils/SelectWarpXUtils/WarpXUtil.H"
#include "../../Utils/SelectWarpXUtils/WarpXConst.H"
#include "../../Utils/CodeUtils/MathFunctions.H"
#include "../../Utils/CodeUtils/CodeUtil.H"

#include <AMReX_ParmParse.H>
#include <AMReX_Parser.H>

amrex::Real c_PointChargeContainer::V0 = 0.;
amrex::Real c_PointChargeContainer::Et = 0.;
std::mutex c_PointChargeContainer::mtx;

c_PointChargeContainer::c_PointChargeContainer(const amrex::Geometry& geom,
                           const amrex::DistributionMapping& dm,
                           const amrex::BoxArray& ba)
    : amrex::ParticleContainer<PCA::realPD::NUM, PCA::intPD::NUM,
                               PCA::realPA::NUM, PCA::intPA::NUM>(geom, dm, ba)
{
    Read_PointCharges();
    amrex::Print() << "\nPrinting just after read:\n";
    Print_Container(true);
}


void c_PointChargeContainer::Set_V0(amrex::Real value) 
{
    std::lock_guard<std::mutex> lock(mtx);   
    V0 = value;
}


void c_PointChargeContainer::Set_Et(amrex::Real value) 
{
    std::lock_guard<std::mutex> lock(mtx);   
    Et = value;
}


void c_PointChargeContainer::Read_PointCharges()
{
    amrex::Vector<amrex::Real> charge_units;
    amrex::Vector<amrex::Real> occupations;

    amrex::ParmParse pp_main("pc");

    amrex::Vector<amrex::Real> vec_offset(AMREX_SPACEDIM, 0.0);
    queryArrWithParser(pp_main, "offset", vec_offset, 0, AMREX_SPACEDIM);

    amrex::Vector<amrex::Real> vec_scaling(AMREX_SPACEDIM, 1.0);
    queryArrWithParser(pp_main, "scaling", vec_scaling, 0, AMREX_SPACEDIM);

    amrex::Real V0_val=0;
    pp_main.query("V0", V0_val);
    Set_V0(V0_val);
    
    amrex::Real Et_val=0;
    pp_main.query("Et", Et_val);
    Set_Et(Et_val);

    pp_main.query("flag_vary_occupation", flag_vary_occupation);

    amrex::Real default_charge_unit=1.;
    pp_main.query("default_charge_unit", default_charge_unit);

    amrex::Real default_occupation=1.;
    pp_main.query("default_occupation", default_occupation);

    int num = 0;
    pp_main.get("num", num);

    amrex::Vector<amrex::Vector<amrex::Real>> particles(num, amrex::Vector<amrex::Real>(AMREX_SPACEDIM, 0.0));

    for (int s = 0; s < num; ++s)
    {
        amrex::ParmParse pp_pc("pc_" + std::to_string(s + 1));

        amrex::Vector<amrex::Real> pos(AMREX_SPACEDIM);
        getArrWithParser(pp_pc, "location", pos, 0, AMREX_SPACEDIM);

        for (int k = 0; k < AMREX_SPACEDIM; ++k)
        {
            particles[s][k] = pos[k] * vec_scaling[k] + vec_offset[k];
        }

        amrex::Real charge_unit = default_charge_unit;
        pp_pc.query("charge_unit", charge_unit);
        charge_units.push_back(charge_unit);

        amrex::Real occupation = default_occupation;
        pp_pc.query("occupation", occupation);
        occupations.push_back(occupation);
   
    }

    Define(particles, charge_units, occupations);
}

void c_PointChargeContainer::Define(const amrex::Vector<amrex::Vector<amrex::Real>>& particles, 
                                    const amrex::Vector<amrex::Real>& charge_units,
                                    const amrex::Vector<amrex::Real>& occupations)
{

    if (ParallelDescriptor::IOProcessor())
    {

        for (int p=0; p<particles.size(); ++p)
        {
            ParticleType new_particle;
            new_particle.id() = ParticleType::NextID();
            new_particle.cpu() = ParallelDescriptor::MyProc();

            for (int j = 0; j < AMREX_SPACEDIM; ++j)
            {
                new_particle.pos(j) = particles[p][j];
            }

            std::array<int, PCA::intPA::NUM> int_attribs = {};
            int_attribs[PCA::intPA::id] = new_particle.id();

            std::array<ParticleReal, PCA::realPA::NUM> real_attribs = {};
            real_attribs[PCA::realPA::charge_unit]= charge_units[p]; 
            real_attribs[PCA::realPA::occupation] = occupations[p]; 
            real_attribs[PCA::realPA::potential]  = 0.0; 

            std::pair<int, int> key {0, 0}; // {grid_index, tile_index}
            int lev = 0;
            auto& particle_tile = GetParticles(lev)[key];

            particle_tile.push_back(new_particle);
            particle_tile.push_back_int(int_attribs);
            particle_tile.push_back_real(real_attribs);

        }
    }

    Redistribute(); 

}


void c_PointChargeContainer::Print_Container(bool print_positions) 
{
    amrex::Print() << "flag_vary_occupation: " << flag_vary_occupation << "\n";
    amrex::Print() << "V0: " << V0 << "\n";
    amrex::Print() << "Et: " << Et << "\n";
    amrex::Print() << "Point Charges: \n";
    int lev = 0;

    for (c_PointChargeContainer::ParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        auto np = pti.numParticles();
        amrex::Print() << "number of charges: " << np << "\n";
        const auto& particles = pti.GetArrayOfStructs();
        const auto* p_par = particles().data();
        const auto& soa_real = pti.get_realPA();
        
        for(int p=0; p<np; ++p) 
        {    
            amrex::Real charge_unit = soa_real[PCA::realPA::charge_unit][p];
            amrex::Real occupation  = soa_real[PCA::realPA::occupation][p];
            amrex::Real potential   = soa_real[PCA::realPA::potential][p];

            amrex::Print() << "ID: " << p
                           << ", charge_unit: " << charge_unit 
                           << ", occupation: " << occupation
                           << ", potential: " << potential;
            if (print_positions)
            {
                amrex::Print() << ", Position: (";
                for (int d = 0; d < AMREX_SPACEDIM; ++d)
                {
                    amrex::Print() << p_par[p].pos(d);
                    if (d < AMREX_SPACEDIM - 1)
                    {
                        amrex::Print() << ", ";
                    }
                }
                amrex::Print() << ")";
            }
            amrex::Print() << "\n";
        }
    }
}


void c_PointChargeContainer::Compute_Occupation()
{
    auto& rCode  = c_Code::GetInstance();
    auto& rGprop = rCode.get_GeometryProperties();

    const auto& plo = rGprop.geom.ProbLoArray();
    const auto& dx = rGprop.geom.CellSizeArray();

    int lev = 0;
    for (c_PointChargeContainer::ParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        auto np = pti.numParticles();
        const auto& particles = pti.GetArrayOfStructs();
        const auto* p_par = particles().data();

        const auto& soa_real = pti.get_realPA();

        auto get_V0 = Get_V0();
        auto get_Et = Get_Et();

        auto& par_potential = pti.get_potential();
        auto* p_potential = par_potential.data();

        auto& par_occupation = pti.get_occupation();
        auto* p_occupation = par_occupation.data();

        amrex::Real unit_charge = PhysConst::q_e;

        amrex::ParallelFor (np, [=] AMREX_GPU_DEVICE (int p)
        {
            p_occupation[p] = 
            MathFunctions::Sigmoid ( - (p_potential[p] - get_V0()) / get_Et() );
        });
    }

}
