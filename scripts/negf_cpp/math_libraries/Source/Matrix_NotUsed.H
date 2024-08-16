#include <AMReX.H>
#include <AMReX_GpuComplex.H>
#include <AMReX_REAL.H>
#include<AMReX_Print.H>
#include<AMReX_TableData.H>
#include <AMReX_GpuUtility.H>

using namespace amrex;
using MatrixDType = amrex::GpuComplex<amrex::Real>;
using Matrix1D = TableData<MatrixDType, 1>;
using Matrix2D = TableData<MatrixDType, 2>;

//template<typename T>
//class TD;

enum class MemSpace : int {Arena, Host, Device, Managed};

class Matrix {

public:    
    Matrix () {}
    Matrix (int height, int width, MemSpace mem_space=MemSpace::Arena) 
    {
        resize(height, width, mem_space);
    }

    void resize(int height, int width, MemSpace mem_space)
    {
        _mem_space = mem_space;
        #ifdef AMREX_USE_GPU
        if(_mem_space == MemSpace::Arena) {
            d_tab.resize({0, 0}, {height, width}, The_Arena());
        } 
        else if(_mem_space == MemSpace::Host) {
            h_tab.resize({0, 0}, {height, width}, The_Pinned_Arena());
        } 
        else if(_mem_space == MemSpace::Device) {
            d_tab.resize({0, 0}, {height, width}, The_Device_Arena());
        }
        else if(_mem_space == MemSpace::Managed) {
            d_tab.resize({0, 0}, {height, width}, The_Managed_Arena());
        } 
        #else
        h_tab.resize({0, 0}, {height, width});
        #endif
    }        

    const auto& hi() const {
        #ifdef AMREX_USE_GPU
        return d_tab.hi();
        #else
        return h_tab.hi();
        #endif
    }

    MemSpace mem_space() const {
        return _mem_space;
    }

    const auto& host_table() const {
        return h_tab.table();
    }

    const auto& host_const_table() const {
        return h_tab.const_table();
    }

    #ifdef AMREX_USE_GPU
    void copy_to_host() {
        const auto& hi = d_tab.hi();
        h_tab.resize({0, 0}, {hi[0], hi[1]}, The_Pinned_Arena());

        d_tab.copy(h_tab);
        Gpu::streamSynchronize();  
    }
    void copy_to_device() {
        const auto& hi = h_tab.hi();
        d_tab.resize({0, 0}, {hi[0], hi[1]}, The_Device_Arena());

        h_tab.copy(d_tab);
        Gpu::streamSynchronize();  
    }
    const auto& device_const_table() const {
        return d_tab.const_table();
    }
    const auto& device_table() const {
        return d_tab.table();
    }
    void clear_device() {
        d_tab.clear();
    }
    #endif

    void clear_host() {
        h_tab.clear();
    }
    void clear() {
        #ifdef AMREX_USE_GPU
        h_tab.clear();
        d_tab.clear();
        #else
        h_tab.clear();
        #endif
    }
private:

    MemSpace _mem_space;
    #ifdef AMREX_USE_GPU
    Matrix2D h_tab;
    Matrix2D d_tab;
    #else
    Matrix2D h_tab;
    #endif
};
