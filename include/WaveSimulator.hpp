#pragma once

#include "data_structures.hpp"
#include "routines.hpp"
#include <fstream>
#include <vector>

template <class ExecSpace> 
class WaveSimulator {
  private:
    using MemSpace = typename ExecSpace::accessible_space;
    Fields<MemSpace> _fields_old;
    Fields<MemSpace> _fields_new;

    size_t _nt;
    size_t _nz;
    size_t _nx;
    float_type _dt;
    float_type _dz;
    float_type _dx;
    float_type _srcz;
    float_type _srcx;
    std::vector<float_type> _wavelet;

  public:
    WaveSimulator() = default;
    ~WaveSimulator() = default;

    // Set-methods
    void set_time_step(float_type dt)    { _dt = dt; }
    void set_spatial_step(float_type dx) { _dx = dx; }
    void set_depth_step(float_type dz)   { _dz = dz; }
    void set_number_of_time_steps(size_t nt)  { _nt = nt; }
    void set_dimensions(size_t nz, size_t nx) { _nz = nz; _nx = nx; }
    void set_source_position_x(float_type x)  { _srcx = x; }
    void set_source_position_z(float_type z)  { _srcz = z; }

    // Get-methods
    float_type get_time_step() const    { return _dt; }
    float_type get_spatial_step() const { return _dx; }
    float_type get_depth_step() const   { return _dz; }
    size_t get_number_of_time_steps() const    { return _nt; }
    size_t get_number_of_spatial_steps() const { return _nx; }
    size_t get_number_of_depth_steps() const   { return _nz; }
    size_t get_source_position_x() const { return _srcx; }
    size_t get_source_position_z() const { return _srcz; }
    Fields<MemSpace> const* get_field() const { return &_fields_new; }

    // Other public methods
    void make_ricker(float_type fpeak);
    void check_CFL_condition( Models<MemSpace> const& models ) const;
    bool wrong_initialization_of_simulator_encountered() { return false; } // NOT IMPLEMENTED YET; SO JUST LET PROGRAM PROCEED UNCHECKED!
    void allocate_internal_data_structures();
    void clean_internal_data_structures();

    // main algorithm
    int run( Models<MemSpace> const& models ){
        if ( wrong_initialization_of_simulator_encountered() ) return 1;
        for (size_t i(0); i < _nt; ++i) {
            if (i % 250 == 0) std::cout << "time-step: " << i << std::endl;
            if ( inject_source( _fields_old, _wavelet[i], (int)(_srcz/_dz), (int)(_srcx/_dx), ExecSpace() ) ) return 1;
            if ( velocity_update( _fields_new, _fields_old, models, _dt, _dz, _dx, ExecSpace() ) ) return 1;
            if ( stresses_update( _fields_new, _fields_old, models, _dt, _dz, _dx, ExecSpace() ) ) return 1;
            std::swap(_fields_new,_fields_old);
        }
        return 0;
    }
};

template <class ExecSpace> void WaveSimulator<ExecSpace>::check_CFL_condition( Models<MemSpace> const& models ) const {
    float_type *vs;
    ppt::MemSpaceHost::allocate(&vs, _nz * _nx);
    MemSpace::copyToHost( vs, models.Vs->get_ptr(), _nz * _nx );
    float_type vmin = 1e10;
    for ( int i=0; i < _nz * _nx; ++i ) 
        if ( vs[i] < vmin ) 
            vmin = vs[i];
    ppt::MemSpaceHost::release(vs);
    float_type cfl = ( (_dt/_dx)*std::sqrt(2)*vmin );
    std::cout << "CFL: " << cfl << std::endl;
    if ( cfl > 1 ) throw std::runtime_error("CFL STABILITY CONDITION NOT SATISFIED!");
}

/**
 * @brief Create a Ricker-wavelet:
 * R(t) = (1 - 2 pi^2 fpeak^2 t^2) * exp( -1 * pi^2 fpeak^2 t^2 )
 */
template <class ExecSpace> void WaveSimulator<ExecSpace>::make_ricker(float_type fpeak)
{
    _wavelet.resize(_nt);
    // compute Ricker-wavelet with user defined peak frequency
    _wavelet[0] = 1.0;
    for (int it = 1; it < _nt; ++it)
    {
        float_type t         = it * _dt;
        float_type term      = M_PI * M_PI * fpeak * fpeak * t * t;
        _wavelet[it]         = (1 - 2 * term) * exp(-term);
    }
}

template <class ExecSpace> void WaveSimulator<ExecSpace>::allocate_internal_data_structures()
{
    _fields_new.V     = new ScalarField<MemSpace>( _nz, _nx );
    _fields_new.S     = new ScalarField<MemSpace>( _nz, _nx );
    _fields_new.T     = new ScalarField<MemSpace>( _nz, _nx );

    _fields_old.V     = new ScalarField<MemSpace>( _nz, _nx );
    _fields_old.S     = new ScalarField<MemSpace>( _nz, _nx );
    _fields_old.T     = new ScalarField<MemSpace>( _nz, _nx );
}

template <class ExecSpace> void WaveSimulator<ExecSpace>::clean_internal_data_structures()
{
    delete _fields_new.V; _fields_new.V = nullptr;
    delete _fields_new.S; _fields_new.S = nullptr;
    delete _fields_new.T; _fields_new.T = nullptr;

    delete _fields_old.V; _fields_old.V = nullptr;
    delete _fields_old.S; _fields_old.S = nullptr;
    delete _fields_old.T; _fields_old.T = nullptr;
}


template<typename MemSpace>
void store_to_binary(const char *filename, ScalarField<MemSpace> const* const field) {
    std::fstream file;
    file.open(filename, std::ios_base::out | std::ios_base::binary);

    if (!file.is_open()) throw std::runtime_error("Unable to open file!");

    // allocate array on host, copy the data into, and then store to binary file
    float_type *data_host;
    ppt::MemSpaceHost::allocate(&data_host, field->get_nElems());
    MemSpace::copyToHost(data_host, field->get_ptr(), field->get_nElems());

    // write data to binary file
    file.write((char *)data_host, sizeof(float_type) * field->get_nElems());

    // deallocate the host-array and close file
    ppt::MemSpaceHost::release(data_host);
    file.close();
}