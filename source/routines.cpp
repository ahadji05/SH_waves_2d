
#include "routines.hpp"
#include "ppt/execution/ExecutionSpacesInc.hpp"

using namespace ppt;

template<>
int inject_source( Fields<MemSpaceHost> &fields, float_type amplitude, int iz, int ix, ExecutionSpaceSerial ){
    int const nz = fields.V->get_nz();
    int const nx = fields.V->get_nx();
    if ( iz >= nz || ix >= nx ) return 1;
    fields.V->get_ptr()[ iz * nx + ix ] = amplitude;
    return 0;
}

template<>
int velocity_update( Fields<MemSpaceHost> &fields_new, Fields<MemSpaceHost> const& fields_old, Models<MemSpaceHost> const& model,
    float_type dt, float_type dz, float_type dx, ExecutionSpaceSerial ){
    
    int const nz = fields_new.V->get_nz();
    int const nx = fields_new.V->get_nx();

    float_type const* const Vold = fields_old.V->get_ptr();
    float_type const* const Sold = fields_old.S->get_ptr();
    float_type const* const Told = fields_old.T->get_ptr();
    float_type const* const L    = model.L->get_ptr();

    float_type const dt_div_dz = dt/dz;
    float_type const dt_div_dx = dt/dx;

    float_type * const Vnew = fields_new.V->get_ptr();

    for (int iz=0; iz < nz - 1; ++iz)
        for (int ix=0; ix < nx - 1; ++ix){
            int idx = iz * nx + ix;
            Vnew[idx] = Vold[idx] + dt_div_dz * L[idx] * (Sold[(iz+1)*nx + ix] - Sold[idx]) + dt_div_dx * L[idx] * (Told[idx+1] - Told[idx]);
        }

    return 0;
}

template<>
int stresses_update( Fields<MemSpaceHost> &fields_new, Fields<MemSpaceHost> const& fields_old, Models<MemSpaceHost> const& model, 
    float_type dt, float_type dz, float_type dx, ExecutionSpaceSerial ){

    int const nz = fields_new.V->get_nz();
    int const nx = fields_new.V->get_nx();

    float_type const* const Vnew = fields_new.V->get_ptr();
    float_type const* const Sold = fields_old.S->get_ptr();
    float_type const* const Told = fields_old.T->get_ptr();

    float_type const* const M    = model.M->get_ptr();

    float_type const dt_div_dz = dt/dz;
    float_type const dt_div_dx = dt/dx;

    float_type * const Snew = fields_new.S->get_ptr();
    float_type * const Tnew = fields_new.T->get_ptr();

    for (int iz=1; iz < nz; ++iz)
        for (int ix=1; ix < nx; ++ix){
            int idx = iz * nx + ix;
            Snew[idx] = Sold[idx] + dt_div_dz * M[idx] * (Vnew[idx] - Vnew[(iz-1)*nx+ix]);
            Tnew[idx] = Told[idx] + dt_div_dx * M[idx] * (Vnew[idx] - Vnew[idx-1]);
        }

    return 0;
}
