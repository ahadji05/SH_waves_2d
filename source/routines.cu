
#include "routines.hpp"
#include "ppt/execution/ExecutionSpacesInc.hpp"

using namespace ppt;

#define BLOCKDIM_X 64
#define BLOCKDIM_Z 1

__global__ void inject_source_kernel( 
    float_type *const p, 
    float_type amplitude,
    int iz, 
    int ix, 
    int nz, 
    int nx )
{
    int idx = blockDim.x * blockIdx.x + threadIdx.x;
    if ( idx > 0 ) return;
    if ( iz >= nz || ix >= nx ) return;

    p[ iz * nx + ix ] = amplitude;
}

template<>
int inject_source( Fields<MemSpaceCuda> &fields, float_type amplitude, int iz, int ix, ExecutionSpaceCuda ){
    
    int const nz = fields.V->get_nz();
    int const nx = fields.V->get_nx();

    float_type *const p = fields.V->get_ptr();
    
    inject_source_kernel<<<1,1>>>( 
        p, amplitude, iz, ix, nz, nx 
    );

    return 0;
}

__global__ void velocity_update_kernel(
    float_type * const Vnew, 
    float_type const* const Vold, 
    float_type const* const Sold, 
    float_type const* const Told, 
    float_type const* const L, 
    int nz, 
    int nx,
    float_type dt_div_dz, 
    float_type dt_div_dx )
{
    int ix = blockDim.x * blockIdx.x + threadIdx.x;
    int iz = blockDim.y * blockIdx.y + threadIdx.y;

    if ( ix >= nx - 1 ) return;
    if ( iz >= nz - 1 ) return;

    int idx = iz * nx + ix;
    Vnew[idx] = Vold[idx] + dt_div_dz * L[idx] * (Sold[(iz+1)*nx + ix] - Sold[idx]) + dt_div_dx * L[idx] * (Told[idx+1] - Told[idx]);
}

template<>
int velocity_update( Fields<MemSpaceCuda> &fields_new, Fields<MemSpaceCuda> const& fields_old, Models<MemSpaceCuda> const& model,
    float_type dt, float_type dz, float_type dx, ExecutionSpaceCuda ){
    int const nz = fields_new.V->get_nz();
    int const nx = fields_new.V->get_nx();

    float_type const* const Vold = fields_old.V->get_ptr();
    float_type const* const Sold = fields_old.S->get_ptr();
    float_type const* const Told = fields_old.T->get_ptr();
    float_type const* const L    = model.L->get_ptr();

    float_type const dt_div_dz = dt/dz;
    float_type const dt_div_dx = dt/dx;
    
    float_type * const Vnew = fields_new.V->get_ptr();

    dim3 nThreads(BLOCKDIM_X, BLOCKDIM_Z, 1);
    size_t nBlock_x = nx % BLOCKDIM_X == 0 ? size_t(nx / BLOCKDIM_X) : size_t(1 + nx / BLOCKDIM_X);
    size_t nBlock_z = nz % BLOCKDIM_Z == 0 ? size_t(nz / BLOCKDIM_Z) : size_t(1 + nz / BLOCKDIM_Z);
    dim3 nBlocks(nBlock_x, nBlock_z, 1);

    velocity_update_kernel<<<nBlocks,nThreads>>>(
        Vnew, Vold, Sold, Told, L, nz, nx, dt_div_dz, dt_div_dx 
    );

    return 0;
}


__global__ void stresses_update_kernel(
    float_type * const Snew, 
    float_type * const Tnew, 
    float_type const* const Vnew, 
    float_type const* const Sold, 
    float_type const* const Told, 
    float_type const* const M, 
    int nz, 
    int nx, 
    float_type dt_div_dz, 
    float_type dt_div_dx )
{
    int ix = blockDim.x * blockIdx.x + threadIdx.x;
    int iz = blockDim.y * blockIdx.y + threadIdx.y;

    if ( ix == 0 || ix >= nx ) return;
    if ( iz == 0 || iz >= nz ) return;
    
    int idx = iz * nx + ix;
    Snew[idx] = Sold[idx] + dt_div_dz * M[idx] * (Vnew[idx] - Vnew[(iz-1)*nx+ix]);
    Tnew[idx] = Told[idx] + dt_div_dx * M[idx] * (Vnew[idx] - Vnew[idx-1]);
}

template<>
int stresses_update( Fields<MemSpaceCuda> &fields_new, Fields<MemSpaceCuda> const& fields_old, Models<MemSpaceCuda> const& model,
    float_type dt, float_type dz, float_type dx, ExecutionSpaceCuda ){
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

    dim3 nThreads(BLOCKDIM_X, BLOCKDIM_Z, 1);
    size_t nBlock_x = nx % BLOCKDIM_X == 0 ? size_t(nx / BLOCKDIM_X) : size_t(1 + nx / BLOCKDIM_X);
    size_t nBlock_z = nz % BLOCKDIM_Z == 0 ? size_t(nz / BLOCKDIM_Z) : size_t(1 + nz / BLOCKDIM_Z);
    dim3 nBlocks(nBlock_x, nBlock_z, 1);

    stresses_update_kernel<<<nBlocks,nThreads>>>(
        Snew, Tnew, Vnew, Sold, Told, M, nz, nx, dt_div_dz, dt_div_dx
    );

    return 0;
}
