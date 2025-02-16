
#include "ppt/types.hpp"
#include "WaveSimulator.hpp"

using namespace ppt;

// valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes ./main.exe
// compute-sanitizer --tool memcheck --leak-check=full ./main.exe

#if !defined(PPT_ENABLE_OPENMP_BACKEND) && !defined(PPT_ENABLE_CUDA_BACKEND)
bool backend_specified = false;
static_assert( backend_specified, "NO BACKEND HAS BEEN SPECIFIED!!!" );
#endif

void create_model( Models<memo_space> &model, int nz, int nx );
void delete_model( Models<memo_space> &model );

int main( int argc, char *argv[] ){

    // Modelling parameters
    int nz = 400;
    int nx = 1000;
    int nt = 1501;
    float_type dz = 1.5;
    float_type dx = 1.5;
    float_type dt = 0.0001;
    float_type source_position_z = 400;
    float_type source_position_x = 750;
    float_type ricker_fpeak = 30;

    // Create model
    Models<memo_space> model;
    create_model( model, nz, nx );

    // Create Simulator
    WaveSimulator<exec_space> Sim;
    Sim.set_number_of_time_steps(nt);
    Sim.set_dimensions(nz, nx);
    Sim.set_time_step(dt);
    Sim.set_spatial_step(dx);
    Sim.set_depth_step(dz);
    Sim.set_source_position_z(source_position_z);
    Sim.set_source_position_x(source_position_x);
    Sim.make_ricker(ricker_fpeak);

    Sim.allocate_internal_data_structures();

    try { Sim.check_CFL_condition( model ); }
    catch(const std::exception& e){ std::cerr << e.what() << '\n'; return EXIT_FAILURE; }

    Sim.run( model );

    store_to_binary( "V.bin", Sim.get_field()->V );
    store_to_binary( "S.bin", Sim.get_field()->S );
    store_to_binary( "T.bin", Sim.get_field()->T );
    store_to_binary( "Vs.bin",  model.Vs );
    store_to_binary( "L.bin",  model.L );
    store_to_binary( "M.bin",  model.M );

    Sim.clean_internal_data_structures();

    delete_model( model );

    return 0;
}

void create_model( Models<memo_space> &model, int nz, int nx ){

    // Create model on host 
    float_type *vs, *m, *l;
    MemSpaceHost::allocate( &vs, nz*nx );
    MemSpaceHost::allocate( &l, nz*nx );
    MemSpaceHost::allocate( &m, nz*nx );
    for ( int iz=0; iz < nz; ++iz )
        for ( int ix=0; ix < nx; ++ix ){
            int idx = iz*nx+ix;
            float_type Vs  = 2000; // S-wave velocity
            float_type rho = 2400; // Density
            if ( iz == (int)nz/2 )   rho /= 2; // add a reflector in the middle of the model
            if ( iz == (int)2*nz/3 ) rho *= 2; // add another one a few depth levels deeper
            vs[idx] = Vs;
            l[idx]  = 1.0 / rho;
            m[idx]  = Vs * Vs * rho;
        }
    
    // Copy model from host to the memory space that is compatible with the execution space
    model.Vs = new ScalarField<memo_space>( nz, nx );
    model.M  = new ScalarField<memo_space>( nz, nx );
    model.L  = new ScalarField<memo_space>( nz, nx );
    memo_space::copyFromHost( model.Vs->get_ptr(), vs, nz*nx );
    memo_space::copyFromHost( model.L->get_ptr(),  l, nz*nx );
    memo_space::copyFromHost( model.M->get_ptr(),  m, nz*nx );

    MemSpaceHost::release( vs );
    MemSpaceHost::release( l );
    MemSpaceHost::release( m );
}

void delete_model( Models<memo_space> &model ){
    delete model.Vs;
    delete model.L;
    delete model.M;
}
