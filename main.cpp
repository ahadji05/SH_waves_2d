
#include "ppt/types.hpp"
#include "WaveSimulator.hpp"

using namespace ppt;

// valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes ./main.exe

int main(){

    // Modelling parameters
    int nz = 400;
    int nx = 1000;
    int nt = 2251;
    float_type dz = 2.5;
    float_type dx = 2.5;
    float_type dt = 0.0002;
    float_type source_position_z = 50;
    float_type source_position_x = 1250;
    float_type ricker_fpeak = 20;

    // Create model
    Models<MemSpaceHost> model;
    model.Vs = new ScalarField<MemSpaceHost>( nz, nx );
    model.M  = new ScalarField<MemSpaceHost>( nz, nx );
    model.L  = new ScalarField<MemSpaceHost>( nz, nx );
    for ( int iz=0; iz < nz; ++iz )
        for ( int ix=0; ix < nx; ++ix ){
            int idx = iz*nx+ix;
            float_type Vs  = 2000; // S-wave velocity
            float_type rho = 2400; // Density
            if ( iz == (int)nz/2 )   rho /= 2; // add a reflector in the middle of the model
            if ( iz == (int)2*nz/3 ) rho *= 2; // add another one a few depth levels deeper
            model.Vs->get_ptr()[idx] = Vs;
            model.L->get_ptr()[idx]  = 1.0 / rho;
            model.M->get_ptr()[idx]  = Vs * Vs * rho;
        }

    // Create Simulator
    WaveSimulator<ExecutionSpaceSerial> Sim;
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

    delete model.Vs;
    delete model.L;
    delete model.M;

    return EXIT_SUCCESS;
}
