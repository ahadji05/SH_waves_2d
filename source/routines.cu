
#include "routines.hpp"
#include "ppt/execution/ExecutionSpacesInc.hpp"

using namespace ppt;

template<>
int velocity_update( Fields<MemSpaceCuda> &fields_new, Fields<MemSpaceCuda> const& fields_old, Models<MemSpaceCuda> const& model, ExecutionSpaceCuda ){
    throw std::runtime_error("NOT IMPLEMENTED YET");
}

template<>
int stresses_update( Fields<MemSpaceCuda> &fields_new, Fields<MemSpaceCuda> const& fields_old, Models<MemSpaceCuda> const& model, ExecutionSpaceCuda ){
    throw std::runtime_error("NOT IMPLEMENTED YET");
}
