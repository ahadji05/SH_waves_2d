#pragma once

#include "data_structures.hpp"

template<class ExecSpace, class MemSpace>
int inject_source( Fields<MemSpace> &fields, float_type amplitude, int iz, int ix, ExecSpace tag );

template<class ExecSpace, class MemSpace>
int velocity_update( Fields<MemSpace> &fields_new, Fields<MemSpace> const& fields_old, Models<MemSpace> const& model, 
    float_type dt, float_type dz, float_type dx, ExecSpace tag );

template<class ExecSpace, class MemSpace>
int stresses_update( Fields<MemSpace> &fields_new, Fields<MemSpace> const& fields_old, Models<MemSpace> const& model, 
    float_type dt, float_type dz, float_type dx, ExecSpace tag );
