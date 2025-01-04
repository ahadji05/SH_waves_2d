#pragma once

#include "ppt/containers/ScalarField.hpp"

template<typename MemSpace>
struct Models {
    ScalarField<MemSpace> *Vs = nullptr; // Shear wave velocity
    ScalarField<MemSpace> *L = nullptr;  // lightness ( inverse density )
    ScalarField<MemSpace> *M = nullptr;  // shear modulus = Vs * Vs * density
};

template<typename MemSpace>
struct Fields {
    ScalarField<MemSpace> *V = nullptr;
    ScalarField<MemSpace> *S = nullptr;
    ScalarField<MemSpace> *T = nullptr;
};

