#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct BiasedXYZ : GradientBase<GradientBiased, Cartesian> {
    template <class TTraits, class TParameter>
    inline double compute(const int direction, const int k);
};

template <class TTraits, class TParameter>
inline double BiasedXYZ::compute(const int direction, const int k) {
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    /*return 0.5 *
           (TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction]) -
            TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]));*/

    /*return 0.25 *
           (-TParameter::template get<Lattice>(
                data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + direction] * Stencil::Q + direction]) +
            4 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction]) -
            3 * TParameter::template get<Lattice>(k));*/

    double gradientsum = 0;

    for (int idx = 1; idx < Stencil::Q; idx++) {
        gradientsum +=
            Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] *
            (-TParameter::template get<Lattice>(
                 data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + idx] * Stencil::Q + idx]) +
             4 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + idx]) -
             3 * TParameter::template get<Lattice>(k)
             );
    }

    return 0.25 / (Stencil::Cs2 * Lattice::DT) * gradientsum;
}
