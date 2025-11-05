#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct BiasedQ : GradientBase<GradientBiased, AllDirections> {
    template <class TTraits, class TParameter>
    inline double compute(const int direction, const int k);
};

template <class TTraits, class TParameter>
inline double BiasedQ::compute(const int direction, const int k) {
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    /*return 0.5 *
           (TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction]) -
            TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]));*/

    return 0.5 *
           (-TParameter::template get<Lattice>(
                data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + direction] * Stencil::Q + direction]) +
            4 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction]) -
            3 * TParameter::template get<Lattice>(k));
}
