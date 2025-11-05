#pragma once
#include <cmath>

#include "../Service.hh"
#include "GradientBase.hh"

struct LaplacianCentralWetting : GradientBase<Laplacian, One> {
    template <class TTraits, class TParameter>
    inline double compute(int direction, int k);
};

template <class TTraits, class TParameter>
inline double LaplacianCentralWetting::compute(int direction, int k) {
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;
    constexpr int N = TTraits::NumberOfComponents;
    using DataType = Data_Base<Lattice, Stencil>;
    
    DataType& data = DataType::getInstance();
    
    double laplaciansum = 0;
    const static auto& param = TParameter::template get<Lattice>();
    
    // int boundaryID = Geometry<Lattice>::getBoundaryType(k);
    // double prefactor = getPrefactorForBoundary(boundaryID);
    
    // for (int idx = 1; idx < Stencil::Q; idx++) {
    //     int neighbor = data.getNeighbor(k, idx);
    //     if (!this->isBoundary<Lattice>(neighbor)) {
    //         laplaciansum += Stencil::Weights[idx] * 2 * (param[neighbor] - param[k]);
    //     } else {
    //         int neighborIndex = data.getNeighbor(k, idx);
    //         // print the neighbor index and boundary ID
    //         int boundaryID = Geometry<Lattice>::getBoundaryType(neighborIndex);
    //         double prefactor = getPrefactorForBoundary(boundaryID);
    //         // std::cout << "Neighbor Index: " << neighborIndex << "prefactor" << prefactor << std::endl;
    //         // we should print the boundary ID and prefactor
    //         // if (boundaryID != 0) {
    //         //         std::cout << "Boundary ID: " << boundaryID << ", Prefactor: " << prefactor << std::endl;
    //         //     }
    //         const int& normalq =
    //             TTraits::Stencil::QMap
    //                 .find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(neighbor).NormalDirection)
    //                 ->second;
    //         double factor = 0.5;

    //         if constexpr (N <= 2) {
    //             double csolid = param[data.getNeighbor(neighbor, normalq)];
    //             laplaciansum += Stencil::Weights[idx] * 2 *
    //                 ((csolid - factor * prefactor * (csolid - pow(csolid, 2))) - param[k]);
    //         } else {
    //             double csolidsum = 0;
    //             double csum = 0;
    //             auto& opcurrent = OrderParameter<TParameter::instance>::template get<Lattice>();
    //             const double& csolidcurrent = opcurrent[data.getNeighbor(neighbor, normalq)];
    //             for (int component = 0; component < N - 1; component++) {
    //                 auto& op = getInstance<OrderParameter, N - 1, Lattice>(component);
    //                 const double& csolid = op[data.getNeighbor(neighbor, normalq)];
    //                 const double& c = op[k];
    //                 csolidsum += csolid;
    //                 csum += c;
    //                 laplaciansum += Stencil::Weights[idx] * 2 *
    //                     ((-factor * prefactor * (csolid * csolidcurrent)));
    //             }
    //             double csolidN = 1.0 - csolidsum;
    //             laplaciansum += Stencil::Weights[idx] * 2 *
    //                 ((csolidcurrent - factor * prefactor * (csolidN * csolidcurrent)) - opcurrent[k]);
    //         }
    //     }
    // }
        for (int idx = 1; idx < Stencil::Q; ++idx) {
        int neighbor = data.getNeighbor(k, idx);

        if (!this->isBoundary<Lattice>(neighbor)) {
            laplaciansum += Stencil::Weights[idx] * 2 * (param[neighbor] - param[k]);
        } else {
            // Neighbor is solid: fluid node k "feels" wetting
            
            int boundaryID = Geometry<Lattice>::getBoundaryType(neighbor);
            double prefactor = getPrefactorForBoundary(boundaryID);

            // Sample a node inside the solid along normal
            // const int normalq = TTraits::Stencil::QMap.find(
            //     BoundaryLabels<TTraits::Lattice::NDIM>::template get<Lattice>(neighbor).NormalDirection
            // )->second;

            const int& normalq =
                TTraits::Stencil::QMap
                    .find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(
                              data.getNeighbor(k, idx))
                              .NormalDirection)
                    ->second;
            double csolid = param[data.getNeighbor(neighbor, normalq)];

            // Update Laplacian at fluid node
            laplaciansum += Stencil::Weights[idx] * 2 *
                            ((csolid - 0.5 * prefactor * (csolid - csolid*csolid)) - param[k]);
        }
    }

    return 1.0 / (Stencil::Cs2 * Lattice::DT * Lattice::DT) * laplaciansum;
}