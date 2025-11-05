#pragma once
#include "../Geometry.hh"
#include "../Template.hh"
#include <execinfo.h>
#include <iostream>
#include <cstdlib>
std::unordered_map<int, double> mBoundaryPrefactorMap;

void print_stacktrace() {
    void *array[10];
    size_t size = backtrace(array, 10);
    char **strings = backtrace_symbols(array, size);

    std::cerr << "Stack trace:\n";
    for (size_t i = 0; i < size; ++i)
        std::cerr << strings[i] << std::endl;

    free(strings);
}


template <template <class> class TGradientType, class TDirections = Cartesian>
struct GradientBase {
    std::vector<double> mvPrefactor = {0.0};
    const double& mPrefactor = mvPrefactor[0];
    std::unordered_map<int, double> mBoundaryPrefactorMap;
    GradientBase& operator=(const GradientBase& other) {
        mBoundaryID = other.mBoundaryID;
        preset_warning = other.preset_warning;
        mvPrefactor = other.mvPrefactor;
        return *this;
    }

    template <class TTraits, class TParameter>
    inline double compute(const int direction, const int k);

    template <class TObj>
    using GradientType = TGradientType<TObj>;

    template <class TStencil>
    inline static constexpr int getNumberOfDirections() {
        if constexpr (std::is_same_v<TDirections, Cartesian>)
            return TStencil::D;
        else if constexpr (std::is_same_v<TDirections, AllDirections>)
            return TStencil::Q;
        else if constexpr (std::is_same_v<TDirections, One>)
            return 1;
        else
            return TStencil::D;
    }

    template <class TLattice>
    inline bool isBoundary(int k) {
        for (int i : mBoundaryID) {
            // TMP: Default BoundaryID warning
            if (Geometry<TLattice>::getBoundaryType(k) == i) {
                if (preset_warning) {
#pragma omp critical
                    //print the function call stack
                    print_stacktrace();
                   std::cout << "Warning: Using default BoundaryID " << i << " for gradient calculation. "
                             << "Consider setting a specific BoundaryID using setBoundaryID()." << std::endl;
                    preset_warning = false;  // Only warn once
                }
                return true;
            }
        }
        return false;
    }

    // inline void setBoundaryID(int id, bool preset = false) {
    //     mBoundaryID = {id};
    //     preset_warning = preset;
    // };
    // inline void setBoundaryID(const std::vector<int>& id, bool preset = false) {
    //     mBoundaryID = id;
    //     preset_warning = preset;
    // };
    inline void setBoundaryID(int id, bool preset = false) {
    mBoundaryID = {id};
    preset_warning = preset;
    std::cout << "mBoundaryID set to: ";
    for (auto v : mBoundaryID) std::cout << v << " ";
    std::cout << std::endl;
    }
    inline void setBoundaryID(const std::vector<int>& id, bool preset = false) {
        mBoundaryID = id;
        preset_warning = preset;
        std::cout << "mBoundaryID set to: ";
        for (auto v : mBoundaryID) std::cout << v << " ";
        std::cout << std::endl;
    }


    inline std::vector<int> getBoundaryIDs() const { return mBoundaryID; }
    std::vector<int> mBoundaryID = {1};
    bool preset_warning = true;

    inline void setPrefactor(double prefactor) { mvPrefactor[0] = prefactor; }
    inline void setPrefactor(const std::vector<double>& prefactor) {
        if (!prefactor.empty()) mvPrefactor[0] = prefactor[0];
    }
    inline void setPrefactor(const std::vector<int>& boundaryIDs, const std::vector<double>& prefactors) {
        mBoundaryPrefactorMap.clear();
        for (size_t i = 0; i < boundaryIDs.size() && i < prefactors.size(); ++i)
            mBoundaryPrefactorMap[boundaryIDs[i]] = prefactors[i];
    }
    inline double getPrefactorForBoundary(int boundaryID) const {
    auto it = mBoundaryPrefactorMap.find(boundaryID);
    if (it != mBoundaryPrefactorMap.end())
        return it->second;
    return mPrefactor; // fallback to default
}
    inline void setInterfaceDistance(double (*distance)(int k, int idx)) {}

    inline void setInterfaceVal(double value) {}
    inline void setPrefactorForBoundary(int boundaryID, double value) {
    mBoundaryPrefactorMap[boundaryID] = value; // insert or overwrite
}

};
