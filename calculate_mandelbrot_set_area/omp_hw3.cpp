#include "MandelbrotSet.hpp"

#include "omp.h"

#include <cmath>
#include <complex>
#include <cstddef>
#include <iostream>
#include <iomanip>
// #include <mutex>
// #include <queue>
// #include <shared_mutex>
// #include <type_traits>
// #include <tuple>
// #include <unordered_map>
// #include <vector>

namespace
{
    using std::size_t;

    // using RealNumber = long double;
    using RealNumber = double;
    using ComplexNumber = std::complex<RealNumber>;
}

namespace
{
    using namespace Fpdc2019;
    using MS = MandelbrotSet<ComplexNumber>;

    template <typename RealNumber, typename Duration>
    void outputResults(RealNumber const result, Duration const duration)
    {
        using std::cout;
        using std::setprecision;

        size_t numOfThreads {static_cast<size_t>(omp_get_max_threads())};
        RealNumber const error {std::abs(MS::MandelbrotSet::k_area - result)};

        cout << "Number of threads = ";
        cout << numOfThreads << "\n";
        cout << "The real area is approximately ";
        cout << setprecision(16) << MS::MandelbrotSet::k_area << "\n";
        cout << "The calculated area is ";
        cout << setprecision(16) << result << "\n";
        cout << "The error is ";
        cout << std::scientific << setprecision(3) << error << "\n";
        cout << "Wall clock time = ";
        cout << std::fixed << setprecision(8) << duration << "\n";
    }

    template
    <
        bool usingUpperHalfOnly=true,
        std::size_t numOfRows=k_defaultNumOfRows,
        std::size_t numOfColumns=k_defaultNumOfColumns
    >
    // constexpr RealNumber computeMandelbrotSetArea()
    RealNumber computeMandelbrotSetArea
    (
        std::size_t const maxIterationCount=
        IterationStopCriterion<ComplexNumber>::k_defaultMaxIterationCount
    )
    {
        MandelbrotSetGrid
        <
            ComplexNumber,
            usingUpperHalfOnly,
            numOfRows,
            numOfColumns
        > mandelbrotSetGrid {};

        IterationStopCriterion<ComplexNumber> stopCriterion {maxIterationCount};

        mandelbrotSetGrid.iterateAllUntil(stopCriterion);

        return mandelbrotSetGrid.getApproximatedArea();
    }
}

int main(int argc, char* argv[])
{
    static_cast<void>(argc);
    static_cast<void>(argv);
    std::ios::sync_with_stdio(false);

    auto const startTime {omp_get_wtime()};
    // constexpr RealNumber computedMandelbrotSetArea {computeMandelbrotSetArea()};
    const RealNumber computedMandelbrotSetArea
    {
        // computeMandelbrotSetArea
        // <
        //     true,
        //     5300,
        //     1000
        // >(10000)
        computeMandelbrotSetArea
        <
            true,
            900,
            100
        >(950)
    };
    auto const duration {omp_get_wtime() - startTime};

    outputResults(computedMandelbrotSetArea, duration);

    return 0;
}