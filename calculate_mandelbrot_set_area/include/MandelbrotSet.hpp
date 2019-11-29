#ifndef FPDC_2019__MANDELBROT_SET_HPP
#define FPDC_2019__MANDELBROT_SET_HPP

#include "omp.h"

#include <array>
// #include <atomic>
#include <complex>
#include <cstddef>
// #include <functional>
#include <limits>
// #include <optional>
// #include <stdexcept>
#include <type_traits>
#include <utility>
// #include <variant>

namespace Fpdc2019
{
    template <typename ComplexNumber>
    struct MandelbrotSet
    {
        private:
            using _RealNumber = typename ComplexNumber::value_type;

        public:
            MandelbrotSet() = delete;

            static inline constexpr _RealNumber k_area {1.5065918849};
            static inline constexpr _RealNumber k_realPartLowerBound {-2.0};
            static inline constexpr _RealNumber k_realPartUpperBound {0.7};
            static inline constexpr _RealNumber k_imagPartLowerBound {-1.2};
            static inline constexpr _RealNumber
            k_imagPartUpperBound {-k_imagPartLowerBound};

            static inline constexpr _RealNumber k_initValueForIteration {0.0};

            static constexpr ComplexNumber iteratedFunction
            (ComplexNumber const var, ComplexNumber const parameter)
            {return var * var + parameter;}
    };

    // template <typename ComplexNumber, typename StopCriterion>
    // constexpr bool willDiverge
    // (ComplexNumber const parameter, StopCriterion criterion);

    template <typename ComplexNumber>
    class IterationStopCriterion
    {
        private:
            std::size_t _remainingIterationCount;
            bool _shouldStop;

        public:
            static inline constexpr std::size_t
            k_defaultMaxIterationCount {100};

            constexpr IterationStopCriterion()
            :
                _remainingIterationCount {k_defaultMaxIterationCount},
                _shouldStop {false}
            {}

            constexpr IterationStopCriterion
            (std::size_t const maxIterationCount)
            :
                _remainingIterationCount {maxIterationCount},
                _shouldStop {false}
            {}

            constexpr void update
            (ComplexNumber const newValue, ComplexNumber const previousValue);

            constexpr bool shouldStop() const {return _shouldStop;}
    };

    using Coordinates = std::pair<std::size_t, std::size_t>;

    // template <typename... Ts>
    // using GridElement = std::variant<std::monostate, Ts...>;

    template
    <
        std::size_t numOfRows,
        std::size_t numOfColumns,
        typename _GridElement
    >
    class Grid
    {
        private:
            // using _GridElement = GridElement<Ts...>;
            using _Row = std::array<_GridElement, numOfColumns>;
            using _BaseGrid = std::array<_Row, numOfRows>;

            template <typename Function>
            static constexpr bool _isInitFunctionOfTwoIndices_v
            {
                // std::disjunction_v
                // <
                //     std::is_invocable_r
                //     <std::variant<Ts...>, Function, std::size_t, std::size_t>,
                //     std::is_invocable_r
                //     <Ts, Function, std::size_t, std::size_t>...
                // >
                []()
                {
                    if constexpr
                    (std::is_invocable_v<Function, std::size_t, std::size_t>)
                    {
                        return std::is_assignable_v
                        <
                            _GridElement,
                            std::invoke_result_t
                            <Function, std::size_t, std::size_t>
                        >;
                    }
                    else
                    {
                        return false;
                    }
                }()
            };

            template <typename Function>
            static constexpr bool _isInitFunctionOfCoordinates_v
            {
                // std::disjunction_v
                // <
                //     std::is_invocable_r
                //     <std::variant<Ts...>, Function, Coordinates>,
                //     std::is_invocable_r
                //     <Ts, Function, Coordinates>...
                // >
                []()
                {
                    if constexpr (std::is_invocable_v<Function, Coordinates>)
                    {
                        return std::is_assignable_v
                        <
                            _GridElement,
                            std::invoke_result_t<Function, Coordinates>
                        >;
                    }
                    else
                    {
                        return false;
                    }
                }()
            };

            template <typename T>
            static constexpr bool _isDefaultValue_v
            {
                // (
                //     std::disjunction_v<std::is_same<T, Ts>...> ||
                //     std::is_same_v<T, std::variant<Ts...>>
                // )
                // &&
                // (
                //     !_isInitFunctionOfTwoIndices_v<T> &&
                //     !_isInitFunctionOfCoordinates_v<T>
                // )
                std::conjunction_v
                <
                    std::bool_constant<!_isInitFunctionOfTwoIndices_v<T>>,
                    std::bool_constant<!_isInitFunctionOfCoordinates_v<T>>,
                    std::is_assignable<_GridElement, T>
                >
            };

            _BaseGrid _baseGrid;

            template
            <typename T, std::enable_if_t<_isDefaultValue_v<T>, bool> = true>
            static constexpr _BaseGrid _initBaseGrid(T const defaultValue)
            {
                _Row row {};
                // row.fill(defaultValue);
                for (std::size_t columnI {0}; columnI < numOfColumns; ++columnI)
                {
                    // using Variant = std::decay_t<decltype(row.at(columnI))>;
                    // row.at(columnI) = Variant {defaultValue};
                    row.at(columnI) = _GridElement {defaultValue};
                }

                _BaseGrid grid;
                // grid.fill(row);
                for (std::size_t rowI {0}; rowI < numOfColumns; ++rowI)
                {
                    grid.at(rowI) = row;
                }

                return grid;
            }

            template
            <
                typename Function,
                std::enable_if_t
                <_isInitFunctionOfTwoIndices_v<Function>, bool> = true
            >
            static constexpr _BaseGrid _initBaseGrid(Function const initFunc)
            {
                _BaseGrid grid {};

                for (std::size_t rowI {0}; rowI < numOfRows; ++rowI)
                {
                    for
                    (std::size_t columnI {0}; columnI < numOfColumns; ++columnI)
                    {
                        grid.at(rowI).at(columnI) = initFunc(rowI, columnI);
                    }
                }

                return grid;
            }

            template
            <
                typename Function,
                std::enable_if_t
                <_isInitFunctionOfCoordinates_v<Function>, bool> = true
            >
            static constexpr _BaseGrid _initBaseGrid(Function const initFunc)
            {
                auto initFuncOfTwoIndices
                {
                    [initFunc]
                    (std::size_t const rowI, std::size_t const columnI)
                    {
                        return initFunc(Coordinates {rowI, columnI});
                    }
                };
                return _initBaseGrid(initFuncOfTwoIndices);
            }

        public:
            using GridElement_t = _GridElement;
            static inline constexpr std::size_t k_numOfRows {numOfRows};
            static inline constexpr std::size_t k_numOfColumns {numOfColumns};

            constexpr Grid()
            : _baseGrid {}
            {}

            template <typename T>
            constexpr Grid(T const defaultValueOrInitFunc)
            : _baseGrid {_initBaseGrid(defaultValueOrInitFunc)}
            {}

            constexpr _GridElement& get
            (std::size_t const rowI, std::size_t const columnI) &
            {return _baseGrid.at(rowI).at(columnI);}

            constexpr _GridElement&& get
            (std::size_t const rowI, std::size_t const columnI) &&
            {return _baseGrid.at(rowI).at(columnI);}

            constexpr _GridElement const& get
            (std::size_t const rowI, std::size_t const columnI) const&
            {return _baseGrid.at(rowI).at(columnI);}

            constexpr _GridElement const&& get
            (std::size_t const rowI, std::size_t const columnI) const&&
            {return _baseGrid.at(rowI).at(columnI);}

            constexpr decltype(auto) get
            (Coordinates const coords) &
            {return get(coords.first, coords.second);}

            constexpr decltype(auto) get
            (Coordinates const coords) &&
            {return get(coords.first, coords.second);}

            constexpr decltype(auto) get
            (Coordinates const coords) const&
            {return get(coords.first, coords.second);}

            constexpr decltype(auto) get
            (Coordinates const coords) const&&
            {return get(coords.first, coords.second);}

            template <typename T>
            constexpr decltype(auto) get
            (std::size_t const rowI, std::size_t const columnI) &
            {return std::get<T>(_baseGrid.at(rowI).at(columnI));}

            template <typename T>
            constexpr decltype(auto) get
            (std::size_t const rowI, std::size_t const columnI) &&
            {return std::get<T>(_baseGrid.at(rowI).at(columnI));}

            template <typename T>
            constexpr decltype(auto) get
            (std::size_t const rowI, std::size_t const columnI) const&
            {return std::get<T>(_baseGrid.at(rowI).at(columnI));}

            template <typename T>
            constexpr decltype(auto) get
            (std::size_t const rowI, std::size_t const columnI) const&&
            {return std::get<T>(_baseGrid.at(rowI).at(columnI));}

            template <typename T>
            constexpr decltype(auto) get
            (Coordinates const coords) &
            {return get<T>(coords.first, coords.second);}

            template <typename T>
            constexpr decltype(auto) get
            (Coordinates const coords) &&
            {return get<T>(coords.first, coords.second);}

            template <typename T>
            constexpr decltype(auto) get
            (Coordinates const coords) const&
            {return get<T>(coords.first, coords.second);}

            template <typename T>
            constexpr decltype(auto) get
            (Coordinates const coords) const&&
            {return get<T>(coords.first, coords.second);}

            template <typename T>
            constexpr void set
            (
                T const value,
                std::size_t const rowI,
                std::size_t const columnI
            ) {get(rowI, columnI) = value;}

            template <typename T>
            constexpr void set
            (
                T const value,
                Coordinates const coords
            ) {set<T>(value, coords.first, coords.second);}
    };


    inline constexpr std::size_t k_defaultNumOfRows {1000};
    inline constexpr std::size_t k_defaultNumOfColumns {1000};

    template <typename ComplexNumber>
    inline constexpr ComplexNumber k_complexInfinity
    {std::numeric_limits<typename ComplexNumber::value_type>::infinity()};

    template
    <
        typename ComplexNumber,
        bool usingUpperHalfOnly=true,
        std::size_t numOfRows=k_defaultNumOfRows,
        std::size_t numOfColumns=k_defaultNumOfColumns
    >
    class MandelbrotSetGrid
    {
        private:
            using _MS = MandelbrotSet<ComplexNumber>;
            using _RealNumber = typename ComplexNumber::value_type;

            Grid
            <
                numOfRows,
                numOfColumns,
                ComplexNumber
            > _valueAfterIterationsGrid;

            // Grid
            // <
            //     numOfRows,
            //     numOfColumns,
            //     bool
            // > notDivergeGrid;

        public:
            static inline constexpr std::size_t k_numOfRows {numOfRows};
            static inline constexpr std::size_t k_numOfColumns {numOfColumns};
            static inline constexpr _RealNumber k_gridLeftBound
            {_MS::k_realPartLowerBound};
            static inline constexpr _RealNumber k_gridRightBound
            {_MS::k_realPartUpperBound};
            static inline constexpr _RealNumber k_gridBottomBound
            {usingUpperHalfOnly ? 0.0 : _MS::k_imagPartLowerBound};
            static inline constexpr _RealNumber k_gridTopBound
            {_MS::k_imagPartUpperBound};
            static inline constexpr _RealNumber k_gridArea
            {
                (k_gridRightBound - k_gridLeftBound) *
                (k_gridTopBound - k_gridBottomBound)
            };

            // TODO
            static inline constexpr _RealNumber
            k_sampleRealAxisStepOffset {0.0};
            static inline constexpr _RealNumber
            k_sampleImagAxisStepOffset {1.0};

            static constexpr ComplexNumber gridCoordsToComplexNum
            (std::size_t const rowI, std::size_t const columnI);

            static constexpr ComplexNumber gridCoordsToComplexNum
            (Coordinates const coords)
            {return gridCoordsToComplexNum(coords.first, coords.second);}

            static constexpr Coordinates complexNumToGridCoords
            (ComplexNumber const complexNumber);

            static constexpr bool gridCoordsOutOfRange
            (std::size_t const rowI, std::size_t const columnI)
            {return rowI >= numOfRows || columnI >= numOfColumns;}

            static constexpr bool gridCoordsOutOfRange
            (Coordinates const coords)
            {return gridCoordsOutOfRange(coords.first, coords.second);}

            static constexpr bool complexNumOutOfGridRange
            (ComplexNumber const complexNumber);

            static constexpr bool complexNumOutOfMandelbrotSetRange
            (ComplexNumber const complexNumber);

            constexpr MandelbrotSetGrid()
            // : _valueAfterIterationsGrid {_MS::k_initValueForIteration}
            : _valueAfterIterationsGrid {}
            {}

            constexpr void iterateOneMoreTime
            (std::size_t const rowI, std::size_t const columnI);

            constexpr void iterateOneMoreTime
            (Coordinates const coords)
            {iterateOneMoreTime(coords.first, coords.second);}

            template <typename StopCriterion>
            constexpr void iterateUntil
            (
                std::size_t const rowI,
                std::size_t const columnI,
                StopCriterion criterion
            );

            template <typename StopCriterion>
            constexpr void iterateUntil
            (Coordinates const coords, StopCriterion criterion)
            {iterateUntil(coords.first, coords.second, criterion);}

            template <typename StopCriterion>
            constexpr void iterateAllUntil(StopCriterion criterion);

            constexpr _RealNumber getApproximatedArea()
            {
                std::size_t count {0};

                // for (std::size_t rowI {0}; rowI < numOfRows; ++rowI)
                // {
                //     for
                //     (std::size_t columnI {0}; columnI < numOfColumns; ++columnI)
                //     {
                //         auto& element
                //         {_valueAfterIterationsGrid.get(rowI, columnI)};

                //         if
                //         (
                //             // std::holds_alternative<ComplexNumber>(element) &&
                //             std::get<ComplexNumber>(element) != k_complexInfinity<ComplexNumber>
                //         )
                //         {
                //             ++count;
                //         }
                //     }
                // }

                #pragma omp parallel for collapse(2) reduction(+: count)
                for (std::size_t rowI = 0; rowI < numOfRows; ++rowI)
                {
                    for
                    (std::size_t columnI = 0; columnI < numOfColumns; ++columnI)
                    {
                        auto const& value
                        {_valueAfterIterationsGrid.get(rowI, columnI)};

                        if (value != k_complexInfinity<ComplexNumber>)
                        {
                            ++count;
                        }
                    }
                }

                _RealNumber const proportion
                {
                    static_cast<_RealNumber>(count) /
                    static_cast<_RealNumber>(numOfRows * numOfColumns)
                };

                if constexpr (usingUpperHalfOnly)
                {
                    return 2 * k_gridArea * proportion;
                }
                else
                {
                    return k_gridArea * proportion;
                }

            }
    };
}

#include "MandelbrotSet.tpp"

#endif