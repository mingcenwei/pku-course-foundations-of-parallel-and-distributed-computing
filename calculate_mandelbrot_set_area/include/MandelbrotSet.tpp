#ifndef FPDC_2019__MANDELBROT_SET_TPP
#define FPDC_2019__MANDELBROT_SET_TPP

#include "MandelbrotSet.hpp"

#include "omp.h"

#include <cstddef>
#include <complex>
#include <stdexcept>
#include <type_traits>

namespace Fpdc2019
{
    // template <typename ComplexNumber, typename StopCriterion>
    // constexpr bool willDiverge
    // (ComplexNumber const parameter, StopCriterion criterion)
    // {
    //     if constexpr (std::imag(parameter) < 0)
    //     {
    //         return willDiverge(-std::conj(parameter), criterion);
    //     }
    //     else
    //     {
    //         return _willDiverge_impl(parameter, criterion);
    //     }
    // }

    // template <typename ComplexNumber, typename StopCriterion>
    // constexpr bool _willDiverge_impl
    // (
    //     ComplexNumber const parameter,
    //     StopCriterion criterion,
    //     ComplexNumber const previousValue=
    //         MandelbrotSet<ComplexNumber>::k_initValueForIteration
    // )
    // {
    //     static_assert
    //     (
    //         std::imag(parameter) >= 0,
    //         "The imaginary part of the parameter c is less than zero"
    //     );

    //     using MS = MandelbrotSet<ComplexNumber>;

    //     if constexpr
    //     (
    //         std::real(previousValue) < MS::k_realPartLowerBound ||
    //         std::real(previousValue) > MS::k_realPartUpperBound ||
    //         std::imag(previousValue) > MS::imagPartUpperBound
    //     )
    //     {
    //         return true;
    //     }
    //     else
    //     {
    //         ComplexNumber const newValue
    //         {iteratedFunction(previousValue, parameter)};

    //         criterion.update(newValue, previousValue);

    //         if constexpr (criterion.shouldStop())
    //         {
    //             return false;
    //         }
    //         else
    //         {
    //             return _willDiverge_impl(parameter, criterion, newValue);
    //         }
    //     }
    // }

    template <typename ComplexNumber>
    constexpr void IterationStopCriterion<ComplexNumber>::update
    (ComplexNumber const newValue, ComplexNumber const previousValue)
    {
        // switch (_remainingIterationCount)
        // {
        //     case 0:
        //         _shouldStop = true;
        //         break;
        //     case 1:
        //         --_remainingIterationCount;
        //         _shouldStop = true;
        //         break;
        //     default:
        //         if (false /* TODO */)
        //         {}
        //         else
        //         {
        //             static_cast<void>(newValue);
        //             static_cast<void>(previousValue);
        //         }
        //         break;
        // }

        // if (_remainingIterationCount == 0)
        // {
        //     _shouldStop = true;
        //     return;
        // }

        --_remainingIterationCount;

        if (_remainingIterationCount == 0)
        {
            _shouldStop = true;
        }
        else if (false /* TODO */)
        {}
        else
        {
            static_cast<void>(newValue);
            static_cast<void>(previousValue);
        }

        return;
    }

    template
    <
        typename ComplexNumber,
        bool usingUpperHalfOnly,
        std::size_t numOfRows,
        std::size_t numOfColumns
    >
    constexpr ComplexNumber MandelbrotSetGrid
    <ComplexNumber, usingUpperHalfOnly, numOfRows, numOfColumns>::
    gridCoordsToComplexNum(std::size_t const rowI, std::size_t const columnI)
    {
        if (rowI > numOfRows || columnI > numOfColumns)
        {
            throw std::out_of_range {"Grid coordinates out of range"};
        }

        constexpr auto horizontalStep
        {(k_gridRightBound - k_gridLeftBound) / numOfColumns};

        constexpr auto verticalStep
        {(k_gridTopBound - k_gridBottomBound) / numOfRows};

        auto const numericColumnI
        {
            static_cast<decltype(horizontalStep)>(columnI) + k_sampleRealAxisStepOffset
        };

        auto const numericRowI
        {
            static_cast<decltype(verticalStep)>(rowI) + k_sampleImagAxisStepOffset
        };

        auto const realAxis {k_gridLeftBound + numericColumnI * horizontalStep};

        auto const imagAxis {k_gridBottomBound + numericRowI * verticalStep};

        return ComplexNumber {realAxis, imagAxis};
    }

    template <typename RealNumber>
    constexpr std::ptrdiff_t _constexpr_floor(RealNumber const num)
    {
        std::ptrdiff_t const numTruncated {static_cast<std::ptrdiff_t>(num)};
        if
        (
            num == static_cast<RealNumber>(numTruncated) ||
            num > 0
        )
        {
            return numTruncated;
        }
        else
        {
            return (numTruncated - 1);
        }
    }

    template
    <
        typename ComplexNumber,
        bool usingUpperHalfOnly,
        std::size_t numOfRows,
        std::size_t numOfColumns
    >
    constexpr Coordinates MandelbrotSetGrid
    <ComplexNumber, usingUpperHalfOnly, numOfRows, numOfColumns>::
    complexNumToGridCoords(ComplexNumber const complexNumber)
    {
        if (complexNumOutOfGridRange(complexNumber))
        {
            throw std::range_error {"Complex number out of range"};
        }

        constexpr auto horizontalStep
        {(k_gridRightBound - k_gridLeftBound) / numOfColumns};

        constexpr auto verticalStep
        {(k_gridTopBound - k_gridBottomBound) / numOfRows};

        auto const realAxisOffset
        {std::real(complexNumber) - k_gridLeftBound};

        auto const imagAxisOffset
        {std::imag(complexNumber) - k_gridBottomBound};

        auto const numericColumnI {realAxisOffset / horizontalStep};
        auto const numericRowI {imagAxisOffset / verticalStep};

        auto const rowI
        {static_cast<std::size_t>(_constexpr_floor(numericRowI))};
        auto const columnI
        {static_cast<std::size_t>(_constexpr_floor(numericColumnI))};

        return Coordinates {rowI, columnI};
    }

    // template
    // <
    //     typename ComplexNumber,
    //     std::size_t numOfRows,
    //     std::size_t numOfColumns
    // >
    // constexpr bool MandelbrotSetGrid<ComplexNumber, numOfRows, numOfColumns>::
    // gridCoordsOutOfRange(std::size_t const rowI, std::size_t const columnI)
    // {
    //     return rowI < 0 ||
    //         rowI >= numOfRows ||
    //         columnI < 0 ||
    //         columnI >= numOfColumns;
    // }

    template
    <
        typename ComplexNumber,
        bool usingUpperHalfOnly,
        std::size_t numOfRows,
        std::size_t numOfColumns
    >
    constexpr bool MandelbrotSetGrid
    <ComplexNumber, usingUpperHalfOnly, numOfRows, numOfColumns>::
    complexNumOutOfGridRange(ComplexNumber const complexNumber)
    {
        auto const realPart {std::real(complexNumber)};
        auto const imagPart {std::imag(complexNumber)};

        return realPart < k_gridLeftBound ||
            realPart > k_gridRightBound ||
            imagPart < k_gridBottomBound ||
            imagPart > k_gridTopBound;
    }

    template
    <
        typename ComplexNumber,
        bool usingUpperHalfOnly,
        std::size_t numOfRows,
        std::size_t numOfColumns
    >
    constexpr bool MandelbrotSetGrid
    <ComplexNumber, usingUpperHalfOnly, numOfRows, numOfColumns>::
    complexNumOutOfMandelbrotSetRange(ComplexNumber const complexNumber)
    {
        auto const realPart {std::real(complexNumber)};
        auto const imagPart {std::imag(complexNumber)};

        return realPart < -2.0 ||
            realPart > 2.0 ||
            imagPart < -2.0 ||
            imagPart > 2.0;
    }

    template
    <
        typename ComplexNumber,
        bool usingUpperHalfOnly,
        std::size_t numOfRows,
        std::size_t numOfColumns
    >
    constexpr void MandelbrotSetGrid
    <ComplexNumber, usingUpperHalfOnly, numOfRows, numOfColumns>::
    iterateOneMoreTime(std::size_t const rowI, std::size_t const columnI)
    {
        // auto& element {_valueAfterIterationsGrid.get(rowI, columnI)};
        // auto const parameter {gridCoordsToComplexNum(rowI, columnI)};

        // if
        // (
        //     // std::holds_alternative<ComplexNumber>(element) &&
        //     std::get<ComplexNumber>(element) != k_complexInfinity<ComplexNumber>
        // )
        // {
        //     auto& value {std::get<ComplexNumber>(element)};
        //     value = _MS::iteratedFunction(value, parameter);

        //     if (complexNumOutOfRange(value))
        //     {
        //         value = k_complexInfinity<ComplexNumber>;
        //     }
        // }

        auto& value {_valueAfterIterationsGrid.get(rowI, columnI)};

        if (value != k_complexInfinity<ComplexNumber>)
        {
            auto const parameter {gridCoordsToComplexNum(rowI, columnI)};
            value = _MS::iteratedFunction(value, parameter);

            if (complexNumOutOfMandelbrotSetRange(value))
            {
                value = k_complexInfinity<ComplexNumber>;
            }
        }

        return;
    }

    template
    <
        typename ComplexNumber,
        bool usingUpperHalfOnly,
        std::size_t numOfRows,
        std::size_t numOfColumns
    >
    template <typename StopCriterion>
    constexpr void MandelbrotSetGrid
    <ComplexNumber, usingUpperHalfOnly, numOfRows, numOfColumns>::
    iterateUntil
    (
        std::size_t const rowI,
        std::size_t const columnI,
        StopCriterion criterion
    )
    {
        // while (!criterion.shouldStop())
        // {
        //     auto& element {_valueAfterIterationsGrid.get(rowI, columnI)};
        //     auto const parameter {gridCoordsToComplexNum(rowI, columnI)};

        //     if
        //     (
        //         // std::holds_alternative<ComplexNumber>(element) &&
        //         std::get<ComplexNumber>(element) != k_complexInfinity<ComplexNumber>
        //     )
        //     {
        //         auto& value {std::get<ComplexNumber>(element)};
        //         auto const previousValue {value};
        //         value = _MS::iteratedFunction(value, parameter);

        //         if (complexNumOutOfRange(value))
        //         {
        //             value = k_complexInfinity<ComplexNumber>;
        //             break;
        //         }

        //         criterion.update(value, previousValue);
        //     }
        //     else
        //     {
        //         break;
        //     }

        // }

        auto& value {_valueAfterIterationsGrid.get(rowI, columnI)};
        auto const parameter {gridCoordsToComplexNum(rowI, columnI)};

        if (k_complexInfinity<ComplexNumber> == value)
        {
            return;
        }
        else
        {
            while (!criterion.shouldStop())
            {
                auto const previousValue {value};
                value = _MS::iteratedFunction(value, parameter);

                if (complexNumOutOfMandelbrotSetRange(value))
                {
                    value = k_complexInfinity<ComplexNumber>;
                    break;
                }
                else
                {
                    criterion.update(value, previousValue);
                }
            }

            return;
        }
    }

    template
    <
        typename ComplexNumber,
        bool usingUpperHalfOnly,
        std::size_t numOfRows,
        std::size_t numOfColumns
    >
    template <typename StopCriterion>
    constexpr void MandelbrotSetGrid
    <ComplexNumber, usingUpperHalfOnly, numOfRows, numOfColumns>::
    iterateAllUntil(StopCriterion criterion)
    {
        #pragma omp parallel for collapse(2)
        for (std::size_t rowI = 0; rowI < numOfRows; ++rowI)
        {
            for (std::size_t columnI = 0; columnI < numOfColumns; ++columnI)
            {
                iterateUntil(rowI, columnI, criterion);
            }
        }

        return;
    }
}

#endif