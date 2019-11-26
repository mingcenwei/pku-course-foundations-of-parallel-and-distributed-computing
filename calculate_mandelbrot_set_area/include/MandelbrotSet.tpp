#ifndef FPDC_2019__MANDELBROT_SET_TPP
#define FPDC_2019__MANDELBROT_SET_TPP

#include "MandelbrotSet.hpp"

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
        --_remainingIterationCount;

        if (_remainingIterationCount <= 0)
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
    }

    template
    <
        typename ComplexNumber,
        std::size_t numOfRows,
        std::size_t numOfColumns
    >
    constexpr ComplexNumber
    MandelbrotSetGrid<ComplexNumber, numOfRows, numOfColumns>::
    gridCoordsToComplexNum(std::size_t const rowI, std::size_t const columnI)
    {
        if (gridCoordsOutOfRange(rowI, columnI))
        {
            throw std::out_of_range {"Grid coordinates out of range"};
        }

        constexpr auto horizontalStep
        {
            (_MS::k_realPartUpperBound - _MS::k_realPartLowerBound) / numOfColumns
        };

        constexpr auto verticalStep
        {_MS::k_imagPartUpperBound / numOfRows};

        auto const numericColumnI
        {static_cast<decltype(horizontalStep)>(columnI) + 0.5};

        auto const numericRowI
        {static_cast<decltype(verticalStep)>(rowI) + 0.5};

        auto const realAxis
        {_MS::k_realPartLowerBound + numericColumnI * horizontalStep};

        auto const imagAxis {numericRowI * verticalStep};

        return ComplexNumber {realAxis, imagAxis};
    }

    template <typename RealNumber>
    constexpr std::ptrdiff_t _constexpr_floor(RealNumber const num)
    {
        std::ptrdiff_t const numTruncated {static_cast<std::ptrdiff_t>(num)};
        if (num == static_cast<RealNumber>(numTruncated))
        {
            return numTruncated;
        }
        else if (num > 0)
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
        std::size_t numOfRows,
        std::size_t numOfColumns
    >
    constexpr Coordinates
    MandelbrotSetGrid<ComplexNumber, numOfRows, numOfColumns>::
    complexNumToGridCoords(ComplexNumber const complexNumber)
    {
        if (complexNumOutOfRange(complexNumber))
        {
            throw std::range_error {"Complex number out of range"};
        }

        constexpr auto horizontalStep
        {
            (_MS::k_realPartUpperBound - _MS::k_realPartLowerBound) / numOfColumns
        };

        constexpr auto verticalStep {_MS::imagPartUpperBound / numOfRows};

        auto const realAxisOffset
        {std::real(complexNumber) - _MS::k_realPartLowerBound};

        auto const imagAxisOffset {std::imag(complexNumber)};

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
        std::size_t numOfRows,
        std::size_t numOfColumns
    >
    constexpr bool MandelbrotSetGrid<ComplexNumber, numOfRows, numOfColumns>::
    complexNumOutOfRange(ComplexNumber const complexNumber)
    {
        auto const realPart {std::real(complexNumber)};
        auto const imagPart {std::imag(complexNumber)};

        return realPart < _MS::k_realPartLowerBound ||
            realPart > _MS::k_realPartUpperBound ||
            imagPart < _MS::k_imagPartLowerBound ||
            imagPart > _MS::k_imagPartUpperBound;
    }

    template
    <
        typename ComplexNumber,
        std::size_t numOfRows,
        std::size_t numOfColumns
    >
    constexpr void MandelbrotSetGrid<ComplexNumber, numOfRows, numOfColumns>::
    iterateOneMoreTime(std::size_t const rowI, std::size_t const columnI)
    {
        auto& element {_valueAfterIterationsGrid.get(rowI, columnI)};
        auto const parameter {gridCoordsToComplexNum(rowI, columnI)};

        if
        (
            // std::holds_alternative<ComplexNumber>(element) &&
            std::get<ComplexNumber>(element) != k_complexInfinity<ComplexNumber>
        )
        {
            auto& value {std::get<ComplexNumber>(element)};
            value = _MS::iteratedFunction(value, parameter);

            if (complexNumOutOfRange(value))
            {
                value = k_complexInfinity<ComplexNumber>;
            }
        }

        return;
    }

    template
    <
        typename ComplexNumber,
        std::size_t numOfRows,
        std::size_t numOfColumns
    >
    template <typename StopCriterion>
    constexpr void MandelbrotSetGrid<ComplexNumber, numOfRows, numOfColumns>::
    iterateUntil
    (
        std::size_t const rowI,
        std::size_t const columnI,
        StopCriterion criterion
    )
    {
        while (!criterion.shouldStop())
        {
            auto& element {_valueAfterIterationsGrid.get(rowI, columnI)};
            auto const parameter {gridCoordsToComplexNum(rowI, columnI)};

            if
            (
                // std::holds_alternative<ComplexNumber>(element) &&
                std::get<ComplexNumber>(element) != k_complexInfinity<ComplexNumber>
            )
            {
                auto& value {std::get<ComplexNumber>(element)};
                auto const previousValue {value};
                value = _MS::iteratedFunction(value, parameter);

                if (complexNumOutOfRange(value))
                {
                    value = k_complexInfinity<ComplexNumber>;
                    break;
                }

                criterion.update(value, previousValue);
            }
            else
            {
                break;
            }

        }

        return;
    }

    template
    <
        typename ComplexNumber,
        std::size_t numOfRows,
        std::size_t numOfColumns
    >
    template <typename StopCriterion>
    constexpr void MandelbrotSetGrid<ComplexNumber, numOfRows, numOfColumns>::
    iterateAllUntil(StopCriterion criterion)
    {
        for (std::size_t rowI {0}; rowI < numOfRows; ++rowI)
        {
            for (std::size_t columnI {0}; columnI < numOfColumns; ++columnI)
            {
                iterateUntil(rowI, columnI, criterion);
            }
        }

        return;
    }
}

#endif