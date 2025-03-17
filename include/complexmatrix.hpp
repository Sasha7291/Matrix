#pragma once

#include "matrix.hpp"

#include <complex>


namespace Matrix
{

    template<NumberType T, const std::size_t M, const std::size_t N>
    class ComplexMatrix : public Matrix<std::complex<T>, M, N>
    {

    public:
        ComplexMatrix() noexcept : Matrix<T, M, N>() {}
        ComplexMatrix(const std::complex<T> value) noexcept : Matrix<std::complex<T>, M, N>(value) {}
        ComplexMatrix(std::span<std::complex<T>> array)
        try : Matrix<std::complex<T>, M, N>(array) {}
        catch (const MatrixException &exception) { throw exception; }
        ~ComplexMatrix() noexcept = default;

        void conjugate();
        [[nodiscard]] ComplexMatrix conjugated() const;
        void conjugatedTranspose();
        [[nodiscard]] ComplexMatrix conjugatedTransposed() const;

    };

    template<NumberType T, const std::size_t M, const std::size_t N>
    void ComplexMatrix<T, M, N>::conjugate()
    {
        for (auto m = 0ull; m < M; ++m)
            for (auto n = 0ull; n < N; ++n)
                this->operator()(m, n).conj();
    }

    template<NumberType T, const std::size_t M, const std::size_t N>
    ComplexMatrix<T, M, N> ComplexMatrix<T, M, N>::conjugated() const
    {
        ComplexMatrix result(*this);
        result.conjugate();
        return result;
    }

    template<NumberType T, const std::size_t M, const std::size_t N>
    void ComplexMatrix<T, M, N>::conjugatedTranspose()
    {
        *this = this->transposed();
        conjugate();
    }

    template<NumberType T, const std::size_t M, const std::size_t N>
    ComplexMatrix<T, M, N> ComplexMatrix<T, M, N>::conjugatedTransposed() const
    {
        ComplexMatrix result(*this);
        result.conjugatedTranspose();
        return result;
    }

}
