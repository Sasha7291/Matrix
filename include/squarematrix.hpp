#pragma once

#include <cmath>

#include "matrix.hpp"


namespace static_matrix
{

    template<NumberType T, const std::size_t N>
    class SquareMatrix : public Matrix<T, N, N>
    {

    public:
        SquareMatrix() noexcept : Matrix<T, N, N>() {}
        SquareMatrix(const T value) noexcept : Matrix<T, N, N>(value) {}
        SquareMatrix(std::span<T> array) try : Matrix<T, N, N>(array) {} catch (const Exception &exception) { throw exception; }
        ~SquareMatrix() noexcept = default;

        SquareMatrix(const Matrix<T, N, N> &other) noexcept;
        SquareMatrix(Matrix<T, N, N> &&other) noexcept;
        SquareMatrix &operator=(const Matrix<T, N, N> &other) noexcept;
        SquareMatrix &operator=(Matrix<T, N, N> &&other) noexcept;

        [[nodiscard]] T determinant() const;
        void inverte();
        [[nodiscard]] SquareMatrix inverted() const;
        [[nodiscard]] Matrix<T, N, N> toMatrix() const;
        [[nodiscard]] T trace() const;
        [[nodiscard]] SquareMatrix<T, N - 1> minor(const std::size_t m, const std::size_t n) const;

    };

    template<NumberType T, const std::size_t N>
    SquareMatrix<T, N>::SquareMatrix(const Matrix<T, N, N> &other) noexcept
    {
        this->data_ = other.data();
    }

    template<NumberType T, const std::size_t N>
    SquareMatrix<T, N>::SquareMatrix(Matrix<T, N, N> &&other) noexcept
    {
        this->data_ = std::move(other.data());
    }

    template<NumberType T, const std::size_t N>
    SquareMatrix<T, N> &SquareMatrix<T, N>::operator=(const Matrix<T, N, N> &other) noexcept
    {
        this->data_ = other.data();
        return *this;
    }

    template<NumberType T, const std::size_t N>
    SquareMatrix<T, N> &SquareMatrix<T, N>::operator=(Matrix<T, N, N> &&other) noexcept
    {
        this->data_ = std::move(other.data());
        return *this;
    }

    template<NumberType T, const std::size_t N>
    T SquareMatrix<T, N>::determinant() const
    {
        if constexpr (N == 1ull)
        {
            return this->at(0, 0);
        }
        else if constexpr (N == 2ull)
        {
            return this->at(0, 0) * this->at(1, 1) - this->at(0, 1) * this->at(1, 0);
        }
        else
        {
            T result = static_cast<T>(0);

            for (auto n = 0ull; n < N; ++n)
                result += std::pow(-1, n) * this->at(0, n) * minor(0, n).determinant();

            return result;
        }
    }

    template<NumberType T, const std::size_t N>
    void SquareMatrix<T, N>::inverte()
    {
        SquareMatrix temp = this->transposed();

        double factor = 1.0 / determinant();
        for (auto m = 0ull; m < N; ++m)
            for (auto n = 0ull; n < N; ++n)
                this->operator()(m, n) = std::pow(-1, m + n) * factor * temp.minor(m, n).determinant();
    }

    template<NumberType T, const std::size_t N>
    SquareMatrix<T, N> SquareMatrix<T, N>::inverted() const
    {
        SquareMatrix result(*this);
        result.inverte();
        return result;
    }

    template<NumberType T, const std::size_t N>
    Matrix<T, N, N> SquareMatrix<T, N>::toMatrix() const
    {
        Matrix<T, N, N> result;
        result.data() = this->data_;
        return result;
    }

    template<NumberType T, const std::size_t N>
    SquareMatrix<T, N - 1> SquareMatrix<T, N>::minor(const std::size_t m, const std::size_t n) const
    {
        SquareMatrix<T, N - 1> result;

        auto row_offset = 0ull;
        auto column_offset = 0ull;
        for (auto i = 0ull; i < N - 1; ++i)
        {
            if (i == m)
                ++row_offset;

            for (auto j = 0ull; j < N - 1; ++j)
            {
                if (j == n)
                    ++column_offset;

                result(i, j) = this->at(i + row_offset, j + column_offset);
            }

            column_offset = 0;
        }

        return result;
    }

    template<NumberType T, const std::size_t N>
    T SquareMatrix<T, N>::trace() const
    {
        T result;

        for (auto n = 0; n < N; ++n)
            result += this->at(n, n);

        return result;
    }

}
