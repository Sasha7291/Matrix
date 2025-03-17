#pragma once

#include "squarematrix.hpp"


class LinearSystemSolverException : std::runtime_error
{
public:
    LinearSystemSolverException(const std::string &message) : std::runtime_error(message) {}
};


class LinearSystemSolver
{
public:
    LinearSystemSolver() noexcept = default;
    ~LinearSystemSolver() noexcept = default;

    LinearSystemSolver(const LinearSystemSolver &) = delete;
    LinearSystemSolver(LinearSystemSolver &&) = delete;
    LinearSystemSolver &operator=(const LinearSystemSolver &) = delete;
    LinearSystemSolver &operator=(LinearSystemSolver &&) = delete;

    template<NumberType T, const std::size_t N>
    [[nodiscard]] Matrix::Matrix<T, N, 1> operator()(
        const Matrix::SquareMatrix<T, N> &A,
        const Matrix::Matrix<T, N, 1> &B
    ) const;

};

template<NumberType T, const std::size_t N>
Matrix::Matrix<T, N, 1> LinearSystemSolver::operator()(
    const Matrix::SquareMatrix<T, N> &A,
    const Matrix::Matrix<T, N, 1> &B
) const
{
    if (A.determinant() == 0)
        throw LinearSystemSolverException("det(A) == 0, A is irreversible");

    return A.inverted().toMatrix() * B;
}
