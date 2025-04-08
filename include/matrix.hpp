#pragma once

#include <algorithm>
#include <array>
#include <concepts>
#include <iostream>
#include <span>


template<class T>
concept NumberType = std::integral<T> || std::floating_point<T>;

namespace static_matrix
{
	
class Exception : public std::runtime_error
{
public:
	Exception(const std::string &message) : std::runtime_error(message) {};
};
	
template<NumberType T> using Element = T;
template<NumberType T, const std::size_t N> using Row = std::array<Element<T>, N>;
template<NumberType T, const std::size_t M, const std::size_t N> using Data = std::array<Row<T, N>, M>;

template<NumberType T, const std::size_t M, const std::size_t N>
class Matrix;
template<NumberType T, const std::size_t M, const std::size_t N>
std::ostream& operator<<(std::ostream &os, const Matrix<T, M, N> &matrix);

template<NumberType T, const std::size_t M, const std::size_t N>
class Matrix
{

public:
	template<NumberType K, const std::size_t S, const std::size_t P>
	friend std::ostream& operator<<(std::ostream &os, const Matrix<K, S, P> &matrix);

	Matrix() noexcept = default;
	Matrix(const T value) noexcept;
	Matrix(std::span<T> array);
	virtual ~Matrix() noexcept = default;

	Matrix(const Matrix<T, M, N> &other) noexcept;
	Matrix(Matrix<T, M, N> &&other) noexcept;
	Matrix &operator=(const Matrix<T, M, N> &other) noexcept;
	Matrix &operator=(Matrix<T, M, N> &&other) noexcept;

	[[nodiscard]] Matrix operator+(const Matrix<T, M, N> &other) noexcept;
	[[nodiscard]] Matrix operator-(const Matrix<T, M, N> &other) noexcept;
	[[nodiscard]] Matrix operator*(const T value) noexcept;
	[[nodiscard]] T &operator()(const std::size_t m, const std::size_t n) noexcept;
	[[nodiscard]] const T &operator()(const std::size_t m, const std::size_t n) const noexcept;

	[[nodiscard]] T at(const std::size_t m, const std::size_t n) const;
	[[nodiscard]] Row<T, M> column(const std::size_t n) const;
	[[nodiscard]] Data<T, M, N> &data();
	[[nodiscard]] const Data<T, M, N> &data() const;
	void fill(const T value);
	inline void ones() { fill(static_cast<T>(1)); }
	[[nodiscard]] Row<T, N> &row(const std::size_t m);
	[[nodiscard]] const Row<T, N> &row(const std::size_t m) const;
	[[nodiscard]] Matrix<T, N, M> transposed() const noexcept;
	inline void zeros() { fill(static_cast<T>(0)); }

protected:
	Data<T, M, N> data_;

};

template<NumberType T, const std::size_t M, const std::size_t N>
std::ostream& operator<<(std::ostream &os, const Matrix<T, M, N> &matrix)
{
	os << "Matrix(" << M << ", " << N << ")\n{\n";

	for (auto m = 0; m < M; ++m)
	{
		os << "\t";
		for (auto n = 0; n < N; ++n)
			os << matrix(m, n) << " ";
		os << "\n";
	}

	return os << "}" << std::endl;
}

template<NumberType T, const std::size_t M, const std::size_t R, const std::size_t N>
Matrix<T, M, N> operator*(const Matrix<T, M, R> &mat1, const Matrix<T, R, N> &mat2)
{
	Matrix<T, M, N> result;

	for (auto m = 0ull; m < M; ++m)
		for (auto n = 0ull; n < N; ++n)
			for (auto r = 0ull; r < R; ++r)
				result(m, n) += mat1(m, r) * mat2(r, n);

	return result;
}

template<NumberType T, const std::size_t M, const std::size_t N>
Matrix<T, M, N>::Matrix(const T value) noexcept
{
	fill(value);
}

template<NumberType T, const std::size_t M, const std::size_t N>
Matrix<T, M, N>::Matrix(std::span<T> array)
{
	if (array.size() < N * M)
		throw Exception("Initialize array is smaller than nessesary");

	for (auto m = 0ull; m < M; ++m)
		for (auto n = 0ull; n < N; ++n)
			data_[m][n] = array[n + m * N];
}

template<NumberType T, const std::size_t M, const std::size_t N>
Matrix<T, M, N>::Matrix(const Matrix<T, M, N> &other) noexcept
{
	data_ = other.data_;
}

template<NumberType T, const std::size_t M, const std::size_t N>
Matrix<T, M, N>::Matrix(Matrix<T, M, N> &&other) noexcept
{
	data_ = std::move(other.data_);
}

template<NumberType T, const std::size_t M, const std::size_t N>
Matrix<T, M, N> &Matrix<T, M, N>::operator=(const Matrix<T, M, N> &other) noexcept
{
	data_ = other.data_;
	return *this;
}

template<NumberType T, const std::size_t M, const std::size_t N>
Matrix<T, M, N> &Matrix<T, M, N>::operator=(Matrix<T, M, N> &&other) noexcept
{
	data_ = std::move(other.data_);
	return *this;
}

template<NumberType T, const std::size_t M, const std::size_t N>
Matrix<T, M, N> Matrix<T, M, N>::operator+(const Matrix<T, M, N> &other) noexcept
{
	Matrix<T, M, N> result(other);

	for (auto it1 = result.data_.begin(), it2 = data_.cbegin(); it1 != result.data_.end(); ++it1, ++it2)
		for (auto jt1 = it1->begin(), jt2 = it2->cbegin(); jt1 != it1->end(); ++jt1, ++jt2)
			(*jt1) += (*jt2);

	return result;
}

template<NumberType T, const std::size_t M, const std::size_t N>
Matrix<T, M, N> Matrix<T, M, N>::operator-(const Matrix<T, M, N> &other) noexcept
{
	Matrix<T, M, N> result(other);

	for (auto it1 = result.data_.begin(), it2 = data_.cbegin(); it1 != result.data_.end(); ++it1, ++it2)
		for (auto jt1 = it1->begin(), jt2 = it2->cbegin(); jt1 != it1->end(); ++jt1, ++jt2)
			(*jt1) -= (*it1);

	return result;
}

template<NumberType T, const std::size_t M, const std::size_t N>
Matrix<T, M, N> Matrix<T, M, N>::operator*(const T value) noexcept
{
	Matrix<T, M, N> result(*this);

	for (auto it1 = result.data_.begin(); it1 != result.data_.end(); ++it1)
		for (auto jt1 = it1->begin(); jt1 != it1->end(); ++jt1)
			(*jt1) *= value;

	return result;
}

template<NumberType T, const std::size_t M, const std::size_t N>
T &Matrix<T, M, N>::operator()(const std::size_t m, const std::size_t n) noexcept
{
	return data_[m][n];
}

template<NumberType T, const std::size_t M, const std::size_t N>
const T &Matrix<T, M, N>::operator()(const std::size_t m, const std::size_t n) const noexcept
{
	return const_cast<const T &>(data_[m][n]);
}

template<NumberType T, const std::size_t M, const std::size_t N>
T Matrix<T, M, N>::at(const std::size_t m, const std::size_t n) const
{
	if (m >= M || n >= N)
		throw Exception("Index out of range");

	return data_[m][n];
}

template<NumberType T, const std::size_t M, const std::size_t N>
Row<T, M> Matrix<T, M, N>::column(const std::size_t n) const
{
	if (n >= N)
		throw Exception("Index out of range");

	Row<T, M> result;

	for (auto m = 0ull; m < M; ++m)
		result[m] = data_[m][n];

	return result;
}

template<NumberType T, const std::size_t M, const std::size_t N>
Data<T, M, N> &Matrix<T, M, N>::data()
{
	return data_;
}

template<NumberType T, const std::size_t M, const std::size_t N>
const Data<T, M, N> &Matrix<T, M, N>::data() const
{
	return static_cast<const Data<T, M, N> &>(data_);
}

template<NumberType T, const std::size_t M, const std::size_t N>
inline void Matrix<T, M, N>::fill(const T value)
{
	for (auto &row : data_)
		for (auto &element : row)
			element = value;
}

template<NumberType T, const std::size_t M, const std::size_t N>
Row<T, N> &Matrix<T, M, N>::row(const std::size_t m)
{
	if (m >= M)
		throw Exception("Index out of range");

	return data_[m];
}

template<NumberType T, const std::size_t M, const std::size_t N>
const Row<T, N> &Matrix<T, M, N>::row(const std::size_t m) const
{
	if (m >= M)
		throw Exception("Index out of range");

	return static_cast<const Row<T, N> &>(data_[m]);
}

template<NumberType T, const std::size_t M, const std::size_t N>
Matrix<T, N, M> Matrix<T, M, N>::transposed() const noexcept
{
	Matrix<T, N, M> result;

	for (auto m = 0ull; m < M; ++m)
		for (auto n = 0ull; n < N; ++n)
			result(n, m) = operator()(m, n);

	return result;
}

}
