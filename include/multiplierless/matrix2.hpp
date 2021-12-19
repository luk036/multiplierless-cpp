#pragma once

#include <tuple>    // import std::tie()
#include <utility>  // import std::move
#include "vector2.hpp"

namespace numeric {

    /**
     * @brief matrix2
     *
     */
    template <typename T1, typename T2 = T1> class matrix2 : public vector2<T1, T2> {
      public:
        /**
         * @brief Construct a new matrix2 object
         *
         * @param x
         * @param y
         */
        constexpr matrix2(T1&& x, T2&& y) noexcept : vector2<T1, T2>{std::move(x), std::move(y)} {}

        /**
         * @brief Construct a new matrix2 object
         *
         * @param x
         * @param y
         */
        constexpr matrix2(const T1& x, const T2& y) : vector2<T1, T2>{x, y} {}


        /** @name Arithmetic operators
         *  definie +, -, *, /, +=, -=, *=, /=, etc.
         */
        ///@{

        /**
         * @brief Negate
         *
         * @return matrix2
         */
        constexpr auto operator-() const -> matrix2 { return matrix2(-this->x(), -this->y()); }

        /**
         * @brief Add
         *
         * @tparam U1
         * @tparam U2
         * @param[in] other
         * @return matrix2&
         */
        template <typename U1, typename U2> constexpr auto operator+=(const matrix2<U1, U2>& other)
            -> matrix2& {
            this->_x += other.x();
            this->_y += other.y();
            return *this;
        }

        /**
         * @brief Substract
         *
         * @tparam U1
         * @tparam U2
         * @param[in] other
         * @return matrix2&
         */
        template <typename U1, typename U2>  //
        constexpr auto operator-=(const matrix2<U1, U2>& other) -> matrix2& {
            this->_x -= other.x();
            this->_y -= other.y();
            return *this;
        }

        /**
         * @brief Multiply
         *
         * @tparam R
         * @param[in] alpha
         * @return matrix2&
         */
        template <typename R> constexpr auto operator*=(const R& alpha) -> matrix2& {
            this->_x *= alpha;
            this->_y *= alpha;
            return *this;
        }

        /**
         * @brief Divide
         *
         * @tparam R
         * @param[in] alpha
         * @return matrix2&
         */
        template <typename R> constexpr auto operator/=(const R& alpha) -> matrix2& {
            this->_x /= alpha;
            this->_y /= alpha;
            return *this;
        }

        /**
         * @brief Add
         *
         * @tparam U1
         * @tparam U2
         * @param[in] x
         * @param[in] y
         * @return matrix2
         */
        template <typename U1, typename U2>  //
        friend constexpr auto operator+(matrix2 x, const matrix2<U1, U2>& y) -> matrix2 {
            return x += y;
        }

        /**
         * @brief Substract
         *
         * @tparam U1
         * @tparam U2
         * @param[in] x
         * @param[in] y
         * @return matrix2
         */
        template <typename U1, typename U2>  //
        friend constexpr auto operator-(matrix2 x, const matrix2<U1, U2>& y) -> matrix2 {
            return x -= y;
        }

        /**
         * @brief Multiply by a scalar
         *
         * @tparam R
         * @param[in] x
         * @param[in] alpha
         * @return matrix2
         */
        template <typename R> friend constexpr auto operator*(matrix2 x, const R& alpha)
            -> matrix2 {
            return x *= alpha;
        }

        /**
         * @brief Multiply (by a scalar)
         *
         * @tparam R
         * @param[in] alpha
         * @param[in] x
         * @return matrix2
         */
        template <typename R> friend constexpr auto operator*(const R& alpha, matrix2 x)
            -> matrix2 {
            return x *= alpha;
        }

        /**
         * @brief Divide (by a scalar)
         *
         * @tparam R
         * @param[in] x
         * @param[in] alpha
         * @return matrix2
         */
        template <typename R> friend constexpr auto operator/(matrix2 x, const R& alpha)
            -> matrix2 {
            return x /= alpha;
        }

        /**
         * @brief
         *
         * @tparam U1
         * @tparam U2
         * @param[in] other
         * @return constexpr auto
         */
        template <typename U1, typename U2>  //
        [[nodiscard]] constexpr auto mdot(const vector2<U1, U2>& other) const {
            return vector2<U1, U2>{this->_x.dot(other), this->_y.dot(other)};
        }


        /**
         * @brief
         *
         * @tparam U1
         * @tparam U2
         * @param[in] other
         * @return constexpr auto
         */
        [[nodiscard]] constexpr auto det() const {
            return this->x().x() * this->y().y() - this->x().y() * this->y().x();
        }

        /**
         * @brief
         *
         * @tparam U1
         * @tparam U2
         * @param[in] other
         * @return constexpr auto
         */
        template <typename U1, typename U2>  //
        [[nodiscard]] constexpr auto solve(const vector2<U1, U2>& v) {
            auto adjoint = matrix2{vector2{this->y().y(), -this->x().y()}, vector2{-this->y().x(), this->x().x()}};
            return adjoint.mdot(v) / this->det();
        }

        ///@}

    };
}  // namespace numeric
