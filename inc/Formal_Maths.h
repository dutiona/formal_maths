/*
*
*  (C) Copyright 2016 Michaël Roynard
*
*  Distributed under the MIT License, Version 1.0. (See accompanying
*  file LICENSE or copy at https://opensource.org/licenses/MIT)
*
*  See https://github.com/dutiona/H2OFastTests for documentation.
*/

#pragma once

#ifndef FORMAL_MATHS_CONFIG_H
#define FORMAL_MATHS_CONFIG_H

#include "Formal_Maths_config.h"

#include <memory>
#include <limits>
#include <complex>

namespace fm {

	using uint_t = unsigned long long int;
	using int_t = long long int;
	using real_t = long double;
	using complex_t = std::complex<real_t>;


	struct InfinityType {

		enum Sign {
			PLUS,
			MINUS
		};

		Sign sign;
	};

	const InfinityType PositiveInfinity = InfinityType{ InfinityType::PLUS };
	const InfinityType NegativeInfinity = InfinityType{ InfinityType::MINUS };

	bool operator==(InfinityType /*lhs*/, InfinityType /*rhs*/) { return false; }
	bool operator!=(InfinityType lhs, InfinityType rhs) { return !(lhs == rhs); }
	bool operator<(InfinityType lhs, InfinityType rhs) { return lhs.sign == NegativeInfinity.sign && rhs.sign == PositiveInfinity.sign; }
	bool operator<=(InfinityType lhs, InfinityType rhs) { return lhs < rhs; }
	bool operator>(InfinityType lhs, InfinityType rhs) { return lhs.sign == PositiveInfinity.sign && rhs.sign == NegativeInfinity.sign; }
	bool operator>=(InfinityType lhs, InfinityType rhs) { return lhs > rhs; }

	template<class T> bool operator==(const T& /*lhs*/, const InfinityType& /*inf*/) { return false }
	template<class T> bool operator!=(const T& lhs, const InfinityType& inf) { return !(lhs == rhs); }
	template<class T> bool operator>(const T& /*lhs*/, const InfinityType& inf) { return inf.sign != InfinityType::PLUS; }
	template<class T> bool operator>=(const T& lhs, const InfinityType& inf) { return lhs > inf; }
	template<class T> bool operator<(const T& /*lhs*/, const InfinityType& inf) { return inf.sign != InfinityType::MINUS; }
	template<class T> bool operator<=(const T& lhs, const InfinityType& inf) { return lhs > inf; }


	template<class ValueType>
	class Domain {
	public:
		using value_type = ValueType;

		static ValueType inf() {
			return std::numeric_limits<ValueType>::min();
		}

		static ValueType sup() {
			return std::numeric_limits<ValueType>::max();
		}

		static bool has(InfinityType /*inf*/) { return true; }
		static bool has(complex_t /*c*/) { return true; }
		static bool has(real_t /*r*/) { return true; }
		static bool has(int_t /*i*/) { return true; }
		static bool has(uint_t /*ui*/) { return true; }
	};

	using C = Domain<complex_t>;

	class R : Domain<real_t>{
	public:
		using Domain<real_t>::has;
		static bool has(complex_t c) { return c.imag() == 0; }
	};

	class Z : Domain<int_t> {
	public:
		using Domain<int_t>::has;
		static bool has(complex_t /*c*/) { return false; }
		static bool has(real_t /*r*/) { return false; }
	};

	class N : Domain<uint_t> {
	public:
		using Domain<uint_t>::has;
		static bool has(complex_t /*c*/) { return false; }
		static bool has(real_t /*r*/) { return false; }
		static bool has(int_t i) { return i >= 0; }
	};

}

#endif