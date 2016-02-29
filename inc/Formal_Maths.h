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

#include <complex>
#include <functional>
#include <limits>
#include <memory>

#ifdef UNUSED
#elif defined(__GNUC__)
# define UNUSED(x) UNUSED_ ## x __attribute__((unused))
#elif defined(__LCLINT__)
# define UNUSED(x) /*@unused@*/ x
#elif defined(__cplusplus)
# define UNUSED(x)
#else
# define UNUSED(x) x
#endif

namespace fm {

	using uint_t = unsigned long long int;
	using int_t = long long int;
	using real_t = double;
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

	bool operator==(InfinityType UNUSED(lhs), InfinityType UNUSED(rhs)) { return false; }
	bool operator!=(InfinityType lhs, InfinityType rhs) { return !(lhs == rhs); }
	bool operator<(InfinityType lhs, InfinityType rhs) { return lhs.sign == NegativeInfinity.sign && rhs.sign == PositiveInfinity.sign; }
	bool operator<=(InfinityType lhs, InfinityType rhs) { return lhs < rhs; }
	bool operator>(InfinityType lhs, InfinityType rhs) { return lhs.sign == PositiveInfinity.sign && rhs.sign == NegativeInfinity.sign; }
	bool operator>=(InfinityType lhs, InfinityType rhs) { return lhs > rhs; }

	template<class T> bool operator==(const T& UNUSED(lhs), const InfinityType& UNUSED(rhs)) { return false }
	template<class T> bool operator!=(const T& lhs, const InfinityType& inf) { return !(lhs == rhs); }
	template<class T> bool operator>(const T& UNUSED(lhs), const InfinityType& inf) { return inf.sign != InfinityType::PLUS; }
	template<class T> bool operator>=(const T& lhs, const InfinityType& inf) { return lhs > inf; }
	template<class T> bool operator<(const T& UNUSED(lhs), const InfinityType& inf) { return inf.sign != InfinityType::MINUS; }
	template<class T> bool operator<=(const T& lhs, const InfinityType& inf) { return lhs > inf; }

	/*
	* Domaine
	*/

	template<class ValueType>
	struct Domain {
		using value_type = ValueType;

		static ValueType inf() {
			return std::numeric_limits<ValueType>::min();
		}

		static ValueType sup() {
			return std::numeric_limits<ValueType>::max();
		}

		static bool has(InfinityType UNUSED(inf)) { return true; }
		static bool has(complex_t UNUSED(c)) { return true; }
		static bool has(real_t UNUSED(r)) { return true; }
		static bool has(int_t UNUSED(i)) { return true; }
		static bool has(uint_t UNUSED(ui)) { return true; }
	};

	using C = Domain<complex_t>;

	struct R : Domain<real_t>{
		using Domain<real_t>::has;
		static bool has(complex_t c) { return c.imag() == 0; }
	};

	struct Z : Domain<int_t> {
		using Domain<int_t>::has;
		static bool has(complex_t UNUSED(c)) { return false; }
		static bool has(real_t UNUSED(r)) { return false; }
	};

	struct N : Domain<uint_t> {
		using Domain<uint_t>::has;
		static bool has(complex_t UNUSED(c)) { return false; }
		static bool has(real_t UNUSED(r)) { return false; }
		static bool has(int_t i) { return i >= 0; }
	};


	/*
	* Usual functions impl
	*/

	// Trivial impl
	template<class ValueType> struct Const { ValueType value; };
	template<class ValueType> struct Pow { ValueType pow; };
	template<class ValueType> struct Linear { };
	template<class ValueType> struct FracPow { ValueType fracPow; };
	template<class ValueType> struct Invert { };
	template<class ValueType> struct SquareRoot{ };
	template<class ValueType> struct Cos { };
	template<class ValueType> struct Sin { };

	//Func exec
	template<class ValueType> ValueType exec_func(const Const<ValueType>& f, ValueType x) { return f.value; }
	template<class ValueType> ValueType exec_func(const Linear<ValueType>& UNUSED(f), ValueType x) { return x; }
	template<class ValueType> ValueType exec_func(const Pow<ValueType>& f, ValueType x) { return std::pow(x, f.pow); }
	template<class ValueType> ValueType exec_func(const SquareRoot<ValueType>& UNUSED(f), ValueType x) { return std::sqrt(x); }
	template<class ValueType> ValueType exec_func(const Invert<ValueType>& UNUSED(f), ValueType x) { return 1 / x; }
	template<class ValueType> ValueType exec_func(const FracPow<ValueType>& f, ValueType x) { return 1 / std::pow(x, f.fracPow); }
	template<class ValueType> ValueType exec_func(const Cos<ValueType>& UNUSED(f), ValueType x) { return std::cos(x); }
	template<class ValueType> ValueType exec_func(const Sin<ValueType>& UNUSED(f), ValueType x) { return std::sin(x); }

	//Deritatives
	template<class ValueType> ValueType exec_deritative(const Const<ValueType>& UNUSED(f), ValueType x) { return 0; }
	template<class ValueType> ValueType exec_deritative(const Linear<ValueType>& UNUSED(f), ValueType x) { return 1; }
	template<class ValueType> ValueType exec_deritative(const SquareRoot<ValueType>& UNUSED(f), ValueType x) { return 1 / (2 * std::sqrt(x)); }
	template<class ValueType> ValueType exec_deritative(const Pow<ValueType>& f, ValueType x) { return f.pow * std::pow(x, f.pow - 1); }
	template<class ValueType> ValueType exec_deritative(const Invert<ValueType>& UNUSED(f), ValueType x) { return -1 / (x * x); }
	template<class ValueType> ValueType exec_deritative(const FracPow<ValueType>& f, ValueType x) { return -f.fracPow / std::pow(x, f.fracPow + 1); }
	template<class ValueType> ValueType exec_deritative(const Cos<ValueType>& UNUSED(f), ValueType x) { return -std::sin(x); }
	template<class ValueType> ValueType exec_deritative(const Sin<ValueType>& UNUSED(f), ValueType x) { return std::cos(x); }

	//Primitives
	template<class ValueType> ValueType exec_primitive(const Const<ValueType>& f, ValueType x) { return f.value * x; }
	template<class ValueType> ValueType exec_primitive(const Linear<ValueType>& UNUSED(f), ValueType x) { return x * x * 0.5; }
	template<class ValueType> ValueType exec_primitive(const SquareRoot<ValueType>& UNUSED(f), ValueType x) { return (2 / 3) * std::pow(x, 1.5); }
	template<class ValueType> ValueType exec_primitive(const Pow<ValueType>& f, ValueType x) { return std::pow(x, f.pow); }
	template<class ValueType> ValueType exec_primitive(const Invert<ValueType>& UNUSED(f), ValueType x) { return std::log(x); }
	template<class ValueType> ValueType exec_primitive(const FracPow<ValueType>& f, ValueType x) { return -1 / ((f.fracPow - 1) * std::pow(x, f.fracPow - 1)); }
	template<class ValueType> ValueType exec_primitive(const Cos<ValueType>& UNUSED(f), ValueType x) { return std::sin(x); }
	template<class ValueType> ValueType exec_primitive(const Sin<ValueType>& UNUSED(f), ValueType x) { return -std::cos(x); }


	// Definition for composed func impl
	template<template <class> class FunctionBase, class ValueType>
	struct Function {
		Function(FunctionBase<ValueType> f) :
			func(f)
		{}

		FunctionBase<ValueType> func;
	};

	template<template <class> class FunctionBase, class ValueType, class ... Args>
	Function<FunctionBase, ValueType> make_func(Args&& ... args) {
		return{ FunctionBase<ValueType>{std::forward<Args>(args)...} };
	}

	template<template <class> class FunctionBase, class ValueType>
	ValueType exec_func(const Function<FunctionBase, ValueType>& f, ValueType x) { return exec_func(f.func, x); }

	template<template <class> class FunctionBase, class ValueType>
	ValueType exec_deritative(const Function<FunctionBase, ValueType>& f, ValueType x) { return exec_deritative(f.func, x); }

	template<template <class> class FunctionBase, class ValueType>
	ValueType exec_primitive(const Function<FunctionBase, ValueType>& f, ValueType x) { return exec_primitive(f.func, x); }

	// Product
	template<template <class> class FunctionBaseLeft, template <class> class FunctionBaseRight, class ValueType>
	struct FunctionProduct {
		FunctionProduct(Function<FunctionBaseLeft, ValueType> fl, Function<FunctionBaseRight, ValueType> fr) :
			funcLeft(fl),
			funcRight(fr)
		{}

		Function<FunctionBaseLeft, ValueType> funcLeft;
		Function<FunctionBaseRight, ValueType> funcRight;
	};

	template<template <class> class FunctionBaseLeft, template <class> class FunctionBaseRight, class ValueType>
	auto operator*(Function<FunctionBaseLeft, ValueType> lhs, Function<FunctionBaseRight, ValueType> rhs)
		-> FunctionProduct<FunctionBaseLeft, FunctionBaseRight, ValueType> {
		return{ lhs, rhs };
	}
	
	template<template <class> class FunctionBaseLeft, template <class> class FunctionBaseRight, class ValueType>
	ValueType exec_func(const FunctionProduct<FunctionBaseLeft, FunctionBaseRight, ValueType>& f, ValueType x) {
		return exec_func(f.funcLeft, x) * exec_func(f.funcRight, x);
	}

	template<template <class> class FunctionBaseLeft, template <class> class FunctionBaseRight, class ValueType>
	ValueType exec_deritative(const FunctionProduct<FunctionBaseLeft, FunctionBaseRight, ValueType>& f, ValueType x) {
		return exec_deritative(f.funcLeft, x) * exec_func(f.funcRight, x) + exec_func(f.funcLeft, x) * exec_deritative(f.funcRight, x);
	}

	template<template <class> class FunctionBaseLeft, template <class> class FunctionBaseRight, class ValueType>
	ValueType exec_primitive(const FunctionProduct<FunctionBaseLeft, FunctionBaseRight, ValueType>& f, ValueType x) {
		static_assert(false, "Primitive(u*v) cannot be calculated :(");
	}
	



}

#endif