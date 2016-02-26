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

		static bool has(InfinityType /*inf*/) { return true; }
		static bool has(complex_t /*c*/) { return true; }
		static bool has(real_t /*r*/) { return true; }
		static bool has(int_t /*i*/) { return true; }
		static bool has(uint_t /*ui*/) { return true; }
	};

	using C = Domain<complex_t>;

	struct R : Domain<real_t>{
		using Domain<real_t>::has;
		static bool has(complex_t c) { return c.imag() == 0; }
	};

	struct Z : Domain<int_t> {
		using Domain<int_t>::has;
		static bool has(complex_t /*c*/) { return false; }
		static bool has(real_t /*r*/) { return false; }
	};

	struct N : Domain<uint_t> {
		using Domain<uint_t>::has;
		static bool has(complex_t /*c*/) { return false; }
		static bool has(real_t /*r*/) { return false; }
		static bool has(int_t i) { return i >= 0; }
	};


	/*
	* Usual functions impl
	*/

	template<class Func, class Ret, class ... Args> Ret exec_func(Func, Args ... args);
	template<class Func, class Ret, class ... Args> Ret exec_deritative(Func, Args ... args);
	template<class Func, class Ret, class ... Args> Ret exec_primitive(Func, Args ... args);

	template<class ValueType>
	struct Const { ValueType value; };

	template<class ValueType>
	struct Pow { ValueType pow; };

	template<class ValueType>
	struct Linear { };

	template<class ValueType>
	struct FracPow { ValueType fracPow; };

	template<class ValueType>
	struct Invert { };

	template<class ValueType>
	struct SquareRoot{ };

	template<class ValueType>
	struct Cos { };

	template<class ValueType>
	struct Sin { };
	
	template<template <class> class FunctionBase, class ValueType>
	struct Function {
		Function(FunctionBase<ValueType> f) : func(f) {}
		FunctionBase<ValueType> func;
	};

	//Func exec

	template<class ValueType>
	ValueType exec_func(Const<ValueType> f, ValueType x) { return f.value; }

	template<class ValueType>
	ValueType exec_func(Linear<ValueType> f, ValueType x) { return x; }

	template<class ValueType>
	ValueType exec_func(Pow<ValueType> f, ValueType x) { return std::pow(x, f.pow); }

	template<class ValueType>
	ValueType exec_func(SquareRoot<ValueType> f, ValueType x) { return std::sqrt(x); }

	template<class ValueType>
	ValueType exec_func(Invert<ValueType> f, ValueType x) { return 1 / x; }

	template<class ValueType>
	ValueType exec_func(FracPow<ValueType> f, ValueType x) { return 1 / std::pow(x, f.fracPow); }

	template<class ValueType>
	ValueType exec_func(Cos<ValueType> f, ValueType x) { return std::cos(x); }

	template<class ValueType>
	ValueType exec_func(Sin<ValueType> f, ValueType x) { return std::sin(x); }

	template<template <class> class FunctionBase, class ValueType>
	ValueType exec_func(Function<FunctionBase, ValueType> f, ValueType x) { return exec_func(f.func, x); }

	//Deritatives

	template<class ValueType>
	ValueType exec_deritative(Const<ValueType> f, ValueType x) { return 0; }

	template<class ValueType>
	ValueType exec_deritative(Linear<ValueType> f, ValueType x) { return 1; }

	template<class ValueType>
	ValueType exec_deritative(SquareRoot<ValueType> f, ValueType x) { return 1 / (2 * std::sqrt(x)); }

	template<class ValueType>
	ValueType exec_deritative(Pow<ValueType> f, ValueType x) { return f.pow * std::pow(x, f.pow - 1); }

	template<class ValueType>
	ValueType exec_deritative(Invert<ValueType> f, ValueType x) { return -1 / (x * x); }

	template<class ValueType>
	ValueType exec_deritative(FracPow<ValueType> f, ValueType x) { return -f.fracPow / std::pow(x, f.fracPow + 1); }

	template<class ValueType>
	ValueType exec_deritative(Cos<ValueType> f, ValueType x) { return -std::sin(x); }

	template<class ValueType>
	ValueType exec_deritative(Sin<ValueType> f, ValueType x) { return std::cos(x); }

	template<template <class> class FunctionBase, class ValueType>
	ValueType exec_deritative(Function<FunctionBase, ValueType> f, ValueType x) { return exec_deritative(f.func, x); }

	//Primitives

	template<class ValueType>
	ValueType exec_primitive(Const<ValueType> f, ValueType x) { return f.value * x; }

	template<class ValueType>
	ValueType exec_primitive(Linear<ValueType> f, ValueType x) { return x * x * 0.5; }

	template<class ValueType>
	ValueType exec_primitive(SquareRoot<ValueType> f, ValueType x) { return (2 / 3) * std::pow(x, 1.5); }

	template<class ValueType>
	ValueType exec_primitive(Pow<ValueType> f, ValueType x) { return std::pow(x, f.pow); }

	template<class ValueType>
	ValueType exec_primitive(Invert<ValueType> f, ValueType x) { return std::log(x); }

	template<class ValueType>
	ValueType exec_primitive(FracPow<ValueType> f, ValueType x) { return -1 / ((f.fracPow - 1) * std::pow(x, f.fracPow - 1)); }

	template<class ValueType>
	ValueType exec_primitive(Cos<ValueType> f, ValueType x) { return std::sin(x); }

	template<class ValueType>
	ValueType exec_primitive(Sin<ValueType> f, ValueType x) { return -std::cos(x); }

	template<template <class> class FunctionBase, class ValueType>
	ValueType exec_primitive(Function<FunctionBase, ValueType> f, ValueType x) { return exec_primitive(f.func, x); }

}

#endif