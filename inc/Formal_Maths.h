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

#include <cassert>
#include <complex>
#include <functional>
#include <limits>
#include <memory>
#include <type_traits>

#ifndef UNUSED
#if defined(__GNUC__)
#	define UNUSED(x) UNUSED_ ## x __attribute__((unused))
#elif defined(__LCLINT__)
#	define UNUSED(x) /*@unused@*/ x
#elif defined(__cplusplus)
#	define UNUSED(x)
#else
#	define UNUSED(x) x
#endif
#endif

namespace fm {

	using uint_t = unsigned long long int;
	using int_t = long long int;
	using real_t = double;
	using complex_t = std::complex<real_t>;


	template<class ValueTypeLeft, class ValueTypeRight> struct Select;

	template<> struct Select<uint_t, uint_t> { using value_type = uint_t; };
	template<> struct Select<uint_t, int_t> { using value_type = int_t; };
	template<> struct Select<uint_t, real_t> { using value_type = real_t; };
	template<> struct Select<uint_t, complex_t> { using value_type = complex_t; };

	template<> struct Select<int_t, uint_t> { using value_type = int_t; };
	template<> struct Select<int_t, int_t> { using value_type = int_t; };
	template<> struct Select<int_t, real_t> { using value_type = real_t; };
	template<> struct Select<int_t, complex_t> { using value_type = complex_t; };

	template<> struct Select<real_t, uint_t> { using value_type = real_t; };
	template<> struct Select<real_t, int_t> { using value_type = real_t; };
	template<> struct Select<real_t, real_t> { using value_type = real_t; };
	template<> struct Select<real_t, complex_t> { using value_type = complex_t; };

	template<> struct Select<complex_t, uint_t> { using value_type = complex_t; };
	template<> struct Select<complex_t, int_t> { using value_type = complex_t; };
	template<> struct Select<complex_t, real_t> { using value_type = complex_t; };
	template<> struct Select<complex_t, complex_t> { using value_type = complex_t; };

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

	// General interface
	template<class ValueType>
	class ComputableHolder;

	template<class ValueType>
	class Computable {
	public:

		using value_type = typename ValueType;

		Computable() {}
		Computable(const Computable&) = delete;
		Computable& operator=(const Computable&) = delete;
		//Computable(Computable&&) = default;
		//Computable&& operator=(Computable&&) = default;

		virtual ~Computable() {}

		ValueType operator()(ValueType value) const { return computeValue(value); }
		ComputableHolder<ValueType> deritative() const { return computeDeritative(); }
		ComputableHolder<ValueType> primitive() const { return computePrimitive(); }

	private:

		virtual ValueType computeValue(ValueType) const = 0;
		virtual ComputableHolder<ValueType> computeDeritative() const = 0;
		virtual ComputableHolder<ValueType> computePrimitive() const = 0;

	};

	template<class ValueType>
	class ComputableHolder {
	public:

		ComputableHolder(){} //TODELETE
		ComputableHolder(const std::shared_ptr<Computable<ValueType>>& computable) :
			computable_(computable)
		{}

		ComputableHolder(const ComputableHolder& comp) : computable_(comp.computable_) {}
		ComputableHolder& operator=(const ComputableHolder& comp) {
			computable_ = comp.computable_;
			return *this;
		}

		ValueType operator()(ValueType value) const { return (*computable_)(value); }
		ComputableHolder<ValueType> deritative() const { return computable_->deritative(); }
		ComputableHolder<ValueType> primitive() const { return computable_->primitive(); }

		Computable<ValueType>& operator*() { return computable_; }
		Computable<ValueType>* operator&() { return computable_.get(); }
		Computable<ValueType>* operator->() { return computable_.get(); }

	private:

		std::shared_ptr<Computable<ValueType>> computable_;

	};

	template<class Function, class ... Args>
	ComputableHolder<typename Function::value_type> make_computable(Args&& ... args) {
		return{ std::make_shared<Function>(std::forward<Args>(args)...) };
	}


	//Constant definition

	template<class ValueType>
	class Const : public Computable<ValueType> {
	public:

		Const(ValueType value) :
			value_(value)
		{}

	private:

		virtual ValueType computeValue(ValueType) const override {
			return value_;
		}

		virtual ComputableHolder<ValueType> computeDeritative() const override {
			return make_computable<Const<ValueType>>(static_cast<ValueType>(0));
		}

		virtual ComputableHolder<ValueType> computePrimitive() const override {
			return make_computable<Const<ValueType>>(value_);
		};

		ValueType value_;

	};

	// Linear

	template<class ValueType>
	class Linear : public Computable<ValueType> {
	public:

		Linear(ValueType value) :
			value_(value)
		{}

	private:

		virtual ValueType computeValue(ValueType x) const override {
			return value_ * x;
		}

		virtual ComputableHolder<ValueType> computeDeritative() const override {
			return make_computable<Const<ValueType>>(value_);
		};

		virtual ComputableHolder<ValueType> computePrimitive() const override {
			return{};// Const<real_t>{ 0.5 } *Linear<real_t>{ 1. } *Linear<real_t>{1.};
		}

		ValueType value_;

	};


	template<class ValueTypeLeft, class ValueTypeRight>
	class Addition : public Computable<typename Select<ValueTypeLeft, ValueTypeRight>::value_type> {
	public:

		using value_type = typename Select<ValueTypeLeft, ValueTypeRight>::value_type;

		Addition(const ComputableHolder<ValueTypeLeft>& fl, const ComputableHolder<ValueTypeRight>& fr) :
			func_left_(fl),
			func_right_(fr)
		{}

	private:

		virtual value_type computeValue(value_type arg) const override {
			return func_left_(arg) + func_right_(arg);
		}

		virtual ComputableHolder<value_type> computeDeritative() const override {
			return func_left_.deritative() + func_right_.deritative();
		}

		virtual ComputableHolder<value_type> computePrimitive() const override {
			return func_left_.primitive() + func_right_.primitive();
		}

		ComputableHolder<ValueTypeLeft> func_left_;
		ComputableHolder<ValueTypeRight> func_right_;

	};

	template<class ValueTypeLeft, class ValueTypeRight>
	auto operator+(const ComputableHolder<ValueTypeLeft>& lhs, const ComputableHolder<ValueTypeRight>& rhs)
		-> ComputableHolder<typename Select<ValueTypeLeft, ValueTypeRight>::value_type>
	{
		return make_computable<Addition<ValueTypeLeft, ValueTypeRight>>(lhs, rhs);
	};





	/*
	// Power

	template<class ValueType>
	class Pow : public Computable<ValueType> {
	public:

	Pow(ValueType value) :
	value_(value)
	{}

	template<class ... Args>
	ValueType compute(Args&& ... args) const;

	ValueType compute(ValueType x) const {
	using std::pow;
	return pow(x, value_);
	}

	Computable<ValueType> computeDeritative() const {
	return make_func<Const, ValueType>(value_) * make_func<Pow, ValueType>(value_ - 1);
	}

	Computable<ValueType> computePrimitive() const {
	return make_func<Pow, ValueType>(value_ + 1) * make_func<Const, ValueType>(1 / (value_ + 1));
	}

	private:

	const ValueType value_;

	};


	template<class ValueType> struct FracPow { ValueType fracPow; };
	template<class ValueType> struct Invert { };
	template<class ValueType> struct SquareRoot{ };
	template<class ValueType> struct Log { };
	template<class ValueType> struct Exp { };
	template<class ValueType> struct Cos { };
	template<class ValueType> struct ACos { };
	template<class ValueType> struct Sin { };
	template<class ValueType> struct Asin { };
	template<class ValueType> struct Tan { };
	template<class ValueType> struct Atan { };
	template<class ValueType> struct Tan2 { };
	template<class ValueType> struct Atan2 { };
	*/







	/*
	// Substraction

	template<template <class> class FunctionBaseLeft, template <class> class FunctionBaseRight, class ValueType>
	class Substraction : Computable<ValueType> {
	public:

	Substraction(FunctionBaseLeft<ValueType> fl, FunctionBaseRight<ValueType> fr) :
	func_left_(fl),
	func_right_(fr)
	{}

	template<class ... Args>
	ValueType operator()(Args&& ... args) const {
	return func_left_.compute(std::forward<Args>(args)...) - func_right_.compute(std::forward<Args>(args)...);
	}

	Function<ValueType> computeDeritative() const {
	return func_left_.computeDeritative() - func_right_.computeDeritative();
	}

	Function<ValueType> computePrimitive() const {
	return func_left_.computePrimitive() - func_right_.computePrimitive();
	}

	private:

	const FunctionBaseLeft<ValueType> func_left_;
	const FunctionBaseRight<ValueType> func_right_;

	};

	template<template <class> class FunctionBaseLeft, template <class> class FunctionBaseRight, class ValueType>
	Computable<ValueType> operator-(FunctionBaseLeft<ValueType> lhs, FunctionBaseRight<ValueType> rhs) {
	return{ Substraction<FunctionBaseLeft, FunctionBaseRight, ValueType>{ lhs, rhs } };
	}

	template<template <class> class FunctionBaseLeft, class ValueType>
	Computable<ValueType> operator-(FunctionBaseLeft<ValueType> lhs) {
	return{ Substraction<FunctionBaseLeft, Const, ValueType>{ make_func<Const, ValueType>(0), lhs } };
	}

	// Product

	template<template <class> class FunctionBaseLeft, template <class> class FunctionBaseRight, class ValueType>
	class Product : Computable<ValueType> {
	public:

	Product(FunctionBaseLeft<ValueType> fl, FunctionBaseRight<ValueType> fr) :
	func_left_(fl),
	func_right_(fr)
	{}

	template<class ... Args>
	ValueType operator()(Args&& ... args) const {
	return func_left_.compute(std::forward<Args>(args)...) * func_right_.compute(std::forward<Args>(args)...);
	}

	Function<ValueType> computeDeritative() const {
	return func_left_.computeDeritative() * func_right_ + func_left_ * func_right_.computeDeritative();
	}

	Function<ValueType> computePrimitive() const {
	static_assert(false, "Primitive(u*v) cannot be calculated :(");
	}

	private:

	const FunctionBaseLeft<ValueType> func_left_;
	const FunctionBaseRight<ValueType> func_right_;
	};

	template<template <class> class FunctionBaseLeft, template <class> class FunctionBaseRight, class ValueType>
	Computable<ValueType> operator*(FunctionBaseLeft<ValueType> lhs, FunctionBaseRight<ValueType> rhs) {
	return{ Product<FunctionBaseLeft, FunctionBaseRight, ValueType>{ lhs, rhs } };
	}
	*/

	// Division

	/*
	template<template <class> class FunctionBaseLeft, template <class> class FunctionBaseRight, class ValueType>
	class Division : Computable<ValueType> {
	public:

	Division(FunctionBaseLeft<ValueTypeLeft> fl, FunctionBaseRight<ValueTypeRight> fr) :
	func_left_(fl),
	func_right_(fr)
	{}

	template<class ... Args>
	ValueType operator()(Args&& ... args) const {
	const auto denominator = func_right_.compute(std::forward<Args>(args)...);
	assert(denomiator != 0);
	return func_left_.compute(std::forward<Args>(args)...) / denominator;
	}

	Function<ValueType> computeDeritative() const {
	return (func_left_.computeDeritative() * func_right_ - func_left_ * func_right_.computeDeritative()) * (1 / func_right_ * func_right_);
	}

	Function<ValueType> computePrimitive() const {
	static_assert(false, "Primitive(u*v) cannot be calculated :(");
	}

	private:

	const FunctionBaseLeft<ValueTypeLeft> func_left_;
	const FunctionBaseRight<ValueTypeRight> func_right_;
	};

	template<template <class> class FunctionBaseLeft, template <class> class FunctionBaseRight, class ValueType>
	Computable<ValueType> operator*(FunctionBaseLeft<ValueType> lhs, FunctionBaseRight<ValueType> rhs) {
	return{ Product<FunctionBaseLeft, FunctionBaseRight, ValueType>{ lhs, rhs } };
	}

	*/











	/*

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

	*/


}

#endif