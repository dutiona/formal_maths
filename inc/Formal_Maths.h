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

#include <stdexcept>
#include <cassert>
#include <complex>
#include <functional>
#include <limits>
#include <memory>
#include <string>
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


	template<class ValueType0, class ValueType1, class ... ValueTypes> struct Select {
		using value_type = Select<ValueType0, Select<ValueType1, ValueTypes...>>;
	};

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

	class ComputeDeritativeError : public std::exception {
	public:
		ComputeDeritativeError(const std::string& message) :
			message_(message)
		{}

		ComputeDeritativeError& operator=(const ComputeDeritativeError&) = delete;

		virtual const char* what() const throw() {
			return message_.c_str();
		}
	private:
		const std::string message_;
	};

	class ComputePrimitiveError : public std::exception {
	public:
		ComputePrimitiveError(const std::string& message) :
			message_(message)
		{}

		ComputePrimitiveError& operator=(const ComputePrimitiveError&) = delete;

		virtual const char* what() const throw() {
			return message_.c_str();
		}
	private:
		const std::string message_;
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

		Computable(const Computable&) = delete;
		Computable& operator=(const Computable&) = delete;
		//Computable(Computable&&) = default;
		//Computable&& operator=(Computable&&) = default;

		virtual ~Computable() {}

		ValueType operator()(ValueType value) const { return computeValue(value); }
		ComputableHolder<ValueType> deritative() const { return computeDeritative(); }
		ComputableHolder<ValueType> primitive() const { return computePrimitive(); }
		virtual bool isConst() const { return false; }

	protected:

		Computable() {}

	private:

		virtual ValueType computeValue(ValueType) const = 0;
		virtual ComputableHolder<ValueType> computeDeritative() const = 0;
		virtual ComputableHolder<ValueType> computePrimitive() const = 0;

	};

	template<class ValueType>
	class ComputableHolder {
	public:

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
		bool isConst() const { return computable_->isConst(); }

	private:

		std::shared_ptr<Computable<ValueType>> computable_;

		template<class ValueTypeLeft, class ValueTypeRight>
		friend bool operator==(const ComputableHolder<ValueTypeLeft>& lhs, const ComputableHolder<ValueTypeRight>& rhs) {
			return *lhs.computable_ == *rhs.computable_;
		}
	};

	template<class Function, class ... Args>
	ComputableHolder<typename Function::value_type> make_computable(Args&& ... args) {
		return{ std::make_shared<Function>(std::forward<Args>(args)...) };
	}

	// Linear

	template<class ValueType>
	class LinearPow : public Computable<ValueType> {
	public:

		LinearPow() :
			const_(static_cast<ValueType>(1)),
			pow_(static_cast<ValueType>(1))
		{}

		LinearPow(ValueType c, ValueType pow) :
			const_(c),
			pow_(pow)
		{}

	private:

		virtual ValueType computeValue(ValueType x) const override {
			return const_ * std::pow(x, pow_);
		}

		virtual ComputableHolder<ValueType> computeDeritative() const override {
			if (pow_ == 0) {
				return make_computable<LinearPow<ValueType>>(static_cast<ValueType>(0), static_cast<ValueType>(0));
			}
			return make_computable<LinearPow<ValueType>>(const_ * (pow_), pow_ - 1);
		};

		virtual ComputableHolder<ValueType> computePrimitive() const override {
			if (pow_ == static_cast<ValueType>(-1)) {
				return make_computable<LinearPow<ValueType>>(const_, 0) * make_computable<Log<ValueType>>();
			}

			return make_computable<LinearPow<ValueType>>(const_ / (pow_ + 1), pow_ + 1);
		}

		virtual bool isConst() const {
			return pow_ == static_cast<ValueType>(0);
		}

		ValueType const_;
		ValueType pow_;

		template<class ValueTypeLeft, class ValueTypeRight>
		friend bool operator==(const LinearPow<ValueTypeLeft>& lhs, const LinearPow<ValueTypeRight>& rhs) {
			return lhs.const_ == rhs.const_ && lhs.pow_ == rhs.pow_;
		}
	};

	template<class ValueType>
	class Log : public Computable<ValueType> {
	public:

		Log() {}

	private:

		virtual ValueType computeValue(ValueType x) const override {
			return std::log(x);
		}

		virtual ComputableHolder<ValueType> computeDeritative() const override {
			return make_computable<LinearPow<ValueType>>(1, -1);
		};

		virtual ComputableHolder<ValueType> computePrimitive() const override {
			return make_computable<LinearPow<ValueType>>(1, 1) * make_computable<Log<ValueType>>() - make_computable<LinearPow<ValueType>>(1, 1);
		}

		virtual bool isConst() const {
			return false;
		}

		template<class ValueTypeLeft, class ValueTypeRight>
		friend bool operator==(const Log<ValueTypeLeft>& lhs, const Log<ValueTypeRight>& rhs) {
			return true;
		}
	};

	template<class ValueType>
	class Exp : public Computable<ValueType> {
	public:

		Exp() {}

	private:

		virtual ValueType computeValue(ValueType x) const override {
			return std::exp(x);
		}

		virtual ComputableHolder<ValueType> computeDeritative() const override {
			return make_computable<Exp<ValueType>>();
		};

		virtual ComputableHolder<ValueType> computePrimitive() const override {
			return make_computable<Exp<ValueType>>();
		}

		virtual bool isConst() const {
			return false;
		}

		template<class ValueTypeLeft, class ValueTypeRight>
		friend bool operator==(const Exp<ValueTypeLeft>& lhs, const Exp<ValueTypeRight>& rhs) {
			return true;
		}
	};

	template<class ValueType>
	class Sin : public Computable<ValueType> {
	public:

		Sin() {}

	private:

		virtual ValueType computeValue(ValueType x) const override {
			return std::sin(x);
		}

		virtual ComputableHolder<ValueType> computeDeritative() const override {
			return make_computable<Cos<ValueType>>();
		};

		virtual ComputableHolder<ValueType> computePrimitive() const override {
			return -make_computable<Cos<ValueType>>();
		}

		virtual bool isConst() const {
			return false;
		}

		template<class ValueTypeLeft, class ValueTypeRight>
		friend bool operator==(const Sin<ValueTypeLeft>& lhs, const Sin<ValueTypeRight>& rhs) {
			return true;
		}
	};

	template<class ValueType>
	class Cos : public Computable<ValueType> {
	public:

		Cos() {}

	private:

		virtual ValueType computeValue(ValueType x) const override {
			return std::cos(x);
		}

		virtual ComputableHolder<ValueType> computeDeritative() const override {
			return -make_computable<Sin<ValueType>>();
		};

		virtual ComputableHolder<ValueType> computePrimitive() const override {
			return make_computable<Sin<ValueType>>();
		}

		virtual bool isConst() const {
			return false;
		}

		template<class ValueTypeLeft, class ValueTypeRight>
		friend bool operator==(const Cos<ValueTypeLeft>& lhs, const Cos<ValueTypeRight>& rhs) {
			return true;
		}
	};

	template<class ValueType>
	class Tan : public Computable<ValueType> {
	public:

		Tan() {}

	private:

		virtual ValueType computeValue(ValueType x) const override {
			return std::tan(x);
		}

		virtual ComputableHolder<ValueType> computeDeritative() const override {
			return make_computable<LinearPow<ValueType>>(1, 0) + make_computable<Tan<ValueType>>() * make_computable<Tan<ValueType>>();
		};

		virtual ComputableHolder<ValueType> computePrimitive() const override {
			return -(make_computable<Log<ValueType>>() ^ make_computable<Cos<ValueType>>(1));
		}

		virtual bool isConst() const {
			return false;
		}

		template<class ValueTypeLeft, class ValueTypeRight>
		friend bool operator==(const Tan<ValueTypeLeft>& lhs, const Tan<ValueTypeRight>& rhs) {
			return true;
		}
	};

	// Addition

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

		virtual bool isConst() const {
			return func_left_.isConst() && func_right_.isConst();
		}

		ComputableHolder<ValueTypeLeft> func_left_;
		ComputableHolder<ValueTypeRight> func_right_;

		template<class ValueTypeLeftCompLeft, class ValueTypeLeftCompRight, class ValueTypeRightCompLeft, class ValueTypeRightCompRight>
		friend bool operator==(const Addition<ValueTypeLeftCompLeft, ValueTypeLeftCompRight>& lhs, const Addition<ValueTypeRightCompLeft, ValueTypeRightCompRight>& rhs) {
			return (lhs.func_left_ == rhs.func_left_ && lhs.func_right_ == rhs.func_right_) ||
				(lhs.func_left_ == rhs.func_right_ && lhs.func_right_ == rhs.func_left_);
		}
	};

	template<class ValueTypeLeft, class ValueTypeRight>
	auto operator+(const ComputableHolder<ValueTypeLeft>& lhs, const ComputableHolder<ValueTypeRight>& rhs)
		-> ComputableHolder<typename Select<ValueTypeLeft, ValueTypeRight>::value_type>
	{
		return make_computable<Addition<ValueTypeLeft, ValueTypeRight>>(lhs, rhs);
	}


	// Substraction

	template<class ValueTypeLeft, class ValueTypeRight>
	class Substraction : public Computable<typename Select<ValueTypeLeft, ValueTypeRight>::value_type> {
	public:

		using value_type = typename Select<ValueTypeLeft, ValueTypeRight>::value_type;

		Substraction(const ComputableHolder<ValueTypeLeft>& fl, const ComputableHolder<ValueTypeRight>& fr) :
			func_left_(fl),
			func_right_(fr)
		{}

	private:

		virtual value_type computeValue(value_type arg) const override {
			return func_left_(arg) - func_right_(arg);
		}

		virtual ComputableHolder<value_type> computeDeritative() const override {
			return func_left_.deritative() - func_right_.deritative();
		}

		virtual ComputableHolder<value_type> computePrimitive() const override {
			return func_left_.primitive() - func_right_.primitive();
		}

		virtual bool isConst() const {
			return func_left_.isConst() && func_right_.isConst();
		}

		ComputableHolder<ValueTypeLeft> func_left_;
		ComputableHolder<ValueTypeRight> func_right_;

		template<class ValueTypeLeftCompLeft, class ValueTypeLeftCompRight, class ValueTypeRightCompLeft, class ValueTypeRightCompRight>
		friend bool operator==(const Substraction<ValueTypeLeftCompLeft, ValueTypeLeftCompRight>& lhs, const Substraction<ValueTypeRightCompLeft, ValueTypeRightCompRight>& rhs) {
			return (lhs.func_left_ == rhs.func_left_ && lhs.func_right_ == rhs.func_right_) ||
				(lhs.func_left_ == rhs.func_right_ && lhs.func_right_ == rhs.func_left_);
		}
	};

	template<class ValueTypeLeft, class ValueTypeRight>
	auto operator-(const ComputableHolder<ValueTypeLeft>& lhs, const ComputableHolder<ValueTypeRight>& rhs)
		-> ComputableHolder<typename Select<ValueTypeLeft, ValueTypeRight>::value_type>
	{
		return make_computable<Substraction<ValueTypeLeft, ValueTypeRight>>(lhs, rhs);
	}

	template<class ValueType>
	ComputableHolder<ValueType> operator-(const ComputableHolder<ValueType>& rhs) {
		return make_computable<Substraction<ValueType, ValueType>>(make_computable<Const<ValueType>>(0), rhs);
	}


	// Product

	template<class ValueTypeLeft, class ValueTypeRight>
	class Product : public Computable<typename Select<ValueTypeLeft, ValueTypeRight>::value_type> {
	public:

		using value_type = typename Select<ValueTypeLeft, ValueTypeRight>::value_type;

		Product(const ComputableHolder<ValueTypeLeft>& fl, const ComputableHolder<ValueTypeRight>& fr) :
			func_left_(fl),
			func_right_(fr)
		{}

	private:

		virtual value_type computeValue(value_type arg) const override {
			return func_left_(arg) * func_right_(arg);
		}

		virtual ComputableHolder<value_type> computeDeritative() const override {
			return func_left_.deritative() * func_right_ + func_left_ * func_right_.deritative();
		}

		virtual ComputableHolder<value_type> computePrimitive() const override {
			if (func_left_.isConst()) {
				return func_left_ * func_right_.primitive();
			}
			else if (func_right_.isConst()){
				return func_left_.primitive() * func_right_;
			}
			else {
				throw ComputePrimitiveError{ "Primitive u*v can't be computed" };
			}
		}

		virtual bool isConst() const {
			return func_left_.isConst() && func_right_.isConst();
		}

		ComputableHolder<ValueTypeLeft> func_left_;
		ComputableHolder<ValueTypeRight> func_right_;

		template<class ValueTypeLeftCompLeft, class ValueTypeLeftCompRight, class ValueTypeRightCompLeft, class ValueTypeRightCompRight>
		friend bool operator==(const Product<ValueTypeLeftCompLeft, ValueTypeLeftCompRight>& lhs, const Product<ValueTypeRightCompLeft, ValueTypeRightCompRight>& rhs) {
			return (lhs.func_left_ == rhs.func_left_ && lhs.func_right_ == rhs.func_right_) ||
				(lhs.func_left_ == rhs.func_right_ && lhs.func_right_ == rhs.func_left_);
		}
	};

	template<class ValueTypeLeft, class ValueTypeRight>
	auto operator*(const ComputableHolder<ValueTypeLeft>& lhs, const ComputableHolder<ValueTypeRight>& rhs)
		-> ComputableHolder<typename Select<ValueTypeLeft, ValueTypeRight>::value_type>
	{
		return make_computable<Product<ValueTypeLeft, ValueTypeRight>>(lhs, rhs);
	}

	// Division

	template<class ValueTypeLeft, class ValueTypeRight>
	class Division : public Computable<typename Select<ValueTypeLeft, ValueTypeRight>::value_type> {
	public:

		using value_type = typename Select<ValueTypeLeft, ValueTypeRight>::value_type;

		Division(const ComputableHolder<ValueTypeLeft>& fl, const ComputableHolder<ValueTypeRight>& fr) :
			func_left_(fl),
			func_right_(fr)
		{}

	private:

		virtual value_type computeValue(value_type arg) const override {
			const auto denominator = func_right_(arg);
			assert(denominator != 0);
			return func_left_(arg) / denominator;
		}

		virtual ComputableHolder<value_type> computeDeritative() const override {
			return func_left_.deritative() * func_right_ / func_left_ * func_right_.deritative();
		}

		virtual ComputableHolder<value_type> computePrimitive() const override {
			if (func_right_.isConst()){
				return func_left_.primitive() * func_right_;
			}
			else {
				throw ComputePrimitiveError{ "Primitive u/v can't be computed" };
			}
		}

		virtual bool isConst() const {
			return func_left_.isConst() && func_right_.isConst();
		}

		ComputableHolder<ValueTypeLeft> func_left_;
		ComputableHolder<ValueTypeRight> func_right_;

		template<class ValueTypeLeftCompLeft, class ValueTypeLeftCompRight, class ValueTypeRightCompLeft, class ValueTypeRightCompRight>
		friend bool operator==(const Division<ValueTypeLeftCompLeft, ValueTypeLeftCompRight>& lhs, const Division<ValueTypeRightCompLeft, ValueTypeRightCompRight>& rhs) {
			return lhs.func_left_ == rhs.func_left_ && lhs.func_right_ == rhs.func_right_;
		}
	};

	template<class ValueTypeLeft, class ValueTypeRight>
	auto operator/(const ComputableHolder<ValueTypeLeft>& lhs, const ComputableHolder<ValueTypeRight>& rhs)
		-> ComputableHolder<typename Select<ValueTypeLeft, ValueTypeRight>::value_type>
	{
		return make_computable<Division<ValueTypeLeft, ValueTypeRight>>(lhs, rhs);
	}


	// Composition

	template<class ValueTypeLeft, class ValueTypeRight>
	class Composition : public Computable<typename Select<ValueTypeLeft, ValueTypeRight>::value_type> {
	public:

		using value_type = typename Select<ValueTypeLeft, ValueTypeRight>::value_type;

		Composition(const ComputableHolder<ValueTypeLeft>& fin, const ComputableHolder<ValueTypeRight>& fout) :
			func_out_(fl),
			func_in_(fr)
		{}

	private:

		virtual value_type computeValue(value_type arg) const override {
			return func_out_(func_in_(arg));
		}

		virtual ComputableHolder<value_type> computeDeritative() const override {
			return func_in_.deritative() * (func_out_.deritative() ^ func_in_);
		}

		virtual ComputableHolder<value_type> computePrimitive() const override {
			if (func_out_.isConst()) {
				return func_in_.primitive() * func_out_;
			}
			else {
				throw ComputePrimitiveError{ "Primitive u(v) can't be computed" };
			}
		}

		virtual bool isConst() const {
			return func_out_.isConst() && func_in_.isConst();
		}

		ComputableHolder<ValueTypeLeft> func_out_;
		ComputableHolder<ValueTypeRight> func_in_;

		template<class ValueTypeLeftCompLeft, class ValueTypeLeftCompRight, class ValueTypeRightCompLeft, class ValueTypeRightCompRight>
		friend bool operator==(const Composition<ValueTypeLeftCompLeft, ValueTypeLeftCompRight>& lhs, const Composition<ValueTypeRightCompLeft, ValueTypeRightCompRight>& rhs) {
			return lhs.func_out_ == rhs.func_out_ && lhs.func_in_ == rhs.func_in_;
		}
	};

	template<class ValueTypeLeft, class ValueTypeRight>
	auto operator^(const ComputableHolder<ValueTypeLeft>& lhs, const ComputableHolder<ValueTypeRight>& rhs)
		-> ComputableHolder<typename Select<ValueTypeLeft, ValueTypeRight>::value_type>
	{
		return make_computable<Composition<ValueTypeLeft, ValueTypeRight>>(lhs, rhs);
	}


	/*
	template<class ValueType> struct ACos { };
	template<class ValueType> struct Asin { };
	template<class ValueType> struct Atan { };
	template<class ValueType> struct Tan2 { };
	template<class ValueType> struct Atan2 { };
	*/

}

#endif