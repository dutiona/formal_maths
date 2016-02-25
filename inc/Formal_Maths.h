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

namespace fm {

	enum Sign { PLUS, MINUS };

	struct Infinity { Sign sign; };
	bool operator==(Infinity lhs, Infinity rhs) {
		return lhs.sign == rhs.sign;
	}

	struct real_t {
		Infinity inf;
		nullptr_t null;
		long double val;

		real_t() : null(nullptr) {}
		real_t(Infinity inf) : inf(inf) {}
		real_t(long double val) : val(val) {}
	};

	bool operator==(const real_t& lhs, const real_t& rhs) {
		return lhs.inf == rhs.inf || lhs.val == rhs.val;
	}

	struct int_t {
		Infinity inf;
		nullptr_t null;
		long long int val;

		int_t() : null(nullptr) {}
		int_t(Infinity inf) : inf(inf) {}
		int_t(long long int val) : val(val) {}
	};

	bool operator>=(const int_t& i, int_t nb) {
		return !i.null && i.val >= nb;
	}

	struct uint_t {
		Infinity inf;
		nullptr_t null;
		unsigned long long int val;

		uint_t() : null(nullptr) {}
		uint_t(Infinity inf) : inf(inf) {}
		uint_t(unsigned long long int val) : val(val) {}
	};

	template<class Derived>
	struct Domain_traits;

	template<template <typename T> class Derived, typename T>
	struct Domain {

		using value_type = typename T;

		static value_type getPlusInf() {
			value_type r;
			r.inf.sign = PLUS;
			return r;
		}
		
		static value_type getMinusInf() {
			value_type r;
			r.inf.sign = MINUS;
			return r;
		}
		
		static bool has(real_t r) { return Derived<T>::has(r); }
		static bool has(int_t i) { return Derived<T>::has(i); }
		static bool has(uint_t ui) { return Derived<T>::has(ui); }
	};

	template<class ValueType>
	struct Domain_R_Private : Domain<Domain_R_Private, ValueType> {

		using value_type = typename ValueType;

		static bool has(real_t /*r*/) { return true; }
		static bool has(int_t /*i*/) { return true; }
		static bool has(uint_t /*ui*/) { return true; }
	};

	template<class ValueType>
	struct Domain_N_Private : Domain<Domain_N_Private, ValueType> {

		using value_type = typename ValueType;

		static bool has(real_t /*r*/) { return false; }
		static bool has(int_t i) { return i >= 0; }
		static bool has(uint_t ui) { return true; }
	};

	using Domain_R = Domain_R_Private<real_t>;
	using Domain_N = Domain_N_Private<uint_t>;

}

#endif