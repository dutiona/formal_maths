#include "Formal_Maths.h"

#include <iostream>
#include <limits>

#pragma warning ( push )
#pragma warning ( disable: 4189 )

using namespace fm;

int main(int /*argc*/, char* /*argv[]*/){

	const auto f1 = ComputableHolder<real_t>{ std::move(std::make_unique<Const<real_t>>(15.)) };
	const auto f2 = ComputableHolder<uint_t>{ std::move(std::make_unique<Linear<uint_t>>(40)) };
	const auto f3 = f1 + f2;
	const auto deriv = f3.deritative();
	const auto prim = f3.primitive();
	const auto a = f3(1.);
	const auto b = deriv(15.);
	const auto c = prim(15.);

	/*
	const auto aa = f1;
	const auto bb = exec_deritative(f1, 13.);
	const auto cc = exec_primitive(f1, 5.);

	const auto f2 = make_func<Sin, real_t>();
	const auto f3 = f1 * f2;

	const auto aaa = exec_func(f3, 40.);
	const auto bbb = exec_deritative(f3, 40.);
	//const auto ccc = exec_primitive(f3, 40.);

	*/
	std::cout << "Press enter to continue...";
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	return EXIT_SUCCESS;
}

#pragma warning ( pop )