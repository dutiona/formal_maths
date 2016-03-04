#include "Formal_Maths.h"

#include <iostream>
#include <limits>

#pragma warning ( push )
#pragma warning ( disable: 4189 )

using namespace fm;

int main(int /*argc*/, char* /*argv[]*/){

	const auto f1 = make_computable<Const<real_t>>(15.);
	const auto f2 = make_computable<Linear<real_t>>(40.);
	const auto f3 = f1 + f2;
	const auto deriv_f1 = f1.deritative();
	const auto deriv_f2 = f2.deritative();
	const auto deriv_f3 = f3.deritative();
	const auto prim_f1 = f1.primitive();
	const auto prim_f2 = f2.primitive();
	const auto prim_f3 = f3.primitive();

	const auto a = f1(1.);
	const auto b = deriv_f1(15.);
	const auto c = prim_f1(15.);

	const auto aa = f2(1.);
	const auto bb = deriv_f2(15.);
	const auto cc = prim_f2(15.);

	const auto aaa = f3(1.);
	const auto bbb = deriv_f3(15.);
	const auto ccc = prim_f3(15.);

	const auto f4 = f1 * f2;
	const auto defiv_f4 = f4.deritative();
	const auto prim_f4 = f4.primitive();

	const auto f5 = f2 / f1;
	const auto defiv_f5 = f5.deritative();
	const auto prim_f5 = f5.primitive();

	std::cout << "Press enter to continue...";
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	return EXIT_SUCCESS;
}

#pragma warning ( pop )