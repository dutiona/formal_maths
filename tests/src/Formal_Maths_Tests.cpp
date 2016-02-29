#include "Formal_Maths.h"

#include <iostream>
#include <limits>

#pragma warning ( push )
#pragma warning ( disable: 4189 )

using namespace fm;

int main(int /*argc*/, char* /*argv[]*/){

	const auto f1 = make_func<Const, real_t>(15.);
	const auto aa = exec_func(f1, 1155555.);
	const auto bb = exec_deritative(f1, 13.);
	const auto cc = exec_primitive(f1, 5.);

	const auto f2 = make_func<Sin, real_t>();
	const auto f3 = f1 * f2;

	const auto aaa = exec_func(f3, 40.);
	const auto bbb = exec_deritative(f3, 40.);
	//const auto ccc = exec_primitive(f3, 40.);


	std::cout << "Press enter to continue...";
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	return EXIT_SUCCESS;
}

#pragma warning ( pop )