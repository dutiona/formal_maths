#include "Formal_Maths.h"

#include <iostream>
#include <limits>


using namespace fm;

int main(int /*argc*/, char* /*argv[]*/){

	complex_t c = { 1, 2 };
	uint_t r = 1;

	std::cout << R::has(r) << std::endl;
	std::cout << N::has(int_t(-1)) << std::endl;


	const auto const_func = Const<real_t>{ 15. };
	const auto a = exec_func(const_func, 1155555.);
	const auto b = exec_deritative(const_func, 13.);
	const auto d = exec_primitive(const_func, 5.);

	const auto f = Function<Const, real_t>(const_func);
	const auto aa = exec_func(f, 1155555.);
	const auto bb = exec_deritative(f, 13.);
	const auto dc = exec_primitive(f, 5.);

	std::cout << "Press enter to continue...";
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	return EXIT_SUCCESS;
}