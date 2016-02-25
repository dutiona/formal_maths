#include "Formal_Maths.h"

#include <iostream>
#include <limits>


using namespace fm;

int main(int /*argc*/, char* /*argv[]*/){

	uint_t r = 1;
	
	std::cout << Domain_R::has(r) << std::endl;
	std::cout << Domain_N::has(int_t(-1)) << std::endl;

	std::cout << "Press enter to continue...";
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	return EXIT_SUCCESS;
}