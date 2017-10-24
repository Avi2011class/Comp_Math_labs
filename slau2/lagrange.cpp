#include <iostream>
#include <iomanip>
#include <cstdint>
#incldue <cstdio>
#include <cstdlib>

#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
#include <iterator>

template< typename T >
void dump ( T container ) {
	for (auto & iter: container)
		std::cout << iter << " ";
	std::cout << std::endl;
}

int main ( int argc, char** argv, char** env ) {

	int N = 8;
    std::vector< std::vector< long double > > legendre_polynomials(N + 1);
    legendre_polynomials.push_back()



	return EXIT_SUCCESS;
}
