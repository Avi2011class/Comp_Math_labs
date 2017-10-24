#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <fstream>
#include <iomanip>
#include <vector>
#include <iterator>
#include <algorithm>
#include <functional>
#include <numeric>

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/float128.hpp>

#include "matrix.hpp"


typedef boost::multiprecision::cpp_dec_float_100 lfloat;

int main( int argc, char** argv, char** env ) {

	lfloat lambda = 0.01;
	const int N = 100;

    std::function< lfloat( int, int ) > K_generator = []( int ) -> lfloat {};



	return 0;
}
