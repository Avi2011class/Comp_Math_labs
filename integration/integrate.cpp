#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <iterator>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <vector>
#include <cmath>
#include <omp.h>

template < typename _T >
_T gaussian_4_integrate( std::function< _T(_T) > f, _T from, _T to, _T recomended_step ) {
	_T mul = 1;
	if ( from > to ) {
		mul = -1;
		std::swap( from, to );
	}

	int N_steps = ceil( ( to - from ) / recomended_step );
    _T step = (to - from) / N_steps;
    _T result = 0;


    #pragma omp parallel for reduction(+ : result)
    for ( int i = 0; i < N_steps; i++ ) {
        _T center = from + i * step + step / 2;
        _T a1 = center - step * 0.339981 / 2;
        _T wa1 = 0.652145;
        _T a2 = center + step * 0.339981 / 2;
        _T wa2 = 0.652145;
        _T b1 = center - step * 0.861136 / 2;
        _T wb1 = 0.347855;
        _T b2 = center + step * 0.861136 / 2;
		_T wb2 = 0.347855;
		result += (f(a1) * wa1 + f(a2) * wa2 + f(b1) * wb1 + f(b2) * wb2) * step / 2;
    }
	return result * mul;
}

double f(double x) {
	double res = log(1 / sin(x));
	return res;
}


int main ( int argc, char** argv, char** env ) {
	std::function< double(double) > ptr = f;

	double delta1 = ( gaussian_4_integrate< double >(ptr, 0, M_PI / 2, 0.001) - gaussian_4_integrate< double >(ptr, 0, M_PI / 2, 0.001 / 2) ) / pow(0.001, 4);

	double delta2 = ( gaussian_4_integrate< double >(ptr, 0, M_PI / 2, 0.0005) - gaussian_4_integrate< double >(ptr, 0, M_PI / 2, 0.0005 / 2) ) / pow(0.0005, 4);

	double delta3 = ( gaussian_4_integrate< double >(ptr, 0, M_PI / 2, 0.00025) - gaussian_4_integrate< double >(ptr, 0, M_PI / 2, 0.00025 / 2) ) / pow(0.00025, 4);

	std::cout << "Как можно видеть по значению err/h^4, требуется регуляризация" << std::endl;
	std::cout << delta1 << " " << delta2 << " " << delta3 << std::endl;

	double step_new = 1e-8;

	long double value = gaussian_4_integrate< long double >(ptr, 1e-10, M_PI / 2, step_new);
	long double err3 = ( value - gaussian_4_integrate< long double >(ptr, 1e-10, M_PI / 2, step_new / 2));

	// 2.5e-9 добавка от отсеченного отрезка [0, 1e-10]
	std::cout << "Ожидаемое зничаение интеграла: " << std::fixed << std::setprecision(12) << value + 2.5e-9<< " примерная ошибка: " << err3 << std::endl;



    return 0;
}
