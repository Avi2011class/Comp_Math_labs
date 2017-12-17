#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cerrno>
#include <iomanip>
#include <fstream>

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/float128.hpp>

#include <boost/math/special_functions.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/detail/bessel_j0.hpp>
#include <boost/math/special_functions/detail/bessel_j1.hpp>
#include <boost/math/special_functions/detail/bessel_jy.hpp>
#include <boost/math/special_functions/detail/bessel_jy_asym.hpp>
#include <boost/math/special_functions/detail/bessel_jy_series.hpp>

// Вычисление функции Бесселя с определенной точностью,
// так как eps ожидается не более 10^-6, а тип поддерживает очень высокую точность,
// то основная ошибка берется из того, что выполняется разложение в ряд тейлора
//
// Ошибка ряда принимается за x^{2k+1} / [ (k+1)!k! 2^{2k+1} ]

// Используется ряд (-1)^k / [k!k!] * (x/2)^(2k)
boost::multiprecision::cpp_dec_float_100 bessel (
						boost::multiprecision::cpp_dec_float_100 x,
						boost::multiprecision::cpp_dec_float_100 eps
) {
	int N_steps = 1;
	// Физическое ограничение на случай непредвиденно малого eps.
	const int Max_N_steps = 1000000;
	if ( eps <= 0 ) {
		std::cerr << "Eps must be greater then zero";
		errno = EINVAL;
		return 0;
	}

	// Ищем число шагов для обеспечения нужной точности
	boost::multiprecision::cpp_dec_float_100 expected_error = abs(x) / 2;
	for ( N_steps = 1; abs( expected_error ) > eps; N_steps++ )
        expected_error = expected_error * x * x / N_steps / (N_steps + 1) / 4;

	if (N_steps > Max_N_steps)
		N_steps = Max_N_steps;

	std::cerr << "---Вычисления будут проводиться до " << N_steps << " члена ряда Тейлора" << std::endl;

    // Taylor series
    boost::multiprecision::cpp_dec_float_100 result = 0;
    boost::multiprecision::cpp_dec_float_100 current = 1;
    for ( int i = 1; i <= N_steps; i++ ) {
        result += current;
        current = -current / i / i * (x / 2) * (x / 2);
    }

	return result;
}


// Функция для вычисления производной функции Бесселя с использованием дифференцирования по
// сетке. Так как для вычислений используются высокоточные типы данных (в задании требуется
// работать с погрешностями порядка 10^(-6), что гораздо больше реально обеспечиваемой), а
// функция в узлах сетки вычисляется с точностью eps, заданной в виде аргумента, то в
// качестве машинной ошибки для выбора оптимального шага дифференцирования h будем использовать
// именно эту погрешность.
//
// h = sqrt(k * eps / M2), где M_2 = max|f''(x)|, k = max|f(x)|. В данном случае примем M=k=1 (по заданию)
//
// Используется следующая формула: f'(x) = [ f(x + h) - f(x) ] / h

boost::multiprecision::cpp_dec_float_100 bessel_derivation (
						boost::multiprecision::cpp_dec_float_100 x,
						boost::multiprecision::cpp_dec_float_100 eps
) {
	if ( eps <= 0 ) {
		std::cerr << "Eps must be greater then zero";
		errno = EINVAL;
		return 0;
	}

	boost::multiprecision::cpp_dec_float_100 h = 1e-100;
	boost::multiprecision::cpp_dec_float_100 result = ( bessel(x + h, eps) - bessel(x, eps) ) / h;

	return result;
}

int main( int argc, char** argv, char** env ) {
	boost::multiprecision::cpp_dec_float_100 X = 100;
	/*
    std::ofstream log_file;
    log_file.open("bessel_log.csv");

    log_file << "x,custom,boost,derivation" << std::endl;
    for ( double i = -10; i <= 10; i += 0.001 )
		log_file
			<< std::fixed << std::setprecision(10) << i << ","
			<< std::fixed << std::setprecision(10) << bessel(i, 1e-40) << ","
			<< std::fixed << std::setprecision(10) <<
							boost::math::cyl_bessel_j(0, i,
										boost::math::policies::policy< boost::math::policies::digits10<40> >()) << ","
			<< std::fixed << std::setprecision(10) << bessel_derivation( i, 0.000001 ) << std::endl;


	for ( double i = 20; i <= 220; i += 10 )
		log_file
			<< std::fixed << std::setprecision(10) << i << ","
			<< std::fixed << std::setprecision(10) << bessel(i, 1e-40) << ","
			<< std::fixed << std::setprecision(10) <<
							boost::math::cyl_bessel_j(0, i,
										boost::math::policies::policy< boost::math::policies::digits10<20> >()) << ","
			<< std::fixed << std::setprecision(10) << bessel_derivation( i, 0.000001 ) << std::endl;

    log_file.close();

    system("python3 Bessel_plots.py");
    std::cout << std::endl << "Отчетность в графике bessel.png и файле bessel_log.csv" << std::endl;
    //*/

    std::cout << std::fixed << std::setprecision(40) << bessel(1, boost::multiprecision::cpp_dec_float_100(1e-6));
	return 0;
}
