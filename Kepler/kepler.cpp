#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <cerrno>
#include <cmath>
#include <functional>
#include <string>

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/math/special_functions/sin_pi.hpp>

#include <boost/math/constants/constants.hpp>

#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>

int main ( int argc, char** argv, char** env ) {

	typedef boost::multiprecision::cpp_dec_float_100 lfloat;
	// Уравнение Кеплера имеет следующий вид: M = E - e * sin( E )
	// и его необзодимо рещить при условии e=0.1, M = 5 * pi / 6

	// Константы из условия
    #define const_e 0.1
    #define const_M (5 * boost::math::constants::pi< lfloat >() / 6) // pi максимально точно для используемого типа

    // Решение методом дихотомии
    // Причем достаточно надежно положить E_max = e + M + 0.1, E_min = - e - M - 0.1
    // Так как это гарантированно локализует корень

	// Желаемая абсолютная точность решения (кол-во десятичных цифр)
    #define const_digits 14

	// Тогда ошибка равна
	#define const_error boost::math::pow< const_digits >( lfloat(0.1) )
	// Функция, нули которой ищем
	std::function< lfloat( lfloat ) > f = [](lfloat E) -> lfloat {
        return E - const_e * boost::math::sin_pi(E / boost::math::constants::pi< lfloat >(), boost::math::policies::policy< boost::math::policies::digits10< const_digits + 3 > >()) - const_M;
	};

	std::function < double( double ) > fd = [](double E) -> double {
		return E + 0.5 * sin(E) + 3;
	};

	// Границы, на которых знак гарантированно отличается
	lfloat E_min = -5;
	lfloat E_max = 5;
	lfloat E_mid, F_mid;
	lfloat F_max = f(E_max);
    lfloat F_min = f(E_min);


	while (abs(E_max - E_min) > const_error * 2) { // Имеем право умножить на 2, так как вернем полусумму
        E_mid = (E_max + E_min) / 2;
        F_mid = f(E_mid);

        //std::string log = (boost::format("f(%1%) = %2%, f(%3%) = %4%, f(%5%) = %6%") % E_min % F_min % E_mid % F_mid % E_max % F_max).str();
        //std::cout << log << std::endl;

        std::cout << E_min << "    " << E_max << std::endl;

        if ( (F_mid <= 0 && F_max >= 0) || (F_mid >= 0 && F_max <= 0) ) {
			E_min = E_mid;
			F_min = F_mid;
        }
        else {
			E_max = E_mid;
			F_max = F_mid;
        }
	}

	// Усредним результат
	lfloat result = (E_max + E_min) / 2;

	std::cout << std::fixed << std::setprecision(const_digits) << result << std::endl;

	return EXIT_SUCCESS;
}

