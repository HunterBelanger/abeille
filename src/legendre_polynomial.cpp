#include <tallies/legendre_polynomial.hpp>

LegendrePolynomials::LegendrePolynomials(std::size_t order)
    : legendre_coeff_(), order_(order){
    // evaluate the coefficient for each and every order
    // from zero to given order    

    legendre_coeff_.reserve(order_ + 1);

    double inv_2 = 1;
    std::size_t low_pow = 1;
    
    for ( std::size_t i = 0; i <= order_; i++ ){
        if (low_pow == 1){
            low_pow = 0;
        } else if ( low_pow == 0 ) {
            low_pow = 1;
        }
        std::size_t high_pow = static_cast<std::size_t>(i/2);
        double sign_change = 1.;
        if ( high_pow % 2 == 1 ) sign_change = -1.;
        std::vector<double> coeff;
        coeff.reserve( static_cast<std::size_t>(order / 2 + 1) );
        for ( std::size_t j = 0; j <= high_pow; j++ ){
            std::size_t k = high_pow -j;
            double denomenator = static_cast<double>(factorial(i-k) * factorial(k) * factorial(i - 2*k));
            double value = inv_2 * sign_change * static_cast<double>(factorial(2*i - 2*k)) / denomenator;
            sign_change *= -1;
            coeff.push_back(value);
        }
        legendre_coeff_.push_back(coeff);
        inv_2 *= 0.5;
    }
}


std::vector<double> LegendrePolynomials::evaluate_legendres(const double x) const {
    // get the powers upto that orders
    std::vector<double> x_powers;
    x_powers.reserve(order_+1);
    double x_pow = 1.; 
    for ( std::size_t i = 0; i <= order_; i++){
        x_powers.push_back(x_pow);
        x_pow *= x;
    }
    // evalute the legendre polynomails of each order
    std::vector<double> legndre_values;
    legndre_values.reserve(order_+1);
    std::size_t low_pow = 1;
    for ( std::size_t i = 0; i <= order_; i++){
        if ( low_pow == 0 ) low_pow = 1;
        else low_pow = 0;

        const std::vector<double> coeff = legendre_coeff_[i];
        std::size_t j = 0;
        double sum = 0.;
        for ( std::size_t k = low_pow; k <= i; k += 2){
            sum += coeff[j] * x_powers[k];
            j++;
        }
        legndre_values.push_back(sum);
    }
    return legndre_values;
}