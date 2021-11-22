#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

Rcpp::NumericVector get_dense_col(const arma::sp_mat &m, const int c) {
    Rcpp::NumericVector dense_col(m.n_rows);

    arma::sp_mat::const_iterator start = m.begin_col(c);
    arma::sp_mat::const_iterator end = m.end_col(c);

    for (arma::sp_mat::const_iterator it = start; it != end; ++it) {
        dense_col[it.row()] = *it;
    }

    return dense_col;
}

// [[Rcpp::export]]
arma::cube simulate_epidemic(
    const Rcpp::IntegerVector unit_population_sizes,
    const int simulation_length,
    const arma::cube R0,
    const double dispersion_factor,
    const double shape_gt,
    const double rate_gt,
    const arma::sp_mat between_units_travel_probs,
    const arma::cube seeding_from_abroad
) {
    int nb_units = between_units_travel_probs.n_rows;
    int nb_variants = seeding_from_abroad.n_slices;

    arma::cube infections(nb_units, simulation_length, nb_variants);
    arma::cube local_infections(nb_units, simulation_length, nb_variants);
    arma::cube imported_infections(seeding_from_abroad);

    // if I don't store the result of Rcpp::sum(infections(j, Rcpp::_)) into a variable declared as double,
    // but directly use the expression Rcpp::sum(infections(j, Rcpp::_)) / population to compute effective_R,
    // it's converted to an integer because C++ performs integer division when you use "/" with 2 integers, so
    // 1 - Rcpp::sum(infections(j, Rcpp::_)) / population is always zero and effective_R doesn't have the right value
    Rcpp::NumericVector total_infections(nb_units);

    for (int i = 0; i < simulation_length; i++) {
        for (int j = 0; j < nb_units; j++) {
            double immunity_effect = 1 - total_infections[j] / unit_population_sizes[j];
            for (int k = 0; k < nb_variants; k++) {
                double effective_R = R0.slice(i)(j, k) * immunity_effect;

                int nb_infections = local_infections(j, i, k) + imported_infections(j, i, k);

                if (nb_infections > 0 & effective_R > 0) {
                    int nb_secondary_infections = Rcpp::sum(Rcpp::rnbinom_mu(nb_infections, dispersion_factor, effective_R));
                    total_infections[j] += nb_secondary_infections;

                    Rcpp::IntegerVector nb_secondary_infections_by_unit(nb_units);
                    Rcpp::NumericVector travel_probs = get_dense_col(between_units_travel_probs, j);
                    R::rmultinom(nb_secondary_infections, travel_probs.begin(), nb_units, nb_secondary_infections_by_unit.begin());
                    
                    for (int l = 0; l < nb_units; l++) {
                        if (nb_secondary_infections_by_unit[l] > 0) {
                            Rcpp::NumericVector generation_times = Rcpp::round(Rcpp::rgamma(nb_secondary_infections_by_unit[l], shape_gt, 1 / rate_gt), 0);
                            Rcpp::NumericVector secondary_infection_times = i + generation_times;

                            for (int m = 0; m < nb_secondary_infections_by_unit[l]; m++) {
                                if (secondary_infection_times[m] < simulation_length) {
                                    infections(j, secondary_infection_times[m], k) = infections(j, secondary_infection_times[m], k) + 1;

                                    if (l == j) {
                                        local_infections(l, secondary_infection_times[m], k) = local_infections(l, secondary_infection_times[m], k) + 1;
                                    } else {
                                        imported_infections(l, secondary_infection_times[m], k) = imported_infections(l, secondary_infection_times[m], k) + 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return infections;
}