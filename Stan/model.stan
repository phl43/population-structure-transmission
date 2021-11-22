data {
    int<lower=1> N; // length of data
    int<lower=1> population; // size of population
    int<lower=1> seeding_length; // number of days for seed
    int<lower=1> prediction_length; // length of prediction
    int<lower=0> infections[N];
    vector[seeding_length + prediction_length] X; // features matrix
    real<lower=0> mu; // mean of generation time
    real<lower=0> sigma; // coefficient of variation of generation time
}

transformed data {
    int total_length = seeding_length + prediction_length;
    vector[total_length] times;
    real cv = sigma / mu;
    real alpha = 1 / cv^2;
    real beta = 1 / (cv^2 * mu);
    vector[total_length] gi;
    vector[total_length] gi_rev;

    // compute generation time distribution
    for (i in 1:total_length) {
        times[i] = i;
        gi[i] = gamma_cdf(i, alpha, beta) - gamma_cdf(i - 1, alpha, beta);
    }

    // reverse generation time for convolution
    for(i in 1:total_length) {
        gi_rev[i] = gi[total_length - i + 1];
    }
}

parameters {
    real<lower=0> seeding_rate;
    real<lower=0> seed;
    real<lower=0> R0;
    real gamma;
    real<lower=0> phi;
}

transformed parameters {
    vector[total_length] Rt = R0 * exp(-X * gamma);
    vector[total_length] Rt_adj = append_row(rep_vector(R0, seeding_length), rep_vector(.0, prediction_length));
    vector[total_length] expected_infections = append_row(rep_vector(seed, seeding_length), rep_vector(0, prediction_length));
    vector[total_length] cumulative_infections = cumulative_sum(expected_infections);

    for (i in (seeding_length + 1):total_length) {
        real convolution = dot_product(head(expected_infections, i - 1), tail(gi_rev, i - 1));
        Rt_adj[i] = (population - cumulative_infections[i - 1]) / population * Rt[i];

        expected_infections[i] = Rt_adj[i] * convolution;

        cumulative_infections[i] = cumulative_infections[i - 1] + expected_infections[i];
    }
}

model {
    seeding_rate ~ uniform(0, 100);
    seed ~ exponential(seeding_rate);
    R0 ~ uniform(0, 5);
    gamma ~ normal(0, 1);
    phi ~ normal(0,5);
    for (i in 1:N) {
        infections[i] ~ neg_binomial_2(expected_infections[i], phi);
    }
}
