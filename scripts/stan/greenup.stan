// i = obs
// j = cell
// y_ij ~ N(mu_ij, sigma)
// GDD could be centered, doesn't have to be
// beta is effect of GDD on greenup for each cell
// mu_ij = alpha_j + beta_j * GDD_ij 
// alpha_j ~ N(mu_alpha, sigma_alpha)
// beta_j ~ N(mu_beta_j, sigma_beta)
// lat centered across entire study area
// gamma is effect of GDD on greenup at mean lat
// theta is effect of latitude on GDD 
// mu_beta_j = gamma + theta * lat_j

data {
int<lower=0> N;   // number of data points
int<lower=0> M;   // number cells
vector[N] y;      // greenup date
vector[N] gdd;    // GDD
array[N] int<lower=0> Cell;
vector[M] lat;                     // lat for each cell
}

parameters {
real<lower=0> sigma;
real<lower=0> sigma_alpha;
real<lower=0> sigma_beta;
real mu_alpha;
real theta;
real gamma;
vector[M] alpha;
// vector[M] beta;
vector[M] beta_raw;
}

transformed parameters {
vector[N] mu;
vector[M] mu_beta;
vector[M] beta;

// gamma = effect of GDD at mean lat 
// theta = effect of lat on GDD 
mu_beta = gamma + theta .* lat;
// implies beta ~ normal(mu_beta, sigma_beta); non-centered parameterization
beta = beta_raw * sigma_beta + mu_beta;

// intercept and slope for each species/site
mu = alpha[Cell] + beta[Cell] .* gdd;
}

model {
// priors
sigma ~ normal(0, 20);
sigma_alpha ~ normal(0, 10);
sigma_beta ~ std_normal();
mu_alpha ~ normal(136, 50); //centering this prior on mean greenup value
beta_raw ~ std_normal();
gamma ~ std_normal();
theta ~ std_normal();

alpha ~ normal(mu_alpha, sigma_alpha);
// see non-centered parameterization above:
// beta ~ normal(mu_beta, sigma_beta);

y ~ normal(mu, sigma);
}
