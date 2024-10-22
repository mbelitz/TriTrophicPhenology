// i = obs
// j = station
// k = species
// y_ijk ~ N(mu_ijk, sigma)
// GDD could be centered, doesn't have to be
// beta is effect of GDD on fledge for each species/station
// mu_ijk = alpha_jk + beta_jk * GDD_ijk 
// alpha_jk ~ N(mu_alpha, sigma_alpha)
// beta_jk ~ N(mu_beta_jk, sigma_beta)
// lat is centered within each species
// gamma is effect of GDD on fledge at mean lat for each species
// theta is effect of latitude on GDD effect
// mu_beta_jk = gamma_k + theta_k * lat_jk
// kappa is effect of GDD on fledge at mean species PC1
// phi is effect of PC1 on species-level GDD effect
// gamma_k ~ N(mu_gamma_k, sigma_gamma)
// mu_gamma_k = kappa + phi * PC1[k]
// theta_k ~ N(mu_theta, sigma_theta)

data {
int<lower=0> N;                        // number of data points
int<lower=0> M;                      // number species/stations
int<lower=0> O;                      // number species
vector[N] y;                  // fledge date
vector[N] gdd;                  // GDD
array[N] int<lower=0> SpSt;
array[M] int<lower=0> Sp;          // species id for each species/station
vector[M] lat;                     // lat for each species/station
vector[O] PC1;                     // PC1 for bird traits
}

parameters {
real<lower=0> sigma;
real<lower=0> sigma_alpha;
real<lower=0> sigma_beta;
real<lower=0> sigma_gamma;
real<lower=0> sigma_theta;
real mu_alpha;
real mu_theta;
vector[M] alpha;
// vector[M] beta;
vector[M] beta_raw;
// vector[O] gamma;
vector[O] gamma_raw;
// vector[O] theta;
vector[O] theta_raw;
real kappa;
real phi;
}

transformed parameters {
vector[N] mu;
vector[M] mu_beta;
vector[M] beta;
vector[O] mu_gamma;
vector[O] gamma;
vector[O] theta;

// one intercept and slope
mu_gamma = kappa + phi * PC1;
// implies gamma ~ normal(mu_gamma, sigma_gamma); non-centered parameterization
gamma = gamma_raw * sigma_gamma + mu_gamma;

// implies theta ~ normal(mu_theta, sigma_theta)
theta = theta_raw * sigma_theta + mu_theta;

// intercept and slope for each species
mu_beta = gamma[Sp] + theta[Sp] .* lat;
// implies beta ~ normal(mu_beta, sigma_beta); non-centered parameterization
beta = beta_raw * sigma_beta + mu_beta;

// intercept and slope for each species/site
mu = alpha[SpSt] + beta[SpSt] .* gdd;
}

model {
// priors
sigma ~ normal(0, 20);
sigma_alpha ~ normal(0, 10);
sigma_beta ~ std_normal();
sigma_gamma ~ std_normal();
sigma_theta ~ std_normal();
mu_alpha ~ normal(200, 50);
mu_theta ~ std_normal(); // equivalent to ~ normal(0, 1) but faster
kappa ~ std_normal();
phi ~ std_normal(); 
beta_raw ~ std_normal();
gamma_raw ~ std_normal();
theta_raw ~ std_normal();

alpha ~ normal(mu_alpha, sigma_alpha);
// see non-centered parameterization above:
// beta ~ normal(mu_beta, sigma_beta);
// gamma ~ normal(mu_gamma, sigma_gamma);
// theta ~ normal(mu_theta, sigma_theta);

y ~ normal(mu, sigma);
}
