// i = obs
// j = cell
// k = ows code
// y_ijk ~ N(mu_ijk, sigma)
// GDD could be centered, doesn't have to be
// beta is effect of GDD on emergence for each code/cell
// mu_ijk = alpha_jk + beta_jk * GDD_ijk + beta2Obsijk * Obsijk
// alpha_jk ~ N(mu_alpha, sigma_alpha)
// beta_jk ~ N(mu_beta_jk, sigma_beta)
// lat is centered within each code
// gamma is effect of GDD on fledge at mean lat for each code
// theta is effect of latitude on GDD effect
// mu_beta_jk = gamma_k + theta_k * lat_jk
// theta_k ~ N(mu_theta, sigma_theta)

data {
int<lower=0> N;   // number of data points
int<lower=0> M;   // number code/cells
int<lower=0> O;   // number of OWS
vector[N] y;      // emergence date
vector[N] gdd;    // GDD
vector[N] obs;    // # of incidental observations per code/cell
array[N] int<lower=0> CellCode;
array[M] int<lower=0> code;         // id for each code/cell
vector[M] lat;                     // lat for each code/cell
}

parameters {
real<lower=0> sigma;
real<lower=0> sigma_alpha;
real<lower=0> sigma_beta;
real mu_alpha;
real mu_theta;
real mu_gamma;
real betaObs;
vector[M] alpha;
// vector[M] beta;
vector[M] beta_raw;
// vector[O] gamma;
vector<lower=0>[2] sigma_gt;
cholesky_factor_corr[2] L_Rho_gt;             // cholesky factor of corr matrix
matrix[2, O] z_gt;                          // z-scores
}

transformed parameters {
vector[N] mu;
vector[M] mu_beta;
vector[M] beta;
vector[O] gamma;
vector[O] theta;
matrix[O, 2] gt;                 // gamma, and theta
matrix[2, 2] Rho_gt;               // corr matrix

// cholesky factor of covariance matrix multiplied by z score
// implies gt ~ MVN(0, sigma)
gt = (diag_pre_multiply(sigma_gt, L_Rho_gt) * z_gt)';
// implies Rho = L_Rho * L_Rho';
Rho_gt = multiply_lower_tri_self_transpose(L_Rho_gt);

gamma = mu_gamma + gt[,1];
theta = mu_theta + gt[,2];

// intercept and slope for each code
mu_beta = gamma[code] + theta[code] .* lat;
// implies beta ~ normal(mu_beta, sigma_beta); non-centered parameterization
beta = beta_raw * sigma_beta + mu_beta;

// intercept and slope for each species/site
mu = alpha[CellCode] + beta[CellCode] .* gdd + betaObs .* obs;
}

model {
// priors
sigma ~ normal(0, 20);
sigma_alpha ~ normal(0, 20);
sigma_beta ~ normal(0,0.1);
mu_alpha ~ normal(125, 50); // mean emergence date
mu_theta ~ normal(0,0.01); 
mu_gamma ~ normal(0,0.25); 
beta_raw ~ std_normal();// equivalent to ~ normal(0, 1) but faster
betaObs ~ normal(0,0.2);
sigma_gt ~ normal(0, 0.2);        // std dev gamma, theta

alpha ~ normal(mu_alpha, sigma_alpha);
// see non-centered parameterization above:
// beta ~ normal(mu_beta, sigma_beta);
// gamma ~ normal(mu_gamma, sigma_gamma);
// theta ~ normal(mu_theta, sigma_theta);

to_vector(z_gt) ~ std_normal();
L_Rho_gt ~ lkj_corr_cholesky(6);

y ~ normal(mu, sigma);
}
