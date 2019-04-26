data {
  int<lower=2> Na;
  int<lower=2> Nb;
  int<lower=2> N;
  vector<lower=0, upper=1>[N] Wa_w;
  int<lower=1,upper=N> Wa_v[N];
  int<lower=1> Wa_u[Na+1]; 
  vector<lower=0, upper=1>[N] Wb_w;
  int<lower=1,upper=N> Wb_v[N];
  int<lower=1> Wb_u[Nb+1]; 
  int<lower=0> a_tot[Na];
  int<lower=0> b_tot[Nb];
  int<lower=0> a_count[Na];
  int<lower=0> b_count[Nb];
  int<lower=2> Nb_groups;
  int<lower=1> b_groups[N];
  
}
parameters {
  vector[N] theta; // intersection parameters, on the a side
  vector[N] psi; // intersection parameters, on the b side
  real theta_mu;
  real<lower=0> theta_sigma;
  vector[Nb_groups] psi_mu; // used for vaalipiiri, to correct effecs from candidates
  real<lower=0> psi_mu_sigma;
  real<lower=0> psi_sigma;
  real a;
}
transformed parameters {
  vector<lower=0, upper=1>[Na] theta_cum; // weighted means of theta over intersections
  vector<lower=0, upper=1>[Nb] psi_cum;   // weighted means of psi over intersections
  theta_cum = csr_matrix_times_vector(Na, N, Wa_w, Wa_v, Wa_u, 
                                      inv_logit(theta + theta_mu));
  psi_cum   = csr_matrix_times_vector(Nb, N, Wb_w, Wb_v, Wb_u, 
                                      inv_logit(psi + a*theta + psi_mu[b_groups]));
}
model {
  a_count ~ binomial(a_tot, theta_cum);
  b_count ~ binomial(b_tot, psi_cum);
  theta ~ normal(0, theta_sigma);
  psi ~ normal(0, psi_sigma);
  psi_mu ~ normal(0, psi_mu_sigma);
}
