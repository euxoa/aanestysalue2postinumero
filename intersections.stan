data {
  int<lower=2> Na;
  int<lower=2> Nb;
  int<lower=2> N;
  int<lower=1> parent_a[N];
  int<lower=1> parent_b[N];
  vector<lower=0>[N] wa; // These sum to one over a given i or a
  vector<lower=0>[N] wb; // These sum to one over a given j or b
  vector<lower=0, upper=1>[N] Wa_w;
  int<lower=1,upper=N> Wa_v[N];
  int<lower=1> Wa_u[Na+1]; 
  int<lower=0> a_tot[Na];
  int<lower=0> b_tot[Nb];
  int<lower=0> a_count[Na];
  int<lower=0> b_count[Nb];
}
transformed data {
  int<lower=0, upper=N> indN[N];
  for (i in 1:N) indN[i] = i;
}
parameters {
  vector[N] theta; // intersection parameters, on the a side
  //vector<lower=0, upper=1>[N] psi; // intersection parameters, on the b side
  real theta_mu;
  real<lower=0> theta_sigma;
}
transformed parameters {
  vector<lower=0, upper=1>[Na] theta_cum; // weighted means of theta over intersections
  //vector<lower=0, upper=1>[Nb] psi_cum; // weighted means of psi over intersections
  //theta_cum = rep_vector(0, Na);
  //psi_cum = rep_vector(0, Nb);
  theta_cum = csr_matrix_times_vector(Na, N, Wa_w, Wa_v, Wa_u, inv_logit(theta));
  //for (k in 1:N) { // -> sparse matrix: theta_cum = W * inv_logit(theta)
  //  theta_cum[parent_a[k]] = theta_cum[parent_a[k]] + wa[k] * inv_logit(theta[k]); 
  //  psi_cum[parent_b[k]] = psi_cum[parent_b[k]] + wb[k] * psi[k]; 
  // }
}
model {
  a_count ~ binomial(a_tot, theta_cum);
 //psi_cum ~ beta(b_count + 1, b_tot - b_count + 1);
  theta ~ normal(theta_mu, theta_sigma);
}
