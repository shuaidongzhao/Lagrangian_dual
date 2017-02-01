# Parameters
param choice_size;
param cost { 1..choice_size } >= 0;
param obs_counter { 1..choice_size } > 0;
param obs_sum;
param approximated_likelihood_threshold;
# Variables
var Subchain {i in 1..choice_size} binary;
var Inequality_multiplier_lamda >= 0;
var Equality_multiplier_mu;
# Objective
minimize total_cost : 
 Equality_multiplier_mu + Inequality_multiplier_lamda * ((sum {i in 1..choice_size} obs_counter[i] * log10(obs_counter[i]) ) - obs_sum - approximated_likelihood_threshold ) + (if Inequality_multiplier_lamda = 0 then 0 else obs_sum * Inequality_multiplier_lamda * log10(Inequality_multiplier_lamda)) + sum {j in 1..choice_size} (Inequality_multiplier_lamda * obs_counter[j] * (if sum {k in 1..choice_size} cost[k] * Subchain[k] - Equality_multiplier_mu = 0 then 0 else log10( sum {k in 1..choice_size} cost[k] * Subchain[k] - Equality_multiplier_mu) ));
# Contraints
subject to choice_sum :
sum {i in 1..choice_size} Subchain[i] == 1;
subject to equality_constr :
sum {i in 1..choice_size} cost[i] * Subchain[i] - Equality_multiplier_mu <= 0;