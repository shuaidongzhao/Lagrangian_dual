set SUBCHAIN_CHOICE ;
# Parameters
param cost { SUBCHAIN_CHOICE } >= 0;
param obs_counter { SUBCHAIN_CHOICE } > 0;
param obs_sum;
param approximated_likelihood_threshold;
# Variables
var Subchain {i in SUBCHAIN_CHOICE} binary;
var Inequality_multiplier_lamda >= 0;
var Equality_multiplier_mu;
# Objective
minimize total_cost : 
 Equality_multiplier_mu + Inequality_multiplier_lamda * ((sum {i in SUBCHAIN_CHOICE} obs_counter[i] * log10(obs_counter[i]) ) - obs_sum - approximated_likelihood_threshold ) + (obs_sum * Inequality_multiplier_lamda * log10 (Inequality_multiplier_lamda) ) + sum {j in SUBCHAIN_CHOICE} Inequality_multiplier_lamda * obs_counter[j] * log10 ( sum {k in SUBCHAIN_CHOICE} cost[k] * Subchain[k] - Equality_multiplier_mu) ;
 
# Contraints
subject to choice_sum :
sum {i in SUBCHAIN_CHOICE} Subchain[i] == 1;
subject to equality_constr :
sum {i in SUBCHAIN_CHOICE} cost[i] * Subchain[i] - Equality_multiplier_mu <= 0;