/* Lower Bound model to show Lagrangaean relaxation  */


int build_limit = ...; 

import int nbCities;
range cities 1..nbCities; 

int send[cities] = ...;
int request[cities] = ...; 

int ship_cost[cities, cities]  = ...; 

var int Build[cities] in 0..1; 

var int Ship[cities,cities] in 0..maxint; 

var float lagrangian_obj; 

import float mult[warehouses] ; 

constraint Supply_Constraint[cities];
constraint Limit_Constraint [cities]; 

minimize lagrangian_obj 

subject to {
 
   lagrangian_obj = sum(i in cities, j in cities) ship_cost[i,j] * Ship[i,j] + 
                     sum(j in cities) mult[j] * (request[j] - sum(i in cities) Ship[i,j]); 
  
   forall(i in cities) 
Supply_Constraint[i]:   sum(j in cities) Ship[i,j] <= send[i] * Build[i]; 
          
   sum(i in cities)    Build[i]  <= build_limit;     
      
};


