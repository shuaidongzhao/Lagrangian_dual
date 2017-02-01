/* LP relaxation model */ 

int build_limit = ...; 

import int nbCities;
range cities 1..nbCities; 

int send[cities] = ...;
int request[cities] = ...; 

int ship_cost[cities, cities]  = ...; 

var float Build[cities] in 0..1; 

var float Ship[cities,cities] in 0..maxint;  

var float shipping_obj in 0..maxint; 

constraint Supply_Constraint[cities];
constraint Demand_Constraint [cities]; 
constraint Limit_Constraint [cities]; 

 minimize shipping_obj 

subject to {

  shipping_obj = sum(i in cities, j in cities) ship_cost[i,j] * Ship[i,j];  
    
   forall(i in cities) 
   sum(j in cities) Ship[i,j] <= send[i] * Build[i]; 
      
   forall(j in cities) 
   sum(i in cities) Ship[i,j] >= request[j]; 
      
   sum(i in cities)    Build[i]  <= build_limit;     
      
};


