/* upperbound model to show Lagrangaean relaxation  */ 


int build_limit = ...; 

import int nbCities;
range cities 1..nbCities; 

int send[cities] = ...;
int request[cities] = ...; 

int ship_cost[cities, cities]  = ...; 

var int Ship[cities,cities] in 0..maxint; 

import int SBuild[cities]; 

var int shipping_obj in 0..maxint; 

constraint Supply_Constraint[cities];
constraint Demand_Constraint [cities]; 
constraint Limit_Constraint [cities]; 

 minimize shipping_obj 

subject to {

  shipping_obj = sum(i in cities, j in cities) ship_cost[i,j] * Ship[i,j];    

   forall(i in cities) 
  Supply_Constraint[i]:   sum(j in cities) Ship[i,j] <= send[i] * SBuild[i]; 
      
   forall(j in cities) 
Demand_Constraint[j]:   sum(i in cities) Ship[i,j] >= request[j]; 
      
   sum(i in cities)   SBuild[i]  <= build_limit;    
      
};


