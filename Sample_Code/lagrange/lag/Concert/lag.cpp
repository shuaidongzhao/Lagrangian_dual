/************************************************************************************
 * The model shows a Lagrangian relaxation for a location-transportation problem.   *
 * The original MIP is decomposed into two problems in order to deduce a multiplier *
 * for a particular constraint based on Lagrange relaxation.						*
 * Written by: Sanjay Ramanujan					  									*
 * Version: 1.0, November 1st, 2002, Tested with CPLEX v8.0							*					
 * For consistency, this is built using an equivalent AMPL model at:				*
 * http://www.ampl.com/NEW/loop2.html												* 
 ************************************************************************************/
#include <random>
#include <stdint.h>
#include <algorithm>
#include <vector>
#include <ilcplex/ilocplex.h>
#include <boost/math/distributions/non_central_chi_squared.hpp>
#define smallNumber 0.000001
#include "ampl/ampl.h"
ILOSTLBEGIN

typedef IloArray<IloNumArray> TwoDMatrix;
typedef IloArray<IloNumVarArray> NumVarMatrix;

// typedef __int32 int32_t;
// typedef unsigned __int32 uint32_t;
typedef std::mt19937 MyRNG;  // the Mersenne Twister with a popular choice of parameters
          // populate somehow
uint32_t seed_val;
MyRNG rng;                   // e.g. keep one global instance (per thread)

struct ConnectedVehicleActivityData
{

};

struct IndividualSpaceTimeBicone
{

	vector<vector<string>> all_possible_td_time_key_activity_index_string_as_subchain_choice_2D_array;
	vector<vector<string>> all_feasible_td_time_key_activity_index_string_as_subchain_choice_2D_array;
	vector<float> time_cost_for_each_feasible_td_sub_chain_choice_array;

	//vector<vector<ConnectedVehicleActivityData>> feasible_subchain_choices_from_daily_hitorical_subhcain_array;
	vector<int> feasible_subchain_choices_index_array;
	vector<int> sub_gradient_subchain_choice_vector;
	double sub_gradient_opt_inequality_multi_lamda;
	double sub_gradient_opt_equality_multi_mu;
	//vector<double> sub_gradient_solution_subchain_choice_possibility_vector;
	vector<int> directly_solution_subchain_choice_vector;
	vector<double> directly_solution_subchain_choice_possibility_vector;
	double directly_opt_inequality_multi_lamda;
	double directly_opt_equality_multi_mu;

	vector<float> time_cost_for_each_feasible_td_sub_chain_choice;

	string vehicle_date_id_string;
};

struct IndividualDataDrivenInputData
{
	int sum_subchain_num_observations;
	vector<int> subchain_choice_index_to_counter;
	vector<int> subchain_choice_array;
	double estimated_likelihood_level;

	double get_estimated_likelihood_level_for_individual_from_activity_type()
	{
		this->estimated_likelihood_level = 0;
		float alpha = 0.05;
		float confidence_level = 1 - alpha;

		int degree_of_freedom = 0;
		int N = 0;

		N = this->sum_subchain_num_observations;// num_unique_obs = this->subchain_choice_index_to_counter.size();

		

		double sum_first_term_approximation_likelihood_level = 0.0;

		for (int i = 0; i < this->subchain_choice_index_to_counter.size(); i++)//while (iter != this->td_individual_subchain_activity_type_sequence_string_counter_map.end())
		{
			int N_i = this->subchain_choice_index_to_counter[i];
			
			double portion = double(N_i) / double(N);
			sum_first_term_approximation_likelihood_level += ((double)N_i)* log(portion); //e base, constant value
		}

		degree_of_freedom = N - 1;
		boost::math::chi_squared chi_2(degree_of_freedom);
		double quantile_of_chi_2 = boost::math::quantile(chi_2, confidence_level);

		this->estimated_likelihood_level = sum_first_term_approximation_likelihood_level - quantile_of_chi_2;
		return this->estimated_likelihood_level;
	}
};

static
void displayResults(IloModel, IloCplex&, IloBoolVarArray, NumVarMatrix, IloInt);

static
double solveRelaxed(IloModel, IloBoolVarArray, NumVarMatrix, IloInt, IloNumArray);

static
double solve_LB_relaxed_problem_cplex(
IndividualSpaceTimeBicone &space_time_bicone,
IndividualDataDrivenInputData &individualdatadriveninputdata,
IloCplex &output_cplex,
double approximated_likelihood_threshold);

double solve_LB_relaxed_problem_ampl(
	IndividualSpaceTimeBicone &space_time_bicone,
	IndividualDataDrivenInputData &individualdatadriveninputdata,
	double approximated_likelihood_threshold);

static
double approximate_likelihood_level(
IndividualSpaceTimeBicone &space_time_bicone,
IndividualDataDrivenInputData &individualdatadriveninputdata);

static
void initialization(
IndividualSpaceTimeBicone &space_time_bicone,
IndividualDataDrivenInputData &individualdatadriveninputdata);

int
main() {
	IndividualSpaceTimeBicone space_time_bicone;
	IndividualDataDrivenInputData individualdatadriveninputdata;
	initialization(
		 space_time_bicone,
		 individualdatadriveninputdata);




	double approximated_likelihood_threshold = individualdatadriveninputdata.estimated_likelihood_level;


	solve_LB_relaxed_problem_ampl(
		space_time_bicone,
		individualdatadriveninputdata,
		approximated_likelihood_threshold);

// 		IloCplex output_cplex;
// 	solve_LB_relaxed_problem_cplex
// 		(space_time_bicone,
// 		individualdatadriveninputdata,
// 		output_cplex,
// 		approximated_likelihood_threshold);
		
#ifdef LINUX
	system("read");
#else
	system("PAUSE");
#endif	

	return 1;
}

int langrange_dual(int argc, char **argv)
{

 IloEnv env;

 try {

	 IloInt i=0, j=0, k=0;
	 IloInt nbCities=0;
	 IloInt build_limit=0;


	 IloNumArray send(env), request(env);

	 TwoDMatrix ship_cost(env);



///////////////// DATA FILE READING ////////////////////////////////
	  const char* filename  = "C:\\Users\\sdzhao\\Dropbox\\Papers\\CV activity chain\\C++_implementation\\Lagrangian_dual\\Sample_Code\\lagrange\\lag\\Concert\\data\\lagdata.dat";
      if (argc > 1)
         filename = argv[1];
      ifstream file(filename);
      if (!file) {
         cerr << "ERROR: could not open file '" << filename
              << "' for reading" << endl;
         cerr << "usage:   " << argv[0] << " <file>" << endl;
         throw(-1);
      }

	  	file >> build_limit >> send >> request >> ship_cost ;
		nbCities = send.getSize();

		env.out() << "Total number of cities: " << nbCities << endl;
		env.out() <<"build_limit value: " << build_limit << endl;


		IloBool consistentData = (request.getSize() == nbCities && ship_cost.getSize() == nbCities);
		if (!consistentData) {
		 cerr << "ERROR: data file '"
              << filename << "' contains inconsistent data" << endl;
         throw(-1);
		}



/////////////////// DECISION VARIABLES WITH NAMES  /////////////////////////////
		IloBoolVarArray  Build(env, nbCities);
		NumVarMatrix Ship(env, nbCities);
		IloNumArray mult(env, nbCities);

		for(k=0; k < nbCities; k++) {
			mult[k] = 0.0;
		}

		for(i=0; i <nbCities; i++) {
			Ship[i] = IloNumVarArray(env, nbCities, 0, CPX_INFBOUND, ILOINT);
		}

		for(i=0; i < nbCities; i++) {
			for(j=0; j < nbCities; j++) {
				Ship[i][j] = IloNumVar(env, 0, CPX_INFBOUND, ILOINT); }
		}


		IloModel model(env);

		IloCplex cplex(env);

//		ofstream fout("mylog.log");
//	    cplex.setOut(fout);


		cplex.setOut(env.getNullStream());

		cplex.setWarning(env.getNullStream());

		cplex.setError(env.getNullStream());


		char *buffer = NULL;
		buffer = new char[200];

		for(i=0; i< nbCities; i++) {
		 sprintf(buffer, "Build(%d)",i);
		 Build[i].setName(buffer);
		}

		for(i=0; i< nbCities; i++) {
			for(j=0; j< nbCities; j++) {
		 sprintf(buffer, "Ship(%d,%d)",i,j);
		 Ship[i][j].setName(buffer);
			}
		}

		delete[] buffer;

		IloInt  iter_limit;
		cout << "Enter number of iterations desired: " << endl;
		cin >> iter_limit;


////////////DEVELOP GENERIC MODEL //////////////////////////

		//shipping objective function
		IloExpr shipobj(env);
		for(i=0; i< nbCities; i++) {
			shipobj += IloScalProd(Ship[i],ship_cost[i]);
		}

		IloObjective shipping_obj(env);
		shipping_obj.setName("shipping_obj");
		shipping_obj = IloAdd(model,IloMinimize(env, shipobj));


		// supply_constraint
		//IloRangeArray supply_constr(env);
		IloConstraintArray supply_constr(env);
		for(i=0; i< nbCities; i++) {
			supply_constr.add((IloSum(Ship[i]) <= Build[i] * send[i]));
			//model.add((IloSum(Ship[i]) <= Build[i] * send[i]));
		}

		model.add(supply_constr);

		// limit_constraint
		IloRange limit_constr(env, -IloInfinity, IloSum(Build), build_limit, "limit_constr");
		model.add(limit_constr);



////////////////////// SOLVE THE RELAXED MODEL NOW //////////////////////////////////////

		double LB = solveRelaxed(model, Build, Ship, nbCities, request);

		env.out() << endl << "LP Relaxation value: " << LB << endl << endl;

		model.remove(shipping_obj);

		IloModel mlowerbound(env);

		mlowerbound.add(model);

		IloCplex clb(env);

		clb.extract(mlowerbound);

		clb.setOut(env.getNullStream());

		clb.setWarning(env.getNullStream());

		clb.setError(env.getNullStream());



		IloModel mupperbound(env);

		shipping_obj = IloAdd(mupperbound,IloMinimize(env, shipobj));

		for(j=0; j< nbCities; j++) {
			IloExpr expD(env);
			for(i=0; i< nbCities; i++) {
				expD += Ship[i][j];
				}
				mupperbound.add(expD >= request[j]);
				expD.end();
			}

		IloCplex cub(env);

		cub.setOut(env.getNullStream());

		cub.setWarning(env.getNullStream());

		cub.setError(env.getNullStream());


		cub.extract(mupperbound);



///////////////// DEFINE LAGRANGE VARIABLES  ////////////////////////////////
		IloInt same		  = 0,
			   same_limit = 3;


		 IloNum scale = 1.0,
				norm  = 0.0,
				step  = 0.0,
				UB	  = 0.0;


		IloNumArray LBlog(env, iter_limit);
		for(k=0; k< iter_limit; k++) {
		LBlog[k] = 0.0;
		}

		LBlog[0] = LB;


		IloNumArray slack(env, nbCities);
			for (i=0; i< nbCities; i++) {
				slack[i] = 0.0;
			}


		IloNumArray temp(env, nbCities);
		for(i=0; i< nbCities; i++) {
		temp[i] = IloMax(ship_cost[i]);
		UB = IloSum(temp);
		}
		temp.end();


		IloNumArray UBlog(env, iter_limit);
		UBlog[0] = UB;


		IloNumArray scalelog(env, iter_limit);
		IloNumArray steplog(env, iter_limit);

		IloNum Lagrangian = 0.0;



////////// BEGIN LAGRANGE ITERATIONS HERE /////////////////////////////////////

		for(k=0; k<iter_limit; k++) {

			env.out() << "  " << endl;

			env.out() << "  ITERATION  " << k+1 << endl;

		IloExpr lagrobj1(env), lagrobj2(env), lagrobj3(env), lagrobj(env);
		for(i=0; i< nbCities; i++) {
			lagrobj1 += IloScalProd(Ship[i],ship_cost[i]);
		}
		for(j=0; j< nbCities; j++) {
			lagrobj2 += mult[j] * request[j];
		}
		for(j=0; j< nbCities; j++) {
			lagrobj3 += IloScalProd(mult, Ship[j]);
		}


		lagrobj += lagrobj1 + lagrobj2 - lagrobj3;

		IloObjective lagrange_obj(env);

		lagrange_obj.setName("lagrange_obj");

	    lagrange_obj = IloAdd(mlowerbound,IloMinimize(env, lagrobj));

		lagrobj.end();


		 if (clb.solve()) {

			Lagrangian = clb.getObjValue();

			env.out() << "lower bound model obj value " << Lagrangian << endl;

			TwoDMatrix tempShipMatrix(env, nbCities);
			for(i=0; i< nbCities; i++) {
				tempShipMatrix[i] = IloNumArray(env, nbCities);
			}

			for(i=0; i< nbCities; i++) {
				clb.getValues(tempShipMatrix[i], Ship[i]);
			}


		IloNumArray tempSumArray(env, nbCities);

		for(i=0; i< nbCities; i++) {
			for(j=0; j< nbCities; j++) {
			tempSumArray[i] += tempShipMatrix[j][i]; }
			}


		for(i=0; i< nbCities; i++) {
		  slack[i] = tempSumArray[i] - request[i];
			}

		    tempShipMatrix.end();
			tempSumArray.end();


			if (Lagrangian > LB + smallNumber) {
				LB = clb.getObjValue();
				same = 0;
				}
			else {	same = same + 1;  }


		 }


			if (same == same_limit) {
				scale = scale/2;
				same  = 0;
			}


			IloNumArray normtemp(env, nbCities);

			for(j=0; j< nbCities; j++) {
				normtemp[j] = IloPower(slack[j], 2); }

			norm = IloSum(normtemp);
			normtemp.end();


			step = scale * ( (UB - Lagrangian) / norm );

			norm = 0.0;


			IloNumArray SBuild(env, nbCities);

			IloNum tolerance = clb.getParam(IloCplex::EpInt);

			for(j=0; j< nbCities; j++) {

			if(clb.getValue(Build[j]) > 1 - tolerance) {

					SBuild[j] = clb.getValue(Build[j]);
				}
			}

			IloRangeArray ub_supply_constr(env);

			for(i=0; i< nbCities; i++) {
				ub_supply_constr.add(IloSum(Ship[i]) <= SBuild[i] * send[i]);
			}

			mupperbound.add(ub_supply_constr);



		if ( (IloScalProd(send, SBuild) ) >=
				( IloSum(request) - 1.0/IloPower(10,8) ) ) {

			if (cub.solve()) {

			cub.solve();

			env.out() << "upper bound model obj value " << cub.getObjValue() << endl;

			if(cub.getObjValue() <= UB) {
				UB = cub.getObjValue();
			}
			else {	UB = UB; }

			}

		}



	LBlog[k] 	=  LB;
	UBlog[k]	=  UB;
	scalelog[k]	=  scale;
	steplog[k]  =  step;



	for(j=0; j< nbCities; j++) {

		if(mult[j] - (step * slack[j]) > 0) {
			mult[j] = mult[j] - (step  * slack[j]);
		}

		else {	mult[j] = 0; }
	}



			//remove for next set of runs
				mupperbound.remove(ub_supply_constr);
				ub_supply_constr.end();
				mlowerbound.remove(lagrange_obj);
				lagrange_obj.end();



				if (k == (iter_limit-1)) {

						cout << " " << endl << endl;

						env.out() << " Results " << endl << endl;

					for(i=0; i< iter_limit; i++) {
						cout << "LBlog[" << i << "]" << LBlog[i] << endl; }

						cout << " " << endl;

						for(i=0; i< iter_limit; i++) {
							cout << "UBlog[" << i << "]" << UBlog[i] << endl; }

						cout << " " << endl;

						for(i=0; i< iter_limit; i++)  {
							cout << "scalelog[" << i << "]" << scalelog[i] << endl; }

						cout << " " << endl;

						for(i=0; i< iter_limit; i++) {
							cout << "steplog[" << i << "]" << steplog[i] << endl; }

						cout << " " << endl;


				}


				env.out() << " " << endl;


}	// end of for(k=0; k< iter_limit; k++)


 }
 catch(IloException &e) {
	env.out() << "ERROR: " << e << endl;
 }
 catch(...){
	env.out() << "Unknown exception" << endl;
 }

 env.end();

#ifdef LINUX
 system("read");
#else
 system("PAUSE");
#endif	

return 0;
}


static
void displayResults(IloModel mdl,
					IloCplex& cplex,
					IloBoolVarArray Build,
					NumVarMatrix Ship,
					IloInt nbCities)
{
		IloEnv env = mdl.getEnv();
		IloInt j=0, k=0;
		env.out() << "Optimal value: " << cplex.getObjValue() << endl << endl;

		env.out() << " --------------------------- " << endl;

		IloNum tolerance = cplex.getParam(IloCplex::EpInt);

		for(k=0; k< nbCities; k++) {
			if(cplex.getValue(Build[k]) > 1 - tolerance)
		env.out() << "Build["<< k <<"] " <<" = " << cplex.getValue(Build[k]) << endl;
			}

		env.out() << " --------------------------- " << endl;

		for(k=0; k< nbCities; k++)
			for(int j=0; j< nbCities; j++) 		{
			if(cplex.getValue(Ship[k][j]) >= 1 - tolerance)
		env.out() << "Ship["<< k << "]" << "[" << j <<"] " <<" = " << cplex.getValue(Ship[k][j]) << endl;
			}

        env.out() << "Time: " << cplex.getTime() << endl << endl;

        env.end();
}


static
double solveRelaxed(IloModel mdl,
				  IloBoolVarArray bvar,
				  NumVarMatrix nvar,
				  IloInt nbCities,
				  IloNumArray req)

{
	IloEnv env = mdl.getEnv();
	IloModel relax(env);
	relax.add(mdl);

	relax.add(IloConversion(env, bvar, ILOFLOAT));
	
	//env.out() << endl << req.get << endl << endl;

	for(int i=0; i<nbCities; i++) {
	relax.add(IloConversion(env, nvar[i], ILOFLOAT)); }

	for(int j=0; j< nbCities; j++) {
		IloExpr expR(env);
			for(int i=0; i< nbCities; i++) {
				expR += nvar[i][j];
					}
	relax.add( expR >= req[j]);
	expR.end();
			}

	IloCplex cplex(env);

	cplex.setOut(env.getNullStream());

	cplex.extract(relax);

	cplex.solve();
	return cplex.getObjValue();

	env.end();
}

static 
double solve_LB_relaxed_problem_cplex(
IndividualSpaceTimeBicone &space_time_bicone,
IndividualDataDrivenInputData &individualdatadriveninputdata,
IloCplex &output_cplex,
double approximated_likelihood_threshold)
{
	IloEnv env;
	double result = DBL_MAX;
	int total_num_feasible_subchain_choices = space_time_bicone.feasible_subchain_choices_index_array.size();
	try
	{
		IloInt i, j;
		IloModel mod(env);
		//----------------------decision variable----------------------------------------
		string c = "equality_multiplier_mu";
		IloNumVar equality_multiplier_mu = IloNumVar(env, -IloInfinity, IloInfinity, IloNumVar::Type::Float, c.c_str());// (env, total_num_feasible_subchain_choices);
		c = "inequality_multiplier_lamda";
		IloNumVar inequality_multiplier_lamda = IloNumVar(env, 0, IloInfinity, IloNumVar::Type::Float, c.c_str());;

		IloBoolVarArray subchain_choice(env, total_num_feasible_subchain_choices);
		for (i = 0; i < total_num_feasible_subchain_choices; i++)
		{
			string c = "subchain_choice[" + to_string(i) + "]";
			subchain_choice[i] = IloBoolVar(env, c.c_str());
		}

		//----------------------Constraints----------------------------------------

		IloExpr equality_cons_expr(env);
		IloExpr inequality_cons_expr(env);
		for (i = 0; i < total_num_feasible_subchain_choices; i++)
		{
			float cost_value = space_time_bicone.time_cost_for_each_feasible_td_sub_chain_choice[i];// space_time_bicone.time_cost_for_each_feasible_td_sub_chain_choice_array[i];

			inequality_cons_expr += cost_value * subchain_choice[i];
			equality_cons_expr += subchain_choice[i];
		}
		mod.add(equality_cons_expr == 1);
		mod.add(inequality_cons_expr - equality_multiplier_mu <= 0);
		//mod.add(inequality_multiplier_lamda >= 0);

		//----------------------Objective func----------------------------------------

		IloExpr obj_expr(env);

		obj_expr += equality_multiplier_mu;



		int N = individualdatadriveninputdata.sum_subchain_num_observations;
		for (i = 0; i < total_num_feasible_subchain_choices; i++)
		{
			//string activity_type_sequence_string = space_time_bicone.feasible_subchain_choices_activity_type_sequence_from_daily_hitorical_subhcain_array[i];
		
			int feasible_subchain_index = space_time_bicone.feasible_subchain_choices_index_array[i];

			int N_i = individualdatadriveninputdata.subchain_choice_index_to_counter[feasible_subchain_index];// individualdatadriveninputdata.td_individual_subchain_activity_type_sequence_string_counter_map[activity_type_sequence_string];

			obj_expr += inequality_multiplier_lamda *(((double)N_i)* log(N_i) - N - approximated_likelihood_threshold);
			IloExpr inner_summation(env);
			for (j = 0; j < total_num_feasible_subchain_choices; j++)
			{
				float cost_value = space_time_bicone.time_cost_for_each_feasible_td_sub_chain_choice[j];// space_time_bicone.time_cost_for_each_feasible_td_sub_chain_choice_array[i];
				inner_summation += cost_value* subchain_choice[j];
			}

			obj_expr += inequality_multiplier_lamda * N_i * IloLog10(inner_summation - equality_multiplier_mu);

		}

		obj_expr += N*inequality_multiplier_lamda*IloLog10(inequality_multiplier_lamda);

		mod.add(IloMinimize(env, obj_expr));
		IloCplex cplex(mod);
	

		cplex.setParam(IloCplex::Param::Preprocessing::Presolve, IloFalse);

		cplex.setParam(IloCplex::Param::Threads, 1);

		cplex.setParam(IloCplex::Param::MIP::Strategy::Search, IloCplex::Traditional);
		if (cplex.solve())
		{
			cplex.solve();
			result = cplex.getObjValue();

			IloNum tolerance = cplex.getParam(IloCplex::EpInt);

			space_time_bicone.directly_opt_inequality_multi_lamda = cplex.getValue(inequality_multiplier_lamda);
			space_time_bicone.directly_opt_equality_multi_mu = cplex.getValue(equality_multiplier_mu);

			space_time_bicone.directly_solution_subchain_choice_vector.clear();
			space_time_bicone.directly_solution_subchain_choice_possibility_vector.clear();

			for (int k = 0; k < total_num_feasible_subchain_choices; k++)
			{
				if (cplex.getValue(subchain_choice[k]) >= 1 - tolerance)
				{
					int feasible_subchain_index = space_time_bicone.feasible_subchain_choices_index_array[k];

					int N_i = individualdatadriveninputdata.subchain_choice_index_to_counter[feasible_subchain_index];

					float cost_value = space_time_bicone.time_cost_for_each_feasible_td_sub_chain_choice[k];
					space_time_bicone.directly_solution_subchain_choice_vector.push_back(cplex.getValue(subchain_choice[k]));

					double possibility = (space_time_bicone.directly_opt_inequality_multi_lamda * N_i) / (cost_value * cplex.getValue(subchain_choice[k]) - space_time_bicone.directly_opt_equality_multi_mu);
					space_time_bicone.directly_solution_subchain_choice_possibility_vector.push_back(possibility);
				}
			}

		}
		output_cplex = cplex;
	}
	catch (IloException& ex) {
		cerr << "IloException Error: " << ex << endl;
	}
	catch (...) {
		cerr << "Unknown Error" << endl;
	}
	return result;

	env.end();
};






static
void initialization(
IndividualSpaceTimeBicone &space_time_bicone,
IndividualDataDrivenInputData &individualdatadriveninputdata)
{
	seed_val = 5000;
	rng.seed(seed_val);

	std::uniform_int_distribution<uint32_t> uint_dist;         // by default range [0, MAX]
	std::uniform_int_distribution<uint32_t> uint_dist1_10(1, 10); // range [0,10]
	//std::normal_distribution<double> normal_dist(mean, stddeviation);  // N(mean, stddeviation)

	int choice_size = 30;
	individualdatadriveninputdata.subchain_choice_index_to_counter.resize(choice_size);

	//generate N_i;
	int N = 0;
	for (int i = 0; i < choice_size; i++)
	{
		int ctr = uint_dist1_10(rng);
		individualdatadriveninputdata.subchain_choice_index_to_counter[i] = ctr;
		N += ctr;
	}
	individualdatadriveninputdata.sum_subchain_num_observations = N;
	std::uniform_int_distribution<uint32_t> uint_dist0_29(0, 30); 
	int feasible_subchain_choice_size = uint_dist0_29(rng);
	std::uniform_real_distribution<float> float_dist20(0.01, 20);
	space_time_bicone.feasible_subchain_choices_index_array.resize(feasible_subchain_choice_size);
	space_time_bicone.time_cost_for_each_feasible_td_sub_chain_choice.resize(feasible_subchain_choice_size);
	for (int i = 0; i < feasible_subchain_choice_size; i++)
	{
		int rnd_subchain_index;
		while (true)
		{
			rnd_subchain_index = uint_dist0_29(rng);
			if (find(space_time_bicone.feasible_subchain_choices_index_array.begin(), space_time_bicone.feasible_subchain_choices_index_array.end(), rnd_subchain_index) == space_time_bicone.feasible_subchain_choices_index_array.end())
				break;
		}
		space_time_bicone.feasible_subchain_choices_index_array[i] = rnd_subchain_index;

		space_time_bicone.time_cost_for_each_feasible_td_sub_chain_choice[i] = float_dist20(rng);
	}
	
	


	individualdatadriveninputdata.get_estimated_likelihood_level_for_individual_from_activity_type();

	cout << "Finish Initialization" << endl;
	//space_time_bicone.feasible_subchain_choices_index_array;
};

double solve_LB_relaxed_problem_ampl(
	IndividualSpaceTimeBicone &space_time_bicone,
	IndividualDataDrivenInputData &individualdatadriveninputdata,
	double approximated_likelihood_threshold)
{
	ampl::AMPL ampl;
	string fn_dir = "C:/Users/sdzhao/Dropbox/Papers/CV activity chain/C++_implementation/Lagrangian_dual";
	ampl.read(fn_dir + "/activity_estimation.mod");

	int total_num_feasible_subchain_choices = space_time_bicone.feasible_subchain_choices_index_array.size();
	//vector<double> feasible_subchain_choices_index_array(space_time_bicone.feasible_subchain_choices_index_array.begin(), space_time_bicone.feasible_subchain_choices_index_array.end());

	vector<double> cost_array(space_time_bicone.time_cost_for_each_feasible_td_sub_chain_choice.begin(), space_time_bicone.time_cost_for_each_feasible_td_sub_chain_choice.end());

	ampl::DataFrame df(1, "SUBCHAIN_CHOICE");
	
	//df.setColumn("SUBCHAIN_CHOICE", feasible_subchain_choices_index_array.data(), total_num_feasible_subchain_choices);
	//df.addColumn("cost", cost_array.data());
	//ampl.setData(df);
	//choice_size
	ampl::Parameter choice_size = ampl.getParameter("choice_size");
	choice_size.set(total_num_feasible_subchain_choices);

	ampl::Parameter cost = ampl.getParameter("cost");
	cost.setValues(cost_array.data(), cost_array.size());

	ampl::Parameter approx_likelihood_threshold = ampl.getParameter("approximated_likelihood_threshold");
	approx_likelihood_threshold.set(approximated_likelihood_threshold);

	ampl::Parameter obs_sum = ampl.getParameter("obs_sum");
	obs_sum.set(individualdatadriveninputdata.sum_subchain_num_observations);

	//obs_counter
	vector<double> observation_ctr;
	int N = individualdatadriveninputdata.sum_subchain_num_observations;
	for (int i = 0; i < total_num_feasible_subchain_choices; i++)
	{
		//string activity_type_sequence_string = space_time_bicone.feasible_subchain_choices_activity_type_sequence_from_daily_hitorical_subhcain_array[i];

		int feasible_subchain_index = space_time_bicone.feasible_subchain_choices_index_array[i];

		int N_i = individualdatadriveninputdata.subchain_choice_index_to_counter[feasible_subchain_index]; 
		observation_ctr.push_back(double(N_i));
	}
	ampl::Parameter obs_counter = ampl.getParameter("obs_counter");
	obs_counter.setValues(observation_ctr.data(), observation_ctr.size());

	ampl.solve();

	std::cout << "Objective: " << ampl.getObjective("total_cost").value() << "\n";

	ampl::DataFrame results = ampl.getVariable("Subchain").getValues();
	std::cout << results.toString() << "\n";
// 	for (int i = 0; i < results.getNumRows(); i++)
// 	{
// 		for (int j = 0; j < results.getNumRows(); j++)
// 		{
// 			//std::cout << results.getRowByIndex(j)[0].c_str() << "\n";
// 			// std::cout << results.getRowByIndex(j)[1].dbl() << "\n";
// 		}
// 	}
	for (int i = 0; i < space_time_bicone.feasible_subchain_choices_index_array.size(); i++)
	{
		cout << space_time_bicone.feasible_subchain_choices_index_array[i] << endl;
	}

	return ampl.getObjective("total_cost").value();
};