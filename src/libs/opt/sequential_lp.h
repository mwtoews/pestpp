#pragma once
#ifndef SEQUENTIAL_LP_H
#define SEQUENTIAL_LP_H

#include "pest.h"
#include "Jacobian_1to1.h"
#include "ClpSimplex.hpp"
#include "RunManagerAbstract.h"
#include "TerminationController.h"
#include "covariance.h"
#include "FileManager.h"
#include "OutputFileWriter.h"
#include "Transformable.h"
#include "ModelRunPP.h"

class sequentialLP
{
	enum ConstraintSense {less_than,greater_than,equal_to};
public:
	sequentialLP(Pest &_pest_scenario, RunManagerAbstract* _run_mgr_ptr, 
		         TerminationController* _termination_ctl_ptr, Covariance &_parcov, 
				 FileManager &_file_mgr, OutputFileWriter* _out_wtr_ptr);
	void initialize_and_check();
	void solve();
	
	//ModelRun get_optimum_run() { return optimum_run; }

private:
	string obj_func_str;
	map<string, double> obj_func_coef_map;
	Observations obj_func_obs;
	ObservationInfo obj_func_info;
	double* ctl_ord_obj_func_coefs;
	int slp_iter;
	map<string, ConstraintSense> constraint_sense_map;
	map <string, string> constraint_sense_name;
	vector<string> ctl_ord_dec_var_names;
	vector<string> ctl_ord_constraint_names;
	double* dec_var_lb;
	double* dec_var_ub;
	double* constraint_lb;
	double* constraint_ub;
	PriorInformation* null_prior = new PriorInformation();
	//ModelRun current_run;
	//ModelRun optimum_run;
	ObjectiveFunc obj_func;
	Parameters decision_vars;
	Observations constraints_obs;
	Observations constraints_sim;
	Pest pest_scenario;
	Pest opt_scenario;
	RunManagerAbstract* run_mgr_ptr;
	TerminationController* termination_ctl_ptr;
	Covariance parcov;
	FileManager file_mgr;
	OutputFileWriter* out_wtr_ptr;
	ClpSimplex solve_lp_problem(Jacobian_1to1 &jco);

	//void initialize_obj_function_components();
	void initial_report();
	void constraint_report();
	void decision_var_report();
	void update(ClpSimplex &model);
	void update_decision_vars(ClpSimplex &model);
	void update_constraints(ClpSimplex &model);
	void separate_scenarios();
	void make_runs(Jacobian_1to1 &jco);
	CoinPackedMatrix jacobian_to_coinpackedmatrix(Jacobian_1to1 &jco);
	void build_constraint_bound_arrays();
	void throw_sequentialLP_error(string message);
	vector<double> get_constraint_residual_vec();
	void build_obj_func_coef_array();

};




#endif