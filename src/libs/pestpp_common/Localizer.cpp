#include <random>
#include <iomanip>
#include <unordered_set>
#include <iterator>
#include "Ensemble.h"
#include "RestartController.h"
#include "utilities.h"
#include "ParamTransformSeq.h"
#include "ObjectiveFunc.h"
#include "RedSVD-h.h"
#include "covariance.h"
#include "PerformanceLog.h"
#include "system_variables.h"
#include "Localizer.h"

bool Localizer::initialize(PerformanceLog *performance_log)
{
	stringstream ss;
	how == How::OBSERVATIONS; //set this for the case with no localization
	string filename = pest_scenario_ptr->get_pestpp_options().get_ies_localizer();
	if (filename.size() == 0)
	{
		use = false;
		return false;
	}
	use = true;
	
	mat.from_file(filename);
	
	string how_str = pest_scenario_ptr->get_pestpp_options().get_ies_localize_how();
	if (how_str[0] == 'P')
	{
		how = How::PARAMETERS;
	}
	else if (how_str[0] == 'O')
	{
		how = How::OBSERVATIONS;
	}
	else
	{
		throw runtime_error("Localizer.initialize(): 'ies_localize_how' must start with 'P' (pars) or 'O' (obs) not " + how_str[0]);
	}

	//error checking and building up container of names
	vector<string> names = pest_scenario_ptr->get_ctl_ordered_adj_par_names();
	set<string> par_names(names.begin(), names.end());
	names = pest_scenario_ptr->get_ctl_ordered_nz_obs_names();
	set<string> obs_names(names.begin(), names.end());

	map<string, vector<string>> pargp_map;
	ParameterGroupInfo *pi = pest_scenario_ptr->get_base_group_info_ptr();

	for (auto &pg : pest_scenario_ptr->get_ctl_ordered_par_group_names())
	{

		names.clear();
		for (auto &p : par_names)
			if (pi->get_group_name(p) == pg)
				names.push_back(p);
		pargp_map[pg] = names;
	}

	map<string, vector<string>> obgnme_map;
	for (auto &og : pest_scenario_ptr->get_ctl_ordered_obs_group_names())
	{
		ObservationInfo *oi = pest_scenario_ptr->get_observation_info_ptr();
		names.clear();
		for (auto &o : obs_names)
			if (oi->get_group(o) == og)
				names.push_back(o);
		obgnme_map[og] = names;
	}

	vector<string> missing, dups;
	vector<vector<string>> obs_map;
	set<string> dup_check;
	//for (auto &o : mat.get_row_names())
	string o;
	vector<string> row_names = mat.get_row_names();
	for (int i=0;i<mat.nrow();i++)
	{
		o = row_names[i];
		if (obs_names.find(o) != obs_names.end())
		{
			obs2row_map[o] = i;
			obs_map.push_back(vector<string>{o});
			if (dup_check.find(o) != dup_check.end())
				dups.push_back(o);
			dup_check.emplace(o);
		}
		else if (obgnme_map.find(o) != obgnme_map.end())
		{
			obs_map.push_back(obgnme_map[o]);
			if (obgnme_map[o].size() == 0)
				throw runtime_error("Localizer::initialize() error: listed observation group '" + o + "' has no non-zero weight observations");
			for (auto &oo : obgnme_map[o])
			{
				obs2row_map[oo] = i;
				if (dup_check.find(oo) != dup_check.end())
					dups.push_back(oo);
				dup_check.emplace(oo);

			}
		}
		else
			missing.push_back(o);
	}
	if (missing.size() > 0)
	{
		ss << " the following rows in " << filename << " were not found in the non-zero-weight observation names or observation group names: ";
		for (auto &m : missing)
			ss << m << ',';
		performance_log->log_event("error:" + ss.str());
		throw runtime_error(ss.str());
	}
	if (dups.size() > 0)
	{
		ss << " the following observations were listed more than once (possibly through an obs group): ";
		for (auto & d : dups)
			ss << d << ',';
		performance_log->log_event("error:" + ss.str());
		throw runtime_error(ss.str());
	}

	vector<vector<string>> par_map;
	dup_check.clear();
	string p;
	vector<string> col_names = mat.get_col_names();
	//for (auto &p : mat.get_col_names())
	for (int i=0;i<mat.ncol();++i)
	{
		p = col_names[i];
		if (par_names.find(p) != par_names.end())
		{
			par2col_map[p] = i;
			par_map.push_back(vector<string>{p});
			if (dup_check.find(p) != dup_check.end())
				dups.push_back(p);
			dup_check.emplace(p);
		}
		else if (pargp_map.find(p) != pargp_map.end())
		{
			par_map.push_back(pargp_map[p]);
			if (pargp_map[p].size() == 0)
				throw runtime_error("Localizer::initialize() error: listed parameter group '" + p + "' has no adjustable parameters");
			for (auto &pp : pargp_map[p])
			{
				par2col_map[pp] = i;
				if (dup_check.find(pp) != dup_check.end())
					dups.push_back(pp);
				dup_check.emplace(pp);

			}
		}
		else
			missing.push_back(p);
	}
	if (missing.size() > 0)
	{
		ss << " the following cols in " << filename << " were not found in the active parameter names or parameter group names: ";
		for (auto &m : missing)
			ss << m << ',';
		performance_log->log_event("error:" + ss.str());
		throw runtime_error(ss.str());
	}
	if (dups.size() > 0)
	{
		ss << " the following parameters were listed more than once (possibly through a par group): ";
		for (auto & d : dups)
			ss << d << ',';
		performance_log->log_event("error:" + ss.str());
		throw runtime_error(ss.str());
	}


	//map all the nz locations in the matrix
	map<int, vector<int>> idx_map;
	vector<string> vobs, vpar;

	//int i, j;
	if (how == How::PARAMETERS)
	{
		vector<string> col_names = mat.get_col_names();
		for (int k = 0; k < mat.e_ptr()->outerSize(); ++k)
		{
			for (Eigen::SparseMatrix<double>::InnerIterator it(*mat.e_ptr(), k); it; ++it)
			{
				if (idx_map.find(it.col()) == idx_map.end())
					idx_map[it.col()] = vector<int>{ it.row() };
				else
					idx_map[it.col()].push_back(it.row());
				//std::cout << "(" << it.row() << ","; // row index
				//std::cout << it.col() << ")\t"; // col index (here it is equal to k)

			}
		}
		//populate the localizer map
		for (auto &idx : idx_map)
		{
			vpar = par_map[idx.first];
			vobs.clear();
			for (auto &i : idx.second)
			{
				vobs.insert(vobs.end(), obs_map[i].begin(), obs_map[i].end());
			}
			pair<vector<string>, vector<string>> p(vobs, vpar);
			//localizer_map.push_back(p);
			localizer_map[col_names[idx.first]] = p;
		}

	}
	else
	{
		vector<string> row_names = mat.get_row_names();
		for (int k = 0; k < mat.e_ptr()->outerSize(); ++k)
		{
			for (Eigen::SparseMatrix<double>::InnerIterator it(*mat.e_ptr(), k); it; ++it)
			{
				if (idx_map.find(it.row()) == idx_map.end())
					idx_map[it.row()] = vector<int>{ it.col() };
				else
					idx_map[it.row()].push_back(it.col());
				//std::cout << "(" << it.row() << ","; // row index
				//std::cout << it.col() << ")\t"; // col index (here it is equal to k)

			}
		}
		//populate the localizer map
		for (auto &idx : idx_map)
		{
			vobs = obs_map[idx.first];
			vpar.clear();
			for (auto &i : idx.second)
			{
				vpar.insert(vpar.end(), par_map[i].begin(), par_map[i].end());
			}
			pair<vector<string>, vector<string>> p(vobs, vpar);
			//localizer_map.push_back(p);
			localizer_map[row_names[idx.first]] = p;
		}

	}

	
	//cout << "done" << endl;
	return true;
}

Eigen::MatrixXd Localizer::get_localizing_obs_hadamard_matrix(int num_reals, string col_name, vector<string> &obs_names)
{

	vector<double> values;
	vector<string> mat_cols = mat.get_col_names();
	vector<string>::iterator it = find(mat_cols.begin(), mat_cols.end(), col_name);
	if (it == mat_cols.end())
		throw runtime_error("Localizer::get_localizing_obs_hadamard_matrix() error: col_name not found: " + col_name);
	int idx = it - mat_cols.begin();
	Eigen::VectorXd mat_vec = mat.e_ptr()->col(idx);
	int col_idx;
	Eigen::MatrixXd loc(obs_names.size(), num_reals);
	for (int i=0;i<obs_names.size();i++)
	{
		col_idx = obs2row_map[obs_names[i]];
		loc.row(i).setConstant(mat_vec[col_idx]);

	}
	return loc;

}


Eigen::MatrixXd Localizer::get_localizing_par_hadamard_matrix(int num_reals, string row_name, vector<string> &par_names)
{

	vector<double> values;
	vector<string> mat_rows = mat.get_row_names();
	vector<string>::iterator it = find(mat_rows.begin(), mat_rows.end(), row_name);
	if (it == mat_rows.end())
		throw runtime_error("Localizer::get_localizing_par_hadamard_matrix() error: col_name not found: " + row_name);
	int idx = it - mat_rows.begin();
	Eigen::VectorXd mat_vec = mat.e_ptr()->row(idx);
	int col_idx;
	Eigen::MatrixXd loc(par_names.size(), num_reals);
	for (int i = 0; i<par_names.size(); i++)
	{
		col_idx = par2col_map[par_names[i]];
		loc.row(i).setConstant(mat_vec[col_idx]);
	}
	return loc;

}