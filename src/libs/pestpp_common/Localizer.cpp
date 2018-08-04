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
	string filename = pest_scenario_ptr->get_pestpp_options().get_ies_localizer();
	if (filename.size() == 0)
		return false;
	/*string ext = filename.substr(filename.size() - 3, 3);
	pest_utils::upper_ip(ext);
	if ((ext == "JCB") || (ext == "JCO"))
	{
		performance_log->log_event("loading localizer from binary file");
		mat.from_binary(filename);
	}
	else if (ext == "MAT")
	{
		performance_log->log_event("loading localizer from binary file");
		mat.from_ascii(filename);
	}
	else if (ext == "CSV")
	{
		performance_log->log_event("loading localizer from csv file");
		mat.from_csv(filename);
	}
	else
	{
		ss << "unrecognnized localizer extension '" << ext << "', should be JCB, JCO, or MAT";
		performance_log->log_event("error: "+ss.str());
		throw runtime_error(ss.str());
	}*/
	mat.from_file(filename);


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
	for (auto &p : mat.get_col_names())
	{
		if (par_names.find(p) != par_names.end())
		{
			par_map.push_back(vector<string>{p});
			if (dup_check.find(p) != dup_check.end())
				dups.push_back(p);
			dup_check.emplace(p);
		}
		else if (pargp_map.find(p) != pargp_map.end())
		{
			par_map.push_back(pargp_map[p]);

			for (auto &pp : pargp_map[p])
			{
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
	int i, j;
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
	vector<string> vobs, vpar;
	vector<string> col_names = mat.get_col_names();

	for (auto &idx : idx_map)
	{
		vpar = par_map[idx.first];
		vobs.clear();
		for (auto &i : idx.second)
		{
			vobs.insert(vobs.end(), obs_map[i].begin(), obs_map[i].end());
		}
		pair<vector<string>, vector<string>> p(vobs,vpar);
		//localizer_map.push_back(p);
		localizer_map[col_names[idx.first]] = p;
	}

	//cout << "done" << endl;
}

Eigen::MatrixXd Localizer::get_localizing_hadamard_matrix(int num_reals, string col_name, vector<string> &obs_names)
{

	vector<double> values;
	vector<string> mat_cols = mat.get_col_names();
	vector<string>::iterator it = find(mat_cols.begin(), mat_cols.end(), col_name);
	if (it == mat_cols.end())
		throw runtime_error("Localizer::get_localizing_vector() error: row_name not found: " + col_name);
	int idx = it - mat_cols.begin();
	Eigen::VectorXd mat_vec = mat.e_ptr()->col(idx);
	vector<double> full_vec;
	int col_idx;
	Eigen::MatrixXd loc(obs_names.size(), num_reals);
	//loc.setOnes();
	//for (auto &name : obs_names)
	for (int i=0;i<obs_names.size();i++)
	{
		col_idx = obs2row_map[obs_names[i]];
		loc.row(i).setConstant(mat_vec[col_idx]);

	}
	//cout << loc << endl;
	return loc;

}
