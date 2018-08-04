#ifndef LOCALIZER_H_
#define LOCALIZER_H_

#include <map>
#include <random>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "FileManager.h"
#include "ObjectiveFunc.h"
#include "OutputFileWriter.h"
#include "PerformanceLog.h"
#include "RunStorage.h"
#include "covariance.h"
#include "RunManagerAbstract.h"
#include "PerformanceLog.h"



class Localizer
{
public:
	Localizer() { ; }
	Localizer(Pest *_pest_scenario_ptr) { pest_scenario_ptr = _pest_scenario_ptr; }
	bool initialize(PerformanceLog *performance_log);
	const map<string,pair<vector<string>, vector<string>>> get_localizer_map() { return localizer_map; }
	void set_pest_scenario(Pest *_pest_scenario_ptr) { pest_scenario_ptr = _pest_scenario_ptr; }
	Eigen::MatrixXd get_localizing_hadamard_matrix(int num_reals,string row_name,vector<string> &obs_names);

private:
	Pest * pest_scenario_ptr;
	Mat mat;
	map<string,pair<vector<string>, vector<string>>> localizer_map;
	map<string, int> obs2row_map;
};

#endif
