#include "tools.h"

using Eigen::VectorXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */

VectorXd rmse = VectorXd::Zero(4);
//   rmse << 0.0,0.0,0.0,0.0;
  for (int i=0; i<estimations.size(); ++i){
     VectorXd error;
     error = estimations[i] - ground_truth[i];
     
     error  = error.array()*error.array();

     rmse += error;
  }

    rmse = rmse/estimations.size();
    rmse = rmse.array().sqrt();
    return rmse;
}