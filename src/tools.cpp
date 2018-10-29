#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

    VectorXd rmse(4);
    rmse << 0,0,0,0;

      // TODO: YOUR CODE HERE

    // check the validity of the following inputs:
    //  * the estimation vector size should not be zero
    //  * the estimation vector size should equal ground truth vector size
    // ... your code here
    if (estimations.size() == 0 || estimations.size() != ground_truth.size()) {
        return rmse;
    }

    //accumulate squared residuals
    std::vector<VectorXd> residuals(estimations.size());
    for(int i=0; i < estimations.size(); ++i){
          VectorXd residual = estimations[i] - ground_truth[i];
          residual = residual.array()*residual.array();
          rmse += residual;
    }

    //calculate the mean
    rmse = rmse/estimations.size();

    //calculate the squared root
    rmse = rmse.array().sqrt();

    //return the result
    return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

  MatrixXd Hj(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//TODO: YOUR CODE HERE
  float squared_pos = px*px + py*py;
  float squared_pos_sqrt = sqrt(squared_pos);
  float squared_pos_pow3_2 = squared_pos * squared_pos_sqrt;

	//check division by zero
	if (fabs(squared_pos) < 0.001) {
	    cout << "Error: division by zero!" << endl;
      return Hj;
	}
	else {

	    float dro1 = px / squared_pos_sqrt;
	    float dro2 = py / squared_pos_sqrt;
	    float dro3 = 0;
	    float dro4 = 0;
	    float dfi1 = -py/squared_pos;
	    float dfi2 = px/squared_pos;
	    float dfi3 = 0;
	    float dfi4 = 0;
	    float drod1 = py*(vx*py-vy*px)/squared_pos_pow3_2;
	    float drod2 = px*(vy*px-vx*py)/squared_pos_pow3_2;
	    float drod3 = px/squared_pos_sqrt;
	    float drod4 = py/squared_pos_sqrt;

	    Hj << dro1, dro2, dro3, dro4,
	          dfi1, dfi2, dfi3, dfi4,
	          drod1, drod2, drod3, drod4;
	}

	return Hj;
}
