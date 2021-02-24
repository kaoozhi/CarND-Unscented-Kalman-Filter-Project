#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd::Identity(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = M_PI/4.;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */

  
  is_initialized_ = false;
  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // predicted sigma points vector
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);


  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;
  
  // set augmented sigma points weights
  weights_ = VectorXd::Ones(2*n_aug_+1)*0.5/(lambda_ + n_aug_);
  weights_(0) = lambda_/(lambda_ + n_aug_);

 // Initialize state convariance matrix
  P_(0,0) = 1;
  P_(1,1) = 1;
  P_(2,2) = 10;
  P_(3,3) = 10;
  P_(4,4) = 10;

  // Initialize nis score
  nis_lidar_score_ = 0;
  nis_radar_score_ = 0;

}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */

    if (!is_initialized_) {
    /*
     Initialize the state x_ with the first measurement.
     */
    time_us_ = meas_package.timestamp_;

    // first measurement
    cout << "UKF: " << endl;

    // x_.fill(1.0);
    
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      // Convert radar from polar to cartesian coordinates 
      //         and initialize state.
        double rho = meas_package.raw_measurements_(0);
        double theta = meas_package.raw_measurements_(1);
        double rho_dot = meas_package.raw_measurements_(2);

        x_ << rho*cos(theta), 
        rho*sin(theta), 
        0, 
        0,
        0;

    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      // Initialize state using lidar measurements
        x_ << meas_package.raw_measurements_(0), 
        meas_package.raw_measurements_(1), 
        0, 
        0,
        0;
    }
    is_initialized_ = true;
    return;
  }

    double dt = (meas_package.timestamp_ - time_us_)/1000000.0; // time diff in sec      
    time_us_ = meas_package.timestamp_;

    Prediction(dt);
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar measurement update
       UpdateRadar(meas_package);
       // Nis monitoring
       cout<<"Radar NIS>7.8 percentage:"<<nis_radar_score_/double(nis_radar_.size())*100<<"%"<<endl;


    } 
    else {
    // Lidar measurement update
       UpdateLidar(meas_package);
    // Nis monitoring
       cout<<"Lidar NIS>7.8 percentage:"<<nis_lidar_score_/double(nis_lidar_.size())*100<<"%"<<endl;

    }
    

    // cout << "time = " << time_us_ << endl;

}

void UKF::Prediction(double dt) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */

  /**
   * Generation of augmented sigma points
   */
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_+1);
  GenerateSigmaPoints(&Xsig_aug);

  /**
   * Sigma points prediction
   */
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);
  for (int i=0; i<2*n_aug_+1; ++i)
  { 
    double v_k = Xsig_aug(2,i);
    double yaw_k = Xsig_aug(3,i);
    double yaw_dot_k = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);
    
    VectorXd F_k(n_x_);
    VectorXd W_k(n_x_);
    if (yaw_dot_k>0.001)
    {
      F_k<< v_k/yaw_dot_k*(sin(yaw_k+yaw_dot_k*dt)-sin(yaw_k)),
      v_k/yaw_dot_k*(-cos(yaw_k+yaw_dot_k*dt)+cos(yaw_k)),
      0,
      yaw_dot_k*dt,
      0;
    }
    else
    {
      F_k<< v_k*cos(yaw_k)*dt,
      v_k*sin(yaw_k)*dt,
      0,          
      yaw_dot_k*dt,
      0;
    }

    W_k<<0.5*dt*dt*cos(yaw_k)*nu_a,
    0.5*dt*dt*sin(yaw_k)*nu_a,
    dt*nu_a,
    0.5*dt*dt*nu_yawdd,
    dt*nu_yawdd;
    Xsig_pred.col(i) = Xsig_aug.col(i).head(n_x_)+F_k+W_k;
  }
  /**
   * States mean and covariance prediction
   */
  VectorXd x = VectorXd::Zero(n_x_);

  // create covariance matrix for prediction
  MatrixXd P = MatrixXd::Zero(n_x_, n_x_);

  // predict state mean
  for (int i=0; i<(2*n_aug_+1); ++i)
  {
      x = x+Xsig_pred.col(i)*weights_(i);
  }
  
  for (int i=0; i<(2*n_aug_+1); ++i)
  {
      VectorXd x_diff(n_x_);
      x_diff = Xsig_pred.col(i) - x;
      while (x_diff(3)> M_PI) {x_diff(3)-=2.*M_PI;}
      while (x_diff(3)<-M_PI) {x_diff(3)+=2.*M_PI;}
      P = P + weights_(i)*x_diff* x_diff.transpose();

  }
  x_=x;
  P_=P;
  Xsig_pred_ =Xsig_pred;
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  VectorXd z(2);
  z << meas_package.raw_measurements_(0),
  meas_package.raw_measurements_(1);

  MatrixXd R(2, 2);
  R << std_laspx_*std_laspx_, 0,
        0, std_laspy_*std_laspy_;
  MatrixXd H(2, 5);
  H<< 1, 0, 0, 0, 0,
      0, 1, 0, 0, 0;

  VectorXd z_pred = H * x_;
  VectorXd z_diff = z - z_pred;
  MatrixXd S = H * P_ * H.transpose() + R;
  MatrixXd K = P_ * H.transpose() * S.inverse();
  // Update state and covariance matrix
  x_ = x_ + (K * z_diff);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H) * P_;


  // Caculate NIS
  double nis;
  nis = z_diff.transpose()*S.inverse()*z_diff;

  if (nis>7.8) nis_lidar_score_++; 
  nis_lidar_.push_back(nis);
  
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
  int n_z = 3;

  /**
   * Predict radar measurement
   */

  // Measurement sigma points generation
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  // mean predicted measurement
  VectorXd z_pred = VectorXd::Zero(n_z);
  
  // measurement covariance matrix S
  MatrixXd S = MatrixXd::Zero(n_z,n_z);
  
  // measurement noise matrix R
  MatrixXd R = MatrixXd::Identity(n_z,n_z);
  R(0,0) = std_radr_*std_radr_;
  R(1,1) = std_radphi_* std_radphi_;
  R(2,2) = std_radrd_ * std_radrd_;

  for (int i=0; i<(2*n_aug_+1); ++i)
  {
      double px =  Xsig_pred_(0,i);
      double py =  Xsig_pred_(1,i);
      double v =  Xsig_pred_(2,i);
      double yaw =  Xsig_pred_(3,i);
    //   double psi_dot=  Xsig_pred(4,i);
      double c1 = sqrt(px*px+py*py);
      if (c1<0.001) {c1 = 0.001;}
      // VectorXd Zsig_k(3);
      Zsig.col(i)<< c1,
                    atan2(py, px),
                    (px*cos(yaw)*v + py*sin(yaw)*v)/c1;
  }

  for (int i=0; i<(2*n_aug_+1); ++i)
  {
      z_pred += Zsig.col(i)*weights_(i);
  }
  
  for (int i=0; i<(2*n_aug_+1); ++i)
  {
      VectorXd z_diff(n_z);
      z_diff = Zsig.col(i) - z_pred;
      // while(z_diff(1) > M_PI) z_diff(1) -= 2*M_PI;
      // while(z_diff(1) < -M_PI) z_diff(1) += 2*M_PI;
      S = S + weights_(i)*z_diff*z_diff.transpose();
  }
  
  S += R;
  while (S(1,1) > M_PI) {S(1,1) -= 2.*M_PI;}
  while (S(1,1) < -M_PI) {S(1,1) += 2.*M_PI;}


  /**
   * Update state and covariance matrix
   */

  VectorXd z(3);
  z << meas_package.raw_measurements_(0),
  meas_package.raw_measurements_(1),
  meas_package.raw_measurements_(2);

  MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);
  // calculate cross correlation matrix
  for (int i=0; i<(2*n_aug_ +1); ++i)
  {
    VectorXd x_diff(n_x_);
    x_diff = Xsig_pred_.col(i) - x_ ;
    if (x_diff(3)> M_PI) {x_diff(3)-=2.*M_PI;}
    else if (x_diff(3)<-M_PI) {x_diff(3)+=2.*M_PI;}

    VectorXd z_diff(n_z);
    z_diff = Zsig.col(i) - z_pred;
    if (z_diff(1)> M_PI) {z_diff(1)-=2.*M_PI;}
    else if (z_diff(1)<-M_PI) {z_diff(1)+=2.*M_PI;}
    Tc += weights_(i)*x_diff*z_diff.transpose();
  }
  // calculate Kalman gain;
  MatrixXd K = MatrixXd(n_x_, n_z);
  K = Tc*S.inverse();
  
  // update state mean and covariance matrix
  VectorXd z_err(n_z);
  z_err = (z - z_pred);
  if(z_err(1)> M_PI) {z_err(1) -= 2.*M_PI;}
  else if(z_err(1)< -M_PI) {z_err(1) += 2.*M_PI;}
  
  x_ = x_ + K*z_err;
  P_ = P_ - K*S*K.transpose();

  // Caculate NIS
  double nis;
  nis = z_err.transpose()*S.inverse()*z_err;

  if (nis>7.8) nis_radar_score_++; 
  nis_radar_.push_back(nis);

}


void UKF::GenerateSigmaPoints(Eigen::MatrixXd* Xsig_out){

  // Generate augmented sigmapoints
  // create augmented mean vector
  VectorXd x_aug = VectorXd::Zero(n_aug_);

  // create augmented state covariance
  MatrixXd P_aug = MatrixXd::Identity(n_aug_, n_aug_);
  
  // create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
 
  // create augmented mean state
  x_aug.head(n_x_) = x_;

  // create augmented covariance matrix
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) = std_a_*std_a_;
  P_aug(n_x_+1, n_x_+1) = std_yawdd_*std_yawdd_;
  // create square root matrix
  
  MatrixXd A_aug = P_aug.llt().matrixL();

  // create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  
  for (int i = 0; i<n_aug_; ++i)
  {
     Xsig_aug.col(i+1) = x_aug + sqrt(lambda_+n_aug_)*A_aug.col(i);
     Xsig_aug.col(i+n_aug_+1) = x_aug - sqrt(lambda_+n_aug_)*A_aug.col(i);

  }
  *Xsig_out = Xsig_aug;

}
