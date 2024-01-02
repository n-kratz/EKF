#include "ekf_models.hpp"
#include <tf/tf.h>
#include "utilities.h"

/**
   TODO
   Fill in the value of the process covariance matrix. The rows/columns of WMWt are
   in the following order [POS_X POS_Y POS_Z ROT_R ROT_P ROT_Y ].
   \param[out] WMWt Covariance matrix of the system.
   \param state_in    The current state estimate
   \param v           The input linear velocity
   \param w           The input angular velocity
   \param dt          Delta time
*/
void sys_evaluate_WMWt( double WMWt[6][6], const State& state_in, double v, double w, double dt ){

  for( int r=0; r<6; r++ )
    for( int c=0; c<6; c++ )
      WMWt[r][c] = 0.0;

  // TODO fill in the matrix WMWt
  double x = state_in.x[0];
  double y = state_in.x[1];
  double z = state_in.x[2];
  double r = state_in.x[3];
  double p = state_in.x[4];
  double yaw = state_in.x[5];

  double a1 = 0.01;
  double a2 = 0.05;
  double a3 = 0.01;
  double a4 = 0.05;

  WMWt[0][0] = (dt*dt)*pow(cos(p),2.0)*pow(cos(yaw),2.0)*(a1*(v*v)+a2*(w*w));
  WMWt[0][1] = (dt*dt)*cos(p)*cos(yaw)*(cos(r)*sin(yaw)+cos(yaw)*sin(p)*sin(r))*(a1*(v*v)+a2*(w*w));
  WMWt[0][2] = (dt*dt)*cos(p)*cos(yaw)*(sin(r)*sin(yaw)-cos(r)*cos(yaw)*sin(p))*(a1*(v*v)+a2*(w*w));
  WMWt[1][0] = (dt*dt)*cos(p)*cos(yaw)*(cos(r)*sin(yaw)+cos(yaw)*sin(p)*sin(r))*(a1*(v*v)+a2*(w*w));
  WMWt[1][1] = (dt*dt)*pow(cos(r)*sin(yaw)+cos(yaw)*sin(p)*sin(r),2.0)*(a1*(v*v)+a2*(w*w));
  WMWt[1][2] = (dt*dt)*(cos(r)*sin(yaw)+cos(yaw)*sin(p)*sin(r))*(sin(r)*sin(yaw)-cos(r)*cos(yaw)*sin(p))*(a1*(v*v)+a2*(w*w));
  WMWt[2][0] = (dt*dt)*cos(p)*cos(yaw)*(sin(r)*sin(yaw)-cos(r)*cos(yaw)*sin(p))*(a1*(v*v)+a2*(w*w));
  WMWt[2][1] = (dt*dt)*(cos(r)*sin(yaw)+cos(yaw)*sin(p)*sin(r))*(sin(r)*sin(yaw)-cos(r)*cos(yaw)*sin(p))*(a1*(v*v)+a2*(w*w));
  WMWt[2][2] = (dt*dt)*pow(sin(r)*sin(yaw)-cos(r)*cos(yaw)*sin(p),2.0)*(a1*(v*v)+a2*(w*w));
  WMWt[3][3] = .001;
  WMWt[4][4] = .001;
  WMWt[5][5] = (dt*dt)*(a3*(v*v)+a4*(w*w));
}

/**
   TODO
   Fill in the value of the measurement covariance matrix. The rows/columns of C
   are in the following order [POS_X POS_Y POS_Z ROT_R ROT_P ROT_Y ]
   \param[out] R Covariance matrix of the sensors.
   \param state_in    The current state estimate
*/
void meas_evaluate_R( double R[6][6], const State& state ){

  for( int r=0; r<6; r++ )
    for( int c=0; c<6; c++ )
      R[r][c] = 0.0;

  // TODO fill in the matrix R
  R[0][0] = 0.1;
  R[1][1] = 0.1;
  R[2][2] = 0.1;
  R[3][3] = 0.001;
  R[4][4] = 0.001;
  R[5][5] = 0.001;
}


/**
   TODO
   Evaluate the system function.
   Compute the process model.
   This function returns the prediction of the next state based on the 
   current state estimate and the commmand input (linear/angular velocities).
   \param state_in    The current state estimate
   \param v           The input linear velocity
   \param w           The input angular velocity
   \param dt          Delta time
*/
State sys_evaluate_g( const State& state_in, double v, double w, double dt ){

  State state_out;
  double x = state_in.x[0];
  double y = state_in.x[1];
  double z = state_in.x[2];
  double r = state_in.x[3];
  double p = state_in.x[4];
  double yaw = state_in.x[5];

  // TODO Given state_in and v and w and dt (time increment) determine the prior
  // estimate state_out

  state_out.x[0] = x+dt*v*cos(p)*cos(yaw);
  state_out.x[1] = y+dt*v*(cos(r)*sin(yaw)+cos(yaw)*sin(p)*sin(r));
  state_out.x[2] = z+dt*v*(sin(r)*sin(yaw)-cos(r)*cos(yaw)*sin(p));
  state_out.x[3] = r;
  state_out.x[4] = p;
  state_out.x[5] = yaw+dt*w;
  return state_out;
}

/**
   TODO
   Evaluate the system Jacobian.
   This function evaluates the Jacobian of the system functions g (see 
   sys_evaluate_g). The entry G[i][j] represents ( d g_i / d s_j )
   \param[out] G      The 6x6 Jacobian of the function g
   \param state       The state of the robot
   \param v           The input linear velocity
   \param w           The input angular velocity
   \param dt          Delta time
*/
void sys_evaluate_G( double G[6][6], const State& state_in, double v, double w, double dt ){
    double x = state_in.x[0];
    double y = state_in.x[1];
    double z = state_in.x[2];
    double r = state_in.x[3];
    double p = state_in.x[4];
    double yaw = state_in.x[5];

  for( int r=0; r<6; r++ )
    for( int c=0; c<6; c++ )
      G[r][c] = 0.0;
  
  // TODO
  // Given state, v and w, compute the system Jacobian G
  G[0][0] = 1.0;
  G[0][4] = -dt*v*cos(yaw)*sin(p);
  G[0][5] = -dt*v*cos(p)*sin(yaw);
  G[1][1] = 1.0;
  G[1][3] = -dt*v*(sin(r)*sin(yaw)-cos(r)*cos(yaw)*sin(p));
  G[1][4] = dt*v*cos(p)*cos(yaw)*sin(r);
  G[1][5] = dt*v*(cos(r)*cos(yaw)-sin(p)*sin(r)*sin(yaw));
  G[2][2] = 1.0;
  G[2][3] = dt*v*(cos(r)*sin(yaw)+cos(yaw)*sin(p)*sin(r));
  G[2][4] = -dt*v*cos(p)*cos(r)*cos(yaw);
  G[2][5] = dt*v*(cos(yaw)*sin(r)+cos(r)*sin(p)*sin(yaw));
  G[3][3] = 1.0;
  G[4][4] = 1.0;
  G[5][5] = 1.0;
}

/**
   TODO
   Evaluate the GPS observation function.
   This function returns the expected satellite fix given the state of the robot
   \param state The state estimate
   \return      A satellite navigation fix (only the latitute, longitude
                and altitude members are used)
*/
sensor_msgs::NavSatFix meas_evaluate_gps( const State& state ){

  sensor_msgs::NavSatFix nsf;

  // TODO
  // Given prior estimate state, determine the expected GPS measurement nsf
  double x = state.x[0];
  double y = state.x[1];
  double z = state.x[2];

  nsf.latitude = (M_PI*3.5859457E+1+cos(M_PI*(3.7E+1/9.0E+1))*(z-9.0/5.0E+1)*2.822140697197316E-5+sin(M_PI*(3.7E+1/9.0E+1))*(y+6.3E+1/5.0E+2)*2.822140697197316E-5)/M_PI;
  nsf.longitude = -(M_PI*1.0823684E+2-((cos(M_PI*(3.7E+1/9.0E+1))*(y+6.3E+1/5.0E+2)-sin(M_PI*(3.7E+1/9.0E+1))*(z-9.0/5.0E+1))*2.822140697197316E-5)/cos(M_PI*1.992192055555556E-1))/M_PI;
  nsf.altitude = 0.999*z + 14;

  return nsf;
}

/**
   TODO
   Evaluate the IMU observation function.
   This function computes the expected imu orientation given the state of the 
   robot.
   \param state_in The current state estimate
   \return         A inertial navigation unit measurement (only the orientation
                   member is used).
*/
sensor_msgs::RPY meas_evaluate_imu( const State& state ){
  sensor_msgs::RPY rpy;

  // TODO
  // Given the prior estimate state, determine the expected RPY measurement rpy
  double r = state.x[3];
  double p = state.x[4];
  double yaw = state.x[5];

  rpy.roll = r;
  rpy.pitch = p;
  rpy.yaw = yaw;
  return rpy;
}

/** 
    TODO
    Observation Jacobian of the GPS
    This function returns the 3x3 observation Jacobian of the GPS. Essentially,
    this is the Jacobian of your meas_evaluate_gps function.
    \param[out] Hgps The 3x3 GPS Jacobian.
    \param[in]  state The state of the robot
*/
void meas_evaluate_Hgps( double Hgps[3][3], const State& state ){
    for( int r=0; r<3; r++ )
      for( int c=0; c<3; c++ )
        Hgps[r][c] = 0.0;

  // TODO
  // Fill the Jacobian matrix Hgps of the GPS observations

    Hgps[0][1] = (sin(M_PI*(3.7E+1/9.0E+1))*2.822140697197316E-5)/M_PI;
    Hgps[0][2] = (cos(M_PI*(3.7E+1/9.0E+1))*2.822140697197316E-5)/M_PI;
    Hgps[1][1] = (cos(M_PI*(3.7E+1/9.0E+1))*2.822140697197316E-5)/(M_PI*cos(M_PI*1.992192055555556E-1));
    Hgps[1][2] = (sin(M_PI*(3.7E+1/9.0E+1))*(-2.822140697197316E-5))/(M_PI*cos(M_PI*1.992192055555556E-1));
    Hgps[2][2] = 0.999;

}

/** 
    Observation Jacobian of the IMU
    This function returns the 3x3 observation Jacobian of the IMU. Essentially,
    this is the Jacobian of your meas_evaluate_imu function.
    \param[out] Himu The 3x3 IMU Jacobian.
    \param[in]  state The state of the robot
*/
void meas_evaluate_Himu( double Himu[3][3], const State& state ){
    for( int r=0; r<3; r++ )
      for( int c=0; c<3; c++ )
        Himu[r][c] = 0.0;

  // TODO
  // Fill the Jacobian matrix Himu of the IMU observations
    Himu[0][0] = 1;
    Himu[1][1] = 1;
    Himu[2][2] = 1;
}

