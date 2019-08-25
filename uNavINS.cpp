/*
uNavINS.cpp

Original Author:
Adhika Lie
2012-10-08
University of Minnesota
Aerospace Engineering and Mechanics
Copyright 2011 Regents of the University of Minnesota. All rights reserved.

Updated to be a class, use Eigen, and compile as an Arduino library.
Added methods to get gyro and accel bias. Added initialization to
estimated angles rather than assuming IMU is level. Added method to get psi,
rather than just heading, and ground track.
Brian R Taylor
brian.taylor@bolderflight.com
2017-12-20
Bolder Flight Systems
Copyright 2017 Bolder Flight Systems

Permission is hereby granted, free of charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include "uNavINS.h"

void uNavINS::update(double TOW,double vn,double ve,double vd,double lat,double lon,double alt,float p,float q,float r,float ax,float ay,float az,float hx,float hy, float hz) {
  //TODO proceed to init in another function?
  if (!initialized) {
    // initial attitude and heading
    theta = asinf(-ax/G);
    phi = asinf(ay/(G*cosf(theta)));
    // magnetic heading correction due to roll and pitch angle
    Bxc = hx*cosf(theta) + (hy*sinf(phi) + hz*cosf(phi))*sinf(theta);
    Byc = hy*cosf(phi) - hz*sinf(phi);
    // finding initial heading
    if (-Byc > 0) {
      psi = M_PI/2.0f - atanf(Bxc/-Byc);
    } else {
      psi= 3.0f*M_PI/2.0f - atanf(Bxc/-Byc);
    }

psi = 0;

    psi = constrainAngle180(psi);
    psi_initial = psi;
    // euler to quaternion
    quat(0) = cosf(psi/2.0f)*cosf(theta/2.0f)*cosf(phi/2.0f) + sinf(psi/2.0f)*sinf(theta/2.0f)*sinf(phi/2.0f);
    quat(1) = cosf(psi/2.0f)*cosf(theta/2.0f)*sinf(phi/2.0f) - sinf(psi/2.0f)*sinf(theta/2.0f)*cosf(phi/2.0f);
    quat(2) = cosf(psi/2.0f)*sinf(theta/2.0f)*cosf(phi/2.0f) + sinf(psi/2.0f)*cosf(theta/2.0f)*sinf(phi/2.0f);
    quat(3) = sinf(psi/2.0f)*cosf(theta/2.0f)*cosf(phi/2.0f) - cosf(psi/2.0f)*sinf(theta/2.0f)*sinf(phi/2.0f);
    // Assemble the matrices
    // ... gravity
    grav(2,0) = G;
    // ... H
    H.block(0,0,5,5) = Eigen::Matrix<float,5,5>::Identity();
    // // ... Rw
    // Rw.block(0,0,3,3) = powf(SIG_W_A,2.0f)*Eigen::Matrix<float,3,3>::Identity();
    // Rw.block(3,3,3,3) = powf(SIG_W_G,2.0f)*Eigen::Matrix<float,3,3>::Identity();
    // Rw.block(6,6,3,3) = 2.0f*powf(SIG_A_D,2.0f)/TAU_A*Eigen::Matrix<float,3,3>::Identity();
    // Rw.block(9,9,3,3) = 2.0f*powf(SIG_G_D,2.0f)/TAU_G*Eigen::Matrix<float,3,3>::Identity();
    // // ... P
    // P.block(0,0,3,3) = powf(P_P_INIT,2.0f)*Eigen::Matrix<float,3,3>::Identity();
    // P.block(3,3,3,3) = powf(P_V_INIT,2.0f)*Eigen::Matrix<float,3,3>::Identity();
    // P.block(6,6,2,2) = powf(P_A_INIT,2.0f)*Eigen::Matrix<float,2,2>::Identity();
    // P(8,8) = powf(P_HDG_INIT,2.0f);
    // P.block(9,9,3,3) = powf(P_AB_INIT,2.0f)*Eigen::Matrix<float,3,3>::Identity();
    // P.block(12,12,3,3) = powf(P_GB_INIT,2.0f)*Eigen::Matrix<float,3,3>::Identity();
    // // ... R
    // R.block(0,0,2,2) = powf(SIG_GPS_P_NE,2.0f)*Eigen::Matrix<float,2,2>::Identity();
    // R(2,2) = powf(SIG_GPS_P_D,2.0f);
    // R.block(3,3,3,3) = powf(SIG_GPS_V,2.0f)*Eigen::Matrix<float,3,3>::Identity();
    // ... Rw
    Rw(0,0) = SIG_W_A*SIG_W_A;          Rw(1,1) = SIG_W_A*SIG_W_A;            Rw(2,2) = SIG_W_A*SIG_W_A;
    Rw(3,3) = SIG_W_G*SIG_W_G;          Rw(4,4) = SIG_W_G*SIG_W_G;            Rw(5,5) = SIG_W_G*SIG_W_G;
    Rw(6,6) = 2.0f*SIG_A_D*SIG_A_D/TAU_A; Rw(7,7) = 2.0f*SIG_A_D*SIG_A_D/TAU_A;   Rw(8,8) = 2.0f*SIG_A_D*SIG_A_D/TAU_A;
    Rw(9,9) = 2.0f*SIG_G_D*SIG_G_D/TAU_G; Rw(10,10) = 2.0f*SIG_G_D*SIG_G_D/TAU_G; Rw(11,11) = 2.0f*SIG_G_D*SIG_G_D/TAU_G;
    // ... P
    P(0,0) = P_P_INIT*P_P_INIT;       P(1,1) = P_P_INIT*P_P_INIT;       P(2,2) = P_P_INIT*P_P_INIT;
    P(3,3) = P_V_INIT*P_V_INIT;       P(4,4) = P_V_INIT*P_V_INIT;       P(5,5) = P_V_INIT*P_V_INIT;
    P(6,6) = P_A_INIT*P_A_INIT;       P(7,7) = P_A_INIT*P_A_INIT;       P(8,8) = P_HDG_INIT*P_HDG_INIT;
    P(9,9) = P_AB_INIT*P_AB_INIT;     P(10,10) = P_AB_INIT*P_AB_INIT;   P(11,11) = P_AB_INIT*P_AB_INIT;
    P(12,12) = P_GB_INIT*P_GB_INIT;   P(13,13) = P_GB_INIT*P_GB_INIT;   P(14,14) = P_GB_INIT*P_GB_INIT;
    // ... R
    R(0,0) = SIG_GPS_P_NE*SIG_GPS_P_NE; R(1,1) = SIG_GPS_P_NE*SIG_GPS_P_NE; R(2,2) = SIG_GPS_P_D*SIG_GPS_P_D;
    R(3,3) = SIG_GPS_V*SIG_GPS_V;       R(4,4) = SIG_GPS_V*SIG_GPS_V;       R(5,5) = SIG_GPS_V*SIG_GPS_V;
    // .. then initialize states with GPS Data
    lat_ins = lat;
    lon_ins = lon;
    alt_ins = alt;
    vn_ins = vn;
    ve_ins = ve;
    vd_ins = vd;
    // specific force -> unbiased measurement
    f_b(0,0) = ax;
    f_b(1,0) = ay;
    f_b(2,0) = az;
    /* initialize the time */
    _t = TOW;
    // initialized flag
    initialized = true;
  } else {
    // get the change in time
    _dt = (TOW-_t);
    _t = TOW;
    // _dt = (float)_t/1000.0;
    // _t = 0;

    if ( _dt < 0) {
      printf("Warning !! dt<0...\n\n",_dt);
      return;
    }

#ifdef DEBUG
    printf("%12.10lf,%12.10lf - %12.10lf,%12.10lf  -  %12.10lf,%12.10lf\n",_t,_dt,lat,lat_ins,lon,lon_ins);
#endif  

    // Keep in memory the velocity and ll before update
    lla_ins(0,0) = lat_ins; 
    lla_ins(1,0) = lon_ins;
    lla_ins(2,0) = alt_ins;
    V_ins(0,0) = vn_ins; 
    V_ins(1,0) = ve_ins;
    V_ins(2,0) = vd_ins;

// ----------------------- PROPAGATION = PREDICTION -----------------------------------

    // AHRS Transformations
    // compute this guy using the PREVIOUS quat... why? -> needed for Jacobian = linearization. ok
    C_N2B = quat2dcm(quat);
    C_B2N = C_N2B.transpose();
    
    // Time update of the Attitude p,q,r
    // = time integration of the quaternion
    dq(0) = 1.0f; 
    dq(1) = 0.5f*om_ib(0,0)*_dt;
    dq(2) = 0.5f*om_ib(1,0)*_dt;
    dq(3) = 0.5f*om_ib(2,0)*_dt;
    printf("qd: %f,%f,%f,%f  \n ",dq(0),dq(1,0),dq(2,0),dq(3,0));
    printf("Q: %f,%f,%f,%f \n  ",quat(0,0),quat(1,0),quat(2,0),quat(3,0));
    quat = qmult(quat,dq);
    
    quat.normalize();
    // Avoid quaternion flips sign
    if (quat(0) < 0) {
      quat = -1.0f*quat;
    }
    //quat can now be used to retreive updated Euler angles

    // obtain euler angles from quaternion
    theta = asinf(-2.0f*(quat(1,0)*quat(3,0)-quat(0,0)*quat(2,0)));
    phi = atan2f(2.0f*(quat(0,0)*quat(1,0)+quat(2,0)*quat(3,0)),1.0f-2.0f*(quat(1,0)*quat(1,0)+quat(2,0)*quat(2,0)));
    psi = atan2f(2.0f*(quat(1,0)*quat(2,0)+quat(0,0)*quat(3,0)),1.0f-2.0f*(quat(2,0)*quat(2,0)+quat(3,0)*quat(3,0)));



#ifdef DEBUG    
    printf("Q: %f,%f,%f,%f \n",quat(0,0),quat(1,0),quat(2,0),quat(3,0));
    printf("angles: %f,%f,%f\n",theta,phi,psi);
#endif
    // Acceleration time update:
    // -> assuming ddt = 0, so do nothing. 
    // f_b is the acceleration from previous time step in BODY FRAME !!!
    // -> C_B2N changes frame between PREVIOUS body frame and earth frame.
    // NOTE TODO: we could actually rotate it, according to our new Euler angles...
    // ... It might be a good idea to assume that the acceleration 
    //    (i.e. the forces) remain ~constant in the body axes !
    // -> 'just' need to compute an update of C_B2N... see the impact on performances?

    // Velocity time Update
    dx = C_B2N*f_b + grav;   //DG: deleting the gravity vector from f_b !!

printf("dx/fb: %f,%f,%f,%f,%f,%f,%f,%f,%f \n ",dx(0,0),dx(1,0),dx(2,0),f_b(0,0),f_b(1,0),f_b(2,0),ax,ay,az);

    vn_ins += _dt*dx(0,0);
    ve_ins += _dt*dx(1,0);
    vd_ins += _dt*dx(2,0);

#ifdef DEBUG
  printf("dx/fb: %f,%f,%f,%f,%f,%f \n ",dx(0,0),dx(1,0),dx(2,0),f_b(0,0),f_b(1,0),f_b(2,0));
#endif      

    // Position time Update
    dxd = llarate(V_ins,lla_ins); //RIGHT! take the velocity before update to remain consistent
    lat_ins += _dt*dxd(0,0);
    lon_ins += _dt*dxd(1,0);
    alt_ins += _dt*dxd(2,0);

    // ----------------------------Jacobian-------------------------------------------
    // We did the time update of the states, and we linearize the Transformations
    // to obtain the Jacobian that we feed in Kalman stuff
    Fs.setZero();
    // ... pos2gs
    Fs.block(0,3,3,3) = Eigen::Matrix<float,3,3>::Identity();
    // ... gs2pos
    Fs(5,2) = -2.0f*G/EARTH_RADIUS;
    // ... gs2att
    Fs.block(3,6,3,3) = -2.0f*C_B2N*sk(f_b);
    // ... gs2acc
    Fs.block(3,9,3,3) = -C_B2N;
    // ... att2att
    Fs.block(6,6,3,3) = -sk(om_ib);
    // ... att2gyr
    Fs.block(6,12,3,3) = -0.5f*Eigen::Matrix<float,3,3>::Identity();
    // ... Accel Markov Bias
    Fs.block(9,9,3,3) = -1.0f/TAU_A*Eigen::Matrix<float,3,3>::Identity();
    Fs.block(12,12,3,3) = -1.0f/TAU_G*Eigen::Matrix<float,3,3>::Identity();

    // State Transition Matrix
    PHI = Eigen::Matrix<float,15,15>::Identity()+Fs*_dt;
    
//DG: I don't understand here: 
// Fs is the Jacobian, so I agree that [x_t+1] = [x_t] + Fs * [x_t] = PHI * [x_t]
// but the covariance is supposed to be obtained with Fs, not with PHI???

    // Process Noise
    Gs.setZero();
    Gs.block(3,0,3,3) = -C_B2N;
    Gs.block(6,3,3,3) = -0.5f*Eigen::Matrix<float,3,3>::Identity();
    Gs.block(9,6,6,6) = Eigen::Matrix<float,6,6>::Identity();

    // Discrete Process Noise
    Q = PHI*_dt*Gs*Rw*Gs.transpose();
    Q = 0.5f*(Q+Q.transpose());

    // Covariance Time Update
    P = PHI*P*PHI.transpose()+Q;
    P = 0.5f*(P+P.transpose());

// ----------------------- CORRECTION = KALMAN -----------------------------------

    // if ((TOW - previousTOW) > 0) {
      previousTOW = TOW;
      lla_gps(0,0) = lat;
      lla_gps(1,0) = lon;
      lla_gps(2,0) = alt;
      V_gps(0,0) = vn;
      V_gps(1,0) = ve;
      V_gps(2,0) = vd;
      lla_ins(0,0) = lat_ins;
      lla_ins(1,0) = lon_ins;
      lla_ins(2,0) = alt_ins;
      V_ins(0,0) = vn_ins;
      V_ins(1,0) = ve_ins;
      V_ins(2,0) = vd_ins;

#ifdef DEBUG
  printf("Vs: %lf,%lf,%lf,%lf,%lf,%lf  ",V_gps(0,0),V_gps(1,0),V_gps(2,0),V_ins(0,0),V_ins(1,0),V_ins(2,0));
#endif      
      // Position, converted to NED
      pos_ecef_ins = lla2ecef(lla_ins);
      pos_ned_ins = ecef2ned(pos_ecef_ins,lla_ins);
      pos_ecef_gps = lla2ecef(lla_gps);
      pos_ned_gps = ecef2ned(pos_ecef_gps,lla_ins);
      // Create measurement Y = actual meas - projected meas (not the opposite??)
      y(0,0) = (float)(pos_ned_gps(0,0) - pos_ned_ins(0,0));
      y(1,0) = (float)(pos_ned_gps(1,0) - pos_ned_ins(1,0));
      y(2,0) = (float)(pos_ned_gps(2,0) - pos_ned_ins(2,0));
      y(3,0) = (float)(V_gps(0,0) - V_ins(0,0));
      y(4,0) = (float)(V_gps(1,0) - V_ins(1,0));
      y(5,0) = (float)(V_gps(2,0) - V_ins(2,0));

    #ifdef DEBUG
      printf("meas: %f,%f,%f,%f,%f,%f  ",y(0,0),y(1,0),y(2,0),y(3,0),y(4,0),y(5,0));
    #endif

      // Kalman gain
      K = P*H.transpose()*(H*P*H.transpose() + R).inverse();
      // Covariance update (a posteriori)
      P = (Eigen::Matrix<float,15,15>::Identity()-K*H)*P*(Eigen::Matrix<float,15,15>::Identity()-K*H).transpose() + K*R*K.transpose();
//wtf? the article says P = (I-KH)*Pprev, not P = (I-KH)*Pprev*(I-KH)^T + K*R+K^-1

      // State update
      Dx = K*y;

#ifdef DEBUG
      printf("S: %f,%f,%f,%f,%f,%f,%f   \n",Dx(0,0),Dx(1,0),Dx(2,0),Dx(3,0),Dx(4,0),Dx(5,0),Dx(6,0),Dx(7,0));
#endif

      //Updating all other variables related to states:  

      denom = (1.0 - (ECC2 * pow(sin(lla_ins(0,0)),2.0)));
      denom = sqrt(denom*denom);
      Re = EARTH_RADIUS / sqrt(denom);
      Rn = EARTH_RADIUS*(1.0-ECC2) / denom*sqrt(denom);
      
//DG: I don't understand why 

      alt_ins = alt_ins - Dx(2,0);
      lat_ins = lat_ins + Dx(0,0) / (Re + alt_ins);
      lon_ins = lon_ins + Dx(1,0) / (Rn + alt_ins) / cos(lat_ins);

      vn_ins = vn_ins + Dx(3,0);
      ve_ins = ve_ins + Dx(4,0);
      vd_ins = vd_ins + Dx(5,0);
      
      // Attitude correction
      //DX(6:8) are phi theta psi, but we want to uptade quat as well...
      dq(0,0) = 1.0f;
      dq(1,0) = Dx(6,0);
      dq(2,0) = Dx(7,0);
      dq(3,0) = Dx(8,0);
      quat = qmult(quat,dq);
      quat.normalize();

// #ifdef DEBUG
      printf("TPS : %f,%f,%f    ",Dx(6,0),Dx(7,0),Dx(8,0));
// #endif

      // obtain euler angles from quaternion
      theta = asinf(-2.0f*(quat(1,0)*quat(3,0)-quat(0,0)*quat(2,0)));
      phi = atan2f(2.0f*(quat(0,0)*quat(1,0)+quat(2,0)*quat(3,0)),1.0f-2.0f*(quat(1,0)*quat(1,0)+quat(2,0)*quat(2,0)));
      psi = atan2f(2.0f*(quat(1,0)*quat(2,0)+quat(0,0)*quat(3,0)),1.0f-2.0f*(quat(2,0)*quat(2,0)+quat(3,0)*quat(3,0)));

//#ifdef DEBUG
      printf("TPS : %f,%f,%f    ",theta, phi, psi);
//#endif
      // These are the biases on acceleration and gyro measurements, that we track as states.
      abx = abx + Dx(9,0);
      aby = aby + Dx(10,0);
      abz = abz + Dx(11,0);
      gbx = gbx + Dx(12,0);
      gby = gby + Dx(13,0);
      gbz = gbz + Dx(14,0);

#ifdef DEBUG
      printf("update: %f,%f,%f\n",theta, phi, psi);
#endif
    }
    // Get the new Specific forces and Rotation Rate,
    // use in the next time update
    f_b(0,0) = ax; //- abx;
    f_b(1,0) = ay; //- aby;
    f_b(2,0) = az; //- abz;

#ifdef DEBUG
      printf("f_bEND: %f,%f,%f\n",f_b(0,0),f_b(1,0),f_b(2,0));
#endif

    //DG TODO: why do you use p,q,r here??? you have it in your states! 6->8
    om_ib(0,0) = p; //- gbx;
    om_ib(1,0) = q; //- gby;
    om_ib(2,0) = r; //- gbz;
  // }

  printf("\n\n");
}

// returns the pitch angle, rad
float uNavINS::getPitch_rad() {
  return theta;
}

// returns the roll angle, rad
float uNavINS::getRoll_rad() {
  return phi;
}

// returns the yaw angle, rad
float uNavINS::getYaw_rad() {
  return constrainAngle180(psi-psi_initial);
}

// returns the heading angle, rad
float uNavINS::getHeading_rad() {
  return constrainAngle360(psi);
}

// returns the INS latitude, rad
double uNavINS::getLatitude_rad() {
  return lat_ins;
}

// returns the INS longitude, rad
double uNavINS::getLongitude_rad() {
  return lon_ins;
}

// returns the INS altitude, m
double uNavINS::getAltitude_m() {
  return alt_ins;
}

// returns the INS north velocity, m/s
double uNavINS::getVelNorth_ms() {
  return vn_ins;
}

// returns the INS east velocity, m/s
double uNavINS::getVelEast_ms() {
  return ve_ins;
}

// returns the INS down velocity, m/s
double uNavINS::getVelDown_ms() {
  return vd_ins;
}

// returns the INS ground track, rad
float uNavINS::getGroundTrack_rad() {
  return atan2f((float)ve_ins,(float)vn_ins);
}

// returns the gyro bias estimate in the x direction, rad/s
float uNavINS::getGyroBiasX_rads() {
  return gbx;
}

// returns the gyro bias estimate in the y direction, rad/s
float uNavINS::getGyroBiasY_rads() {
  return gby;
}

// returns the gyro bias estimate in the z direction, rad/s
float uNavINS::getGyroBiasZ_rads() {
  return gbz;
}

// returns the accel bias estimate in the x direction, m/s/s
float uNavINS::getAccelBiasX_mss() {
  return abx;
}

// returns the accel bias estimate in the y direction, m/s/s
float uNavINS::getAccelBiasY_mss() {
  return aby;
}

// returns the accel bias estimate in the z direction, m/s/s
float uNavINS::getAccelBiasZ_mss() {
  return abz;
}

// This function gives a skew symmetric matrix from a given vector w
Eigen::Matrix<float,3,3> uNavINS::sk(Eigen::Matrix<float,3,1> w) {
  Eigen::Matrix<float,3,3> C;
  C(0,0) = 0.0f;    C(0,1) = -w(2,0); C(0,2) = w(1,0);
  C(1,0) = w(2,0);  C(1,1) = 0.0f;    C(1,2) = -w(0,0);
  C(2,0) = -w(1,0); C(2,1) = w(0,0);  C(2,2) = 0.0f;
  return C;
}

// This function calculates the rate of change of latitude, longitude, and altitude.
Eigen::Matrix<double,3,1> uNavINS::llarate(Eigen::Matrix<double,3,1> V,Eigen::Matrix<double,3,1> lla) {
  double Rew, Rns, denom;
  Eigen::Matrix<double,3,1> lla_dot;

  denom = (1.0 - (ECC2 * pow(sin(lla(0,0)),2.0)));
  denom = sqrt(denom*denom);

  Rew = EARTH_RADIUS / sqrt(denom);
  Rns = EARTH_RADIUS*(1.0-ECC2) / denom*sqrt(denom);

  lla_dot(0,0) = V(0,0)/(Rns + lla(2,0));
  lla_dot(1,0) = V(1,0)/((Rew + lla(2,0))*cos(lla(0,0)));
  lla_dot(2,0) = -V(2,0);

  return lla_dot;
}

// This function calculates the ECEF Coordinate given the Latitude, Longitude and Altitude.
Eigen::Matrix<double,3,1> uNavINS::lla2ecef(Eigen::Matrix<double,3,1> lla) {
  double Rew, denom;
  Eigen::Matrix<double,3,1> ecef;

  denom = (1.0 - (ECC2 * pow(sin(lla(0,0)),2.0)));
  denom = sqrt(denom*denom);

  Rew = EARTH_RADIUS / sqrt(denom);

  ecef(0,0) = (Rew + lla(2,0)) * cos(lla(0,0)) * cos(lla(1,0));
  ecef(1,0) = (Rew + lla(2,0)) * cos(lla(0,0)) * sin(lla(1,0));
  ecef(2,0) = (Rew * (1.0 - ECC2) + lla(2,0)) * sin(lla(0,0));

  return ecef;
}

// This function converts a vector in ecef to ned coordinate centered at pos_ref.
Eigen::Matrix<double,3,1> uNavINS::ecef2ned(Eigen::Matrix<double,3,1> ecef,Eigen::Matrix<double,3,1> pos_ref) {
  Eigen::Matrix<double,3,1> ned;
  ned(2,0)=-cos(pos_ref(0,0))*cos(pos_ref(1,0))*ecef(0,0)-cos(pos_ref(0,0))*sin(pos_ref(1,0))*ecef(1,0)-sin(pos_ref(0,0))*ecef(2,0);
  ned(1,0)=-sin(pos_ref(1,0))*ecef(0,0) + cos(pos_ref(1,0))*ecef(1,0);
  ned(0,0)=-sin(pos_ref(0,0))*cos(pos_ref(1,0))*ecef(0,0)-sin(pos_ref(0,0))*sin(pos_ref(1,0))*ecef(1,0)+cos(pos_ref(0,0))*ecef(2,0);
  return ned;
}

// quaternion to dcm
Eigen::Matrix<float,3,3> uNavINS::quat2dcm(Eigen::Matrix<float,4,1> q) {
  Eigen::Matrix<float,3,3> C_N2B;
  C_N2B(0,0) = 2.0f*powf(q(0,0),2.0f)-1.0f + 2.0f*powf(q(1,0),2.0f);
  C_N2B(1,1) = 2.0f*powf(q(0,0),2.0f)-1.0f + 2.0f*powf(q(2,0),2.0f);
  C_N2B(2,2) = 2.0f*powf(q(0,0),2.0f)-1.0f + 2.0f*powf(q(3,0),2.0f);

  C_N2B(0,1) = 2.0f*q(1,0)*q(2,0) + 2.0f*q(0,0)*q(3,0);
  C_N2B(0,2) = 2.0f*q(1,0)*q(3,0) - 2.0f*q(0,0)*q(2,0);

  C_N2B(1,0) = 2.0f*q(1,0)*q(2,0) - 2.0f*q(0,0)*q(3,0);
  C_N2B(1,2) = 2.0f*q(2,0)*q(3,0) + 2.0f*q(0,0)*q(1,0);

  C_N2B(2,0) = 2.0f*q(1,0)*q(3,0) + 2.0f*q(0,0)*q(2,0);
  C_N2B(2,1) = 2.0f*q(2,0)*q(3,0) - 2.0f*q(0,0)*q(1,0);
  return C_N2B;
}

// quaternion multiplication
Eigen::Matrix<float,4,1> uNavINS::qmult(Eigen::Matrix<float,4,1> p, Eigen::Matrix<float,4,1> q) {
  Eigen::Matrix<float,4,1> r;
  r(0,0) = p(0,0)*q(0,0) - (p(1,0)*q(1,0) + p(2,0)*q(2,0) + p(3,0)*q(3,0));
  r(1,0) = p(0,0)*q(1,0) + q(0,0)*p(1,0) + p(2,0)*q(3,0) - p(3,0)*q(2,0);
  r(2,0) = p(0,0)*q(2,0) + q(0,0)*p(2,0) + p(3,0)*q(1,0) - p(1,0)*q(3,0);
  r(3,0) = p(0,0)*q(3,0) + q(0,0)*p(3,0) + p(1,0)*q(2,0) - p(2,0)*q(1,0);
  return r;
}

// bound yaw angle between -180 and 180
float uNavINS::constrainAngle180(float dta) {
  if(dta >  M_PI) dta -= (M_PI*2.0f);
  if(dta < -M_PI) dta += (M_PI*2.0f);
  return dta;
}

// bound heading angle between 0 and 360
float uNavINS::constrainAngle360(float dta){
  dta = fmod(dta,2.0f*M_PI);
  if (dta < 0)
    dta += 2.0f*M_PI;
  return dta;
}
