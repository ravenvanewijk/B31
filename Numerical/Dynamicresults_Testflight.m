
%u_tas = flightdata.Dadc1_tas.data();
%alpha = deg2rad(flightdata.vane_AOA.data());
theta = deg2rad(flightdata.Ahrs1_Pitch.data());
q = deg2rad(flightdata.Ahrs1_bPitchRate.data());
elev = deg2rad(flightdata.delta_e.data());
phi = flightdata.Ahrs1_Roll.data()/180*pi;
Rollrate= flightdata.Ahrs1_bRollRate.data()/180*pi;
yawrate = flightdata.Ahrs1_bYawRate.data()/180*pi;
rud_def = flightdata.delta_r.data()/180*pi;
ail_def = flightdata.delta_a.data()/180*pi;
timevector = 0:0.1:4212;

plot(timevector,[  theta, q, elev, phi, Rollrate, yawrate,rud_def,ail_def])
legend('\theta','q','elev','phi','roll rate','yaw rate','rudder','ail_def')