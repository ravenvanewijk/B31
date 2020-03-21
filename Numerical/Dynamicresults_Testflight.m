%Data from the flight test dynamic measurements is collected here

% Phugoid 00:42
% u [m/s], aplha [rad], theta [rad], q [rad/s]
% time, Dadc1_tas, vane_AOA, Ahrs1_Pitch, Ahrs1_bPitchRate

t1 = 42*60;
I1 = find(flightdata.time.data==t1);

u_tas1 = flightdata.Dadc1_tas.data((I1+272):(I1+2000));
alpha1 = flightdata.vane_AOA.data((I1+272):(I1+2000));
theta1 = flightdata.Ahrs1_Pitch.data((I1+272):(I1+2000));
q1 = flightdata.Ahrs1_bPitchRate.data((I1+272):(I1+2000));

elev_dte = flightdata.delta_e.data((I1+272):(I1+2000));

%plot([u_tas1,alpha1,theta1,q1,elev_dte])
%legend('uTAS [m/s]','\alpha [rad]','\theta [rad]','q [rad/s]','elev')

% Short period 00:45
% u [m/s], aplha [rad], theta [rad], q [rad/s]
% time, Dadc1_tas, vane_AOA, Ahrs1_Pitch, Ahrs1_bPitchRate

t2 = 42*60;
I2 = find(flightdata.time.data==t2);

u_tas2 = flightdata.Dadc1_tas.data((I2+272):(I2+722));
alpha2 = flightdata.vane_AOA.data((I2+272):(I2+722));
theta2 = flightdata.Ahrs1_Pitch.data((I2+272):(I2+722));
q2 = flightdata.Ahrs1_bPitchRate.data((I2+272):(I2+722));
elev_dte2 = flightdata.delta_e.data((I2+272):(I2+722));

time2 = 0:0.1:45;
plot(time2,[alpha2,q2,elev_dte2]);

legend('\alpha [rad]','q [rad/s]','elev')


% Dutch roll 00:46
% Aperiodic roll 00:48
% Spiral 00:52