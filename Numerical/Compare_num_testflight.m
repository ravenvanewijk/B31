%% Short period
% %stationary values:
% %V0=171.9655, alpha0= 5.4218/180*pi, theta0=2.0054/180*pi, q0 = 0.1651/180*pi;     %Initial values
% %hp0 = 1714.3, mass = 6398.919kg
% %Interesting are: alpha and q
% 
% t1 = 0:0.1:20;
% x0_1 = [0, 0.1, 0, 0];
% elev = deg2rad(flightdata.delta_e.data((I1+360):(I1+560)));
% y1 = lsim(sys1, elev, t1);
% 
% ti1 = 45*60;
% I1 = find(flightdata.time.data==ti1);
% 
% u_tas1 = flightdata.Dadc1_tas.data((I1+360):(I1+560));
% alpha1 = deg2rad(flightdata.vane_AOA.data((I1+360):(I1+560)));
% theta1 = deg2rad(flightdata.Ahrs1_Pitch.data((I1+360):(I1+560)));
% q1 = deg2rad(flightdata.Ahrs1_bPitchRate.data((I1+360):(I1+560)));
% u0_1 = flightdata.Dadc1_tas.data(I1+310);
% hp0_1 = flightdata.Dadc1_alt.data(I1+310)*0.3048;
% alpha0_1 = flightdata.vane_AOA.data(I1+310);
% theta0_1 = flightdata.Ahrs1_Pitch.data(I1+310);
% q0_1 = flightdata.Ahrs1_bPitchRate.data(I1+310);
% 
% time1 = 0:0.1:20;
% 
% % subplot(2,1,1)
% % plot(t1, y1(:,2),time1,[alpha1-alpha0,elev])
% % legend('numerical','test','\delta_e [rad]')
% % title('Angle of attack')
% % grid()
% % xlabel('time [s]')
% % ylabel('\alpha [rad]')
% %subplot(2,1,2)
% plot(t1, y1(:,4)+0.004541,time1,[q1])
% legend('numerical','test')
% title('Pitch rate')
% grid()
% xlabel('time [s]')
% ylabel(' q [rad/s]')
%% Phugoid
% % Stationary values:
% % V0 = 171.4245, alpha_0 = 5.8008/180*pi, theta_0 = 4.5117/180*pi, q0 = 0.0935
% % hp0 = 1701.7, m = 6421.26kg
% % Interesting are: u and theta
% 
% ti2 = 42*60;
% I2 = find(flightdata.time.data==ti2);
% elev2 = deg2rad(flightdata.delta_e.data((I2+272):(I2+1772)));
% t2 = 0:0.1:150;      %Elevator input in rad   %Initial values
% y2 = lsim(sys1, elev2, t2);
% 
% u0_2 = flightdata.Dadc1_tas.data(I2+272);
% alpha0_2 = flightdata.vane_AOA.data(I2+272);
% theta0_2 = flightdata.Ahrs1_Pitch.data(I2+272);
% q0_2 = deg2rad(flightdata.Ahrs1_bPitchRate.data(I2+272));
% hp0_2 = flightdata.Dadc1_alt.data(I2+272)*0.3048;
% 
% u_tas2 = flightdata.Dadc1_tas.data((I2+272):(I2+1772));
% alpha2 = deg2rad(flightdata.vane_AOA.data((I2+272):(I2+1772)));
% theta2 = deg2rad(flightdata.Ahrs1_Pitch.data((I2+272):(I2+1772)));
% q2 = deg2rad(flightdata.Ahrs1_bPitchRate.data((I2+272):(I2+1772)));
% time2 = 0:0.1:150;
% 
% subplot(4,1,1)
% plot(t2, y2(:,1)+u0_2,time2,u_tas2)
% legend('numerical','test')
% title('Velocity')
% grid minor
% xlabel('time [s]')
% ylabel('u [m/s]')
% subplot(4,1,2)
% plot(t2, y2(:,3),time2,(theta2-th0))
% legend('numerical','test')
% title('Pitch')
% grid minor
% xlabel('time [s]')
% ylabel('\theta [rad]')
% subplot(4,1,3)
% plot(t2, y2(:,4)+q0_2,time2,q2)
% legend('numerical','test')
% title('Pitch rate')
% grid minor
% xlabel('time [s]')
% ylabel('q [rad/s]')
% subplot(4,1,4)
% plot(time2,elev2)
% title('Elevator deflection')
% grid minor
% xlabel('time [s]')
% ylabel('\delta_e [rad]')
%% Aperiodic roll
% % Stationary values:
% % V0 = 175.3652, hp0 = 1634, alpha0 = 7.1213/180*pi, th0 = 4.3314/180*pi
% % m = 6378.41kg
% 
% ti3 = 48*60;
% I3 = find(flightdata.time.data==ti3);
% 
% 
% ail_def1 = zeros(2,601);
% ail_def1(1,:) = flightdata.delta_a.data(I3+340:I3+940)/180*pi;
% 
% t3 = 0:0.1:60;
% y3 = lsim(sys2,ail_def1,t3);
% 
% %phi, p
% %roll angle, roll rate
% %Ahrs1_Roll, Ahrs1_bRollRate
% 
% V0_3 = flightdata.Dadc1_tas.data(I3+340);
% hp0_3 = flightdata.Dadc1_alt.data(I3+340)*0.3048;
% alpha0_3 = flightdata.vane_AOA.data(I3+340);
% theta0_3= flightdata.Ahrs1_Pitch.data(I3+340);
% 
% phi1 = flightdata.Ahrs1_Roll.data(I3+340:I3+940)/180*pi;
% rollrate1 = flightdata.Ahrs1_bRollRate.data(I3+340:I3+940)/180*pi;
% yawrate1 = flightdata.Ahrs1_bYawRate.data(I3+340:I3+940)/180*pi;
% 
% time3 = 0:0.1:60;
% 
% %  t3,y3(:,2:3)
% 
% subplot(2,1,1)
% plot(t3,-y3(:,3)-0.0133,time3,[rollrate1])
% legend('numerical','test')
% grid minor
% title('Roll rate')
% xlabel('time [s]')
% ylabel('p [rad/s]')
% subplot(2,1,2)
% plot(time3,ail_def1(1,:))
% grid minor
% title('Aileron deflection')
% xlabel('time [s]')
% ylabel('\delta_a [rad]')

%% Spiral
% % Stationary values
% % hp0 = 2508.3, alpha0 = 5.8124/180*pi, theta0 = 4.7420/180*pi, 
% % V0 = 177.0704, m =  6341.62kg
% 
% ail_def2 = zeros(2,1501);
% ail_def2(1,:) = flightdata.delta_a.data(I4+80:I4+1580)/180*pi;
% t4 = 0:0.1:150;
% y4 = lsim(sys2,ail_def2,t4);
% 
% ti4 = 52*60;
% I4 = find(flightdata.time.data==ti4);
% 
% V0_4 = flightdata.Dadc1_tas.data(I4+80);
% hp0_4 = flightdata.Dadc1_alt.data(I4+80)*0.3048;
% alpha0_4 = flightdata.vane_AOA.data(I4+80);
% theta0_4 = flightdata.Ahrs1_Pitch.data(I4+80);
% 
% phi2 = flightdata.Ahrs1_Roll.data(I4+80:I4+1580)/180*pi;
% rollrate2= flightdata.Ahrs1_bRollRate.data(I4+80:I4+1580)/180*pi;
% yawrate2 = flightdata.Ahrs1_bYawRate.data(I4+80:I4+1580)/180*pi;
% 
% time4= 0:0.1:150;
% 
% %t4,y4(:,2:3),
% 
% % subplot(4,1,1)
% % plot( t4,y4(:,2),time4, [phi2])
% % title('Roll angle')
% % legend('numerical','test')
% % xlabel('time [s]')
% % ylabel('\phi [rad]')
% % grid minor
% subplot(3,1,1)
% plot( t4,-y4(:,3)/20-0.003695,time4, [rollrate2])
% title('Roll rate')
% legend('numerical','test')
% xlabel('time [s]')
% ylabel('p [rad/s]')
% grid minor
% subplot(3,1,2)
% plot( t4,y4(:,4)/20-0.01012,time4, [yawrate2])
% title('Yaw rate')
% legend('numerical','test')
% xlabel('time [s]')
% ylabel('r [rad/s]')
% grid minor
% subplot(3,1,3)
% plot(time4, (ail_def2(1,:)))
% title('Aileron deflection ')
% xlabel('time [s]')
% ylabel('\delta_a [rad]')
% grid
%% Dutch roll
% Stationary values
% V0 = 168.8398, hp0 = 1736.7, alpha0 = 5.4285/180*pi, th0 = 3.8658/180*pi
% m = 6392.77kg
% 
% ti5 = 46*60;
% I5 = find(flightdata.time.data==ti5);
% 
% t5 = 0:0.1:25;
% rud_def1 = zeros(2,251);
% rud_def1(2,:) = flightdata.delta_r.data((I5+324):(I5+574))/180*pi;
% y5 = lsim(sys2,rud_def1,t5);
% 
% V0_5 = flightdata.Dadc1_tas.data(I5+324);
% hp0_5 = flightdata.Dadc1_alt.data(I5+324)*0.3048;
% alpha0_5 = flightdata.vane_AOA.data(I5+324);
% theta0_5 = flightdata.Ahrs1_Pitch.data(I5+324);
% 
% phi3 = flightdata.Ahrs1_Roll.data((I5+324):(I5+574))/180*pi;
% rollrate3= flightdata.Ahrs1_bRollRate.data((I5+324):(I5+574))/180*pi;
% yawrate3 = flightdata.Ahrs1_bYawRate.data((I5+324):(I5+574))/180*pi;
% time5= 0:0.1:25;
% 
% subplot(3,1,1)
% plot(t5, -y5(:,3)-0.1378,time5, [rollrate3])
% title('Roll rate')
% legend('numerical','test')
% xlabel('time [s]')
% ylabel('p [rad/s]')
% grid minor
% subplot(3,1,2)
% plot(time5, [yawrate3])
% title('Yaw rate')
% legend('numerical','test')
% xlabel('time [s]')
% ylabel('r [rad/s]')
% grid minor
% subplot(3,1,3)
% plot(time5, rud_def1(2,:))
% xlabel('time [s]')
% title('Rudder deflection ')
% ylabel('\delta_r [rad]')
% grid minor











