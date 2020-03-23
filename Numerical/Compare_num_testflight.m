%% Short period
% t1 = 0:0.01:15;
% %x0_1 = [0, 0.1, 0, 0];
% u1 = zeros(1,length(t1));
% u1(t1>=1 & t1<=12) = -0.6761;        %Elevator input in rad
% x0_1 = [174.4,5.595,3.931,0.2602];     %Initial values
% y1 = lsim(sys1, u1, t1, x0_1);
% 
% ti1 = 45*60;
% I1 = find(flightdata.time.data==ti1);
% 
% u_tas1 = flightdata.Dadc1_tas.data((I1+397):(I1+547));
% alpha1 = flightdata.vane_AOA.data((I1+397):(I1+547));
% theta1 = flightdata.Ahrs1_Pitch.data((I1+397):(I1+547));
% q1 = flightdata.Ahrs1_bPitchRate.data((I1+397):(I1+547));
% time1 = 0:0.1:15;
% 
% subplot(2,2,1)
% plot(t1, y1(:,1),time1,u_tas1)
% legend('numerical','test')
% title('u')
% subplot(2,2,2)
% plot(t1, y1(:,2),time1,alpha1)
% legend('numerical','test')
% title('AoA')
% subplot(2,2,3)
% plot(t1, y1(:,3),time1,theta1)
% legend('numerical','test')
% title('Pitch')
% subplot(2,2,4)
% plot(t1, y1(:,4),time1,q1)
% legend('numerical','test')
% title('Pitch rate')

%% Phugoid
t2 = 0:0.01:171.5;
u2 = zeros(1,length(t2));
u2(t2>=1 & t2<=12) = -0.6761;        %Elevator input in rad
x0_2 = [174.4,5.595,3.931,0.2602];     %Initial values
y2 = lsim(sys1, u2, t2, x0_2);

ti2 = 42*60;
I2 = find(flightdata.time.data==ti2);

u_tas2 = flightdata.Dadc1_tas.data((I2+272):(I2+1987));
alpha2 = flightdata.vane_AOA.data((I2+272):(I2+1987));
theta2 = flightdata.Ahrs1_Pitch.data((I2+272):(I2+1987));
q2 = flightdata.Ahrs1_bPitchRate.data((I2+272):(I2+1987));
time2 = 0:0.1:171.5;

subplot(2,2,1)
plot(t2, y2(:,1),time2,u_tas2)
legend('numerical','test')
title('u')
subplot(2,2,2)
plot(t2, y2(:,2),time2,alpha2)
legend('numerical','test')
title('AoA')
subplot(2,2,3)
plot(t2, y2(:,3),time2,theta2)
legend('numerical','test')
title('Pitch')
subplot(2,2,4)
plot(t2, y2(:,4),time2,q2)
legend('numerical','test')
title('Pitch rate')

%% Aperiodic roll
% t3 = 0.:0.01:45;
% %u3 = zeros(1,length(t3));
% %u3(t3>=1 & t3<=12) = -0.6761;
% %x0_3 = [];
% %y3 = lsim(sys2,u3,t3);
% x0 = [0, 0.1, 0, 0];
% y3 = initial(sys2, x0, t3);
% 
% 
% ti3 = 48*60;
% I3 = find(flightdata.time.data==ti3);
% 
% %phi, p
% %roll angle, roll rate
% %Ahrs1_Roll, Ahrs1_bRollRate
% 
% phi1 = flightdata.Ahrs1_Roll.data(I3-322:I3+128)/180*pi;
% rollrate1 = flightdata.Ahrs1_bRollRate.data(I3-322:I3+128)/180*pi;
% ail_def1 = flightdata.delta_a.data(I3-322:I3+128)/180*pi;
% time3 = 0:0.1:45;
% 
% plot( time3, [phi1,rollrate1,ail_def1])
% legend('phi','roll rate','aileron deflection')
% 



