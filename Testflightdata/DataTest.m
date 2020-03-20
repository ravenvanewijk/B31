%%Testflight data viewer%%
load FTISxprt-20200305_flight3                           % Data March 5 flight 3  

%Plots
subplot(2,2,1)
plot(flightdata.time.data ,flightdata.Dadc1_tas.data)   %Velocity Plot

subplot(2,2,2)
plot(flightdata.time.data ,flightdata.Ahrs1_bNormAcc.data)   %Acceleration Plot

subplot(2,2,3)
plot(flightdata.time.data ,flightdata.Ahrs1_bYawRate.data)   %Yaw Rate


%{ Aperiodic rol
%subplot(2,3,3)
%t = 0:0.01:15;
%x0 = [0, 0.1, 0, 0];
%y3 = initial(sys2, x0, t);
%plot(t, y3(:,:))
%title('Aperiodic roll')
%legend('\beta [rad]','\phi [rad]','p [rad/s]','r [rad/s]')
%}
