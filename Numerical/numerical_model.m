%% Variables
% Citation 550 - Linear simulation

% Stationary flight condition

hp0    = 1736.7;     % pressure altitude in the stationary flight condition [m]
V0     = 168.8398;      % true airspeed in the stationary flight condition [m/sec]
alpha0 = 5.4285/180*pi;       	  % angle of attack in the stationary flight condition [rad]
th0    = 3.8658/180*pi;       	  % pitch angle in the stationary flight condition [rad]

% Aircraft mass
m      = 6392.77;         	  % mass [kg]

% aerodynamic properties
e      = 0.873;            % Oswald factor [ ]
CD0    = 0.02205;          % Zero lift drag coefficient [ ]
CLa    = 4.68107;            % Slope of CL-alpha curve [ ]

% Longitudinal stability
Cma    =  -0.60358;      %-0.5626;     % longitudinal stabilty [ ]
Cmde   =  -1.38183;           %-1.1642;     % elevator effectiveness [ ]

% Aircraft geometry

S      = 30.00;	          % wing area [m^2]
Sh     = 0.2*S;           % stabiliser area [m^2]
Sh_S   = Sh/S;	          % [ ]
lh     = 0.71*5.968;      % tail length [m]
c      = 2.0569;	  % mean aerodynamic cord [m]
lh_c   = lh/c;	          % [ ]
b      = 15.911;	  % wing span [m]
bh     = 5.791;	          % stabilser span [m]
A      = b^2/S;           % wing aspect ratio [ ]
Ah     = bh^2/Sh;         % stabilser aspect ratio [ ]
Vh_V   = 1;		  % [ ]
ih     = -2*pi/180;       % stabiliser angle of incidence [rad]

% Constant values concerning atmosphere and gravity

rho0   = 1.2250;          % air density at sea level [kg/m^3] 
lambda = -0.0065;         % temperature gradient in ISA [K/m]
Temp0  = 288.15;          % temperature at sea level in ISA [K]
R      = 287.05;          % specific gas constant [m^2/sec^2K]
g      = 9.81;            % [m/sec^2] (gravity constant)

%rho    = rho0*((1+(lambda*hp0/Temp0)))^(-((g/(lambda*R))+1));   % [kg/m^3]  (air density)
rho = 1.2;
W      = m*g;				                        % [N]       (aircraft weight)

% Constant values concerning aircraft inertia

muc    = m/(rho*S*c);
mub    = m/(rho*S*b);
KX2    = 0.019;
KZ2    = 0.042;
KXZ    = 0.002;
KY2    = 1.25*1.114;

% Aerodynamic constants

Cmac   = 0;                     % Moment coefficient about the aerodynamic centre [ ]
CNwa   = CLa;   		        % Wing normal force slope [ ]
CNha   = 2*pi*Ah/(Ah+2);        % Stabiliser normal force slope [ ]
depsda = 4/(A+2);               % Downwash gradient [ ]

% Lift and drag coefficient

CL = 2*W/(rho*V0^2*S);             % Lift coefficient [ ]
CD = CD0 + (CLa*alpha0)^2/(pi*A*e);  % Drag coefficient [ ]

% Stabiblity derivatives

CX0    = W*sin(th0)/(0.5*rho*V0^2*S);
CXu    = -0.095;
CXa    = +0.47966;  %was negative 
CXadot = +0.08330;
CXq    = -0.28170;
CXde   = -0.03728;

CZ0    = -W*cos(th0)/(0.5*rho*V0^2*S);
CZu    = -0.37616;
CZa    = -5.74340;
CZadot = -0.00350;
CZq    = -5.66290;
CZde   = -0.69612;

Cmu    = +0.06990;
Cmadot = +0.17800;
Cmq    = -8.79415;

CYb    = -0.7500;
CYbdot =  0     ;
CYp    = -0.0304;
CYr    = +0.8495;
CYda   = -0.0400;
CYdr   = +0.2300;

Clb    = -0.10260;
Clp    = -0.71085;
Clr    = +0.23760;
Clda   = -0.23088;
Cldr   = +0.03440;

Cnb    =  +0.1348;
Cnbdot =   0     ;
Cnp    =  -0.0602;
Cnr    =  -0.2061;
Cnda   =  -0.0120;
Cndr   =  -0.0939;

xcg = 0.25*c;

%% Variables obtained from the flight test

% xcg
% Speed standardization to V_eq
% Lift-aoa curve: CLa
% Thrust calculation
% Drag curve: Cd0
% Elevator effectiveness Cmde
% Elevation deflection and force control curves
% Pitch stability Cma



%% Numerical Model 

% flight test data sheet  

% Angle of attack [deg]
% Deflection of elevator trim [deg]
% Force on elevator control wheel [N]
% Engine 1 (left): Fuel mass flow [lbs/hr]
% Engine 2 (right): Fuel mass flow [lbs/hr]
% Engine 1: Inter Turbine Temperature (ITT) [deg C]
% Engine 2: Inter Turbine Temperature (ITT) [deg C]
% Engine 1: Oil pressure [psi]
% Engine 2: Oil pressure [psi]
% Deflection of the control column (Se or DCOC) [deg]
% Engine 1: Fan speed (N1) [%]
% Engine 1: Turbine speed (N2) [%]
% Engine 2: Fan speed (N1) [%]
% Engine 2: Turbine speed (N2) [%]
% Engine 1: Calculated fuel used by fuel mass flow [lbs]
% Engine 2: Calculated fuel used by fuel mass flow [lbs]
% Deflection of aileron (right wing?) [deg]
% Deflection of elevator [deg]
% Deflection of rudder [deg]
% UTC Date DD:MM:YY
% UTC Seconds
% Roll angle [deg]
% Pitch angle [deg]
% 
% GNSS Latitdue [deg]
% GNSS Longitude [deg]
% Body Roll rate [deg/s]
% Body Pitch rate [deg/s]
% Body Yaw rate [deg/s]
% Body Longitudinal acceleration [g]
% Body Lateral acceleration [g]
% Body Normal acceleration [g]
% Along heading acceleration [g]
% Cross heading acceleration [g]
% Vertical acceleration [g]
% Static air temperature [deg C]
% Total air temperature [deg C]
% Pressure Altitude (1013.25 mB) [ft]
% Baro corrected altitude #1 [ft]
% 
% Mach number [-]
% Computed Airspeed [knots]
% True Airspeed [knots]
% Altitude rate [ft/min]
% 
% 
% 
% 
% Time [sec]

%% Symmetric case

C1s = [- 2* muc* (c/V0^2), 0, 0, 0;
     0, (CZadot - 2* muc)* (c/V0), 0, 0;
     0, 0, -(c/V0), 0;
     0, Cmadot* (c/V0), 0, -2* muc* KY2*(c/V0)^2];

C2s = [CXu/V0, CXa, CZ0, CXq*(c/V0);
     CZu/V0, CZa, -CX0, (CZq + 2* muc)*(c/V0);
     0, 0, 0,(c/V0);
     Cmu/V0, Cma, 0, Cmq*(c/V0)];

C3s = [CXde;
    CZde;
    0;
    Cmde];

As = -inv(C1s)*C2s;
Bs = -inv(C1s)*C3s;
Cs = eye(4);
Ds = zeros(4,1);
sys1 = ss(As,Bs,Cs,Ds);


pole(sys1)
%rltool(sys1(1,:))

%% Asymmetric case

Db = (b/V0);

C1a = [(CYbdot- 2* mub)* Db, 0, 0, 0;
     0, -0.5 * Db, 0, 0;
     0, 0, -2*mub*KX2*Db^2, 2*mub*KXZ*Db^2;
     Cnbdot*Db, 0, 2*mub*KXZ*Db^2, -2*mub*KZ2*Db^2];

C2a = [CYb, CL, CYp*0.5*Db, (CYr - 4*mub)*0.5*Db;
     0, 0, 1*0.5*Db, 0;
     Clb, 0, Clp*0.5*Db, Clr*0.5*Db;
     Cnb, 0, Cnp*0.5*Db, Cnr*0.5*Db];

C3a = [CYda, CYdr;
    0, 0;
    Clda, Cldr;
    Cnda, Cndr];

Aa = -inv(C1a)*C2a;
Ba = -inv(C1a)*C3a;
Ca = eye(4);
Da = zeros(4,2);
sys2 = ss(Aa,Ba,Ca,Da);

pole(sys2)
%rltool(sys2(1,1))

%% Plotting initial responses to disturbance input (requirement no.11)

% Symmetric

%% Short period
% % AoA and pitch rate the short period damping the best
% subplot(1,2,1)
% t = 0:0.01:10;
% y1 = -1/180*pi*step(sys1, t); %distubance input on d_e pf -1 deg
% plot(t, y1(:,2))
% xlabel('time [s]')
% ylabel('\alpha [rad]')
% grid()
% title('Angle of attack')
% 
% subplot(1,2,2)
% plot(t,y1(:,4))
% xlabel('time [s]')
% ylabel('q [rad/s]')
% grid()
% title('Pitch rate')


%% Phugoid
% %U and theta show this eigenmotion response the best
% subplot(2,2,1)
% t = 0:0.01:400;
% y2 = -1/180*pi*step(sys1, t); %distubance input on d_e pf -1 deg
% plot(t, y2(:,1))            % speed = 150 m/s
% xlabel('time [s]')
% ylabel('u [m/s]')
% grid()
% title('Velocity')
% 
% subplot(2,2,2)
% t = 0:0.01:400;
% y2 = -1/180*pi*step(sys1, t); %distubance input on d_e pf -1 deg
% plot(t, y2(:,2))            % speed = 150 m/s
% xlabel('time [s]')
% ylabel('\alpha [rad]')
% grid()
% title('Angle of attack')
% 
% subplot(2,2,3)
% t = 0:0.01:400;
% y2 = -1/180*pi*step(sys1, t); %distubance input on d_e pf -1 deg
% plot(t, y2(:,3))            % speed = 150 m/s
% xlabel('time [s]')
% ylabel('\theta [rad]')
% grid()
% title('Pitch')
% 
% subplot(2,2,4)
% t = 0:0.01:400;
% y2 = -1/180*pi*step(sys1, t); %distubance input on d_e pf -1 deg
% plot(t, y2(:,4))
% xlabel('time [s]')
% ylabel('q [rad/s]')
% grid()
% title('Pitch rate')

%% Asymmetric

%% Aperiodic roll
% subplot(1,2,1)
% t = 0:0.01:10;
% x0 = [0, 1/180*pi, 0, 0];   %Initial roll angle
% y3 = initial(sys2, x0, t);
% plot(t, y3(:,:))
% grid()
% title('Aperiodic roll')
% xlabel('time [s]')
% x = 'x';
% s2 = ['$\bar{' x '}$'];
% ylabel(s2, 'Interpreter', 'LaTeX');
% legend('\beta [rad]','\phi [rad]','p [rad/s]','r [rad/s]')
% 
% %% Spiral
% 
% subplot(1,2,2)
% t = 0:0.01:100;
% x0 = [0, 1/180*pi, 0, 0];     %Initial roll angle
% y5 = initial(sys2, x0, t);
% plot(t, y5(:,:))
% grid()
% xlabel('time [s]')
% x = 'x';
% s3 = ['$\bar{' x '}$'];
% ylabel(s3, 'Interpreter', 'LaTeX');
% title('Spiral')
% legend('\beta [rad]','\phi [rad]','p [rad/s]','r [rad/s]')

%% Dutch roll
% 
% t = 0:0.01:15;
% x0 = [1/180*pi, 0, 0, 0];       %Initial sideslip angle
% y4 = initial(sys2, x0, t);
% plot(t, y4(:,:))
% title('Dutch roll')
% grid()
% xlabel('time [s]')
% x = 'x';
% s4 = ['$\bar{' x '}$'];
% ylabel(s4, 'Interpreter', 'LaTeX');
% legend('\beta [rad]','\phi [rad]','p [rad/s]','r [rad/s]')
% 

