%%Calculate W(t) Xcg(t)
load FTISxprt-20200305_flight3 %Initiate data

%Passenger weights in ([kg's;inch;inch])
%p1 = [80; 131; 1];                      %pilot 1
%p2 = [102; 131; 1];                     %pilot 2
%Marta = [60; 214: 1];                        %co-ordinator
%Jari  = [75; 214; 1]     ;                    %observer 1L:		
%Martin = [83; 251; 1]      ;                  %observer 1R:		
%Wessel = [66; 251; 1]     ;                   %observer 2L:		
%Simon =	[89; 288; 1]     ;                    %observer 2R:		
%Niek =	[85; 288; 1]    ;                     %observer 3L:		
%Julian	= [90; 170; 1]  ;                    %observer 3R:		

%% Rampmass calculator
%Payload 
%Position    =   [P1 ; P2 ; 1L ; 1R ; 2L ; 2R ; 3L ; 3R ; emp; Coi];
MassKg       =   [90 ; 102; 75 ; 83 ; 66 ; 89 ; 85 ; 90 ; 00 ; 60 ];
PositionInch =   [131; 131; 214; 214; 251; 251; 288; 288; 170; 170]; 
MassLbs      =   MassKg/ 0.45359237;

MomentInchLbs=   MassLbs .* PositionInch;

PaySummass   =   sum(MassLbs);          
PaySumMoment =   sum(MomentInchLbs);
PayXcgDatum  =   PaySumMoment/PaySummass;  

% OEM
BemptyW      = 9165.0;           % Basic Empty Weight (In pounds)
BemptyWXcg   = 291.65;           % cg from datum in inches 
BemptyWMoment= 2672953.5;        % Moment in Pound-Inch

%Zerofuel 
ZeroFuelMass =   PaySummass + BemptyW;
ZeroFuelMom  =   PaySumMoment + BemptyWMoment;
ZeroFuelXcg  =   ZeroFuelMom / ZeroFuelMass;

%Fuel
Fuelmass     = 4100 ;                               %Blockfuel
FuelMoment   = 11705.5*100;                         %From E2 Appendix assigment

%Rampmass 
RampMass     =  ZeroFuelMass + Fuelmass;
RampMoment   =  ZeroFuelMom  + FuelMoment;
RampXcg      =  RampMoment / RampMass ;


%% In flight Weight and Xcg calculator
%time        = time in seconds * 10
t            = 10            * 10
FuelUsed     = flightdata.lh_engine_FU.data(t) + flightdata.rh_engine_FU.data(t)

FuelLeft     = Fuelmass - FuelUsed


steeringmat = [...
   713100   2500
   741533   2600
   769960   2800
   826906   2900
   855405   3000
   883904   3100
   912480   3200
   941062   3300
   969697   3400
   998340   3500
   1027008  3600
   1055684  3700
   1084387  3800
   1113100  3900
   1141820  4000
   1170550  4100
   1199331  4200
   ]';
    
FlightFuelMoment   = interp1(steeringmat(2,:), steeringmat(1,:), FuelLeft, 'linear', 'extrap');
FlightMoment = ZeroFuelMom + FlightFuelMoment;
FlightMass   = RampMass - FuelUsed;
FlightXcg    = FlightMoment / FlightMass;

%% Other units
C_bar = 2.0569                              %MAC lenght

XcgMAC       = FlightXcg - 261.45           %Inches from the leading MAC
XcgMACmeter  = XcgMAC*0.0254                %Meters from the leadiing edge MAC
XcgMACper    = XcgMACmeter * 100/ C_bar     %percentage of MAC
