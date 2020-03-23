%%Calculate W(t) Xcg(t)
clear
load ReferenceData %Initiate data	

%% Rampmass calculator
%Payload 
%Position    =   [P1 ; P2 ; 1L ; 1R ; 2L ; 2R ; 3L ; 3R ; emp; Coi];
MassKg       =   [95 ; 92 ; 66 ; 61 ; 75 ; 78 ; 86 ; 68 ; 00 ; 74 ];
PositionInch =   [131; 131; 214; 214; 251; 251; 131; 288; 170; 170]; 
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
Fuelmass     = 4050 ;                            %Blockfuel
FuelMoment   = 11705.5*100;                         %From E2 Appendix assigment

%Rampmass 
RampMass     =  ZeroFuelMass + Fuelmass;
RampMoment   =  ZeroFuelMom  + FuelMoment;
RampXcg      =  RampMoment / RampMass ;


%% In flight Weight and Xcg calculator
%time        = time in seconds * 10
t            = (52*60+46)           * 10;
FuelUsed     = flightdata.lh_engine_FU.data(t) + flightdata.rh_engine_FU.data(t)

FuelLeft     = Fuelmass - FuelUsed


steeringmat = [...
   29816    100
   59118    200
   87908    300
   116542   400
   144840   500
   173252   600
   201480   700
   229884   800
   258192   900
   286630   1000
   315018   1100
   343452   1200
   371852   1300
   400323   1400
   428776   1500
   457224   1600
   485656   1700
   514116   1800
   542564   1900
   570990   2000
   599404   2100
   627847   2200
   656282   2300
   684696   2400
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
FlightMassKG = FlightMass * 0.45359237;
FlightXcg    = FlightMoment / FlightMass;

%% Other units
C_bar = 2.0569                              %MAC lenght

XcgMAC       = FlightXcg - 261.45 ;          %Inches from the leading MAC
XcgMACmeter  = XcgMAC*0.0254       ;         %Meters from the leadiing edge MAC
XcgMACper    = XcgMACmeter * 100/ C_bar;     %percentage of MAC
