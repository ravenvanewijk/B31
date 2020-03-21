function [W,Xcg] = MassBalance(t)
% For calculating W(t) and Xcg(t)

IMPtoMET3 = [0.45359237; 0.0254; 0.0254];
IMPtoMET1 = [0.45359237];

%%Initial weights and positions (Weight,Xposition, Yposition) %%
BemptyWimp = 9165.0;           % Basic Empty Weight (In pounds)
BlockFuelimp = 2600;           % Block fuel (In pounds)

BasicEmptyWeight    =   BemptyWimp * IMPtoMET1;
BlockFuel           =   BlockFuelimp * IMPtoMET1;  
                      
%Passenger weights in (kg's)

p1 = [80; 1; 1];                      % pilot 1
p2 = [102;                            % pilot 2
Marta = 60;                   % co-ordinator
Jari  = 75;                  %observer 1L:		
Martin = 83;                 %observer 1R:		
Wessel = 66 ;                %observer 2L:		
Simon =	89    ;               %observer 2R:		
Niek =	85     ;              %observer 3L:		
Julian	= 90    ;             %observer 3R:		

W=30;
Xcg =5

end

