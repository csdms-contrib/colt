function [dX r]=funBAY(X,PAR)

%Determines change in bay depth and width by solving mass balance between
%fluxes of sediment into and out of the bay from marsh edge erosion, tidal
%exchange with the outside sediment source, and sediment deposited onto the
%marsh surface

global C_e_ODE
global Fc_ODE

%%%%%%%%%%
%rename parameters that were loaded into the function
rhos=PAR(1);
P=PAR(2);
B=PAR(3);
ws=PAR(4);
tcr=PAR(5);
Co = PAR(6);
wind=PAR(7);
Ba=PAR(8);
Be=PAR(9);
amp=PAR(10);
RSLR=PAR(11);
Fm2 = PAR(12);
rhos=PAR(13);
lamda=PAR(14);
dist=PAR(15);
dmo = PAR(16);
rhob = PAR(17);
rhom = PAR(18);

%%%%%%%%%%%% dynamic variable
fetch=X(1); %mudflat width
df=X(2); %mudflat depth

%%%%%%%%%%%%%%%%%%%%%%%%%
fac=min(1,df/(2*amp));%Proportion of tide that the bay is flooded, dimensionless
Df=(df+(df-fac*2*amp))/2;%[m] average bay depth over tidal cycle
dm = dmo;%[m] marsh edge depth

tw=  wavetau(fetch,wind,Df); %Calculates wave bed shear stress [Pa] as a function of fetch, wind speed, and the depth of the bay, function from Giulio Mariotti

tau=max((tw-tcr)/tcr,0)*lamda;%Excess shear stress, dimensionless
Cr=rhos*tau/(1+tau); %Reference suspended sediment concentration in the basin [kg/m3]

hb=dm+(df-dm)*(1-exp(-dist*0.1/df)); %[m] scarp height at a fixed distance from the marsh according to semi-empirical shoaling profile
W = waveTRNS(amp,wind,fetch,hb); %[W] Wave power density at the marsh boundary, function from Giulio Mariotti

E=(Be*W/(hb-dm)-Ba*Cr*ws/rhom); %(m2/s) Net flux of sediment eroded from/deposited to the marsh edge, balance between erosion and deposition at the marsh edge

Fc=(Cr-Co)*(fac*2*amp)/P/rhob; %(m2/s) Net flux of sediment lost or gained through tidal exchange with external sediment supply/sink

Fc_ODE(numel(Fc_ODE)+1)=Fc*rhob*fetch; %Save Fc as a mass flux (kg/s) for each iteration of the ODE
C_e_ODE(numel(C_e_ODE)+1) = Cr; %Save C_e (susp. sed. con. at marsh edge, kg/m3) for each iteration of the ODE to use in marsh model
%%%%%%%%%%%%%%%%%%%%%
dX = zeros(2,1);
dX(1)=E; %(m2/s, or m/s if integrated over 1m transect width) Change in bay width due to erosion
dX(2)=-E*(df-dm)/fetch +Fm2/fetch/rhos + Fc +RSLR; %(m/s) Change in bay depth due to mass balance between fluxes into and out of bay