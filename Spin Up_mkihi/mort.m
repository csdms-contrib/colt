function [mortality]=mort(Bpeak)

%function for calculating annual change in abovegroudn biomass and annual
%aboveground mortality from Mudd, S. M., S. M. Howell, and J. T. Morris. 2009. 
%Impactof dynamic feedbacks between sedimentation, sea-level rise, and biomass
%production on near-surface marsh stratigraphy and carbon accumulation.
%Estuarine, Coastal and Shelf Science 82:377–389.
%Pass to evolvemarsh


day_peak_season=240;
Gmin=0;
Gpeak=0.0138;
phi_g=56;
Bmin=0;
dt=365;

M_daily=zeros(365,1);
for day=1:365;
    M_daily(day+1)= (1/(4*pi))*(2*dt*(Gmin+Gpeak)*pi+2*(Bmin-Bpeak)*pi*(cos(2*pi*(day-day_peak_season)/365)-cos(2*pi*(dt-day+day_peak_season)/365))-365*(Gmin-Gpeak)*(sin(2*pi*(day-day_peak_season+phi_g)/365)+sin(2*pi*(dt-day+day_peak_season-phi_g)/365)));
end

mortality=sum(M_daily);


figure
hold on
% plot((1:366),B_daily,'r')
plot((1:366),M_daily,'k')
hold off                    




