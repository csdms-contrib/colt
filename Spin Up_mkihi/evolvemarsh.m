function [marshelevation,organic_autoch,organic_alloch,mineral,Fm_min,Fm_org] = evolvemarsh(marshelevation,msl,C_e,OCb)

global tr
global numiterations
global P
global dt
global ws
global timestep
global BMax
global Dmin
global Dmax
global rhoo
global rhos
global yr

%Use this code to run the evolvemarsh program independently
% tr = 1.5;
% numiterations = 500;
% P = 45000;
% dt = P/numiterations;
% ws = 0.037*10^-3;
% timestep=365*(24/12.5);

L = length(marshelevation);

time_submerged(1:L) = 0;
floodfraction(1:L) = 0;
mineralcycle(1:L) = 0;
depth = 0;
vegtype(1:L) = 0;
bgb(1:L) = 0;
d = 0;

for i = 2:numiterations
    depth(i,1:L)=0.5*tr*sin(2*pi*(i*dt/P))+(msl-marshelevation(1:L));
    ind=find(depth(i,1:L)>0);  % finds indices of points where depth >0
    time_submerged(i,ind)=dt;
    ind2=find(depth(i,1:L)<0);
    time_submerged(i,ind2)=0;
    clear ind
    clear ind2ii
end

% Calculate belowground productivity. 

%%%%%%%% Biomass curves used by Carr model %%%%%%%%%%%%%%%%%%%%%%
%Creates a biomass curve where peak biomass occurs at a depth halfway b/t
%the maximum depth for vegetation and the minimum (here, mean high water level)

dm = msl + tr/2 - marshelevation; %[m] gives the depth of the marsh surface below HWL at any given point

bgb = zeros(1,L); %[g] Belowground Biomass
organic_autoch=zeros(1,L);

for ii = 1:L
    if ii == L
        p=9;
    end
    if dm(ii) > Dmax
        bgb(ii) = 0; %[g] water
        organic_autoch(ii)=0;
    elseif dm(ii) <= Dmin
%         [forest_mortality]=forestdep(forestage); %[g] all forest is assumed to be the same age (for Goodwin we set to 66yr)
        bgb(ii) = 0; %forest_mortality; %[g] forest
        organic_autoch(ii)= bgb(ii); %forest_mortality; %[g]
    else
        bgb(ii)=BMax*(dm(ii)-Dmax)*(dm(ii)-Dmin)/(.25.*(-Dmin-Dmax)*(Dmax-3*Dmin)); %[g] marsh
%         [mortality]=mort(bgb(ii));
        organic_autoch(ii)=bgb(ii); %[g]
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% bgb_time(:)=bgb;
marshwidth=length(find(bgb>0.001)); %[m]

distance = 0; %[m] 

%Now determine mineral deposition
for xx = 1:L
    if(bgb(xx)>0)
        distance=distance+1; %[m]
        C(xx) = C_e*exp(-0.0031.*distance); %[kg/m3] coefficient "-.0031" is a fitted parameter for realistic marsh topography
    else
        distance=1; %[m]
        C_e=0.9*C_e; %[kg/m3] Decrease concentration at the new "marsh edge" by 10% with each subsequent pond formation
        C(xx) = C_e*exp(-0.0031.*distance); %[kg/m3] coefficient "-.0031" is a fitted parameter for realistic marsh topography
    end
end

floodfraction=sum(time_submerged,1)/P;
    
for i = 2:numiterations
    ind=find(depth(i,1:L)>0);  % finds indices of points where depth >0
    mineralcycle(i,ind)=C(1,ind)*ws*dt; %[kg] no depletion for now
    ind2=find(depth(i,1:L)<0);
    mineralcycle(i,ind2)=0; %[kg] no depletion for now
    clear ind
    clear ind2
end

susp_dep=sum(mineralcycle(:,1:L),1)*timestep; %[kg/yr] suspended sediment deposition in an entire year. Summing columns
susp_dep=susp_dep.*1000; %[g] Convert from kg to g.
mineral=susp_dep.*(1-OCb); %[g] Mineral deposition of suspended sediment in a given year is equal determined by the organic content of the bay sediment 
organic_alloch=susp_dep.*OCb; %[g] Organic deposition of suspended sediment is equal determined by the organic content of the bay sediment 

Fm_min = sum(mineral)/1000; %[kg/yr] Flux of mineral sediment from the bay
Fm_org = sum(organic_alloch)/1000; %[kg/yr] Flux of organic sediment from the bay

% Calculate thickness of new sediment (mineral+organic) based off LOI and
% its effect on density
loi=(organic_autoch+organic_alloch)./(mineral+organic_autoch+organic_alloch); %[LOI]
density=1./((loi./rhoo)+((1-loi)./rhos)); %[kg/m3] Bulk density is calculated according to Morris et al. (2016)
density=density.*1000; %[g/m3]

density(isnan(density))=1; %If there is zero deposition, loi calculation will divide by zero and make density nan. Set density equal to one in this case, so that accretion is zero, instead of nan.

accretion(1:L)=(mineral+organic_autoch+organic_alloch)./density(1:L); %[m] accretion in a given year

%loi and density aren't changing through time because all cells are being
%replaced each time interval.

%Update elevation
marshelevation=marshelevation+accretion; %[m]
