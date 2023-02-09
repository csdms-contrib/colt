function spinup

%Function to simulate the formation of a 1000m wide marsh onto a gentle upland slope (m) under a given set
%of conditions for sea level rise (R) and suspended sediment conentration
%(C), to be used as the initial condition for model runs. 

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
global organic_dep_autoch
global B
global mui
global mki
global elevation

clc

m=.001;
R=.001; % (m/yr)
C=50/1000; % (kg/m3)
amp=1.5; % (m)
tr=amp*2; % (m)
P=12.5*3600; % tidal period [s]
numiterations=500;
dt=P/numiterations;
rhos=2000;%bulk densities [kg/m3]
rhoo=85.0;%bulk density of organic matter [kg/m3]
ws=0.05*10^-3;%wsi typically set to .05 %effective settling velocity, m sec-1, Marani 2010 uses 0.2 mm s-1, Mudd 2009 uses .037 mm/s
timestep=365*(24/12.5); %[tidal cycles per year] number to multiply accretion simulated over a tidal cycle by
BMax=2500;% g/m2
Dmin=0; % (m)
Dmax=.7167*2*amp-.483; % (m)
OCb= 0.01; %0.01 is the average value (over time) from simulations in the main model with the same conditions.
B=1000; %Transect width (m)
dur=550; % (yr)

%Decomposition parameters
mui = 0.4; %[m] Depth below which decomposition goes to zero in the marsh
mki = 0.2; %Coefficient of decomposition in the marsh

tempelevation = (1:B).*m+amp+.02-Dmax;
% dur = (tempelevation(B)-amp)/R+50;
% dur=ceil(dur);

elevation = ones(dur,B);
elevation(1,:) = tempelevation;
organic_dep_autoch=zeros(dur,B);
organic_dep_alloch=zeros(dur,B);


close all

figure
xlabel('Distance (m)')
ylabel('Elevation relative to initial SL (m)')
plot(1:B,elevation(1,:),'k--')
hold on

yri=1;

for yr = 2:dur
    yr
    msl=yr*R;
    [tempelevation,temporg_autoch,temporg_alloch,tempmin]=evolvemarsh(elevation(yr-1,:),msl,C,OCb);
    mineral_dep(yr,:) = tempmin; %[g] mineral sediment deposited in a given year
    organic_dep_autoch(yr,:) = temporg_autoch; %[g] belowground plant material deposited in a given year
    organic_dep_alloch(yr,:) = temporg_alloch; %[g] allochthanous organic material deposited in a given year
    elevation(yr,:)=tempelevation;
    x_f = find(elevation(yr,:)>(msl+amp-Dmin),1);
    [compaction]=decompose(1,x_f);
    elevation(yr,:) = elevation(yr,:)-compaction; %Adjust marsh and forest elevation due to compaction from decomposition
    if yr > 25*yri
        plot(1:B,elevation(yr,:),'k-')
        yri = yri+1;
    end
end

elev=elevation;
save(['elev.mat'],'elev')
elev_25=elev;
save(['mineral deposition.mat'],'mineral_dep')
min_25=mineral_dep;
save(['organic deposition.mat'],'organic_dep_autoch','organic_dep_alloch')
orgAT_25=organic_dep_autoch;
orgAL_25=organic_dep_alloch;
save(['MarshStrat_all_3mtr_RSLR' num2str(R*1000) '_CO' num2str(C*1000) '.mat'],'elev_25','min_25','orgAL_25','orgAT_25')
