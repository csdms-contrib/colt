function calcDMequil(distance)

%Function to find the initial condition for marsh depth at which the model
%will experience the lowest rate of change. Outputs an array of marsh depths
%for a given combination of rate of sea level rise and external sediment
%supply. Will produce a different array for each fetch and wind speed
%conditions.

close all

fetch = 5000; %mudeflat width [m]
wind = 6; %reference wind speed [m/s]

filename = ['Fetch' num2str(fetch) '_Wind' num2str(wind)];

if exist(filename) == 0
    mkdir(filename)
end

B = 10000; %basin width [m]
tcr=0.1; %tau cr mudflat [Pa]
lamda=0.0001; % mudflat erodability coeff [-]
rhos=1000;rhom=1000;%bulk densities [kg/m3]
amp=1.4/2; % tidal amplitude [m]
Dmax=.7167*2*amp-.483; %Maximum depth for vegetation below HWL
Dmin=0;
tr = amp*2; %Tidal range
numiterations=500; %DCW
P=12.5*3600*1; % tidal period [s]
dt=P/numiterations; %DCW
ws=0.05*10^-3; %effective settling velocity, m sec-1, Marani 2010 uses 0.2 mm s-1, Mudd 2009 uses .037 mm/s
msl = 0;
BMax=350;% g/m2
nuGp=.0138;
por=1000/2650;
chiref=0.15;
timestep=365*(24/12.5); %DCW


dfs = 1.7:.01:2.69;

for i = 1:numel(dfs)
    df = dfs(i);
    fac=min(1,df/(2*amp));
    Df=(df+(df-fac*2*amp))/2;
    tw=  wavetau(fetch,wind,Df,B);
    tau=max((tw-tcr)/tcr,0)*lamda;
    Cr(i)=rhos*tau/(1+tau);
end

fh1 = figure;
hold on
plot(dfs,Cr)

Ci = Cr*exp(-0.0031*distance);
plot(dfs,Ci,'-r')

legend('SSC @ Marsh Edge (ME)','SSC @ 500m from ME')
xlabel('Bay Depth (m)','FontSize',15)
ylabel('SSC at Marsh Edge (kg/m^3)','FontSize',15)
saveas(fh1,[filename '/SSC at Marsh Edge.fig'])
print('-dpng',fh1,[filename '/SSC at Marsh Edge.png'])

%%%%%
fh2 = figure;
hold on
xlabel('Depth below HWL (m)','FontSize',15)
ylabel('Accretion (mm/yr)','FontSize',15)
dm_equil = zeros(15,numel(Ci)); %Create an array for storing the equilibrium marsh depth for a given bay depth (and associated Ci) and RSLR.
for i = 1:numel(Ci)
    [temp_dm_equil] = accretion(distance,Ci(i),amp,numiterations,tr,dt,P,msl,Dmax,Dmin,BMax,ws,timestep,i);
    dm_equil(1:15,i) = temp_dm_equil;
end
xlim([0 Dmax])
saveas(fh2,[filename '/ACCvDEP.fig'])
print('-dpng',fh2,[filename '/ACCvDEP.png'])
% save([filename '/Equilibrium Marsh Depth.mat'],'dm_equil')


%%%%%%%%Convert from bay depth to Css%%%%%%%%%
load([filename '/Equilibrium Bay Depth'])

db_eq = round(100*db_eq)/100;
dm_eq = zeros(15,15);
size(dm_equil)
for Co = 1:15
    for RSLR = 1:15
        dm_eq(Co,RSLR) = dm_equil(RSLR,round((db_eq(RSLR,Co)-1.69)*100));
    end
end

fh3 = figure;
xlabel('RSLR (mm/yr)','FontSize',15)
ylabel('Co (kg/m^3)','FontSize',15)
hold on
surf(dm_eq)
set(gca,'YLim',[1 15],'XLim',[1 15],'YTick',1:15,'YTickLabel',10:10:150,'XTick',1:15)
colorbar

saveas(fh3,[filename '/@' num2str(distance) 'm DM_phase.fig'])
print('-dpng',fh3,[filename '/@' num2str(distance) 'm DM_phase.png'])
save([filename '/Equilibrium Marsh Depth.mat'],'dm_eq')