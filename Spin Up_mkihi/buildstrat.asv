function buildstrat

%Function to take outputs from  model spinup and build a stratigraphy
%to use for all model runs

figure
hold on

%load in variables from the spinup
load('elev.mat')
load('mineral deposition.mat')
load('organic deposition.mat')

%Create matrices to save the consolidated data for each of elevation, 
%mineral deposition, autochthanous organic deposition, and allochthanous 
%organic deposition
elev_25=zeros(25,1000);
min_25=zeros(25,1000);
orgAT_25=zeros(25,1000);
orgAL_25=zeros(25,1000);

%Loop through to consolidate 500yrs of deposition into 25 20-yr packets
yr1=1;
for i = 1:25
    yr2 = i*12;

    elev_25(i,:) = elev(yr2,:);
    min_25(i,:) = sum(mineral_dep(yr1:yr2,:));
    orgAT_25(i,:) = sum(organic_dep_autoch(yr1:yr2,:));
    orgAL_25(i,:) = sum(organic_dep_alloch(yr1:yr2,:));
    
    yr1=yr2+1;
    
    plot(1:1000,elev_25(i,:),'k-')
end

save('MarshStrat_RSLR2_CO50.mat','elev_25','min_25','orgAL_25','orgAT_25')