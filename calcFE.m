function [FE_org,FE_min] = calcFE(bfoc,bfop)

global elevation
global yr
global organic_dep_autoch
global organic_dep_alloch
global mineral_dep
global rhos

%Function to calculate the flux of organic matter (FE_org) and the flux of
%mineral sediment (FE_min) from the marsh to the bay,using the fetch for
%the current year (bfoc) the fetch for the previous year (bfop) and the
%stragtigraphy of organic and mineral deposition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%Calculate OM eroded from the marsh platform
organic_dep = organic_dep_autoch + organic_dep_alloch; %calculate total organic deposition
pyr = yr-1; %Previous year
E=bfoc-bfop; %Amount of erosion b/t the previous yr and the current yr

x_m1 = ceil(bfop); %first marsh cell to erode
x_m2 = ceil(bfoc); %last marsh cell to erode

if  E <= 0 
    FE_org = 0; %[g] If there is no erosion, no OM is eroded
    FE_min = 0; %[g] If there is no erosion, no OM is eroded
elseif x_m1 == x_m2 %If actively eroding marsh edge has not changed since preious year
    if  elevation(pyr,1) < elevation(1,x_m1) %If depth of erosion is below the lowest marsh deposit
        us=elevation(1,x_m1)-elevation(yr-1,1); %[m] Depth of underlying stratigraphy
        usmass=us*rhos*1000; %[g] Mass of pure mineral sediment underlying marsh at marsh edge
        FE_org = sum(organic_dep(1:pyr,x_m1))*E; %[g] OM eroded is equal to the total amount of OM in the eroding marsh edge (including both initial deposit and OM deposited since the model run began) times the fraction of the marsh edge cell that is eroded
        FE_min = sum(mineral_dep(1:pyr,x_m1))*E+usmass; %[g] MIN eroded is equal to the total amount of OM in the eroding marsh edge (including both initial deposit and OM deposited since the model run began) times the fraction of the marsh edge cell that is eroded
    else %If depth of erosion is less than marsh deposits
        boundyr = find(elevation(:,x_m1) > elevation(1,pyr),1); %Year at which deposit is above the boundary for the depth of erosion
        FE_org = sum(organic_dep(boundyr:pyr,x_m1))*E; %[g] Total mass of OM deposited in the marsh edge cell
        FE_min = sum(mineral_dep(boundyr:pyr,x_m1))*E; %[g] Total mass of MIN deposited in the marsh edge cell
    end
else
    ecells = x_m2 - x_m1; %Number of cells eroded
    Hfrac_ero=ones(1,ecells);
    Hfrac_ero(1) = x_m1 - bfop; % Horizontal fraction of previous marsh edge that is eroded
    Hfrac_ero(end) = bfoc-floor(bfoc); %Horizontal fraction of current marsh edge that is eroded    
    
    FE_org=0;
    FE_min=0;
    i=1;
    for x_m = x_m1:x_m2 %for cells between previous marsh edge to current marsh edge
        if  elevation(pyr,1) < elevation(1,x_m) %If depth of erosion is below the lowest marsh deposit
            us=elevation(1,x_m)-elevation(yr-1,1); %[m] Depth of underlying stratigraphy
            usmass=us*rhos*1000; %[g] Mass of pure mineral sediment underlying marsh at marsh edge
            FE_org = FE_org + sum(organic_dep(1:pyr,x_m))*Hfrac_ero(i); %[g] OM eroded from previous marsh edge cell
            FE_min = FE_min + sum(mineral_dep(1:pyr,x_m))*Hfrac_ero(i) + usmass; %[g] MIN eroded from previous marsh edge cell
        else %If depth of erosion is less than marsh deposit
            boundyr = find(elevation(:,x_m) > elevation(1,pyr),1); %Year at which deposit is above the boundary for the depth of erosion

            FE_org = FE_org + sum(organic_dep(boundyr:pyr,x_m))*Hfrac_ero(i); %[g] OM eroded from previous marsh edge cell
            FE_min = FE_min + sum(mineral_dep(boundyr:pyr,x_m))*Hfrac_ero(i); %[g] MIN eroded from the marsh edge cell
        end    
    end
end