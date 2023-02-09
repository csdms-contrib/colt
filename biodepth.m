function [CC]=biodepth(dddd,tempbgb,x_m,x_f)
global yr
global organic_dep_autoch
global elevation
global B
global mui
global mki
global rhoo

 AA=[];
 for x = x_m:x_f %for all cells in the marsh
     for cohort = yr:-1:2 %starting at the surface and increasing the depth
         topdepth=dddd(cohort,x-(x_m-1)); %calculate the top of each sediment pocket
         bottomdepth=dddd(cohort-1,x-(x_m-1)); %calculate the bottom of each sediment pocket
         if bottomdepth>mui  %if the bottom depth of the sediment pocket is below the rooting depth/max depth for decomp, OM in that layer is zero
             AA(cohort,x)=0;
         elseif bottomdepth<topdepth %if the top depth is less than the bottom depth (due to numerical rounding errors), OM in that layer is zero
             AA(cohort,x)=0;
         else %if soil cohort is above depomp depth, 
             g=0.27*mui; %calculate constant to define exponential decay with depth
             bro=tempbgb(x-(x_m-1))./g; %calculate constant for each x position in the marsh (the total amount of biomass in each cell, divided by a constant
             fun1=@(depth) bro.*exp(-depth/g); %Defines the function that distributes biomass exponentially with depth in soil profile 
             AA(cohort,x)=integral(fun1,topdepth,bottomdepth,'ArrayValued',true); %take the integral of the function from the top of each soil cohort to the bottom of each soil cohort for each marsh cell, set as OM for that layer
         end
     end

 end
 CC=organic_dep_autoch.*0; %create empty matrix to match dimensions
 CC(1:yr,1:length(AA))=AA; %put OM (i.e. matrix AA) into empty matrix CC, which matches the correct dimensions to add to previous autochthonous OM deposition in the function transect
end
        
