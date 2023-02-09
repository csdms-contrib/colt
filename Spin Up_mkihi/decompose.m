function [compaction,Fd] = decompose(x_m,x_f)

global yr
global organic_dep_autoch
global elevation
global B
global mui
global mki
global rhoo

compaction = zeros(1,B);
Fd=0;

%Decompose the marsh sediment
for x = x_m:x_f-1 %Loop through each marsh and upland cell in the domain
    for tempyr = yr:-1:1 %Loop through each pocket of sediment in each cell, starting at the most recently deposited packet of sediment at the surface
        depth = elevation(yr,x) - elevation(tempyr,x); %Depth of sediment pocket below the surface
        if depth > mui %Maximum depth at which decomposition occurs
            decomp(tempyr)=0;
            break
        else
            decomp(tempyr) = organic_dep_autoch(tempyr,x)*(mki * exp(-depth/mui));
            organic_dep_autoch(tempyr,x) = organic_dep_autoch(tempyr,x) - decomp(tempyr);
        end
    end
    compaction(x)=sum(decomp)/1000/rhoo; %Total compact in a given cell is a result of the sum of all decomposition in that cell
    Fd=Fd+sum(decomp); %[kg] Flux of organic matter out of the marsh due to decomposition
    clear decomp
end

%Decompose for forest sediment
% for x = x_f:B
%     for tempyr = yr:-1:1 %Loop through each pocket of sediment in each cell, starting at the most recently deposited packet of sediment at the surface
%         depth = elevation(yr,x) - elevation(tempyr,x); %Depth of sediment pocket below the surface
%         if depth > fui %Maximum depth at which decomposition occurs
%             decomp=0;
%             break
%         else
%             decomp(tempyr) = organic_dep_autoch(tempyr,x)*(fki * exp(-depth/fui));
%             organic_dep_autoch(tempyr,x) = organic_dep_autoch(tempyr,x) - decomp(tempyr);
%         end
%     end
%     compaction(x)=sum(decomp)/1000/rhoo; %Total compact in a given cell is a result of the sum of all decomposition in that cell
%     clear decomp
% end