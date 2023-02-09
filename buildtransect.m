function [B,dfo,elevation] = buildtransect(R,C,slope,mwo,elev_25)

global amp
global wind
global bfo
global filename
global endyear

spinupdur=min(size(elev_25)); %determine length of spinup (this is not the most flexible code if you change the marsh width)
% Dmin=0; %these are not used if there is no ramp
% Dmax=.7167*2*amp-.483; %these are not used if there is no ramp

fh1 = figure; %create figure to visualize initial domain configuration
hold on

filename1 = ['Calibrate init conditions/Fetch' num2str(bfo) '_Wind' num2str(wind)]; %name of file with initial conditions to load (in subsequent lines)

%%%% Determine initial bay depth, such that change in depth will be small %%%%
if exist(filename1) == 0
    warning('Initial conditions have not been calibrated for these paramemters')
    dfo = 2; %set bay depth to 2 m
else %if there is a initial conditions file, load those
    load([filename1 '/Equilibrium Bay Depth.mat'])
    load([filename1 '/Equilibrium Marsh Depth.mat'])    
    if C/10 > size(db_eq,1) || R > size(db_eq,2) %if requested model run is outside the range of initial conditions that were calibrated (10-100 mg/L or SLR 1-15 mm/yr)
        warning('Initial conditions have not been calibrated for these paramemters 1') %give this error
        dfo = 2; %set bay depth equal to 2 m
    elseif C/10 < .5 || R < .5 %if SSC is less than 5 mg/L or SLR is less than 0.5 mm/yr
        warning('Initial conditions have not been calibrated for these paramemters 2') %give this error
        dfo = 2; %set bay depth equal to 2 m
    else
        dfo = db_eq(round(C/10),round(R)); %rounds both SSC and SLR inputs to integer, then uses lookup table to set the bay depth
    end
end

x_m = bfo+1; %marsh edge position is the width of the bay (i.e. the fetch) plus 1
%ramp_width=ceil((Dmax-Dmin)/slope); % only used if adding a slope between
%marsh and forest

maxY = R/1000 * endyear + amp + .5; %gives maximum elevation that the water level will reach by the end of the scenario (assumes a constant rate of SLR) [m]
max_potentialwidth = ceil(maxY/slope);%convert elevation change of water level to the width of forest that will be invaded by water/become marsh [m]
%upland_width = 3500; %If you want, you can set an upland (i.e. forest
%width manually), turn off next line if you want to do this
upland_width=ceil(maxY/slope); %automatically calculates the exact amount of forest needed to allow the marsh to migrate 

if max_potentialwidth > upland_width %if you set forest width manually, this checks to make sure the domain is large enough 
    disp('WARNING! Slope/Sea level rise conditions are such that the model domain is too small for the scenario. Adjust the upland width accordingly.')
end

B = bfo + mwo + upland_width; %Total domain width [m], and also number of cells in domain each with 1-m width
%%%B = bfo + mwo + upland_width+ramp_width; %only use if you add in a ramp

x = 1:B; %X-position of each cell in model domain
elevation = zeros(endyear,B); %create empty matrix to store elevation values
elevation(1:spinupdur,1:x_m-1) = amp-dfo; %Bay depth for spinup duration is at equilibrium
elevation(1:spinupdur,x_m:x_m+mwo-1) = elev_25-(spinupdur*(1/1000)); %Marsh elevation comes from model spin up, adjusted to modern sea level

%Now form the underlying forest stratigraphy for future organic deposits at
%depth. 
modernslope = slope*(1:upland_width)+elevation(spinupdur,x_m+mwo-1); %create the elevation profile of the forest for the spinup

for i = 1:spinupdur %for the length of the spinup
    elevation(i,x_m+mwo:B) = modernslope-((spinupdur-i)*.025);% implement elevation profile of the forest -- it is always the same slope, creates layers that are 2.5 cm thick, one for each year of the spinup
    plot(x,elevation(i,:),'-','Color',[.5 .5 .5]) %plots initial stratigraphy
end

plot(x,elevation(spinupdur,:),'k-','LineWidth',1.5) %plots initial stratigraphy

set(gcf,'units','Inches','position',[1 2 12 5],'PaperPositionMode','auto')
set(gca,'FontSize',12)

text(B*.1,max(elevation(spinupdur,:))*.9,['RSLR = ' num2str(R) ' mm/yr'],'FontSize',12)
text(B*.1,max(elevation(spinupdur,:))*.7,['Wind Speed = ' num2str(wind) ' m/s'],'FontSize',12)
text(B*.1,max(elevation(spinupdur,:))*.8,['C = ' num2str(C) ' kg/m^3'],'FontSize',12)

text(B*.1,max(elevation(1,:))*.3,['d_b = ' num2str(round(100*dfo)/100) ' m'],'FontSize',12)

text(bfo*.4,0,'Fetch')
text(bfo*.4,-.2,[num2str(bfo) 'm'])
text(bfo+mwo*.2,amp+.25,'Platform')
text(bfo+mwo*.2,amp+.05,[num2str(mwo) 'm'])
text(bfo+mwo+upland_width*.3,elevation(spinupdur,round(bfo+mwo+upland_width*.3))+.8,'Upland Slope')
text(bfo+mwo+upland_width*.3,elevation(spinupdur,round(bfo+mwo+upland_width*.3))+.6,[num2str(slope)])

plot([x(x_m) x(x_m+mwo) x(x_m+mwo)],[elevation(spinupdur,x_m) elevation(spinupdur,x_m+mwo) elevation(spinupdur,x_m+mwo)],'xg','MarkerSize',10,'LineWidth',2)

ylabel('Elevation Relative to Initial Sea Level (m)')
xlabel('Distance (m)')

outputfilename = [filename 'Initial Surface'];
print('-dpng',outputfilename)
saveas(fh1,[outputfilename '.fig'])