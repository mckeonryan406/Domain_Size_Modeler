%function[arrDATA] = arrhenius_calc(data_in)

data = load('MI43d2_example_Arrhenius_Calc_StandAlone.txt');

%  This function calculates the arrhenius data for an initially 
%  square profile using the equations from F+K 19XX

% The inputs for this function are a three column matrix of step heating 
% data with the following:
%   col 1 = step number
%   col 2 = step Temp (deg C)
%   col 3 = step duration (hr)
%   col 4 = step F (0-1) not %

% This function returns a two column array with 10000/K and LnDa^2 for each
% step.

time_s = data(:,3)*3600;
TK = data(:,2) + 273.15;
Fcum = data(:,4);
arrDATA =[];
lnDa2 =[];

% trap out Fcum = 1

for i = 1:length(Fcum)
    if Fcum(i)==1
        Fcum(i) = 0.999999;
    end
end


% Square Profile (1st step)

F = Fcum(1);
t = time_s(1);

if F < 0.1
    lnDa2(1) = log((F^2*pi)/(36*t)); % NOTE this is natural log!!  

elseif F > 0.9
    a = 1/((pi^2)*t);
    b = log((pi^2/6)*(1-F));
    lnDa2(1) = log(-a*b);
    
else
    a = 1/((pi^2)*t);
    b = (2*pi)-((pi^2/3)*F);
    c = (2*pi)*(sqrt(1-(pi/3))*F);
    lnDa2(1) = log(a*(b-c));
end

% Rounded Profile -- rest of steps

for i = 2:length(Fcum)
    
    F = Fcum(i);
    Fbefore = Fcum(i-1);
    t = time_s(i);
    
    if F < 0.1
        lnDa2(i) = log((((F^2)-(Fbefore^2))*pi)/(36*t));

    elseif F > 0.9
        lnDa2(i) = log((1/((pi^2)*t))*log((1-Fbefore)/(1-F)));
    
    else
        a = 1/((pi^2)*t);
        b = (-1*(pi^2/3))*(F-Fbefore);
        c = (2*pi)*((sqrt(1-((pi/3)*F)))-(sqrt(1-((pi/3)*Fbefore))));
        lnDa2(i) = log(a*(b-c));
    end
end
    

% calc 10000/T
for i = 1:length(Fcum)
    arrDATA(i,1) = 10000/TK(i);
end

arrDATA(:,2) = lnDa2;



% Plot Arrhenius Data with Step Number Labels
labels = cellstr(num2str(data(:,1)));
plot(arrDATA(:,1),arrDATA(:,2),'ko','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',3,'LineStyle','-','LineWidth',0.5)
text(arrDATA(:,1),arrDATA(:,2),labels,'VerticalAlignment','bottom','HorizontalAlignment','left')
text(11,-34,'Enter the Step Numbers to Regress in the Command Window!','Color','red','FontSize',14)
xlim([10 33])
ylim([-35 -5])
title('Observed Arrhenius Data -- Select Steps to Regress for Kinetics','FontSize',16)
xlabel('10000/K','FontSize',14)
ylabel('ln(D/a2)','FontSize',14)

% Steps to Regress for Kinetics  -- Get from USER INPUT from the command
% window
prompt = 'Enter the step numbers you would like to regress separated by spaces: ';
x = input(prompt, 's')
regress_steps = str2num(x)

% Store Arrhenius Data for steps to regress
regressX = arrDATA(regress_steps,1);
regressY = arrDATA(regress_steps,2);

% run a linear Regression Model to get necessary info
linMod = fitlm(regressX,regressY);

% linear model results
Fit_lnDoa2 = linMod.Coefficients.Estimate(1)
Fit_slope = linMod.Coefficients.Estimate(2)
Fit_Ea = -1 * Fit_slope* 1.987 * 10   % Ea in kcal/mol!!
Fit_R2 = linMod.Rsquared.Ordinary

% build text output
a = num2str(Fit_Ea);
b = num2str(Fit_lnDoa2);
c = num2str(Fit_R2);
aa = strcat('Ea = ',' ',a,' kcal/mol');
bb = strcat('Ln(Do/a2 =',' ',b);
cc = strcat('R squared =',' ',c);

% plot the data and the steps to regress with the results of regession
cla % clear the old plot with all the labels
plot(arrDATA(:,1),arrDATA(:,2),'ko','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',3,'LineStyle','-','LineWidth',0.5)
xlim([10 33])
ylim([-35 -5])
title('Observed Arrhenius Data with Kinetics Regression','FontSize',16)
xlabel('10000/K','FontSize',14)
ylabel('ln(D/a2)','FontSize',14)
hold on
plot(regressX,regressY,'ro','MarkerFaceColor','r','MarkerSize',6)
lsline  % plot a linear regession line through the regress_steps
text(20,-7,aa,'FontSize',14)
text(20,-10,bb,'FontSize',14)
text(20,-13,cc,'FontSize',14)







