% FINAL RESULTS PLOTTER 

% PURPOSE ----------
% Plot the results for a saved model run using the output .csv file.

% INPUTS -----------
% The correctly formatted .csv file from a previous model run.  See below
% for formatting if you need to build the csv from scratch.

%Input CSV File:
input_file = 'modeled_data_OUT.csv';
domain_model_to_plot = csvread(input_file);

% OUTPUTS --------------------------------------
% This script will populate a Matlab figure window with an Arrhenius Plot 
% plot and a Delta plot (Lovera's r/ro), and a CSV with all of the observed and
% modeled data for this result.  Here are specifics:

% Plot Figure:
% Arrhenius Plot -- Here you have observed data in black, modeled kinectics
% for each domain in red, and reference lines for domains of known physical
% size for hematite in grey (1 mm, 10 microns, and 100 nanometers).

% Delta Plot --  Here the modeled data  is in red and the observed data is
% in black with the number of domains, the misfit, and the BIC score
% indicated.

% Input Data File Formating --------------------- 
%   This is the output from Domain_Model_Master.m)
%
%   row 1 = misfit score, BIC (bayesian information criterion), ndomains, Ea, lnDo/a2, slope
%   row 2 = Estimated physical domain size (in microns) for each domain modeled 
%   row 3 = Gas fraction of each domain modeled
%   row 4 = log(Do/a2) for each domain --> for input to QtQT
%   row 5 = closure temperature for each domain (assuming a 10 degC/Ma cooling rate)
%   row 6 - n:      
%       col 1 = step number
%       col 2 = Temp (deg C)
%       col 3 = time (hr)
%       col 4 = Fcum modeled
%       col 5 = 10000/K
%       col 6 = lnDa2 modeled
%       col 7 = delta modeled --> ln(r/ro) in Lovera-speak
%       col 8 = delta OBSERVED --> included for plotting comparison
%       col 9 = lnDa2 OBSERVED --> included for plotting
%       col 10 = Fcum OBSERVED --> included for plotting


%% Delta Plot -- BIC Best Result
%make plot --> delta observed vs. delta modeled

FigHandle = figure('Position', [100, 100, 425, 750]);


% build stairs for OBSERVED data
x=[]; % clear out x and y
y=[];
data_in = domain_model_to_plot;
Fcum = data_in(6:end,10);
r_ro = data_in(6:end,8);
for a = 1:length(Fcum)-1;
        x(2*a-1) = Fcum(a);
        x(2*a) = Fcum(a+1);
        y(2*a-1) = r_ro(a);
        y(2*a) = r_ro(a);
end;
subplot(2,1,1)
figure1 = plot(x,y,'k','LineWidth',2);
hold on

textY = min(y)+1;

% build stairs for MODELED data

Fcum = domain_model_to_plot(6:end,4);
r_ro = domain_model_to_plot(6:end,7);
for a = 1:length(Fcum)-1;
        x(2*a-1) = Fcum(a);
        x(2*a) = Fcum(a+1);
        y(2*a-1) = r_ro(a);
        y(2*a) = r_ro(a);
end;

% get Misfit and BIC scores and nDomains
misfitNOW = num2str(domain_model_to_plot(1,1));
BICnow = num2str(domain_model_to_plot(1,2));
ndomainsNOW = num2str(domain_model_to_plot(1,3));

plot(x,y,'r','LineWidth',2);
xname = 'Sum 3He Released';
yname = 'Delta';
title('Modeled Delta (red) vs. Observed Delta (black)','FontSize',15)
xlabel(xname,'fontsize',12);  %axis labels
ylabel(yname,'fontsize',12);
xlim([0 1]);
text(0.7,textY,strcat({'Misfit = '},{misfitNOW}),'FontSize',12);
text(0.7,textY+1,strcat({'BIC = '},{BICnow}),'FontSize',12);
text(0.7,textY+2,strcat({'Domains = '},{ndomainsNOW}),'FontSize',12);
box on  


%% Arrhenius Plot -- Observed Data with desired result as lines
% calc slope
R = 1.987;
%Ea = domain_model_to_plot(1,4);

%from Farley 2018 GCA
Ea = 40.6 % ~170 kJ/mol converted to kcal/mol

slope = -Ea/R/10; % divided by 10 due to kcal/mol units of Ea

% generate ln(Do/a2) for reference diffusion domain sizes -->  using single
% crystal data of Ken's and your GSA poster (2015) Do estimate.
%Do = 0.1; % cm/s^2

% Updated August 2018 to use Published kinetics
Do = exp(-0.66) % cm/s^2 from Farley 2018 GCA - reported as ln Do = -0.66

dd_size = [0.00001, 0.001, 0.1]; % for 100 nm, 10 um, and 1 mm, numbers are relative to cm's
lnDoa2_ref=[];
for i = 1:length(dd_size)
    lnDoa2_ref(i) = log(Do/dd_size(i)^2);
end

%plot lines for reference domain sizes
X=[0,40];
Y=[0,0];

subplot(2,1,2)
for i = 1:length(dd_size)
    Y(1) = lnDoa2_ref(i);
    Y(2) = X(2)*slope + lnDoa2_ref(i);
    figure3 = plot(X,Y,'Color', [.7,.7,.7], 'LineWidth',1, 'LineStyle','--');
    hold on
end

% calc ln(Do/a2) for each domain
ndomains = domain_model_to_plot(1,3);
lnDoa2 = [];
for i = 1:ndomains
    logDoa2 = domain_model_to_plot(4,i);
    Doa2 = 10^logDoa2;
    lnDoa2(i) = log(Doa2);  % remember in Matlab log() calculates the Natural Log!
end

%plot lines for modled domains
X=[0,40];
Y=[0,0];
for i = 1:ndomains
    Y(1) = lnDoa2(i);
    Y(2) = X(2)*slope + lnDoa2(i);
    plot(X,Y,'Color',[1,0,0],'LineWidth',1);
    hold on
end


% plot observed arrhenius data
obs_data = [];
obs_data(:,1) = domain_model_to_plot(6:end,5);
obs_data(:,2) = domain_model_to_plot(6:end,9);

plot(obs_data(:,1),obs_data(:,2),'ko','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',3,'LineStyle','-','LineWidth',0.5)
xlim([5 35]);
ylim([-35 -5]);
xname = '10000/K';
yname = 'ln(D/a2) ln(S-1)';
title('Modeled (red) vs. Observed (black) with Reference (gray)','FontSize',15)
xlabel(xname,'fontsize',12);  %axis labels
ylabel(yname,'fontsize',12);
% label the reference lines
text(18.5,-34,'1 mm','FontSize',12,'Color',[.5 .5 .5]);
text(22.2,-34,'10 um','FontSize',12,'Color',[.5 .5 .5]);
text(26,-34,'100 nm','FontSize',12,'Color',[.5 .5 .5]);
box on  %draws a border around the plot

% %% Save Outputs
% % build file strings
% plot_file = strcat(filename_prefix,'_plots.pdf');
% data_file = strcat(filename_prefix,'_modeled_data.xls');
% 
% % save files
% xlswrite(data_file, domain_model_out_best)
% saveas(FigHandle,plot_file);

