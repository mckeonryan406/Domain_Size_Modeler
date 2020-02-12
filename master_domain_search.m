function[modeled_data_OUT] = master_domain(data_file, kinetics_file, domains_file)

% Probably want things like nloops to try, domain range... 


% Takes in observed data and models the Fcum, lnDa2, and delta according to the given domain
% structure --> output file described below

% data_in = observed data array
%   col 1 = step #
%   col 2 = Temp (deg C)
%   col 3 = time (hr)
%   col 4 = Fcum observed
%   col 5 = 10000/K
%   col 6 = lnDa2 observed
%   col 7 = ln(r/ro) observed ---  here referred to as delta observed

% kinetics_in = array of observed kinetics that control modeling
%   value 1 = Ea in kcal/mol!
%   value 2 = lnDoa2 ---> important because this ties the reference size
%   value 3 = reference domain size for these kinetics

% domains_in = array of domain data for modeling F
%   col 1 = size of domain relative to ref_a
%   col 2 = gas fraction contained in domain

% OUTPUT 

% generates a text file with the following

%   row 1 = misfit score, BIC (bayesian information criterion) Ea, lnDoa2, slope
%   row 2 = Relative domain size for each domain modeled
%   row 3 = Gas fraction of each domain modeled
%   row 4 - n      
%       col 1 = step number
%       col 2 = Temp (deg C)
%       col 3 = time (hr)
%       col 4 = Fcum modeled
%       col 5 = 10000/K
%       col 6 = lnDa2 modeled
%       col 7 = delta modeled --> ln(r/ro) in Lovera-speak
%       col 8 = delta OBSERVED --> included for easy ploting comparison

% now on to modeling!


% unpack Input Data data 

% unpack observed step heating data (data_in)
data_in = load(data_file);  % load data from file as described above
Temp = data_in(:,2);
time = data_in(:,3);
Fcum_obs = data_in(:,4);
delta_obs = data_in(:,7);

% unpack kinetics data (kinetics_in)
kinetics_in = load(kinetics_file);
Ea = kinetics_in(1);
lnDoa2 = kinetics_in(2);
ref_a = kinetics_in(3);

% unpack domain data (domains_in) -- NOTE for INVERSE modeling the domain
% structure, this variable will be what changes.  FORWARD modeling the gas
% release from a single domain structure, this will not change.

domains_in = load(domains_file);

% Calculate Slope of input Arrehius Data

R = 1.987; % Gas Constant for kcal/mol
slope = -Ea/R/10;  % calculate the slope of a regression line given Ea and lnDoa2 --- this is used for calculating delta (ln(r/ro).

% *****************************************************************
% LOOP FOR SEARCH WILL START HERE with new domain data

% calc F for modeled domains
data = [];
data(:,1) = Temp;
data(:,2) = time;

% Generate New Domain Structure Here


% NOW Calculate the results with new domain stucture

Fcum_domains = [];

for i = 1:length(domains_in)
    
    a = domains_in(i,1);
    Fcum_domains(:,i) = calc_F_1_domain(data,Ea,lnDoa2,a);
end

% calculate cummulative F for all domains with weighted gas fracations and
% combine
Fcum_mod = weight_domains(Fcum_domains,domains_in);

% calc arrhenius and delta (ln(r/ro)) for modeled Fcum

data(:,3) = Fcum_mod; % data is the growing array to store results of model

arrDATA = arrhenius_calc(data);
delta = delta_calc(arrDATA,Ea);

data(:,4) = arrDATA(:,1);
data(:,5) = arrDATA(:,2);
data(:,6) = delta;

% calculate misfit between observed and modeled F values

%   Following Farely and Flowers, 2012 (EPSL - supplemental material) I
%   equaly weight the misfit between the step F and cummulative F for each
%   step.  I normalize the misfit by the number of steps.

% calculate F step release for obs and mod data
Fstep_obs = [];
Fstep_mod = [];

Fcum_obs = transpose(Fcum_obs);
Fstep_obs(1) = Fcum_obs(1);
for i = 2:length(Fcum_obs)
    Fstep_obs(i) = Fcum_obs(i) - Fcum_obs(i-1);
end

Fstep_mod(1) = Fcum_mod(1);
for i = 2:length(Fcum_mod)
    Fstep_mod(i) = Fcum_mod(i) - Fcum_mod(i-1);
end

cum_misfit = abs(Fcum_obs - Fcum_mod)*100;
step_misfit = abs(Fstep_obs - Fstep_mod)*100;

raw_misfit = sum(cum_misfit) + sum(step_misfit);
final_misfit = raw_misfit/length(Fcum_obs)

% Calculate the Bayesian Information Criterion (BIC)

% here BIC is defined by the following equation: 
% BIC = misfit + degrees of freedom * ln(parameters being fit)

nsteps = length(Fcum_mod);
ndomains = length(domains_in);

L = final_misfit;
n = 2 * nsteps; % parameters being fit
K = 2 * ndomains - 1; % degrees of freedom

BIC = L+K*log(n) % remember log(x) in MATLAB is the natural log (ln(x))!


% LOOP for SEARCH will END HERE
%**********************************************************************

% clean up and prepare output

mod_data=[];


if ndomains < 8 % make sure output array is big enough
    cols = 8;
else
    cols = ndomains;
end
modeled_data_OUT = zeros(nsteps+3,cols); % create empty array for output

first_row = [final_misfit,BIC,ndomains,Ea,lnDoa2,slope];
modeled_data_OUT(1,1:6) = first_row;     % add misfit result and kinetic data used for modeling 
modeled_data_OUT(2,1:ndomains) = domains_in(:,1); % add domain size
modeled_data_OUT(3,1:ndomains) = domains_in(:,2); % add gas fraction

for i = 1:length(Fcum_mod)
    mod_data(i,1) = i; % make step number for output
end
mod_data(:,2:7) = data; %block of modeled data results
mod_data(:,8) = data_in(:,7); % put the observed delta in too for easy comparison 
modeled_data_OUT(4:nsteps+3,1:8) = mod_data;

% make plot --> delta observed vs. delta modeled
cla
% build stairs for OBSERVED data
Fcum = data_in(:,4);
r_ro = data_in(:,7);
for a = 1:length(Fcum)-1;
        x(2*a-1) = Fcum(a);
        x(2*a) = Fcum(a+1);
        y(2*a-1) = r_ro(a);
        y(2*a) = r_ro(a);
end;
plot(x,y,'k','LineWidth',2);
hold on

% build stairs for MODELED data

Fcum = mod_data(:,4);
r_ro = mod_data(:,7);
for a = 1:length(Fcum)-1;
        x(2*a-1) = Fcum(a);
        x(2*a) = Fcum(a+1);
        y(2*a-1) = r_ro(a);
        y(2*a) = r_ro(a);
end;
plot(x,y,'r');


xname = 'Sum 3He Released';
yname = 'Delta';
title('Modeled Delta (red) vs. Observed Delta (black)')
xlabel(xname,'fontsize',12);  %axis labels
ylabel(yname,'fontsize',12);
xlim([0 1]);
box on  


save modeled_data_OUT.txt modeled_data_OUT -ascii -tabs







