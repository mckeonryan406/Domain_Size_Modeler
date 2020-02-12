function[modeled_data_OUT] = master_domain(obs_data, obs_kinetics, domains_in)

%% Need to switch this to take in data from obs_data and obs_kinetics

% Takes in observed data and models the Fcum, lnD/a2, and delta according to the given domain
% structure --> output variable described below

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
%   value 3 = reference domain size for these kinetics --> should be 1 for
%   easy math

% domains_in = array of domain data for modeling F
%   col 1 = size of domain relative to ref_a
%   col 2 = gas fraction contained in domain

% OUTPUT 

% generates a variable with the following

%   row 1 = misfit score, BIC (bayesian information criterion), ndomains, Ea, lnDo/a2, slope
%   row 2 = Gas fraction of each domain modeled
%   row 3 = Relative domain size for each domain modeled
%   row 4 - n      
%       col 1 = step number
%       col 2 = Temp (deg C)
%       col 3 = time (hr)
%       col 4 = Fcum modeled
%       col 5 = 10000/K
%       col 6 = lnDa2 modeled
%       col 7 = delta modeled --> ln(r/ro) in Lovera-speak
%       col 8 = delta OBSERVED --> included for plotting comparison
%       col 9 = lnDa2 OBSERVED --> included for plotting 

% now on to modeling!


% unpack Input Data data 

% unpack observed step heating data (obs_data)
step_num = obs_data(:,1);
Temp = obs_data(:,2);
time = obs_data(:,3);
Fcum_obs = obs_data(:,4);
delta_obs = obs_data(:,7);

% unpack kinetics data (obs_kinetics)
Ea = obs_kinetics(1);
lnDoa2 = obs_kinetics(2);
ref_a = obs_kinetics(3);

%get number of domains being modeled
x = size(domains_in);
ndomains = x(1);

% Calculate Slope of input Arrehius Data

R = 1.987; % Gas Constant for kcal/mol
slope = -Ea/R/10;  % calculate the slope of a regression line given Ea and lnDoa2 --- this is used for calculating delta (ln(r/ro).

% calc F for modeled domains
mod_data = [];
mod_data(:,1) = step_num;
mod_data(:,2) = Temp;
mod_data(:,3) = time;
Fcum_domains = [];

for i = 1:ndomains
    a = domains_in(i,1);
    Fcum_domains(:,i) = calc_F_1_domain(mod_data,Ea,lnDoa2,a);
end

% calculate cummulative F for all domains with weighted gas fracations and
% combine
Fcum_mod = weight_domains(Fcum_domains,domains_in);

% calc arrhenius and delta (ln(r/ro)) for modeled Fcum
mod_data(:,4) = Fcum_mod; % data is the growing array to store results of model

arrDATA = arrhenius_calc(mod_data);
delta = delta_calc(arrDATA,Ea,lnDoa2);

mod_data(:,5) = arrDATA(:,1);
mod_data(:,6) = arrDATA(:,2);
mod_data(:,7) = delta;

%% calculate misfit between observed and modeled F values

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
final_misfit = raw_misfit/length(Fcum_obs)*5;

%% Calculate the Bayesian Information Criterion (BIC)

% here BIC is defined by the following equation: 
% BIC = misfit + degrees of freedom * ln(parameters being fit)

nsteps = length(Fcum_mod);
ndomains = length(domains_in);

L = final_misfit;
n = 2 * nsteps; % parameters being fit
K = 2 * ndomains - 1; % degrees of freedom

BIC = L+K*log(n); % remember log(x) in MATLAB is the natural log (ln(x))!

% clean up and prepare output

mod_dataOUT=[];  % array for building output

if ndomains < 10 % make sure output array is big enough
    cols = 10;
else
    cols = ndomains;
end
modeled_data_OUT = zeros(nsteps+3,cols); % create empty array for output

first_row = [final_misfit,BIC,ndomains,Ea,lnDoa2,slope];
modeled_data_OUT(1,1:6) = first_row;     % add misfit result and kinetic data used for modeling 
modeled_data_OUT(2,1:ndomains) = domains_in(:,1); % add domain size
modeled_data_OUT(3,1:ndomains) = domains_in(:,2); % add gas fraction


mod_dataOUT = mod_data; %block of modeled data results
mod_dataOUT(:,8) = obs_data(:,7); % put the observed delta in too for easy comparison 
mod_dataOUT(:,9) = obs_data(:,6); % put in the observed ln(D/a2) as well
mod_dataOUT(:,10) = obs_data(:,4); % and awww heck, let's include the observed Fcum as well.
modeled_data_OUT(4:nsteps+3,1:10) = mod_dataOUT;







