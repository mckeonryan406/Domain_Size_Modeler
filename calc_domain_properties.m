function[domain_model_out_best] = calc_domain_properties(obs_kinetics, domain_model_out_best)

% This function calculates physcial domain size in microns, LOG Do/a2 and closure temperature for each 
% modeled domain.  The cooling rate is set to 10 degC/Myr  


% unpack input data
Ea = obs_kinetics(1);
lnDoa2 = obs_kinetics(2);
ref_a = obs_kinetics(3);
R = 1.987; % gas constant
Farley_Do = exp(-0.66);  % published Do for Hematite from Farley 2018 GCA

% sort modeled domains by size
ndomains = domain_model_out_best(1,3);
domains = [];
domains(1,:) = domain_model_out_best(2,1:ndomains); % size
domains(2,:) = domain_model_out_best(3,1:ndomains); % gas fraction

[y,i] = sort(domains(1,:));
domains_sorted = domains(:,i);

domain_model_out_best(2,1:ndomains) = domains_sorted(1,:);
domain_model_out_best(3,1:ndomains) = domains_sorted(2,:);

% get the relative domain size for each modeled domain

best_domains_size = domain_model_out_best(2,1:ndomains);

% Calc Input Do
Do = exp(lnDoa2)/ref_a^2;

% Put empty rows into output array to catch new domain properties

size_input = size(domain_model_out_best);
ncols = size_input(2);
updater = zeros(2, ncols);
domain_model_out_best = insertrows(domain_model_out_best, updater, 3);

% Calc Do/a2 and physical domain size for each modeled domain

mod_domain_prop = [];
for i = 1:length(best_domains_size)
    Do_domain = Do/best_domains_size(i)^2;  % calc Do/a2 for Tc calc below
    domain_physical_size = sqrt(Farley_Do/Do_domain);  % calc pysical domain size in cm's
    
    mod_domain_prop(1,i) = Do_domain;
    mod_domain_prop(2,i) = domain_physical_size * 10000; % convert physical size to microns
end

% Calc closure temp for each domain

% calc cooling rate --- convert from degC/Myr to degC/sec

cool_rate = 10; % degC/Myr
cool_sec = cool_rate/1e6/365.25/24/3600; % years/days/hours/seconds
Ea_cal_mol = 1000*Ea; % convert to cal/mol

for i = 1:length(best_domains_size)
    Doa2 = mod_domain_prop(1,i);
    guess_temp_start = 400; % starting Tc guess in Kelvin
    guess_temp2 = (Ea_cal_mol/R)/(log((Doa2*guess_temp_start^2*55*R)/(Ea_cal_mol*cool_sec)));  % Iterate through the Tc calcuation
    guess_temp3 = (Ea_cal_mol/R)/(log((Doa2*guess_temp2^2*55*R)/(Ea_cal_mol*cool_sec)));
    Tc = (Ea_cal_mol/R)/(log((Doa2*guess_temp3^2*55*R)/(Ea_cal_mol*cool_sec))) - 273.15;  % convert from Kelvin to deg C
    
    % store results in correct location
    domain_model_out_best(2,i) = mod_domain_prop(2,i);  % overwrite Relative Domain Size with Physcial Size in Microns
    domain_model_out_best(4,i) = log10(Doa2); % calc log(Do/a2) for MDD modeling in QTQt
    domain_model_out_best(5,i) = Tc;
end







