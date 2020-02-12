%% INVERSE DOMAIN SEARCH FOR MDD ANALYSIS

% Created By: Ryan McKeon
% Date: 10 March 2017
% Contact: ryan.e.mckeon@dartmouth.edu

% Updated 13 August 2018 --> Updated hematite diffusion kinetics following
% Farley  2018 GCA "Helium diffusion parameters of hematite from a
% single-diffusion-domain crystal

%% PURPOSE --------------------------------
% This code generates a diffusion domain structure for Multiple Diffusion
% Domain analysis from step heating data (either 4He/3He or Ar/Ar). At 
% present, this code uses a SPHEREICAL diffusion domain only and the
% analytical solution for which is from Farley et al., 2010 G-cubed.  This
% code uses the Bayesian Information Criterion score (BIC) to balance the 
% improvement of a model fit vs. the number of domains used to generate
% said fit.  I.E. is it better to use fewer domains and introduce less 
% complexity to the model?  In the end you can decide, but the results here
% may suprise you.  Define the input file immediatedly below, read on
% further down to learn how this code works and what the structure of the
% inputs and outputs are.

% Keep reading past the inputs!!

%% Main User Inputs:
% INPUT DATA FILE -----------------------
obs_data = load('MI43d2_Example.txt');   % file with observed step heating data

% HOW MANY DOMAINS? 
domain_range = [8] % You can enter a single number of domains to fit your data with
%domain_range = [2,3,4,6,9]; % or you can give it a lot and it will run in a loop.

% NOTE:  It is best to run a wider range of domains numbers to start, then when you are 
% more savvy about the data set, run a longer model on just the specfic
% number of domains you wish to use.  Change "n_model_loops" below to a
% higher number to try more attempts at fitting for a single number of
% domains.

%% KNOBS --------------------------------
% How many times to try each set of domains?
n_model_loops = 5; % number of times to run full model to avoid local minima in search (5 is a good number)

% no real need to change the values below - this controls the search 
n_hone_loops = 30; % total number of honing loops after random search start
adjust_loops = 50; % number of pairs of domains to adjust for gas fraction on each hone loop
size_adjuster = 1.4; % magnitude of size adjustment per step in loop
gas_adjuster = 1.15; % magnitude of gas adjustment per step in loop


%% HOW THE CODE WORKS ---------------------

% The the Readme.m file!!!

%% Inputs and Outputs: ------------------------------------

% Inputs:
% obs_data_file = observed data array - Defined Above!
%   col 1 = step #
%   col 2 = Temp (deg C)
%   col 3 = time (hr)
%   col 4 = Fcum observed

% Outpus: 
% This generates a CSV file called "modeled_data_OUT.csv", which stores the results
% from the model iteration with the lowest BIC score for the whole modeling session.  
% I.E. if you try a wide range of differnt numbers of domains, only one single combination of 
% goodness of fit with allowed complexity will be saved.  This file has the following structure:

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

% Also generates plots -- The running Delta plot shows the progress of each
% iteration as the fit is improved.  At the conclusion of the code, a pdf with
% plots of Delta (model vs. observed), arrhenius data (observed with modeled and
% reference domains of 100 nm, 10 um, and 1 mm radius) and a comparison of best 
% BIC scores for different numbers of domains are saved as
% "domain_model_plots.pdf".

% IMPORTANT ------
% Running this code as a script (by hitting the green play botton above
% rather than as a function from the command window) saves the in-run
% variables to the Workspace.  This is recommended for this code because it
% allows you to quickly plot up the final results for any specific number
% of domains from the previous model run using Final_Results_Plotter.m
% instead of just getting the iteration with the lowest BIC score as this
% code is hardwired to do.

% On to the code!


%% Set up Variables and Running Plot Window

% block annoying warnings from the command window
warning('off','all')
warning

% set up figure window for Arrhenius Data and Evolving Delta plot
FigHandle = figure('Position', [100, 100, 1100, 400]);

% Create (or clear out) Arrays to catch results 
Master_Results_Table = [];
BIC_log = [];


%% Arrhenius Data Calculation ---------------   
% Here we will get the diffusion kinetics that we will need

% calculate data for arrhenius plot
arrDATA = arrhenius_calc(obs_data);

% add arrhenius data to observed data array
obs_data(:,5) = arrDATA(:,1);  % 10000/K
obs_data(:,6) = arrDATA(:,2);  % ln(D/a2)

% Plot Arrhenius Data with Step Number Labels

subplot(1,2,1)
labels = cellstr(num2str(obs_data(:,1)));  % get the step number for labels
plot(arrDATA(:,1),arrDATA(:,2),'ko','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',3,'LineStyle','-','LineWidth',0.5)
text(arrDATA(:,1),arrDATA(:,2),labels,'VerticalAlignment','bottom','HorizontalAlignment','left')
text(11,-32.5,'Choose Step Numbers For Points To Regress','Color','red','FontSize',14)
text(11,-34,'Enter the numbers in the Command Window','Color','red','FontSize',14)
xlim([10 33])
ylim([-35 -5])
title('Observed Arrhenius Data -- Select Steps to Regress for Kinetics','FontSize',16)
xlabel('10000/K','FontSize',14)
ylabel('ln(D/a2)','FontSize',14)

%% Get User Input for Kinetics ----------
% Steps to Regress for Kinetics  -- Get from USER INPUT from the command
% window
fprintf('\n\n')
fprintf('DOMAIN MODELER -- by Ryan McKeon\n\n')
fprintf('Time to calculate the diffusion kinetics.\n\n')
disp('If you know the Activation Energy (in kcal/mol')
fprintf('and LnDo/a2 you want to use, you can enter it below.\n\n')
disp('Otherwise, you use the Arrhenius plot to get your kintics')
disp('by linearly regressing some of your step heating data points.')
disp('If you want to regress your data, using the Magnifier and') 
fprintf('Pan tools in the plot window makes life easier.\n\n')
prompt1 = 'Enter 1 to regress your data OR 2 to manually enter your kinetics:';
check = input(prompt1);

% Check to see if user wants to regress for kinetics or enter them manually
if check == 1   % we are regressing for kinetics
    fprintf('\n\n')
    prompt2 = 'Enter the step numbers to use separated by commas or spaces: ';
    user_steps = input(prompt2, 's')
    regress_steps = str2num(user_steps);

    % Store Arrhenius Data for steps to regress
    regressX = arrDATA(regress_steps,1);
    regressY = arrDATA(regress_steps,2);

    % run a linear Regression Model to get necessary info
    linMod = fitlm(regressX,regressY);

    % linear model results
    Fit_lnDoa2 = linMod.Coefficients.Estimate(1);
    Fit_slope = linMod.Coefficients.Estimate(2);
    Fit_Ea = -1 * Fit_slope* 1.987 * 10;   % Ea in kcal/mol!!
    Fit_R2 = linMod.Rsquared.Ordinary;
    
elseif check == 2  % user is entering Ea and LnDo/a2 manually
    fprintf('\n\n')
    prompt3 = 'ENTER Ea in kcal/mol:';
    Fit_Ea = input(prompt3)
    prompt4 = 'ENTER Ln(Do/a2):';
    Fit_lnDoa2 = input(prompt4)
    Fit_R2 = 0.0;
    
else  % Trap out invalid input!
    fprintf('\n\n')
    disp('Aaaaaaacccckkkk!! You failed me on kinetics input!!!!')
    disp('I am outta here and you aint getting no domains modeled.')
    return
end

%% Build kinetics vector ----------
ref_a = 1  % reference domain size for kinectics from regression --  a value of 1 is arbitrary (but it keeps the math simple)
obs_kinetics = [Fit_Ea Fit_lnDoa2 ref_a];

% calc Delta [Lovera's ln(r/ro)] for observed data with user defined
% kinetics
obs_data(:,7) = delta_calc(arrDATA,Fit_Ea,Fit_lnDoa2);

% build text output
a = num2str(Fit_Ea);
b = num2str(Fit_lnDoa2);
c = num2str(Fit_R2);
aa = strcat({'Ea = '},{a},{' kcal/mol'});
bb = strcat({'Ln(Do/a2) = '},{b});
cc = strcat({'R squared = '},{c});

% plot the data and the steps to regress with the results of regession
subplot(1,2,1,'replace') % replace the previous arrhenius plot to remove the labels

plot(arrDATA(:,1),arrDATA(:,2),'ko','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',3,'LineStyle','-','LineWidth',0.5)
xlim([10 33])
ylim([-35 -5])
title('Observed Arrhenius Data with Kinetics Regression','FontSize',16)
xlabel('10000/K','FontSize',14)
ylabel('ln(D/a2)','FontSize',14)
hold on
if check == 1
    plot(regressX,regressY,'ro','MarkerFaceColor','r','MarkerSize',6)
    lsline  % plot a linear regession line through the regress_steps
    text(24,-7,aa,'FontSize',14)
    text(24,-10,bb,'FontSize',14)
    text(24,-13,cc,'FontSize',14)
end


%% THE BUISNESS -- Create Domains and Hone Size and Gas Content to fit Observed Data


% Master Loop --- Iterates over different numbers of domains to allow the
% model
for iii = 1:length(domain_range)
    ndomains = domain_range(iii);
    disp('Now modeling observed data using');
    disp(ndomains);
    disp('Domains');

    % Domain Loop --- Iterate the full model for a single number of domains 
    % to avoid local minima -- i.e. for a single number of domains allowed, 
    % try to fit the observed data multiple times from scratch.
   
    for ii = 1:n_model_loops

        % Start search with random domain structure
        n_rand_loops = 100;
        BIC_log_rand = zeros(n_rand_loops,1);
        gas_rand = zeros(n_rand_loops,ndomains+1);
        for i = 1:n_rand_loops
            domains_in = create_domain_structure(ndomains);

            % checking gas faction 
            gas_total = sum(domains_in(:,2));
            gas_rand(i,1:ndomains)=domains_in(:,2);
            gas_rand(i,ndomains+1) = gas_total;


            domain_model_out = master_domain(obs_data,obs_kinetics,domains_in);
            BIC_log_rand(i) = domain_model_out(1,2);
            if i == 1
                BIC_best = domain_model_out(1,2); % if it is the first loop, store the BIC
                domain_model_out_best = domain_model_out; % store the domain structure too
            else
                BIC_current = domain_model_out(1,2); % for the rest of the loops compare the BIC of the current loop with the best and store it
                if BIC_current < BIC_best
                    BIC_best = BIC_current;
                    domain_model_out_best = domain_model_out; % store the domain structure too
                end
            end
           % update the Delta plot
           subplot(1,2,2)
           figure1 = make_delta_plot(obs_data, domain_model_out_best);
           drawnow;
           rand_best_domains = domain_model_out_best(2:3,1:ndomains); % store the best rand result
        end
        disp('Model Iteration Number')
        disp(ii)
        disp('')
        disp('Done with Random Domains')
        disp('Now Honing In on Delta fit...')
        disp('Current Best BIC is')
        disp(BIC_best)

        % Now Hone the best attempt from the random search --> use
        % domain_model_out_best

        for z = 1:n_hone_loops

            % first hone in on domain size
            for i = 1:ndomains
                adjuster = 0.01; % how you are modifying the size on the first loop
                for j = 1:adjust_loops
                    domains_in = domain_model_out_best(2:3,1:ndomains);  % grab the current best domain structure 
                    domains_in(1,i) = domains_in(1,i) * adjuster; % change the size of one domain
                    domains_adjusted = transpose(domains_in);  % transpose domains_in for master_domain function

                    % now run the model for the hone iteration and compare the results
                    % to the current best domain structure
                    domain_model_out = master_domain(obs_data,obs_kinetics,domains_adjusted);

                    BIC_current = domain_model_out(1,2);  % get the BIC score of the loop that just ran
                    %BIC_hone_log_size(j,i,z) = BIC_current; % log the hone results

                    % now compare it to the best BIC score
                    if BIC_current < BIC_best
                        BIC_best = BIC_current; % store the current if it is better than the best so far
                        domain_model_out_best = domain_model_out; % store the domain structure too
                        disp(BIC_best)
                    end

                    adjuster = adjuster * size_adjuster; % increment the adjustment

                end
                % update the Delta plot
                subplot(1,2,2,'replace')
                figure1 = make_delta_plot(obs_data, domain_model_out_best);
                drawnow;
            end

            % Now hone in on the Gas Fraction -- obviously doesn't work for
            % 1 domain case.
            
            if ndomains ~= 1

                % Randomly select 2 domains, sum the gas in these domains, then split
                % the gas between the two domains in a systematic way scaling from 0.02
                % to 1.0 of the total gas of the two domains going to one or the other
                % domain.
                gas_hone = zeros(adjust_loops,ndomains+1,ndomains);

                for i = 1:ndomains

                    % Pick two domains to monkey with
                    domains_in = domain_model_out_best(2:3,1:ndomains);  % grab the current best domain structure 
                    test = 0;
                    while test == 0  % check to make sure you grabbed a different domain
                        which_domain = randi(ndomains,1,2);  % randomly select two domains to monkey with
                        if which_domain(1) - which_domain(2) == 0
                            test = 0;
                        elseif which_domain(1) - which_domain(2) ~= 0
                            test = 1;
                        end
                    end

                    % Now do the systematic adjustment of the combined gas of the two
                    % domains

                    adjuster = 0.02; % how you are modifying the gas allocation on the first loop
                    for j = 1:adjust_loops
                        % now sum the gas in the two domains
                        adjust_sum = domains_in(2,which_domain(1)) + domains_in(2,which_domain(2));

                        domains_in(2,which_domain(1)) = adjust_sum * adjuster; % change the gas fraction of one domain
                        domains_in(2,which_domain(2)) = adjust_sum * (1-adjuster);  % change the gas fraction of the other domain

                        %checking gas fraction
                        gas_total = sum(domains_in(2,:));
                        gas_hone(j,1:ndomains,i)=domains_in(2,:);
                        gas_hone(j,ndomains+1,i) = gas_total;

                        %domains_in(2,:)= domains_in(2,:)/sum(domains_in(2,:)); % re-normalize so the sum of the gas fractions = 1, just in case.
                        domains_adjusted = transpose(domains_in);   % transpose domains_in for master_domain function

                        % now run the model for the hone iteration and compare the results
                        % to the current best domain structure
                        domain_model_out = master_domain(obs_data,obs_kinetics,domains_adjusted);

                        BIC_current = domain_model_out(1,2);  % get the BIC score of the loop that just ran
                        %BIC_hone_log_gas(j,i,z) = BIC_current; % log the hone results

                        % now compare it to the best BIC score
                        if BIC_current < BIC_best
                            BIC_best = BIC_current; % store the current if it is better than the best so far
                            domain_model_out_best = domain_model_out; % store the domain structure too
                            disp(BIC_best)
                        end

                        adjuster = adjuster + 0.02; % increment the adjustment

                    end
                end
            end
            % Update the Delta plot
            subplot(1,2,2,'replace')
            figure1 = make_delta_plot(obs_data, domain_model_out_best);
            drawnow;
            
        end

        % Compare the BIC score for each full model iteration, keep the best
        if ii == 1
            best_current_model = domain_model_out_best;  % if it is the first loop, store the result
        else
            BIC_stored = best_current_model(1,2);       % otherwise, compare the current BIC iteration with the stored BIC iteration
            BIC_current_model = domain_model_out_best(1,2);
            if BIC_current_model < BIC_stored
                best_current_model = domain_model_out_best;
            end
        end
    end
    
    % log the results for each full model iteration for a set number of domains 
    BIC_log(iii,1) = ndomains;
    BIC_log(iii,2) = best_current_model(1,2);
    BIC_log(iii,3) = best_current_model(1,1);
    
    % write best results for this set of domains to master results array
    Master_Results_Table(:,:,iii) = best_current_model;
    
    % Compare the BIC score for each full model iteration (i.e. all the attempts for a set number of domains), keep the best
    if iii == 1
        best_domain_model = best_current_model;  % if it is the first loop, store the result
    else
        BIC_stored = best_domain_model(1,2);       % otherwise, compare the current BIC iteration with the stored BIC iteration
        BIC_current_model = best_current_model(1,2);
        if BIC_current_model < BIC_stored
            best_domain_model = best_current_model;
        end
    end
end

domain_model_out_best = best_domain_model;  % put the best model result into the correct variable


%% PREPARE OUTPUTS

% calcuate domain properties for output
domain_model_out_best = calc_domain_properties(obs_kinetics, domain_model_out_best);

% Make arrhenius plot and other output plots with observed data and modeled domains

FigHandle = figure('Position', [100, 100, 500, 450]);
figure2 = make_final_plots(domain_model_out_best, BIC_log);

% Save output figures and model data
csvwrite('modeled_data_OUT.csv', domain_model_out_best);
saveas(FigHandle,'domain_model_plots.pdf');




disp('----- Finished -----')



