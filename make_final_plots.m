function [figure2] = make_final_plots(domain_model_out_best, BIC_log)

%% Delta Plot -- BIC Best Result
%make plot --> delta observed vs. delta modeled
cla
data_in=[]; % clear data_in variable before continuing

% build stairs for OBSERVED data
data_in = domain_model_out_best;
Fcum = data_in(6:end,10);
r_ro = data_in(6:end,8);
for a = 1:length(Fcum)-1;
        x(2*a-1) = Fcum(a);
        x(2*a) = Fcum(a+1);
        y(2*a-1) = r_ro(a);
        y(2*a) = r_ro(a);
end;
subplot(2,2,1)
figure1 = plot(x,y,'k','LineWidth',2);
hold on

textY = min(y)+1;

% build stairs for MODELED data

Fcum = domain_model_out_best(6:length(domain_model_out_best),4);
r_ro = domain_model_out_best(6:length(domain_model_out_best),7);
for a = 1:length(Fcum)-1;
        x(2*a-1) = Fcum(a);
        x(2*a) = Fcum(a+1);
        y(2*a-1) = r_ro(a);
        y(2*a) = r_ro(a);
end;

% get Misfit and BIC scores and nDomains
misfitNOW = num2str(domain_model_out_best(1,1));
BICnow = num2str(domain_model_out_best(1,2));
ndomainsNOW = num2str(domain_model_out_best(1,3));

plot(x,y,'r','LineWidth',2);
xname = 'Sum 3He Released';
yname = 'Delta';
title('Modeled Delta (red) vs. Observed Delta (black) -- BIC Best')
xlabel(xname,'fontsize',12);  %axis labels
ylabel(yname,'fontsize',12);
xlim([0 1]);
text(0.7,textY,strcat({'Misfit = '},{misfitNOW}),'FontSize',12);
text(0.7,textY+1,strcat({'BIC = '},{BICnow}),'FontSize',12);
text(0.7,textY+2,strcat({'Domains = '},{ndomainsNOW}),'FontSize',12);
box on  


%% Arrhenius Plot -- Observed Data with BIC Best result as lines and approximated domain size
% calc slope
R = 1.987;
%Ea = domain_model_out_best(1,4);
Ea_Farley = 40.6; % ~170 kJ/mol converted to kcal/mol
Ea = Ea_Farley;
slope = -Ea/R/10; % divided by 10 due to kcal/mol units of Ea

% generate ln(Do/a2) for reference diffusion domain sizes -->  using single
% crystal data of Ken's and your GSA poster (2015) Do estimate.
%Do = 0.1; % cm/s^2 from GSA Poster 2015
Do_Farley = exp(-0.66); % cm/s^2 from Farley 2018 GCA
Do = Do_Farley;

dd_size = [0.00001, 0.001, 0.1]; % for 100 nm, 10 um, and 1 mm relative to cm's
lnDoa2_ref=[];
for i = 1:length(dd_size)
    lnDoa2_ref(i) = log(Do/dd_size(i)^2);
end

%plot lines for reference domain sizes
X=[0,40];
Y=[0,0];

subplot(2,2,3)
for i = 1:length(dd_size)
    Y(1) = lnDoa2_ref(i);
    Y(2) = X(2)*slope + lnDoa2_ref(i);
    figure2 = plot(X,Y,'Color', [.7,.7,.7], 'LineWidth',1, 'LineStyle','--');
    hold on
end

% calc ln(Do/a2) for each domain
ndomains = domain_model_out_best(1,3);
lnDoa2 = [];
for i = 1:ndomains
    logDoa2 = domain_model_out_best(4,i);
    Doa2 = 10^logDoa2;
    lnDoa2(i) = log(Doa2);  % remember in Matlab log() calculate the Natural Log!
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
obs_data(:,1) = domain_model_out_best(6:end,5);
obs_data(:,2) = domain_model_out_best(6:end,9);

plot(obs_data(:,1),obs_data(:,2),'ko','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',3,'LineStyle','-','LineWidth',0.5)
xlim([5 35]);
ylim([-35 -5]);
xname = '10000/K';
yname = 'ln(D/a2) ln(S-1)';
title('Modeled (red) vs. Observed (black) with Reference (gray) -- BIC Best')
xlabel(xname,'fontsize',12);  %axis labels
ylabel(yname,'fontsize',12); 
box on  %draws a border around the plot

%% plot best BIC score vs. number of domains

subplot(2,2,2)
plot(BIC_log(:,1),BIC_log(:,2),'ko','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',5,'LineStyle','-','LineWidth',0.5)
xname = 'Number of Domains';
yname = 'BIC Score';
title('Goodness of Fit vs. Allowed Complexity')
xlabel(xname,'fontsize',12);  %axis labels
ylabel(yname,'fontsize',12); 

%% plot best misfit score vs. number of domains

subplot(2,2,4)
plot(BIC_log(:,1),BIC_log(:,3),'ko','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',5,'LineStyle','-','LineWidth',0.5)
xname = 'Number of Domains';
yname = 'Misfit Score';
title('Just Goodness of Fit')
xlabel(xname,'fontsize',12);  %axis labels
ylabel(yname,'fontsize',12); 



