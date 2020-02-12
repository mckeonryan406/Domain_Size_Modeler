function [ figure1 ] = make_delta_plot(obs_data, domain_model_out_best)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%make plot --> delta observed vs. delta modeled
cla
% build stairs for OBSERVED data
data_in = obs_data;
Fcum = data_in(:,4);
r_ro = data_in(:,7);
for a = 1:length(Fcum)-1;
        x(2*a-1) = Fcum(a);
        x(2*a) = Fcum(a+1);
        y(2*a-1) = r_ro(a);
        y(2*a) = r_ro(a);
end;
textY = min(y)+1;

figure1 = plot(x,y,'k','LineWidth',2);
hold on

% get Misfit and BIC scores and nDomains
misfitNOW = num2str(domain_model_out_best(1,1));
BICnow = num2str(domain_model_out_best(1,2));
ndomainsNOW = num2str(domain_model_out_best(1,3));

% build stairs for MODELED data
Fcum = domain_model_out_best(4:length(domain_model_out_best),4);
r_ro = domain_model_out_best(4:length(domain_model_out_best),7);
for a = 1:length(Fcum)-1;
        x(2*a-1) = Fcum(a);
        x(2*a) = Fcum(a+1);
        y(2*a-1) = r_ro(a);
        y(2*a) = r_ro(a);
end;
plot(x,y,'r','LineWidth',2);


xname = 'Sum 3He Released';
yname = 'Delta';
title('Modeled Delta (red) vs. Observed Delta (black)','FontSize',16)
xlabel(xname,'fontsize',14);  %axis labels
ylabel(yname,'fontsize',14);
xlim([0 1]);
text(0.7,textY,strcat({'Misfit = '},{misfitNOW}),'FontSize',14);
text(0.7,textY+1,strcat({'BIC = '},{BICnow}),'FontSize',14);
text(0.7,textY+2,strcat({'Domains = '},{ndomainsNOW}),'FontSize',14);
box on  


end

