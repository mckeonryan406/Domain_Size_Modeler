function[Fsum] = weight_domains(Fdata, domains)

% This function takes in F release data for multiple domains and using
% defined gas fractions for each domain it weights the step release and
% returns the sum of F for all domains for each step.

% Inputs:  Fdata = array with Fcum for each domain
%                  col 1 = domain 1 F
%                  col 2 = domain 2 F, etc...
%          domains = array with relative domain size in col 1 and gas
%                    fraction in col 2

bb = size(Fdata);
nsteps = bb(1);
ndomains = bb(2);

weight_F = [];

for i = 1:ndomains
    Fcum = Fdata(:,i);
    weight_F(:,i) = Fcum * domains(i,2);
end
Fsum = [];
for i = 1:nsteps
    Fsum(i) = sum(weight_F(i,:));
end
