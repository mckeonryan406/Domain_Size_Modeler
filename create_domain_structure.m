function[domains_in] = create_domain_structure(ndomains)

% Create Random Domain Structure to start search

% Input --> number of domains to populate with size and gas fraction

% Output --> a random domain structure of ndomains as "domains_in"

% generate array of random numbers
domain_raw = rand(ndomains,3);

% use column 1 for gas fraction
total_gas = sum(domain_raw(:,1));
domain_raw(:,1) = domain_raw(:,1)/total_gas; % set total gas = to 1

% use column 2 for domain size, but adjust it up according to the
% value in column 3...  remember that domains can vary in size by orders of
% magnitude, so hard wire in a lot of variation potential here.

% cludgey if/then for 1 and 2 domain case is necessary due to the for loop
% for 3 and above domain cases... not sure why, but this now works.
% if ndomains == 1
%     if domain_raw(3) < 0.2
%         domain_raw(2) = domain_raw(2) * 1;
%     elseif domain_raw(3) < 0.4
%         domain_raw(2) = domain_raw(2) * 10;
%     elseif domain_raw(3) < 0.6
%         domain_raw(2) = domain_raw(2) * 100;
%     elseif domain_raw(3) < 0.8
%         domain_raw(2) = domain_raw(2) * 1000;
%     else 
%         domain_raw(2) = domain_raw(2) * 10000;
%     end
% elseif ndomains == 2
%     if domain_raw(1,3) < 0.2
%         domain_raw(1,2) = domain_raw(1,2) * 1;
%     elseif domain_raw(1,3) < 0.4
%         domain_raw(1,2) = domain_raw(1,2) * 10;
%     elseif domain_raw(1,3) < 0.6
%         domain_raw(1,2) = domain_raw(1,2) * 100;
%     elseif domain_raw(1,3) < 0.8
%         domain_raw(1,2) = domain_raw(1,2) * 1000;
%     else 
%         domain_raw(1,2) = domain_raw(1,2) * 10000;
%     end
%     if domain_raw(2,3) < 0.2
%         domain_raw(2,2) = domain_raw(2,2) * 1;
%     elseif domain_raw(2,3) < 0.4
%         domain_raw(2,2) = domain_raw(2,2) * 10;
%     elseif domain_raw(2,3) < 0.6
%         domain_raw(2,2) = domain_raw(2,2) * 100;
%     elseif domain_raw(2,3) < 0.8
%         domain_raw(2,2) = domain_raw(2,2) * 1000;
%     else 
%         domain_raw(1,2) = domain_raw(1,2) * 10000;
%     end
% else
%     for i = 1:length(domain_raw)
%         if domain_raw(i,3) < 0.2
%             domain_raw(i,2) = domain_raw(i,2) * 1;
%         elseif domain_raw(i,3) < 0.4
%             domain_raw(i,2) = domain_raw(i,2) * 10;
%         elseif domain_raw(i,3) < 0.6
%             domain_raw(i,2) = domain_raw(i,2) * 100;
%         elseif domain_raw(i,3) < 0.8
%             domain_raw(i,2) = domain_raw(i,2) * 1000;
%         else 
%             domain_raw(i,2) = domain_raw(i,2) * 10000;
%         end
%     end
% end

for i = 1:ndomains
    if domain_raw(i,3) < 0.2
        domain_raw(i,2) = domain_raw(i,2) * 1;
    elseif domain_raw(i,3) < 0.4
        domain_raw(i,2) = domain_raw(i,2) * 10;
    elseif domain_raw(i,3) < 0.6
        domain_raw(i,2) = domain_raw(i,2) * 100;
    elseif domain_raw(i,3) < 0.8
        domain_raw(i,2) = domain_raw(i,2) * 1000;
    else 
        domain_raw(i,2) = domain_raw(i,2) * 10000;
    end
end

%build output
domains_in = [];
domains_in(:,1) = domain_raw(:,2);  % size
domains_in(:,2) = domain_raw(:,1);  % gas fraction



