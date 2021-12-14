%% --- Create a droplet profile to solve for ---


% By Andrew J. Buggee
%%

function inputs = create_droplet_profile(inputs,data_inputs)

profile_type = inputs.model.profile.type;

r_top = inputs.model.profile.r_top;
r_bottom = inputs.model.profile.r_bottom;

% we create a droplet profile for each pixel that we wish to investigate
tau_c = data_inputs.truthTable.modisT17;

% define how many points along our drolet profile we wish to plot
num_query_points = 100;
r_tau = zeros(length(tau_c),num_query_points);
r_z = zeros(length(tau_c),num_query_points);

if strcmp(profile_type,'subadiabatic_aloft')
    
    % if the profile chosen is subadiabatic aloft, then 0<x<1
    x = 1/2;
    
    a0 = r_top^(2*x + 3/x);
    a1 = r_top^(2*x + 3/x) - r_bottom^(2*x + 3/x);
    
    for ii = 1:length(tau_c)
        T = linspace(0,tau_c(ii),num_query_points);
        r_tau(ii,:) = (a0 - a1*(T./tau_c(ii))).^(1/(2*x + 3/x));
    end
    
elseif strcmp(profile_type,'adiabatic')
    
        % if the profile chosen is adiabatic, then x=1
    x = 1;
    
    % boundary conditions for r as a function of tau
    a0 = r_top^(2*x + 3/x);
    a1 = r_top^(2*x + 3/x) - r_bottom^(2*x + 3/x);
    
    % boundary conditions for r as a function of tau
    b0 = r_bottom^3;
    b1 = r_top^3 - r_bottom^3;
    
    for ii = 1:length(tau_c)
        T = linspace(0,tau_c(ii),num_query_points);
        r_tau(ii,:) = (a0 - a1*(T./tau_c(ii))).^(1/(2*x + 3/x));
        r_z(ii,:) = (b0 + b1 * (z-z0)./h).^(1/3); % droplet profile in geometric coordniate system for an adiabatic cloud
    end
    
    
elseif strcmp(profile_type,'subadiabatic_midlevel')
    
    % if the profile chosen is subadiabatic aloft, then is either 3 or -3
    x = 3;
    
    a0 = r_top^(2*x + 3/x);
    a1 = r_top^(2*x + 3/x) - r_bottom^(2*x + 3/x);
    
    for ii = 1:length(tau_c)
        T = linspace(0,tau_c(ii),num_query_points);
        r_tau(ii,:) = (a0 - a1*(T./tau_c(ii))).^(1/(2*x + 3/x));
    end    
    
    
else
    
    error('I dont recognize the droplet profile you want!')
    
end


% gather calculations

inputs.model.profile.x = x;
inputs.model.profile.a0 = a0;
inputs.model.profile.a1 = a1;
inputs.model.profile.r_tau = r_tau;
inputs.model.profile.r_z = r_z;


end