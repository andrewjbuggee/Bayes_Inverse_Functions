%% --- Create a droplet profile to solve for ---


% By Andrew J. Buggee
%%

function inputs = create_droplet_profile(inputs,data_inputs)


% introduce constants
rho_l = 1e6;    % g/m^3 - density of liquid water


profile_type = inputs.model.profile.type;

r_top = inputs.model.profile.r_top;                                         % microns
r_bottom = inputs.model.profile.r_bottom;                                   % microns

% we create a droplet profile for each pixel that we wish to investigate
% lets pull the optical depth values that MODIS calculated for this data
% set. The data set we use also determines what wavelength is used when we
% reference the optical depth
tau_c = data_inputs.truthTable.modisT17;
wl_tau = modisBandsCenter(1);                       % nm - Make sure the number used is the first digit in the above MODIS data set!



% define how many points along our drolet profile we wish to plot
num_query_points = 10;
r_tau = zeros(length(tau_c),num_query_points);
r_z = zeros(length(tau_c),num_query_points);
Nc = zeros(length(tau_c),1);
lwc_z = zeros(length(tau_c),num_query_points); 

% ---- for now.....

% no matter the profile, we are create a droplet profile that extends 1 km,
% from the base of the cloud at 2km above the ground to the top of the
% cloud at 3 km above the ground
% we only define the values at the base of each layer
h = 1000; % meters - cloud thickness
z0 = 2000; % meters - cloud base location
z = linspace(z0,(z0+h),num_query_points); % meters - the locations in geometric space that define the bottom of each layer within the cloud


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
    
    % boundary conditions for r as a function of z
    
    b0 = r_bottom^3;
    b1 = r_top^3 - r_bottom^3;
    
    % Boundary conditions for LWC as a function of z
    
    % Lets step through each pixel according to data_inputs, and compute
    % the profile for a given r_bottom, r_top, and tau_c
    
    for ii = 1:length(tau_c)
        T = linspace(0,tau_c(ii),num_query_points);
        r_tau(ii,:) = (a0 - a1*(T./tau_c(ii))).^(1/(2*x + 3/x));
        r_z(ii,:) = (b0 + b1 * (z-z0)./h).^(x/3);                      % droplet profile in geometric coordniate system for an adiabatic cloud
        
        % ----------------------- ASSUMPTION ----------------------
        % We assume the number concentration is constant with heighy
        % ---------------------------------------------------------
        
        % First lets compute the extinction efficient
        % Lets assume a monodispersed distribution
        yq = interp_mie_computed_tables([linspace(wl_tau, wl_tau,num_query_points)',r_z(ii,:)'],'mono');
        
        Qext = yq(:,5);                                                                          % Extinction efficiency
        Nc(ii) = tau_c(ii)/(pi*trapz((z-z0),Qext'.*(r_z(ii,:)*1e-6).^2));                                     % m^(-3) - number concentration
        
        lwc_z(ii,:) = 4/3 * pi * rho_l * r_z(ii,:).^3 * Nc(ii);
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
inputs.model.profile.Nc = Nc;
inputs.model.profile.lwc_z = lwc_z;
inputs.model.profile.lambda_tau = wl_tau;


end