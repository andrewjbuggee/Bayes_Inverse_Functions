function bayes_inputs = create_bayes_inputs(modisInputs)

% Define the number of pixels to estimate a profile for
bayes_inputs.numPixels2Calculate = 4;


% define the number of iterations for the gauss-newton solver
bayes_inputs.GN_iterations = 5;

% Define the convergence limit. Convergence is defined using the residual,
% which is the true measurement subtracted from the estiamted measurement.
% We take the RMS of the residual using all spectral channels. This is how
% we define the convergence limit. If the residual is the difference
% between the true measurement and the estimated measurement, and the true
% measurement has an uncertainty of 10%, then our estimate measurement
% should be within this uncertainty. For MODIS, the measuremen uncertainty
% for the reflectance is between 3 and 7%. So lets meet in the middle and
% say 5 %
bayes_inputs.convergence_limit = 0.005;

% define the type of model prior pdf
bayes_inputs.model.prior = 'gaussian';


% define the number of model parameters to solve for
bayes_inputs.num_model_parameters = 3;

% Define the spectral channels to use in the gauss-newton solver
% The data from 11-11-2008 at 18:50 measured erroneous values in the 1.6
% micron channel. If using this data, lets ignore this measurement

if strcmp(modisInputs.modisDataFolder(96:end), '/2008_11_11_1850/')==true
    
    bayes_inputs.bands2use = [1:5,7];
    
else 

    bayes_inputs.bands2use = 1:7;  % number of spectral bands to use

end



% -------------------------------------------
% --- Stuff for the Model Parameter Prior ---
% -------------------------------------------

% Using the King and Vaughn (2012) method, we retireve 3 parameters
%   (1) r_top = effective droplet size at the cloud top
%   (2) r_bottom = effective droplet size at the cloud bottom
%   (3) tau_c = cloud optical depth
% a good starting place is to assume the droplet size at cloud top and
% bottom are the same value



    

bayes_inputs.model.param_names = {'Effective Radius at Top of Cloud', 'Effective Radius at Bottom of Cloud',...
    'Cloud Optical Depth'};


% ---------------------------------------
% --- Stuff for the Measurement Prior ---
% ---------------------------------------


bayes_inputs.measurement.prior = 'gaussian';
% covaraince_type can be:
%   (1) 'independent - thus all off diagonal elements are 0
%   (2) 'computed' - uses measured data to compute covaraince
bayes_inputs.measurement.covariance_type = 'independent';

% -----------------------------------------------
% --- Stuff for the Assumed Vertical Profile ---
% -----------------------------------------------

% we have to assume a vertical profile of droplet size with cloud optical
% depth exsists. And we only retrieve the droplet size at the top and
% bottom. This is the method of King and Vaughn (2012)

% the options for the vertical droplet profile are:
%       (a) 'adiabatic' - this assumption forces the liquid water content to
%       be proportionl to z, the altitude.
%       (b) 'subadiabatic_aloft' - this assumption assumes there is
%       increasing entrainment and drying towards the cloud top.
%       (c) 'linear_with_z' - this constraint forces the effective droplet profile
%       to behave linearly with z (re(z)~z). Physically we are forcing subadiabtatic
%       behavior at mid-levels.
%       (d) 'linear_with_tau' - this constraint forces the effective
%       droplet radius to have linearly with optical depth (re(z)~tau).
%       Physically, this too forces subadiabatic behavior at mid-levels.
% x is determined by the choice of droplet profile within the function
% create_droplet_profile.m

bayes_inputs.model.profile.type = 'adiabatic';
bayes_inputs.model.profile.r_top = 10; % microns - value for our model
bayes_inputs.model.profile.r_bottom = 5; % microns - value for our model



% -----------------------------------------------
% ------------- Folder Locations  ---------------
% -----------------------------------------------
bayes_inputs.save_calcs_fileName = ['uvspec_CALCS_4Bayes_',date,'.mat'];




% -----------------------------------------------
% --------------- Define flags  -----------------
% -----------------------------------------------



end