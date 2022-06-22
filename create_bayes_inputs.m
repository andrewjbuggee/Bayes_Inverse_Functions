function bayes_inputs = create_bayes_inputs()

% Define the number of pixels to estimate a profile for
bayes_inputs.numPixels2Calculate = 10;


% define the number of iterations for the gauss-newton solver
bayes_inputs.GN_iterations = 20;

% Define the convergence limit. Convergence is defined using the residual,
% which is the true measurement subtracted from the estiamted measurement.
% We take the RMS of the residual using all spectral channels. This is how
% we define the convergence limit. If the residual is the difference
% between the true measurement and the estimated measurement, and the true
% measurement has an uncertainty of 10%, then our estimate measurement
% should be within this uncertainty. For MODIS, the measuremen uncertainty
% for the reflectance is between 3 and 7%. So lets meet in the middle and
% say 5 %
bayes_inputs.convergence_limit = 0.05;

% define the type of model prior pdf
bayes_inputs.model.prior = 'gaussian';


% define the number of model parameters to solve for
bayes_inputs.num_model_parameters = 3;

% Define the number of spectral channels to use in the gauss-newton solver
bayes_inputs.numBands2use = 7;

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
%   (1) 'subadiabatic_aloft'
%   (2) 'subadiabatic_midlevel' - radius grows linearly with height r ~ z
%   (3) 'adiabatic' - liquid water content grows linearly with height LWZ ~ z
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