
function [GN_output, GN_inputs] = calc_retrieval_gauss_newton_4modis(GN_inputs,modis,modisInputs, pixels2use)

% ----- unpack inputs -----

model_apriori = GN_inputs.model.apriori'; % a priori expected values for the model parameters
model_cov = GN_inputs.model.covariance; % model parameter covariance matrix
measurement_cov = GN_inputs.measurement.covariance; % measurement covaraince matrix
initialGuess = GN_inputs.model.initialGuess';      % Initial guess to start the Gauss-Newton iteration

% Retrieve the convergence limit
convergence_limit = GN_inputs.convergence_limit;

% Create the measurement vectors for each pixel!
% Each column is associated with a specific pixel, and each row represents
% the reflectance measurement at a specific modis band
measurements = create_measurement_vector(modis,GN_inputs, pixels2use); % each column represents one pixel

% % The data from 11-11-2008 at 18:50 measured erroneous values in the 1.6
% % micron channel. If using this data, lets ignore this measurement
% if strcmp(modisInputs.modisDataFolder(96:end), '/2008_11_11_1850/')==true
% 
%     % get rid of the 6th MODIS channel
%     measurements(6,:) = [];
% 
% else
% 
% end


% -----------------------------------------------------------------------
% --------------- Define the spectral response function -----------------
% -----------------------------------------------------------------------

% using the input 'bands2run', this defines the spectral bands that will be
% written into INP files.

% check to see if the MODIS instrument is aboard Terra or Aqua
if strcmp(modisInputs.L1B_filename(1:3), 'MOD')==true
    % Then read in the spectral response functions for the terra instrument
    GN_inputs.spec_response = modis_terra_specResponse_func(GN_inputs.bands2use, GN_inputs.RT.sourceFile_resolution);

elseif strcmp(modisInputs.L1B_filename(1:3), 'MYD')==true
    % Then read in the spectral response functions for the Aqua instrument
    GN_inputs.spec_response = modis_aqua_specResponse_func(GN_inputs.bands2use, GN_inputs.RT.sourceFile_resolution);
end


% ------------------------------------------------------------------------





num_iterations = GN_inputs.GN_iterations; % number of iterations to preform
num_parameters = GN_inputs.num_model_parameters; % number of parameters to solve for

% The number of pixels to solve for using the Gauss-Newton method comes
% from the Gauss-Newton input structure
num_pixels = GN_inputs.numPixels2Calculate;

% ----- define number of spectral bands to use -----

num_bands = length(GN_inputs.bands2use);



% --- Create iterative Gauss-Newton Solver ----

retrieval = zeros(num_parameters,num_iterations+1,num_pixels); % we include the starting point, which is outside the number of iterations
residual = zeros(num_bands,num_iterations,num_pixels); % we include the starting point, which is outside the number of iterations
rms_residual = zeros(num_iterations, num_pixels);      % RMS of the residual across all bands
diff_guess_prior = zeros(num_parameters,num_iterations,num_pixels); % we include the starting point, which is outside the number of iterations
jacobian_diff_guess_prior = zeros(num_bands,num_iterations,num_pixels);      % this is an expression that multiplies the Jacobian with the difference between the current iteration and the a priori
posterior_cov = zeros(num_parameters,num_parameters,num_pixels); % my posterior covariance matrix

for pp = 1:num_pixels
    
    % --- define an initial guess ----
    
    % we've set the effective radius at the top and bottom to be the same
    % value. This is our initial guess. By setting this to be the a priori
    % guess we can compute the inforamtion gain between the bi-spectral
    % approach and the hyperspectal approach
    
    % define the initial guess
    % Here we define it to be the same re retireved for r_top and r_bottom
    retrieval(:,1,pp) = initialGuess(:,pp);
    
    % -----------------------------------------------
    % ----- USING King and Vaughn (2012) Method -----
    % -----------------------------------------------
    % we set the a priori guess to be our initial guess. That way, any
    % change in our posterior tells us the information gained over the
    % bi-spectral method
    
    
    % ---- define which pixel to use -----
    pixel_row = pixels2use.res1km.row(pp);
    pixel_col = pixels2use.res1km.col(pp);
    
    
    
    for ii = 1:num_iterations
        
        disp(['(iteration,pixel): (',num2str(ii),',',num2str(pp),')']) % this tells the user where we are in the process
        
        % at each iteration I need to compute the forward model at my current
        % state vector estimate
        
        
        current_guess = retrieval(:,ii,pp);
        
        % we compute the forward model at our previous estimate of the state vector
        % Therefore, we ask, 'what is the reflectance of a cloud with our
        % current state vector guess?'
        measurement_estimate = compute_forward_model_4modis(modis,current_guess,GN_inputs,pixel_row,pixel_col,modisInputs)';
        
        Jacobian = compute_jacobian_4modis(modis,current_guess,measurement_estimate,GN_inputs,modisInputs, pixel_row,pixel_col, GN_inputs.measurement.variance(:,pp));
        
        
        
        
        residual(:,ii,pp) = measurements(:,pp) - measurement_estimate;
        diff_guess_prior(:,ii,pp) = current_guess - model_apriori(:,pp);
        jacobian_diff_guess_prior(:,ii,pp) = Jacobian*diff_guess_prior(:,ii,pp);
        
        % -----------------------------------------------------------------
        % -----------------------------------------------------------------
        % new_guess using the previous iteration
        %new_guess = current_guess + (model_cov(:,:,pp)^(-1) + jacobian' * measurement_cov^(-1) *jacobian)^(-1) * (jacobian' *  measurement_cov(:,:,pp)^(-1) * residual(:,ii,pp) - model_cov(:,:,pp)^(-1) *diff_guess_prior(:,ii,pp));
        
        % new_guess using the model prior mean value
        new_guess = model_apriori(:,pp) + model_cov(:,:,pp) * Jacobian' * (Jacobian * model_cov(:,:,pp) * Jacobian' + measurement_cov(:,:,pp))^(-1) * (residual(:,ii,pp) + jacobian_diff_guess_prior(:,ii,pp));
        % -----------------------------------------------------------------
        % -----------------------------------------------------------------


        % when using the Hu and Stamnes parameterization, re must be larger
        % than 2.5 and smaller than 60 microns. If the new guess is out of
        % these bounds, fix it!
        if new_guess(1)>60
            disp([newline,'r_top = ',num2str(new_guess(1)),'. Set to 59.5 \mum'])
            new_guess(1) = 59.5; % microns - this may just bump back up to 60, but maybe not. The model prior should help with that
        elseif new_guess(1)<2.5
            disp([newline,'r_top = ',num2str(new_guess(1)),'. Set to 2.6 \mum'])
            new_guess(1) = 2.6; % microns
        end
        
        if new_guess(2)>60
            disp([newline,'r_bottom = ',num2str(new_guess(2)),'. Set to 59.5 \mum'])
            new_guess(2) = 59.5; % microns - this may just bump back up to 60, but maybe not. The model prior should help with that
        elseif new_guess(2)<2.5
            disp([newline,'r_bottom = ',num2str(new_guess(2)),'. Set to 2.6 \mum'])
            new_guess(2) = 2.6; % microns
        end
        
        
        % store the latest guess
        retrieval(:,ii+1,pp) = new_guess;

        % If the residual is below a certain threshold as defined in the
        % GN_inputs strucute, break the for loop. We've converged
        rms_residual(ii,pp) = sqrt(sum(residual(:,ii,pp).^2)/size(residual,1));

        if rms_residual(ii,pp)<convergence_limit
            disp([newline, 'Convergence reached in ', num2str(ii),' iterations.', newline,...
                'RMS = ', num2str(rms_residual(ii,pp))])

            % Clear the rest of the zeros that are place holders for later
            % iterations
            retrieval(:,ii+2:end,pp) = nan;
            rms_residual(ii+1:end,pp) = nan;
            residual(:,ii+1:end,pp) = nan;
            diff_guess_prior(:,ii+1:end,pp) = nan;
            jacobian_diff_guess_prior(:,ii+1:end,pp) = nan;

            break
        end

        
    end
    
    % once convergence has occured, we can compute the posterior covariance
    % matrix
    
    posterior_cov(:,:,pp) = (Jacobian' * measurement_cov(:,:,pp)^(-1) * Jacobian + model_cov(:,:,pp)^(-1))^(-1);

    
    % Compute the retireved Liquid water path with the final profile
    
    z = linspace(GN_inputs.model.cloudTop_height - GN_inputs.model.cloudDepth, GN_inputs.model.cloudTop_height, GN_inputs.model.cloud_layers);        % km - altitude above ground vector

    GN_outputs.LWP(pp) = trapz(z, create_droplet_profile2([GN_outputs.retrieval(1,end), GN_outputs.retrieval(2,end)],...
            GN_outputs.tau_vector, 'altitude', GN_inputs.model.profile.type));

    
end

GN_output.retrieval = retrieval;
GN_output.residual = residual;
GN_output.rms_residual = rms_residual;
GN_output.diff_guess_prior = diff_guess_prior;
GN_output.jacobian_diff_guess_prior = jacobian_diff_guess_prior;
GN_output.posterior_cov = posterior_cov;



end








