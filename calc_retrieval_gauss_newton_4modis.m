
function [output] = calc_retrieval_gauss_newton_4modis(bayes_inputs,modis,data_inputs)

% ----- unpack inputs -----

model_mean = bayes_inputs.model.mean'; % a priori guess
model_cov = bayes_inputs.model.covariance; % model parameter covariance matrix
measurement_cov = bayes_inputs.measurement.covariance; % measurement covaraince matrix

saveCalculations_fileName = bayes_inputs.save_calcs_fileName; % where to save the computed data

%model_mean = bayes_inputs.model.mean; % a priori expected value for the model parameters

measurements = create_measurement_vector(modis,data_inputs); % each column represents one pixel
biSpectral_estimate = [data_inputs.truthTable.estR17,data_inputs.truthTable.estT17];

num_iterations = bayes_inputs.GN_iterations; % number of iterations to preform
num_parameters = bayes_inputs.num_model_parameters; % number of parameters to solve for
num_pixels = data_inputs.inputs.pixels.num_2calculate;

% ----- unpack look-up table information -----

num_bands = bayes_inputs.measurement.num; % number of spectral bands ran in the RT model


% --- Create iterative Gauss-Newton Solver ----


retrieval = zeros(num_parameters,num_iterations+1,num_pixels); % we include the starting point, which is outside the number of iterations
residual = zeros(num_bands,num_iterations,num_pixels); % we include the starting point, which is outside the number of iterations
diff_guess_prior = zeros(num_parameters,num_iterations,num_pixels); % we include the starting point, which is outside the number of iterations
posterior_cov = zeros(num_parameters,num_parameters,num_pixels); % my posterior covariance matrix

for pp = 1:num_pixels
    
    % --- define an initial guess ----
    
    % we've set the effective radius at the top and bottom to be the same
    % value. This is our initial guess. By setting this to be the a priori
    % guess we can compute the inforamtion gain between the bi-spectral
    % approach and the hyperspectal approach
    initial_guess = [biSpectral_estimate(pp,1),biSpectral_estimate(pp,1),biSpectral_estimate(pp,2)];
    retrieval(:,1,pp) = initial_guess;
    
    % ----- USING King and Vaughn (2012) Method -----
    % -----------------------------------------------
    % we set the a priori guess to be our initial guess. That way, any
    % change in our posterior tells us the information gained over the
    % bi-spectral method
    
    model_mean = initial_guess';
    
    % King and Vaughn define the uncertainty in tau_c as 10% the value
    % found by the bi-spectral method
    model_cov(3,3) = 0.1 * initial_guess(3);
    
    % ---- define pixel geometry -----
    pixel_row = data_inputs.pixels2use.res1km.row(pp);
    pixel_col = data_inputs.pixels2use.res1km.col(pp);
    
    
    
    for ii = 1:num_iterations
        
        disp(['(iteration,pixel): (',num2str(ii),',',num2str(pp),')']) % this tells the user where we are in the process
        
        % at each iteration I need to compute the forward model at my current
        % state vector estimate
        
        
        current_guess = retrieval(:,ii,pp);
        
        % we compute the forward model at our previous estimate to the
        % state vector
        measurement_estimate = compute_forward_model_4modis(modis,current_guess,bayes_inputs,pixel_row,pixel_col,data_inputs.inputs.INP_folderName,saveCalculations_fileName)';
        
        jacobian = compute_jacobian_4modis(modis,current_guess,measurement_estimate,bayes_inputs,pixel_row,pixel_col,data_inputs.inputs.INP_folderName,saveCalculations_fileName);
        
        
        
        
        residual(:,ii,pp) = measurements(:,pp) - measurement_estimate;
        diff_guess_prior(:,ii,pp) = model_mean - current_guess;
        
        new_guess = model_mean + model_cov * jacobian' * (jacobian * model_cov * jacobian' + measurement_cov)^(-1) * (residual(:,ii,pp) + jacobian * diff_guess_prior(:,ii,pp));
        
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
        
    end
    
    % once convergence has occured, we can compute the posterior covariance
    % matrix
    
    posterior_cov(:,:,pp) = (jacobian' * measurement_cov^(-1) * jacobian + model_cov^(-1))^(-1);
    
end

output.retrieval = retrieval;
output.residual = residual;
output.diff_guess_prior = diff_guess_prior;
output.posterior_cov = posterior_cov;



end








