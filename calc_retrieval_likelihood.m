% -Compute retrieval likelihood for gaussian model and measurement pdf -


% By Andrew J. Buggee
%%

function retrieval = calc_retrieval_likelihood(inputs,data_inputs,modis,K)

% --- parse through inputs

num_parameters = inputs.num_model_parameters;
measurement_covariance_inv = (inputs.measurement.covaraince)^(-1);
model_covariance_inv = (inputs.model.covaraince)^(-1);
model_mean = inputs.model.mean';

% first, collect the measurements into vector format
% reshape measurement data into Rodgers format for gaussian model and
% measurement pdfs
y = create_measurement_vector(modis,data_inputs);

num_pixels = size(y,2);

retrieval = zeros(num_parameters,num_pixels);



for pp = 1:num_pixels
    
    
    some_K = K{pp}(:,:,1000); % picking a random locaiton in our grid, but this doesn't seem right
    
    retrieval(:,pp) = (some_K' * measurement_covariance_inv * some_K + model_covariance_inv)^(-1) * ...
                    (some_K' * measurement_covariance_inv * y(:,pp) + model_covariance_inv * model_mean);
                
end

    





end