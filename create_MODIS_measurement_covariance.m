%% create measurement prior


% By Andrew J. Buggee
%%

function [GN_inputs] = create_MODIS_measurement_covariance(GN_inputs,modis,modisInputs)


covariance_type = GN_inputs.measurement.covariance_type;



% ---**--- Important Quantity ---**---
GN_inputs.measurement.uncertainty = 0.02; % percentage of measurement uncertainty for reflectance according To "VALIDATION OF MODIS-DERIVED TOP-OF-ATMOSPHERE SPECTRAL RADIANCES BY MEANS OF VICARIOUS CALIBRATION"

% Define the number of spectral channels
n = GN_inputs.numBands2use;


% --------------------------------------------------------
% Create the covariance matrix by taking the cross
% correlation between spectral channels with the modis data
% ---------------------------------------------------------

if strcmp(covariance_type,'computed') == true
    data = cat(3,modis.EV.m250.reflectance,modis.EV.m500.reflectance);
    data = data(:,:,GN_inputs.spectral_bins);
    for bb = 1:length(GN_inputs.spectral_bins)
        for ii = 1:length(modisInputs.pixels2use.res500m.row)
            data2run(ii,bb) = data(modisInputs.pixels2use.res500m.row(ii),modisInputs.pixels2use.res500m.col(ii),bb);
        end
    end
    
    GN_inputs.measurement.covariance = cov(data2run);
    
elseif strcmp(covariance_type,'independent') == true
    % create the covaraince matrix of the model parameters
    % if the covariance matrix is diagonal, then we are assuming each
    % measurement (spectral channel) is independent of one another
    
    for ii = 1:n
        % if each uncertainty represents the standard deviation, the
        % variance is the square of each value
        channel_variance = mean(modis.EV1km.reflectanceUncert(:,:,ii),'all')^2;
        
        % record this as the channel variance
        GN_inputs.measurement.variance(ii) = channel_variance;
        
    end
    
    % Create a diagonal matrix where each entry is the variance of that
    % spectral channel for reflectance measurements
    GN_inputs.measurement.covariance = diag(GN_inputs.measurement.variance);
    
    
    
end




end



