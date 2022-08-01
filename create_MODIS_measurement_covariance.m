%% create measurement prior


% By Andrew J. Buggee
%%

function [GN_inputs] = create_MODIS_measurement_covariance(GN_inputs,modis,modisInputs,pixels2use)


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

    GN_inputs.measurement.variance = zeros(n,length(pixels2use.res1km.linearIndex));
    GN_inputs.measurement.covariance = zeros(n,n,length(pixels2use.res1km.linearIndex));

    % Step through each pixel being used
    for pp = 1:length(pixels2use.res1km.linearIndex)

        % Step through each band
        for ii = 1:n
            % if each uncertainty represents the standard deviation, the
            % variance is the square of each value.
            % the refelctance uncertanties are listed in percentages. So we
            % multiply these percentages with the modis reflectance values to
            % get the uncertainty in reflectance.

            % Grab the row and column
            r = pixels2use.res1km.row(pp);
            c = pixels2use.res1km.col(pp);

            % Lets start by converting the percentage to a decimal
            GN_inputs.measurement.percent_uncertainty = 0.01*modis.EV1km.reflectanceUncert(r,c,ii);

            % Lets assume the percentage given is the standard deviation
            % According to King and Vaughn (2012): 'the values along the main
            % diagonal correspond to the square of the uncertainty estimate for
            % each wavelength channel'

            GN_inputs.measurement.variance(ii,pp) = (modis.EV1km.reflectance(r,c,ii).* GN_inputs.measurement.percent_uncertainty).^2;


        end

        % Create a diagonal matrix where each entry is the variance of that
        % spectral channel for reflectance measurements
        GN_inputs.measurement.covariance(:,:,pp) = diag(GN_inputs.measurement.variance(:,pp));

    end





end




end



