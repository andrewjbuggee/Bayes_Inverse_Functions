%% Plot the VOCALS-Rex in-situ data along with a MODIS retrieved droplet size and optical depth for a single pixel


% By Andrew John Buggee

%%

function plot_vocalsRex_with_MODIS_retrieved_re_and_vertProf_retrieval(vocalsRex, modis, GN_outputs)


% ------------------------------------------------------------------
% ---------------- Compute Liquid Water Path -----------------------
% ------------------------------------------------------------------

LWP_vocals = trapz(vocalsRex.altitude, vocalsRex.lwc);                 % grams of water/m^2

% ----- Compute the CDP uncertainty -----
re_uncertainty = cloud_droplet_probe_uncertainty_estimate(vocalsRex.re);



% I want a subplot with the number concentration and altitude, and the
% effective radius with altitude
nice_blue = [0 0.4470 0.741];
nice_orange = [0.8500, 0.3250, 0.0980];

figure;
errorbar(flipud(vocalsRex.re), vocalsRex.tau', flipud(re_uncertainty), 'horizontal','-o','Color','black', 'MarkerSize',10,...
    'MarkerFaceColor','black','LineWidth',1);
set(gca,'YDir','reverse')
ylabel('$\tau$','interpreter','latex','FontSize',35);
xlabel('$r_{e}$ $$(\mu m)$$','Interpreter','latex')
title('Comparison between in-situ and MODIS retrieved $r_e$', 'Interpreter','latex')
grid on; grid minor; hold on;

% Plot the Gauss-Newton Retreval
for pp = 1:size(GN_outputs.re_profile,2)
    plot(GN_outputs.re_profile(:,pp), GN_outputs.tau_vector(:,pp), 'Color',mySavedColors(pp,'fixed'),'LineStyle',':', 'LineWidth',3)
end

% Fit a curve to the in-situ data to show the capability we are interested
% in devloping

% curve_fit_linewidth = 6;
% curve_fit_color = mySavedColors(1,'fixed');                        % Bright pink
%
% f = fit(tau', fliplr(double(vocalsRex.re))', 'smoothingspline','SmoothingParam',0.9);
% hold on;
% plot(f(tau),tau,'Color',curve_fit_color,'LineStyle',':', 'LineWidth',curve_fit_linewidth);

% Plot the z-space in meters on the right axis
yyaxis right
ylim([0, vocalsRex.altitude(end) - vocalsRex.altitude(1)])
set(gca,'YColor','black')
ylabel('Altitude within cloud $(m)$', 'Interpreter','latex','FontSize',30);
yyaxis left

% Label cloud top and cloud bottom
% Create textbox
annotation('textbox',[0.02,0.865079365079366,0.051,0.077777777777778],...
    'String',{'Cloud','Top'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',22,...
    'FitBoxToText','off');

% Create textbox
annotation('textbox',[0.02,0.096825396825397,0.051,0.077777777777778],...
    'String',{'Cloud','Bottom'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',22,...
    'FitBoxToText','off');



% Plot the modis droplet estimate as a constant vertical line
xl0 = xline(modis.cloud.effRadius17(vocalsRex.modisIndex_minDist_median),':',...
    ['MODIS $$r_{2.1} = $$',num2str(modis.cloud.effRadius17(vocalsRex.modisIndex_minDist_median)), '$$\mu m$$'], 'Fontsize',22,...
    'Interpreter','latex','LineWidth',2,'Color',nice_blue);
xl0.LabelVerticalAlignment = 'bottom';

% Plot the MODIS optical depth estiamte as a constant horizontal line
yl0 = yline(modis.cloud.optThickness17(vocalsRex.modisIndex_minDist_median),':',...
    ['MODIS $$\tau_{2.1} = $$',num2str(modis.cloud.optThickness17(vocalsRex.modisIndex_minDist_median))], 'Fontsize',22,...
    'Interpreter','latex','LineWidth',2,'Color',nice_orange);
yl0.LabelVerticalAlignment = 'top';
yl0.LabelHorizontalAlignment = 'left';





% Let's compute the mean number concentration within this cloud and print
% it on our plot

mean_Nc = mean(vocalsRex.Nc);

dim = [.137 .425 .3 .3];
str = ['$$< N_c >_{in-situ} = \;$$',num2str(round(mean_Nc)),' $$cm^{-3}$$',newline,'$$LWP_{in-situ} = $$',num2str(round(LWP_vocals,1)),' $$g/m^{2}$$'];
annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',20,'FontWeight','bold');
set(gcf,'Position',[0 0 1200 630])

% Create a Legend with only the two black curves
%legend('Vocals Rex In-situ Measurement', 'Desired Retrieval Profile', 'Interpreter','latex', 'Location','best')
legend({'Vocals Rex In-situ Measurement', 'Retrieved Profile - first pixel','Retrieved Profile - median pixel',...
    'Retrieved Profile - last pixel'}, 'Interpreter','latex', 'Location','northwest', 'FontSize', 20)


end