function Fig1_Ernst
% Fig1_Ernst
% Plots Signal Intensity as a function of flip angle for two TRs and a
% range of constrast agent concentrations.
%
% Requires EPG-X function (see cep_doctor for details).
%
% D.Atkinson@ucl.ac.uk
%
% See also cep_doctor ssSPGR
%

% Copyright 2018 University College London.

T10 = 1.400 ; % [s]            Initial T1 of prostate
r1 = 4 ; % [Mm-1 s-1]          Relaxivity of contrast agent
FAs = [0:0.1:50]; % [degrees]  Flip Angles
FAs_rad = d2r(FAs) ; % [radians]

Cs = fliplr([0 0.01 0.025 0.05 0.1]) ; % mM concentrations. Reversed order for clearer legend
Rs = 1/T10 + r1.*Cs ;
T1s = 1000./ Rs ; % ms

TRs = [ 10 5 ] ; % [ms] TR 

hf = figure('DefaultAxesFontSize',12,...
      'DefaultAxesFontWeight', 'bold', ...
      'DefaultAxesLineWidth',2, ...
      'Units','centimeters') ;

lw = [2 2] ; % plot linewidths
col = {[1 0 0], [0 0 1] } ; % line colors


for iTR = 1:length(TRs)
    for iT1 = 1:length(T1s)
        SI = ssSPGR(FAs_rad, TRs(iTR), T1s(iT1)) ; % signal
        lcol = col{iTR} * (length(T1s) -iT1)/(length(T1s)-1) ; % line colour
        plot(FAs, SI, 'LineWidth',lw(iTR), 'Color', lcol, 'DisplayName',[num2str(TRs(iTR),'%2d'), 'ms,', ...
        num2str(Cs(iT1),' %5.3f'), 'mM'])
        hold on
        grid on
    end
end
lgd = legend ;
lgd.FontSize = 12 ;
lgd.FontName = 'FixedWidth' ;
lgd.FontWeight = 'bold' ;

xlabel('Flip Angle (degrees)')
ylabel('Steady-State Signal Intensity')

disp(['T1s: ',num2str(T1s)])
disp(['[CA]: ',num2str(Cs)])


