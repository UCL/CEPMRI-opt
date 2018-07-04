function Fig2_SIT1vsConc
% Fig2_SIT1vsConc
% Plots Signal Intensity and T1 as a function of contrast agent concentration  
%
% Requires functions from EPG-X (run cep_doctor for details)
%
% D.Atkinson@ucl.ac.uk
%
% See also cep_doctor ssSPGR d2r

% Copyright 2018 University College London.

T10 = 1.400 ; % [s]            Initial T1 of prostate
r1 = 4 ; % [Mm-1 s-1]          Relaxivity of contrast agent
FAs = [9 13]; % [degrees]  Flip Angles
TRs = [ 5 10 ] ; % [ms] TR 

FAs_rad = d2r(FAs) ; % [radians]

Cs = [0 : 0.01 : 1] ; % mM concentrations. 
Rs = 1/T10 + r1.*Cs ;
T1s = 1000./ Rs ; % ms


hf = figure('DefaultAxesFontSize',12,...
      'DefaultAxesFontWeight', 'bold', ...
      'DefaultAxesLineWidth',2, ...
      'Units','centimeters') ;

lw = 2 ; % plot linewidth

yyaxis left

for iperm = 1:length(FAs_rad)

    SI = ssSPGR(FAs_rad(iperm), TRs(iperm), T1s) ;

    dSIdC = (SI(2)-SI(1))/(Cs(2)-Cs(1)) ;
    plot(Cs,SI,'LineWidth',lw, 'DisplayName',['FA: ',num2str(FAs(iperm)), ...
        ' TR: ',num2str(TRs(iperm)), ' slope: ',num2str(dSIdC)])
    hold on
%     SI0 = SI(1) ;
%     [mn, loc] = min(abs(SI-SI0*1.5)) ;
%     plot(Cs(loc),SI(loc),'o')
    
end


text( 0.5, 0.08,'\bf\fontsize{12}\leftarrow 5ms, 9\circ')
text( 0.5, 0.116,'\bf\fontsize{12}\leftarrow 10 ms, 13\circ')

grid on
axis([0 min([2 max(Cs)]) 0 0.2])

xlabel('Contrast agent concentration (mM)', 'FontWeight','bold')
ylabel('Steady State Signal Intensity', 'FontWeight','bold')

yyaxis right
plot(Cs,T1s, 'LineWidth',lw, 'DisplayName','T1 vs [CA]')
text(0.5, 390,'\rightarrow','FontWeight','bold', 'FontSize',20)
ylabel('T1 (ms)', 'FontSize', 10, 'FontWeight','bold')


% lgd = legend ;
% lgd.FontSize = 10 ;
% lgd.FontName = 'FixedWidth' ;
% lgd.FontWeight = 'bold' ;



