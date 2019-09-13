function [sq, series] = build_seq(scheme, varargin)
% BUILD_SEQ  Builds a sequence for EPG simulation
%
%  [sq, series] = build_seq(scheme)
%  [sq, series] = build_seq(scheme, param, val, ... )
%
%  sq is a structure array of sequence objects,
%  series is a structure with properties of the overall sequence
%    series also has material properties (T1, T2) which ideally should be
%    in a separate structure.
%    The fields ntotal and nstartup refer to the number of RFs within a
%    single TFE shot. The nstartup correspond to dummies used to allow
%    approach towards steady-state.
%
%  The individual 'schemes' were originally written to model specific
%  scenarios with hard coded T1 etc. In those cases, the 
%  input parameter / value pairs will be overridden. The ouput series
%  structure will contain the used values, except for the rounded SPAIR
%  duration.
%  See discussion below, but the time allocated in the sequence for the 
%  SPAIR duration is rounded to the nearest TR. The inv_dur is not rounded 
%  (though it probably should) so should be passed in as a multiple of a TR.
%
%  scheme can be:
%   'generic_proset' - generic single shot 
%   'generic_SPAIRwsweep' - generic SPAIR with flip angle sweep
%   'generic_SPAIRnosweep'  - without flip angle sweep
%   'generic_inv' - generic with inversion pulse
%   'gel18_noFS', 'gel18_ ...' were 1Dtestmode seqeunces (2D) on TO5 gel 18
%   'Malik1'
%   'MalikSpoiledWPrep'
%   'MalikSSFPWPrep'
%   'ssSPGR'     single shot
%   'msSPGR'     multi-shot
%   'msINV'      multi-shot with inversions on resonance
%
% parameters can be:
%   TR [5], FA [10], T1 [1400], T2 [80], spoil_incr [150], diff [[]], 
%   SPAIR_dur [20.49], inv_dur [TR], sweep ['quad']
%   plot_sequence [false]
%
% For EPG modelling of GRE, the unbalanced gradient area between RFs should
% be the same. If not, there needs to be a modification to the shifting
% operation. If the transverse magentisation is completely spoiled by a
% crusher, then modelling a subsequent delay as one that includes gradients 
% is OK as there is no transverse magnetisation for them to dephase. So, the SPAIR
% will be modelled as a crush and delay. The crush will be implemented by
% setting the F states to zero. This is an approximation as really they
% should be shifted, but there should also be their diffusion effect added.
%
% The inversion has a crusher on the slice select and this means that we
% should consider EPG in another dimension (the first really representing
% the readout). The assumption of an inversion being a 180 followed by a 
% crusher seems to match measurements using the 1Dtestmode.  
%
%
% David Atkinson  D.Atkinson@ucl.ac.uk, with help from Shaihan Malik.
%
% See also QFLIP_SWEEP

% Copyright 2018, University College London.

% defaults - may be overwritten in specific schemes
TR = 5 ;
FA = 10  ;
spoil_incr = 150 ; 
T1 = 1400 ; 
T2 = 80 ;
diff = [] ;
SPAIR_dur = 20.49 ; % this will be rounded to nearest TR in add_seq
inv_dur = TR ; % not rounded 
sweep = 'quad' ; % flip angle sweep type, see qflip_sweep
plot_sequence = false ;

% user settings
for ipv = 1:2:length(varargin)
    param = varargin{ipv} ;
    val = varargin{ipv+1} ;
    switch param
        case {'TR', 'FA', 'spoil_incr', 'T1', 'T2', 'plot_sequence'}
            eval([param ,' = val ;'])
        
        otherwise
            error(['Unknown input param: ',param])
    end
end

series.scheme = scheme ;
series.TR = TR ;
series.FA = FA ;
series.spoil_incr = spoil_incr ; % RF spoiling phase angle incr
series.T1 = T1 ;
series.T2 = T2 ;
series.diff = [] ;
series.SPAIR_dur = SPAIR_dur ;
% series.delay_dur = 4*series.TR ;
series.inv_dur = series.TR ; % approximation
series.sweep = sweep ;

switch scheme
    case 'generic_proset'
        % Generic ProSet for 10s, no flip angle sweep, no dummies,
        % single-shot
        
        series.ntotal = round( 10000 / series.TR) ;
        series.sweep = 'none' ; 
        series.nstartup = 0 ;
        next_tstart = 0 ;
        next_spoil_index = 1 ;
        
        [sq, next_tstart, next_spoil_index] = build_shot(series, next_tstart, next_spoil_index) ;
        
    case 'generic_SPAIRwsweep'
        % Generic SPAIR with sweep
        series.sweep = 'quad' ;
        series.nstartup = 7 ; % taken from actual seq with TR 5.7, FA 10, dyn time 10.3s
        series.ntotal = 38 + series.nstartup ;
        
        nshot = round(10000 / (series.TR * series.ntotal + series.SPAIR_dur) ) ;
                     
        next_tstart = 0 ;
        next_spoil_index = 1 ;
        
        sq = [] ;
        for ishot = 1: nshot
             [sq_add, next_tstart, next_spoil_index] = add_sq('SPAIR', series, next_tstart, next_spoil_index) ;
        
             [sq_shot, next_tstart, next_spoil_index] = build_shot(series, next_tstart, next_spoil_index) ;
        
             sq = cat(2, sq, sq_add, sq_shot) ;
        end
        
    case 'generic_SPAIRnosweep'
        % Generic SPAIR WITHOUT sweep
        series.sweep = 'none' ;
        series.nstartup = 0 ; 
        series.ntotal = 45 + series.nstartup ;
        
        nshot = round(10000 / (series.TR * series.ntotal + series.SPAIR_dur) ) ;
                     
        next_tstart = 0 ;
        next_spoil_index = 1 ;
        
        sq = [] ;
        for ishot = 1: nshot
             [sq_add, next_tstart, next_spoil_index] = add_sq('SPAIR', series, next_tstart, next_spoil_index) ;
        
             [sq_shot, next_tstart, next_spoil_index] = build_shot(series, next_tstart, next_spoil_index) ;
        
             sq = cat(2, sq, sq_add, sq_shot) ;
        end
        
    case 'generic_inv'
        
        % make about 10s long and with about 1.6s between inversions
        series.ntotal = round(1600 / series.TR) ;
        series.sweep = 'none' ; 
        series.nstartup = 0 ;
        next_tstart = 0 ;
        next_spoil_index = 1 ;
        nshot = round(10000 / (series.TR * series.ntotal) ) ;
        series.inv_dur = series.TR ; % use a TR multiple (typical inv time is 10ms)
        
        sq = [] ;
        for ishot = 1: nshot
             [sq_add, next_tstart, next_spoil_index] = add_sq('inv', series, next_tstart, next_spoil_index) ;
        
             [sq_shot, next_tstart, next_spoil_index] = build_shot(series, next_tstart, next_spoil_index) ;
        
             sq = cat(2, sq, sq_add, sq_shot) ;
        end
        
    case 'gel18_noFS'
        % Scan 34, 20 Mar 2018
        series.TR = 10.9593 ;
        series.FA = 12 ;
        series.T1 = 1460 ;
        series.T2 = 150 ;
        series.ntotal = 950 ;
        series.sweep = 'none' ;
        series.nstartup = 0 ;
        next_tstart = 0 ;
        next_spoil_index = 1 ;
        
        [sq, next_tstart, next_spoil_index] = build_shot(series, next_tstart, next_spoil_index) ;
        
    case 'gel18_proset'
        % Scan 46, 20 Mar 2018
        series.TR = 13.2269 ;
        series.FA = 12 ;
        series.T1 = 1460 ; % approx values for gel 18
        series.T2 = 150 ;
        series.ntotal = 950 ;
        series.sweep = 'none' ; 
        series.nstartup = 0 ;
        next_tstart = 0 ;
        next_spoil_index = 1 ;
        
        [sq, next_tstart, next_spoil_index] = build_shot(series, next_tstart, next_spoil_index) ;
    case 'gel18_SPAIR'
        % Scan 50, 20 Mar 2018. 
        % This 2D seq had large gaps between shot and SPAIR - modelled here as
        % a long SPAIR
        series.TR = 10.9593 ;
        series.FA =  12 ;
        series.T1 = 1460 ; % approx values for gel 18
        series.T2 = 150 ;
        series.ntotal = 29 ;
        series.sweep = 'quad' ; 
        series.nstartup = 4 ;
        next_tstart = 0 ;
        next_spoil_index = 1 ;
        nshot = 38 ;
        series.SPAIR_dur = 503.9 - (series.ntotal*series.TR) ; % SPAIR of 21ms plus gap
                           % 503.9 was time between shots
        
        sq = [] ;
        for ishot = 1: nshot
             [sq_add, next_tstart, next_spoil_index] = add_sq('SPAIR', series, next_tstart, next_spoil_index) ;
        
             [sq_shot, next_tstart, next_spoil_index] = build_shot(series, next_tstart, next_spoil_index) ;
        
             sq = cat(2, sq, sq_add, sq_shot) ;
        end
        
    case 'gel18_inv'
        % scan 52, 20 March 2018
        series.TR = 13.17 ;
        series.FA =  12 ;
        series.T1 = 1460 ; % approx values for gel 18
        series.T2 = 150 ;
        series.ntotal = 72 ;
        series.sweep = 'none' ; 
        series.nstartup = 0 ;
        next_tstart = 0 ;
        next_spoil_index = 1 ;
        nshot = 13 ;
        series.inv_dur = series.TR ; % inv dur was really 10ms, but using TR here
        
        sq = [] ;
        for ishot = 1: nshot
             [sq_add, next_tstart, next_spoil_index] = add_sq('inv', series, next_tstart, next_spoil_index) ;
        
             [sq_shot, next_tstart, next_spoil_index] = build_shot(series, next_tstart, next_spoil_index) ;
        
             sq = cat(2, sq, sq_add, sq_shot) ;
        end
        
    case 'Malik1'
        % Gives the same signal exactly as first part of Malik's 
        % test1_steady_state_GRE  (perform a plot([0:5:999],abs(s0)) after running )
        series.TR = 5 ;
        series.FA = 10 ;
        series.spoil_incr = 117 ;
        series.T1 = 779 ;
        series.T2 = 45 ;
        series.ntotal = 200 ; % ntotal is for one shot and excludes prep pulses
        series.nstartup = 0 ; % startup pulses in one shot
        
        next_tstart = 0 ;   % start of entire sequence
        next_spoil_index = 1 ;
        
        [sq, next_tstart, next_spoil_index] = build_shot(series, next_tstart, next_spoil_index) ;
        
    case 'MalikSpoiledWPrep'
        series.TR = 12 ;
        series.FA = 40 ;
        series.spoil_incr = 150 ;
        
        series.nstartup = 0 ;
        
        next_tstart = 0 ;   % start of entire sequence
        next_spoil_index = 1 ;
        
        %%% White matter model (see Gloor 2008)
        f = 0.1166;  %% F=0.132 = f/(1-f) => f=0.1166
        kf = 4.3e-3;
        kb = kf * (1-f)/f;
        R1f = 1/779; % ms^-1
        R1b = 1/779; %<- Gloor fix as T1f
        R2f = 1/45;
        R1obs = 0.5*(R1f + kf + R1b + kb)-0.5*sqrt((R1f + kf + R1b + kb).^2 ...
            -4*(R1f*R1b + R1f*kb + R1b*kf));
        
        %%% Flip angle variation
        ss = [0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 10 10 11 11 12 12 13 13 14 14 15 15]/15;
        ss = [ss fliplr(ss)];
        s0 = zeros([1 16]);
        alphas = [s0 s0 ss ss s0 s0 ss];
        alphas=sin(alphas*pi/2)*series.FA;
        
        series.ntotal = length(alphas) ; %! ntotal is for a shot, not whole seq
        
        series.T1 = 1/R1obs ;
        series.T2 = 1/R2f ;
        series.inv_dur = 0 ;
        
        % add an inversion at start (note currently this code does not
        % spoil)
        sq = [] ;
        [sq_add, next_tstart, next_spoil_index] = ... 
            add_sq('inv', series, next_tstart, next_spoil_index) ;
        sq = cat(2, sq, sq_add) ;
        
        [sq_shot, next_tstart, next_spoil_index] = ...
                build_shot(series, next_tstart, next_spoil_index) ;
        % modify FAs to be as calculated above in alphas (not first as this is the inversion).
        sq = cat(2, sq, sq_shot) ;
        
        for iFA = 1:series.ntotal
          sq(1+iFA).FA = alphas(iFA) ;
        end
        
        
    case 'MalikSSFPWPrep'
        % call 'MalikSpoiledWPrep' and modify the phase cycling
        [sq, series] = build_seq('MalikSpoiledWPrep') ;
        
        for isq = 1:series.ntotal
            if mod(isq,2) == 0
                sq(1+isq).phase =  0 ;
            else
                sq(1+isq).phase =  180 ;
            end
        end
        
        series.kmax = inf ;
        
    case 'ssSPGR' % single-shot
        
        series.nstartup = 4 ; % no. of startup RFs in a shot (no ADC)
        series.ntotal = 1300 ; % no. of RFs in a shot
        series.sweep = 'none' ;
        
        next_tstart = 0 ;   % start of entire sequence
        next_spoil_index = 1 ;
        
        [sq, next_tstart, next_spoil_index] = build_shot(series, next_tstart, next_spoil_index) ;
    
    case 'ProSetSPGR' % single-shot
        series.TR = 6 ;
        series.sweep = 'none' ; % distinguishing feature
        series.nstartup = 4 ; % no. of startup RFs in a shot (no ADC)
        series.ntotal = 1433 ; % no. of RFs in a shot
        
        next_tstart = 0 ;   % start of entire sequence
        next_spoil_index = 1 ;
        
        [sq, next_tstart, next_spoil_index] = build_shot(series, next_tstart, next_spoil_index) ;
        
    case 'SPAIRSPGR'
%         series.TR = 5.8 ;
%         series.FA = 15 ;
        
        series.nstartup = 7 ; %7
        series.ntotal = 38  + series.nstartup ; % ?
        series.spoil_incr = 150 ;
        
        nshot = 35 ;
        sq = [] ;
        
        series.sweep = 'quad' ; 
        next_tstart = 0 ; next_spoil_index = 1 ;
        
        for ishot = 1:nshot
             [sq_add, next_tstart, next_spoil_index] = add_sq('SPAIR', series, next_tstart, next_spoil_index) ;
        
             [sq_shot, next_tstart, next_spoil_index] = build_shot(series, next_tstart, next_spoil_index) ;
        
             sq = cat(2, sq, sq_add, sq_shot) ;
             
        end
        
    case 'SPAIRSPGRnosweep'
%         series.TR = 5.8 ;
%         series.FA = 15 ;
        
        series.nstartup = 7 ; %7
        series.ntotal = 38  + series.nstartup ; % ?
        series.spoil_incr = 150 ;
        
        nshot = 35 ;
        sq = [] ;
        
        series.sweep = 'none' ; 
        next_tstart = 0 ; next_spoil_index = 1 ;
        
        for ishot = 1:nshot
             [sq_add, next_tstart, next_spoil_index] = add_sq('SPAIR', series, next_tstart, next_spoil_index) ;
        
             [sq_shot, next_tstart, next_spoil_index] = build_shot(series, next_tstart, next_spoil_index) ;
        
             sq = cat(2, sq, sq_add, sq_shot) ;
             
        end
        
    case 'splitSPGR' % split single-shot, inserting tests - should make no difference
        
        series.nstartup = 7 ; % no. of startup RFs in a shot (no ADC)
        series.ntotal = 500 ; % no. of RFs in a shot
        
        next_tstart = 0 ;   % start of entire sequence
        next_spoil_index = 1 ;
        
        [sq_shot1, next_tstart, next_spoil_index] = build_shot(series, next_tstart, next_spoil_index) ;
        
        % OK without zero length delay (just concatenated shots).
        % OK with single_TR inserted
        % Not OK with delay of zero.
        
        [sq_add, next_tstart, next_spoil_index] = add_sq('SPAIR', series, next_tstart, next_spoil_index) ;
        
        series.nstartup = 0 ;
        [sq_shot2, next_tstart, next_spoil_index] = build_shot(series, next_tstart, next_spoil_index) ;
        
        sq = cat(2, sq_shot1, sq_add, sq_shot2) ;
                
    case 'msSPGR' % multi-shot. No longer with gaps (for modelling SPAIR)
        series.TR = 10 ;
        series.FA = 15 ;
        series.delay_dur = 20.49 ;% 20.49 ;
        series.nstartup = 7 ; %7
        series.ntotal = 38  + series.nstartup ; % ?
        series.spoil_incr = 150 ;
        series.sweep = 'none' ;
        
        nshot = 35 ;
        sq = [] ;
        
        next_tstart = 0 ; next_spoil_index = 1 ;
        
        for ishot = 1:nshot
            % add pre-pulse (here just delay in this case)
% %             [sq_add, next_tstart, next_spoil_index] = add_sq('delay', series, next_tstart, next_spoil_index) ;
% %             sq_add.crush = true ;
% %             sq = cat(2, sq, sq_add) ;
            
            %!! Should be able to emulate single shot if delay does nothing
            % add shot
            [sq_shot, next_tstart, next_spoil_index] = build_shot(series, next_tstart, next_spoil_index) ;
            sq = cat(2,sq,sq_shot) ;
        end
        
    case 'msINV'  % multi-shot with inversion pulses for contrast manipulation
        series.TR = 7.0 ;
        series.nstartup = 4 ;
        series.ntotal = 225 + 4 ;
        series.sweep  = 'quad' ;
        
        
        nshot = 6 ;
        sq = [] ;
        
        next_tstart = 0 ; next_spoil_index = 1 ;
        
        for ishot = 1:nshot
            % add pre-pulse 
            [sq_add, next_tstart, next_spoil_index] = add_sq('inv', series, next_tstart, next_spoil_index) ;
            sq = cat(2, sq, sq_add) ;
            
            % add shot
            [sq_shot, next_tstart, next_spoil_index] = build_shot(series, next_tstart, next_spoil_index) ;
            sq = cat(2,sq,sq_shot) ;
        end
        
    otherwise
        error(['(build_seq) sequence not implemented: ',scheme])
end

if plot_sequence
    plot_seq(scheme, series,sq) 
end

end


function [sq_add, next_tstart, next_spoil_index] = add_sq(sqtype, series, next_tstart, next_spoil_index) 

switch sqtype
    case 'SPAIR'
        % modelled as a delay and crush
        % duration rounded to nearest TR
        % EPG assumes all gaps have same gradient moment. Whilst this is
        % not true for these SPAIRs, if we assume the crush following the
        % SPAIR removes transverse magnetisation, we can ignore the
        % incorrect gradient moments. T1 effects are still modelled.
        
        nTR = round(series.SPAIR_dur / series.TR) ;
        
        for isq = 1: nTR
            sq_add(isq).dur = series.TR ;
            sq_add(isq).T1 = series.T1 ; sq_add(isq).T2 = series.T2 ; 
            sq_add(isq).tstart = next_tstart ;
            next_tstart = sq_add(isq).tstart + sq_add(isq).dur ;
            sq_add(isq).FA = 0 ;
            sq_add(isq).phase = 0 ;
            sq_add(isq).ADCon = false ;
            sq_add(isq).crush = false ; % except last, set later
        end
            
        sq_add(nTR).crush = true ; % crush at end
        
    case 'delay'
        warning(['Needs checking'])
        sq_add.crush = false ;
        sq_add.dur = series.delay_dur ;
        sq_add.T1 = series.T1; sq_add.T2 = series.T2 ;
        sq_add.tstart = next_tstart ;
        next_tstart = sq_add.tstart + sq_add.dur ;
        sq_add.FA = 0;
        sq_add.phase = 0 ;
        sq_add.ADCon = false ;
   
    case 'single_TR'   % for testing, single RF with rf-spoiling
        sq_add.crush = false ;
        sq_add.dur = series.TR ;
        sq_add.T1 = series.T1; sq_add.T2 = series.T2 ;
        sq_add.tstart = next_tstart ;
        next_tstart = sq_add.tstart + sq_add.dur ;
        sq_add.FA = series.FA ;
        
        np = next_spoil_index  ; % total number of RF spoil pulses from beginning to end of this shot
        p=[1:np];
        phi = cumsum((p-1).* series.spoil_incr) ; % phi will be length np, just use the last ntotal for this shot

        sq_add.phase = phi(end) ;  
        next_spoil_index = next_spoil_index + 1; 
        sq_add.ADCon = true ; 
    
    case 'inv' % inversion
        sq_add.dur = series.inv_dur ;
        sq_add.T1 = series.T1; sq_add.T2 = series.T2 ;
        sq_add.tstart = next_tstart ;
        next_tstart = sq_add.tstart + sq_add.dur ;
        sq_add.FA = 180 ;
        sq_add.phase = 0 ; 
        sq_add.ADCon = false ; 
        sq_add.crush = true ; % Previously (to March 27, 2018) had was false, 
        %  but there is a small crusher on the slice direction. 
        
    otherwise
        error(['sqtype not implemented: ',sqtype])
end


end

function [sq, next_tstart, next_spoil_index] = build_shot(series, next_tstart, next_spoil_index) 

shot_fa = qflip_sweep(series.FA, series.nstartup, series.ntotal, series.sweep) ; % compute flip angle sweep
np = next_spoil_index + (series.ntotal-1) ; % total number of RF spoil pulses from beginning to end of this shot
p=[1:np];
phi = cumsum((p-1).* series.spoil_incr) ; % phi will be length np, just use the last ntotal for this shot

for isq = 1: series.ntotal
    sq(isq).crush = false ; % no crush in a shot
    sq(isq).dur = series.TR ;
    sq(isq).T1 = series.T1 ; sq(isq).T2 = series.T2 ; 
    sq(isq).tstart = next_tstart ;
    next_tstart = sq(isq).tstart + sq(isq).dur ;
    
    sq(isq).FA = shot_fa(isq) ;
    
    sq(isq).phase = phi(np-series.ntotal+isq) ;
    if isq <= series.nstartup
       sq(isq).ADCon = false ;
    else
       sq(isq).ADCon = true ; 
    end
end

next_spoil_index = np + 1;



end
    

function plot_seq(scheme, series, sq)
figure('Name',['(build_seq) ',scheme])
nsq = length(sq) ;

ts = [sq.tstart] ;
FAs = [sq.FA] ;
plot(ts, FAs)
hold on

ADCs = [sq.ADCon] ;
yADC = -1*ones([1 length(ADCs)]) ;
loc = find(ADCs == false) ;
tsa = ts ;
tsa(loc) = [] ; yADC(loc) = [] ;

plot(tsa,yADC,'r.')

tot_dur = sq(nsq).dur + sq(nsq).tstart - sq(1).tstart ;

xlabel(['Time (ms).  dur: ',num2str(tot_dur)])
ylabel('FA (degrees)')

end

        