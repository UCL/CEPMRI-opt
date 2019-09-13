function shot_fa = qflip_sweep(alpha, nstartup, ntotal, sweep_type)
% QFLIP_SWEEP Flip angle sweep (default quadratic) for a GRE shot
%  shot_fa = flip_sweep(alpha, nstartup, ntotal)
%  shot_fa = flip_sweep(alpha, nstartup, ntotal, sweep_type)
% 
% Input
%  alpha - target flip angle
%  nstartup - number of dummy RFs in shot (sweep is up to 3x this number)
%  ntotal   - total number of RFs in shot
%  sweep_type - {'quad'} , 'none', 'quarter' 
%      for 'quarter', first flip angle is alpha/4, remainder are alpha
%
% Output 
%  shot_fa  [1  ntotal]  flip angles for the shot
%
% David Atkinson, D.Atkinson@ucl.ac.uk
%
% See also sq_epg_gre build_seq

% Copyright 2019, University College London.


if nargin < 4
    sweep_type = 'quad' 
end

shot_fa = alpha * ones([1 ntotal]) ; % Flip angles in shot before modification

switch sweep_type
    case 'quad'
        
        var_profs = min( 3*nstartup, ntotal - 1) ; % number of variable flip angles
        
        ind = [1:var_profs+1] ;
        shot_fa = shot_fa .* (1 - (((ind-1)-var_profs)/var_profs).^2 );

    case 'none'
        % keep uniform (no sweep modification)
    case 'quarter'
        % trial 'sweep' to reduce oscillations
        shot_fa(1) = alpha /4 ;
    otherwise
        error(['Unknown flip angle sweep method: ',series.sweep])
end

        
