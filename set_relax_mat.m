function [ES, b] = set_relax_mat(T1, T2, dur, kmax, N, S)
% SET_RELAX_MAT Set Relaxation matrices (no diffusion)
%
%  [ES, b] = set_relax_mat(T1, T2, dur, kmax, N, S) 
%
% Taken from code within EPG_GRE.m written by Shaihan Malik
%
% D.Atkinson@ucl.ac.uk

E1 = exp(-dur/T1);
E2 = exp(-dur/T2);
E = diag([E2 E2 E1]);

%%% regrowth
b = zeros([N 1]);
b(3) = 1-E1;%<--- just applies to Z0

%%% Add in diffusion at this point 
% if exist('diff','var')
%     E = E_diff(E,diff,kmax,N);
% else
    % If no diffusion, E is the same for all EPG orders
    E = spdiags(repmat([E2 E2 E1],[1 kmax+1])',0,N,N);
% end
    

%%% Composite relax-shift
ES=E*S;
ES=sparse(ES);
