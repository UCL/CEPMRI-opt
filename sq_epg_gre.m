function [F0,Fn,Zn,F] = sq_epg_gre(sq, varargin)
% SQ_EPG_GRE Extended phase graph for gradient echo (GRE)
% Adapted from Shaihan Malik's EPG_GRE code in his EPG-X repository.
%
% [F0,Fn,Zn,F] = sq_epg_gre(sq, varargin)
%
% David Atkinson   D.Atkinson@ucl.ac.uk
%
% See also cep_doctor

% Copyright 2018, University College London.


for ipv = 1:2:length(varargin)
   switch varargin{ipv}
       case 'kmax'
           kmax = varargin{ipv+1} ;
       otherwise
           error(['Unknown parameter input: ',varargin{ipv}])
   end
end


nsq = length(sq) ; % no. sequence objects (previously included those ...
                   % with no RF but coimplicated)

np = nsq  ;


% if not defined, assume want max
if ~exist('kmax','var')
    kmax = np - 1;
end

if isinf(kmax)
    % this flags that we don't want any pruning of pathways
    allpathways = true;
    kmax = np - 1; % this is maximum value
else
    allpathways = false;
end

%%% Variable pathways
if allpathways
    kmax_per_pulse = 0:kmax;
else
    kmax_per_pulse = [1:ceil(np/2) (floor(np/2)):-1:1];
    kmax_per_pulse(kmax_per_pulse>kmax)=kmax;
     
    if max(kmax_per_pulse)<kmax
        kmax = max(kmax_per_pulse);
    end
end

%%% Number of states is 6x(kmax +1) -- +1 for the zero order
N=3*(kmax+1);

%%% Build Shift matrix, S
S = EPG_shift_matrices(kmax);
S = sparse(S);

%%% Pre-allocate RF matrix
T = zeros(N,N);
T = sparse(T);

% store the indices of the top 3x3 corner, this helps build_T
i1 = [];
for ii=1:3
    i1 = cat(2,i1,sub2ind(size(T),1:3,ii*ones(1,3)));
end

%%% F matrix (many elements zero, not efficient)
F = zeros([N np]); %%<-- records the state after each RF pulse 

%%% Initial State
FF = zeros([N 1]);
FF(3)=1;   % M0 - could be variable



% now loop over sq objects
ip = 0 ; % RF pulse counter (should be same as isq as now all sequence 
         % objects need to have an RF.

for isq = 1: nsq
    
    %%% Compute ES this sq object for evolution during period
    [ES, b] = set_relax_mat(sq(isq).T1, sq(isq).T2, sq(isq).dur, kmax, N, S) ;
    
    % Calculate RF related variables
    
    
    ip = ip + 1 ;
    
    A = RF_rot( d2r(sq(isq).FA), d2r(sq(isq).phase)) ;
    
    %%% Variable order of EPG, speed up calculation
    kmax_current = kmax_per_pulse(ip);
    kidx = 1:3*(kmax_current+1); %+1 because states start at zero
    
    %%% Replicate A to make large transition matrix
    build_T(A); % I think safe to do this within a loop as it overrides the same elements??
    
    %%% Apply flip and store this: splitting these large matrix
    %%% multiplications into smaller ones might help
    F(kidx,ip)=T(kidx,kidx)*FF(kidx);
    
    
    % update FF
    FF(kidx) = ES(kidx,kidx)*F(kidx,ip)+b(kidx);
    
    
    if sq(isq).crush == true
        mask = zeros(size(FF)) ;
        mask(3:3:end) = 1 ;
        FF = FF .* mask ;
    end
    
    % Deal with complex conjugate after shift
    FF(1)=conj(FF(1)); %<---- F0 comes from F-1 so conjugate 
end

%%% Return signal
F0=F(1,:);

%%% phase demodulate
phases = [sq.phase] ;
F0 = F0(:) .* exp(-1i*d2r(phases(:))) *1i;


%%% Construct Fn and Zn
idx=[fliplr(5:3:size(F,1)) 1 4:3:size(F,1)]; 
kvals = -kmax:kmax;

%%% Now reorder
Fn = F(idx,:);
%%% Conjugate
Fn(kvals<0,:)=conj(Fn(kvals<0,:));

%%% Similar for Zn
Zn = F(3:3:end,:);
 
    %%% NORMAL EPG transition matrix as per Weigel et al JMR 2010 276-285 
    function Tap = RF_rot(a,p)
        Tap = zeros([3 3]);
        Tap(1) = cos(a/2).^2;
        Tap(2) = exp(-2*1i*p)*(sin(a/2)).^2;
        Tap(3) = -0.5*1i*exp(-1i*p)*sin(a);
        Tap(4) = conj(Tap(2));
        Tap(5) = Tap(1);
        Tap(6) = 0.5*1i*exp(1i*p)*sin(a);
        Tap(7) = -1i*exp(1i*p)*sin(a);
        Tap(8) = 1i*exp(-1i*p)*sin(a);
        Tap(9) = cos(a);
    end

   function build_T(AA)
        ksft = 3*(3*(kmax+1)+1);
        for i2=1:9
            T(i1(i2):ksft:end)=AA(i2);
        end
    end


end