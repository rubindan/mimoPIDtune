function [ C,objval ] = auto_mimo_pid( P,w,Smax,Tmax,Qmax,tau,Options )
%AUTO_MIMO_PID Autotune a MIMO PID for a nominal plant by LMIs. 
%  [ S. Boyd, M. Hast, and J. Astrom, "MIMO PID tuning via iterated LMI ]
%  [    restiction", Inter. Jur. of Robust & Nonlinear Control, 2016.   ]
%
%   Usage:
%   [ C,objval] = AUTO_MIMO_PID( P,w,Smax,Tmax,Qmax,tau,Stracture ) 
%   returns the optimized MIMO PID controller C. Can also return objval,  
%   the value of the objective after optimization.
%
%   Inputs:
%   P       	plant in ss, tf, zpk, or frd form
%   w           vector of N frequencies responsibly chosen to cover the 
%               dynamic range 
%   Smax        sensitivity specification upper bound
%   Tmax        low frequency sensitivity specification upper bound
%   Qmax        cost of feedback specification upper bound
%   tau         derivative action time constant 
%   Options     optioanl struct with fields 'Structure', 'Initial', and 
%               'Sign' each with subfields Kp, Ki, and Kd; For 
%               - 'Structure', each subfield is a matrix with 0/1 elements, 
%                 where 0 constrain the corresponding element to be 0.
%               - 'Initial', each subfield is a matrix for initialization.
%               - 'Sign', each subfiled is an integer: 1 for positve, -1
%                 for negative
%
%   Smax, Tmax, and Qmax, are vectors of length N, specifying upper bound
%   at each frequency given in w. Reasonable values for Smax and Tmax are
%   in the range [1.1,1.6]. Reasonable value for Qmax is multipicative of 
%   1/min(svd(dcgain(P))) by a scalar in range [3,10].
%   
%   The time constant tau must be responsibly chosen as e.g. a modest
%   fraction of the desired closed-loop respopnse time.
% 
%   Outputs: 
%   C           optimized MIMO PID controller transfer function.
%   objval      value of the objective after optimization.   
%   
%   Type pid(C(i,j)) to view element ij of C in PID form. 
% 
%Requires YALMIP
% 
%Not yet tested for non-square plants!!!
%
% -----------------------------------------------------------------------
% |     This function is for bakcup only. use AUTO_MIMO_PID instead     |
% |        AUTO_MIMO_PID support both certain and uncertin plants       |
% -----------------------------------------------------------------------
%
%Created: Daniel Rubin, 14-Nov-2017.
%--------------------------------------------------------------------------

epsilon=0.01;

[p,m] = size(P); % number of I/O
N=length(w);
if isa(P,'frd'),
    Pw0=freqresp(P,w(1));
    P0=abs(Pw0).*sign(real(Pw0)); %P0=dcgain(P);
else
    P0=dcgain(P); % accurate solution for non numeric data
end

if nargin<7, Options=struct(); end

t0=0;
if isfield(Options,'Initial'),
    Kp0=Options.Initial.Kp;
    Ki0=Options.Initial.Ki;
    Kd0=Options.Initial.Kd;
else
    Kp0 = zeros(m,p);
    Ki0 = epsilon*pinv(P0);
    Kd0 = zeros(m,p);
end

fprintf('Optimizaing a %gx%g PID controller \n',m,p)
tic

for Iter=1:10,
    
    fprintf('Iteration %g \n',Iter)
    
    t = sdpvar(1);
    if isfield(Options,'Structure'),
        Kp = sparsesdpvar(Options.Structure.Kp);
        Ki = sparsesdpvar(Options.Structure.Ki);
        Kd = sparsesdpvar(Options.Structure.Kd);
    else
        Kp = sdpvar(m,p,'full');
        Ki = sdpvar(m,p,'full');
        Kd = sdpvar(m,p,'full');
    end
    if isfield(Options,'Sign'),
        constraints=[Options.Sign.Kp*Kp(:)>0
                     Options.Sign.Ki*Ki(:)>0
                     Options.Sign.Kd*Kd(:)>0 ];
    else
        constraints=[];
    end
    
    for k=1:N
        
        wk=w(k);
        s=1j*wk;
        Pk = freqresp(P,wk);
        Ck = Kp+(1/s)*Ki+s/(1+tau*s)*Kd;
        Ck0 = Kp0+(1/s)*Ki0+s/(1+tau*s)*Kd0;
        
        LMI0=makeLMI(P0*Ki,P0*Ki0,t*eye(p));
        
        Z = eye(p)+Pk*Ck;      
        Z0 = eye(p)+Pk*Ck0;
        Y1 = (1/Smax(k))*eye(p);
        LMI1 = makeLMI(Z,Z0,Y1);
        
        %Z2 = eye(p)+Pk*Ck;
        %Z20 = eye(p)+Pk*Ck0;
        Y2 = (1/Tmax(k))*Pk*Ck;
        LMI2 = makeLMI(Z,Z0,Y2);
        
        %Z3 = eye(p)+Pk*Ck;
        %Z30 = eye(p)+Pk*Ck0;
        Y3 = (1/Qmax(k))*Ck;
        LMI3 = makeLMI(Z,Z0,Y3);
        
        constraints = [constraints LMI0>0 LMI1>0 LMI2>0 LMI3>0];
        
    end
    
    ops = sdpsettings ('solver','sdpt3');
    solvesdp(constraints,-t,ops);
    
    Kp0 = double(Kp);
    Ki0 = double(Ki);
    Kd0 = double(Kd);
    t=double(t);
    
    if t-t0<0.001, break; end
    t0=t;
end

Kp = Kp0;
Ki = Ki0;
Kd = Kd0;

s = tf('s');
C = Kp+(1/s)*Ki+s/(1+tau*s)*Kd;
objval = norm(inv(P0*Ki)); % spectra norm (largest singular value)

fprintf('Complited after %g iterations! \nObjective value ||(P(0)Ki)^-1||=%g, t=%g.\n',Iter,objval,t);
fprintf('Total time: %g seconds \n',toc);

% NESTED FUNCTIONS
function [LMI] = makeLMI(Z,Z0,Y)
%Constructs the LMI 
LMI = [Z'*Z0+Z0'*Z-Z0'*Z0   Y'
       Y                    eye(p)];
%end of makeLMI function
end 

function [X] = sparsesdpvar(S)
%Define sparse sdpvar with given structure S   
if nnz(S)>0,
    [I,J,~] = find(S);
    X=sparse(I,J,sdpvar(nnz(S),1));
else
    X=S;
end
%end of sparsesdpvar function
end

%end of AUTO_MIMO_PID main function
end 