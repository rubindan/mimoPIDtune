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
%   P       	plant in ss, tf, zpk, or frd form; input an array of plants 
%               to account for uncertaity and variations in parametes
%   w           vector of N frequencies responsibly chosen to cover the 
%               dynamic range 
%   Smax        sensitivity specification upper bound
%   Tmax        low frequency sensitivity specification upper bound
%   Qmax        cost of feedback specification upper bound
%   tau         derivative action time constant 
%   Options     optional struct which specifies additional options listed
%               below
%
%   Additional Options:
%   Structure   substructor with fields Kp, Ki, and Kd; each subfield is a 
%               matrix with 0/1 elements, where 0 constrain the corresponding 
%               element to be 0.
%   Initial,    substructor with fields Kp, Ki, and Kd; each subfield is a 
%               matrix for initialization.
%   Sign        substructor with fields Kp, Ki, and Kd; each subfiled is an 
%               integer: 1 for positve, -1 for negative
%   maxInterp   maximal number of interpolations (def=10)
%
%   For the uncertain or LTV case P is such that P(i,j,k) is the j input to
%   the i output transfer function of the k-th parametic case. P(i,j,1) is
%   assumed to be the nominal case.
%
%   Smax, Tmax, and Qmax, are vectors of length N, specifying upper bound
%   at each frequency given in w. Reasonable values for Smax and Tmax are
%   in the range 1.1--1.6. Reasonable value for Qmax is 3--10 multipicative 
%   of 1/min(svd(dcgain(P))).
%   
%   The time constant tau must be responsibly chosen as e.g. a modest
%   fraction of the desired closed-loop respopnse time.
% 
%   Outputs: 
%   C           optimized MIMO PID controller transfer matrix
%   objval      value of the objective after optimization.   
%   
%   Type pid(C(i,j)) to view element ij of C in PID form. 
% 
%Requires YALMIP and SDPT3
% 
%Not yet tested for non-square plants!!!
%
%Created: Daniel Rubin, 14-Nov-2017
%--------------------------------------------------------------------------


[p,m,l] = size(P);  % p = number of outputs                                
                    % m = number of inputs
                    % l = number of plant cases (uncertainty)
N=length(w);
Pw0=freqresp(P,w(1));
P0=abs(Pw0).*sign(real(Pw0));

if nargin<7, Options=struct(); end
% initialization of Kp, Ki, Kd:
if isfield(Options,'Initial')
    Kp0=Options.Initial.Kp;
    Ki0=Options.Initial.Ki;
    Kd0=Options.Initial.Kd;
else
    epsilon=w(1)/10;
    Kp0 = zeros(m,p);
    Ki0 = epsilon*pinv(P0(:,:,1));
    Kd0 = zeros(m,p);
end
t0=0;

if isfield(Options,'maxInterp')     % specify max number of interpolations
    maxInterp = Options.maxInterp;
else
    maxInterp = 10;
end

fprintf('Optimizaing a %gx%g PID controller \n',m,p)
tic

for Iter=1:maxInterp
    
    fprintf('Iteration %g \n',Iter)
    
    t = sdpvar(1);
    if isfield(Options,'Structure')
        Kp = sparsesdpvar(Options.Structure.Kp);
        Ki = sparsesdpvar(Options.Structure.Ki);
        Kd = sparsesdpvar(Options.Structure.Kd);
    else
        Kp = sdpvar(m,p,'full');
        Ki = sdpvar(m,p,'full');
        Kd = sdpvar(m,p,'full');
    end
    if isfield(Options,'Sign')
        constraints=[Options.Sign.Kp*Kp(:)>0
                     Options.Sign.Ki*Ki(:)>0
                     %Options.Sign.Kd*Kd(:)>0 
                     ];
    else
        constraints=[];
    end
    
    for icase=1:l % repeat for each uncetain plant case
        Pcase=P(:,:,icase);
        P0case=P0(:,:,icase);
        for k=1:N
            
            wk=w(k);
            s=1j*wk;
            
            Pk = freqresp(Pcase,wk);
            Ck = Kp+(1/s)*Ki+s/(1+tau*s)*Kd;
            Ck0 = Kp0+(1/s)*Ki0+s/(1+tau*s)*Kd0;
            
            LMI0=makeLMI(P0case*Ki,P0case*Ki0,t*eye(p));
            
            Z = eye(p)+Pk*Ck;
            Z0 = eye(p)+Pk*Ck0;
            Y1 = (1/Smax(k))*eye(p);
            LMI1 = makeLMI(Z,Z0,Y1);
            
            Y2 = (1/Tmax(k))*Pk*Ck;
            LMI2 = makeLMI(Z,Z0,Y2);
            
            Y3 = (1/Qmax(k))*Ck;
            LMI3 = makeLMI(Z,Z0,Y3);
            
            constraints = [constraints LMI0>0 LMI1>0 LMI2>0 LMI3>0];         
        end
    end
    
    %ops = sdpsettings ('solver','sdpt3');
    ops = sdpsettings ('solver','sedumi');
    %solvesdp(constraints,-t,ops);
    diagnostics = optimize(constraints,-t,ops)
    
    Kp0 = double(Kp);
    Ki0 = double(Ki);
    Kd0 = double(Kd);
    t = double(t);
    
    if abs(t-t0)<0.001
        break
    end
    t0=t;
end

Kp = Kp0;
Ki = Ki0;
Kd = Kd0;

s = tf('s');
C = Kp+(1/s)*Ki+s/(1+tau*s)*Kd;
objval = norm(inv(P0(:,:,1)*Ki)); % spectra norm (largest singular value)

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
if nnz(S)>0
    [I,J,~] = find(S);
    X=sparse(I,J,sdpvar(nnz(S),1));
else
    X=S;
end
%end of sparsesdpvar function
end

%end of AUTO_MIMO_PID main function
end 