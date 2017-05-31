% =========================================================================
% (c) Lars Ruthotto 2013/10/11
% http://www.eos.ubc.ca/about/researcher/L.Ruthotto.html
%
% function [xc,iter] = ActiveSetMethod(G,A,c,b,varargin)
%
% Active set method for convex quadratic programs
%
%
%       minimize_x      q(x) := .5*x'*G*x + x'*c
%
%       subject to      Ax >= b
%
% Input:
%  var1   - description
%
% Output:
%  var2   - description
%
% see also
% =========================================================================

function [xc,iter] = ActiveSetMethod(G,A,c,b,varargin)

if nargin==0,
    help(mfilename); runMinimalExample; return;
end

maxIter = 1000;
tolP    = 1e-14;
x0 = [];
for k=1:2:length(varargin),     % overwrites default parameter
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

if isempty(x0),
    error('Provide valid starting point!');
end

numC = size(A,1);
out = 0;
actC = find(A*x0<=b);
if not(isempty(actC)),
    Wc = actC(1);
else
    Wc = [];
end
xc = x0;

q = @(x) .5*x'*G*x + x'*c + .5*c'*c;

iter = 0;
while 1,
    if out,
        fprintf('i=%d, q(x)=%1.4f, \tnumel(Wc)=%d\n',iter,q(xc),numel(Wc));
    end
    % solve QP
    Aac = A(Wc,:);
    gc  = G*xc+c;
    pc  = solveSchur(G,Aac,gc,zeros(numel(Wc,1)));
    
    if norm(pc)<tolP,
        % compute Lagrange multipliers
        lam = Aac' \ gc;
        if all(lam>=0),
            break;
        else
            [mini,idm] = min(lam,[],1);
%             fprintf('remove  constraint %d\n',Wc(idm));
            Wc(idm)=[]; 
        end
    else
        % compute step length
        ik = (1:numC)'; ik(Wc)=[];
        Aout = A(ik,:);
        t2 = b(ik)-Aout*xc;  
        t1 = Aout*pc;  id = find(t1<0);
        ik = ik(id);
        t2 = t2(id);
        t1 = t1(id);
        
        alpha = min(1,min(t2./t1));
        
        xc = xc + alpha*pc;
        if alpha<1
            idm = ik(find(alpha==t2./t1,1,'first'));
            % add one blocking constraint
            Wc = [Wc;idm];
%             fprintf('add blocking constraint %d, alpha=%f\n',idm,alpha);
        end
    end
    
    
    
    
    % check stopping
    STOP = (iter>maxIter);
    if STOP, break; end;
    iter = iter+1;
end


function [x,lam] = solveSchur(G,A,c,b,x0)
%
% Solve KKT system by Schur complement method
%
%       | G   A^T | |  -p  |   | g |
% KKT = |         | |      | = |   |
%       | A    0  | | lam* |   | h |
%
% h = A*x0-b, g = c+G*x, p=x*-x
%
% now lam*  = (A*inv(G)*A')\(A*G\g -h)
% and

if not(exist('x0','var')) || isempty(x0), x0 = 0*c; end;

g   = c + G*x0;
h   = A*x0-b;

lam = (A*(G\A'))\(A*(G\g)-h);
p   = G\ (A'*lam) - G\g;
x   = p+x0;


function runMinimalExample
G = 2*eye(2);
c = [-2; -5];
A = [1 -2; -1 -2; -1 2; 1 0 ;0  1];
b = [-2;-6;-2;0;0];

x0 = [2;0];

x = feval(mfilename,G,A,c,b,'x0',x0);

fprintf('x=[%f,%f]^T, q(x)=%f\n',x, .5*x'*G*x+x'*c + .5*c'*c)
