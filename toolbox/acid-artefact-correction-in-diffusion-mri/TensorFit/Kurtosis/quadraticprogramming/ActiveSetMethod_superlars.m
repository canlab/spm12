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

function [xc,iter] = ActiveSetMethod_superlars(G,A,c,b,varargin)

if nargin==0,
    help(mfilename);load dummyAS; ActiveSetMethod_superlars(G,A,c,b,'x0',0*c); return;
end
% save dummyAS;

maxIter = 100;
tolP    = 1e-14;
x0      = [];
for k=1:2:length(varargin),     % overwrites default parameter
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

if isempty(x0) || any(A*x0<b),
    error('Provide valid starting point!');
end
out = 0;
actC = find(A*x0<b);
if not(isempty(actC)),
    Wc = actC(1);
else
    Wc = [];
end
xc = x0;

q = @(x) .5*x'*G*x + x'*c + .5*c'*c;
sdiag = @(v) spdiags(v(:),0,numel(v), numel(v));
% normalize constraints
nA = full(sqrt(sum(A.^2,2)));
A  = sdiag(1./nA)*A;
b  = sdiag(1./nA)*b;
    

iter = 0;
while 1,
    if out,
        fprintf('i=%d, q(x)=%1.4e, \tnumel(Wc)=%d\n',iter,q(xc),numel(Wc));
    end
    % gradient
    gc  = G*xc+c;
    
    %update working set
    Aac = A(Wc,:); 
    
    % using nullspace method to solve KKT system
    %
    % determine null space of constraints
    [U,S,V] = svd(Aac);
    Z = V(:,numel(Wc)+1:end);
    Y = V(:,1:numel(Wc));
    
    % project into nullspace of constraints
    if isempty(Wc),
        h = zeros(numel(Wc),1);
    else
        h = -b(Wc);
    end
    py = (Aac*Y)\-h;
    pz = (Z'*G*Z)\(-Z'*G*Y*py - Z'*gc);
    pc  = Y*py + Z * pz;
    lam =  (Aac*Y)\ (Y'*(gc+G*pc));
        
    if norm(pc)<tolP,
        % compute Lagrange multipliers
        if all(lam>=0),
            break;
        else
            [mini,idm] = min(lam,[],1);
            Wc(idm)=[]; 
        end
    else
        % compute step length
        t2 = min(0,(b-A*xc));
        t1 = (A*pc);
        t  = t2./t1;
        t(t1>=0) = 1;
        t(Wc)    = 1;
        alpha = min(1,min(t));
        if alpha<0,
            error('alpha negative');
        end
        
        xc = xc + alpha*pc;
        if alpha<1
            idm = find(alpha==t);
            At = A(alpha==t,:); % blocking constraints
            pp = Z*Z'*At';
            [mx,ii] = max(sum(pp.^2,1));
            if mx>1e-3,
                Wc = [Wc;idm(ii)];
            elseif all(lam)>0,
                break;
            else
                warning('cannot find linear independent constraints');
                break;
            end
            
        end
    end
    % check stopping
    STOP = (iter>maxIter);
    if STOP, break; end;
    iter = iter+1;
end