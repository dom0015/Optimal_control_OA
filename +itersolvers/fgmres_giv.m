function [x,X,nit,err,resvec] = fgmres_giv(A,b,m,tol,maxit,P,x0)
%Flexible GMRES (FGMRES) with restart length m.
%Solves Ax = b to tolerance tol using preconditioner P and initial guess x0
%If P is a function handle, P(A) =approx= I.
%If P is a matrix, P\A =approx= I.
%Optional output X is the whole sequence of FGMRES iterates

%Written by Nick Alger 6/10/2014, released to public domain.
%Reference:
%Saad, "A flexible inner-outer preconditioned GMRES algorithm", SIAM 1993

%Parse input
N = size(b,1);
resvec=ones(1,maxit);

if (isa(A,'numeric'))
    if (size(A,1) == N && size(A,2) == N)
        Afct = @(u) A*u;
    else
        error(['Dimensions of A and b in Ax = b are not consistent: '...
            'A is ', num2str(size(A,1)),'-by-', num2str(size(A,2)),', and '...
            'b is ', num2str(size(b,1)),'-by-', num2str(size(b,2)),'.']);
    end
elseif (isa(A,'function_handle'))
    Afct = A;
else
    error('Forward operator needs to be a matrix or function handle.')
end

if (nargin < 3 || isempty(m))
    m = 50; %Default restart
end

if (nargin < 4 || isempty(tol))
    tol = 1e-6;    
end

if (nargin < 5 || isempty(maxit))
    maxit = min(N,200);
end

if (nargin < 6 || isempty(P))
    Pfct = @(u) u;
elseif (isa(P,'numeric'))
    Pfct = @(u) P\u;
elseif (isa(P,'function_handle'))
    Pfct = P;
else
    error('Preconditioner needs to be a matrix or function handle.')
end

if (nargin < 7 || isempty(x0))
    x0 = zeros(size(b));
end

%Run FGMRES
normb=norm(b);
cs(1:m) = zeros(m,1);
sn(1:m) = zeros(m,1);
e1    = zeros(N,1);
e1(1) = 1.0;
X = [];
nit = 1;
while (nit <= maxit)
    %disp('Starting new FGMRES cycle')
    H = zeros(m+1,m);
    r0 = b - Afct(x0);
    beta = norm(r0);
    V = zeros(N, m+1);
    V(:,1) = (1/beta)*r0;
    s=beta*e1;
    Z = zeros(N,m);
    for i=1:m
        Z(:,i) = Pfct(V(:,i));
        w = Afct(Z(:,i));
        for k=1:i
            H(k,i) = w'*V(:,k);
            w = w - H(k,i)*V(:,k);
        end
        H(i+1, i) = norm(w);
        V(:, i+1) = (1/H(i+1, i))*w;
        
        for k = 1:i-1,                              % apply Givens rotation
            temp     =  cs(k)*H(k,i) + sn(k)*H(k+1,i);
            H(k+1,i) = -sn(k)*H(k,i) + cs(k)*H(k+1,i);
            H(k,i)   = temp;
        end
        [cs(i),sn(i)] = rotmat( H(i,i), H(i+1,i) ); % form i-th rotation matrix
        temp   = cs(i)*s(i);                        % approximate residual norm
        s(i+1) = -sn(i)*s(i);
        s(i)   = temp;
        H(i,i) = cs(i)*H(i,i) + sn(i)*H(i+1,i);
        H(i+1,i) = 0.0;
        err  = abs(s(i+1)) / normb;
        
        resvec(nit)=err;
%         disp(['k=', num2str(nit),  ...
%             ', relative residual= ', num2str(err)])

        if ( err <= tol )
            disp(['FGMRES rel.tolerance ', ...
                num2str(err), ...
                ' iter. ', ...
                num2str(nit)])
            y = H(1:i,1:i) \ s(1:i); 
            x = x0 + Z(:, 1:i)*y;
            return
        end
        
%         if (nargout > 1)
%             X = [X,x];
%         end
        

        nit = nit + 1;
        if (nit > maxit)
            break;
        end
    end
    
    if ( err <= tol ), break, end
    y = H(1:m,1:m) \ s(1:m);
    x = x0 + Z(:, 1:m)*y;                        % update approximation
    r = (( b-Afct(x) ));                              % compute residual
    s(i+1) = norm(r);
    err = s(i+1) / normb;                        % check convergence
    if ( err <= tol ), break, end;
    
    x0 = x;
end