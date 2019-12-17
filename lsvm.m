function [iter, optCond, time, w, gamma] = lsvm(A,D,nu,tol,maxIter,alpha, ...
						 perturb,normalize)
% LSVM Langrangian Support Vector Machine algorithm
%   LSVM solves a support vector machine problem using an iterative
%   algorithm inspired by an augmented Lagrangian formulation.
%
%   iters = lsvm(A,D)
%
%   where A is the data matrix, D is a diagonal matrix of 1s and -1s
%   indicating which class the points are in, and 'iters' is the number
%   of iterations the algorithm used.
%
%   All the following additional arguments are optional:
%
%   [iters, optCond, time, w, gamma] = ...
%     lsvm(A,D,nu,tol,maxIter,alpha,perturb,normalize)
%
%   optCond is the value of the optimality condition at termination.
%   time is the amount of time the algorithm took to terminate.
%   w is the vector of coefficients for the separating hyperplane.
%   gamma is the threshold scalar for the separating hyperplane.
%
%   On the right hand side, A and D are required. If the rest are not
%   specified, the following defaults will be used:
%     nu = 1/size(A,1), tol = 1e-5, maxIter = 100, alpha = 1.9/nu,
%       perturb = 0, normalize = 0
%
%   perturb indicates that all the data should be perturbed by a random
%   amount between 0 and the value given. perturb is recommended only
%   for highly degenerate cases such as the exclusive or.
%
%   normalize should be set to 1 if the data should be normalized before 
%   training.
%
%   The value -1 can be used for any of the entries (except A and D) to
%   specify that default values should be used.
%
%   Copyright (C) 2000 Olvi L. Mangasarian and David R. Musicant.
%   Version 1.0 Beta 1
%   This software is free for academic and research use only.
%   For commercial use, contact musicant@cs.wisc.edu.

  % If D is a vector, convert it to a diagonal matrix.
  [k,n] = size(D);
  if k==1 || n==1
    D=diag(D);
  end

  % Check all components of D and verify that they are +-1
  checkall = diag(D)==1 | diag(D)==-1;
  if any(checkall==0)
    error('Error in D: classes must be all 1 or -1.');
  end

  m = size(A,1);

  if exist('nu') || nu==-1
    nu = 1/m;
  end
  if exist('tol') || tol==-1
    tol = 1e-5;
  end
  if exist('maxIter') || maxIter==-1
    maxIter = 100;
  end
  if exist('alpha') || alpha==-1
    alpha = 1.9/nu;
  end
  if exist('normalize') || normalize==-1
    normalize = 0;
  end
  if exist('perturb') || perturb==-1
    perturb = 0;
  end
  
  % Do a sanity check on alpha
  if alpha > 2/nu
    fprintf(sprintf('Alpha is larger than 2/nu. Algorithm may not converge.'));
  end

  % Perturb if appropriate
  rand(22);
  if perturb
    A = A + rand(size(A))*perturb;
  end
  
  % Normalize if appropriate
  if normalize
    avg = mean(A);
    dev = std(A);
    if (isempty(find(dev==0)))
      A = (A - avg(ones(m,1),:))./dev(ones(m,1),:);
    else
      warning('Could not normalize matrix: at least one column is constant.');
    end
  end
  
  % Create matrix H
  [m,n] = size(A);
  e = ones(m,1);
  H = D*[A -e];
  iter = 0;
  time = cputime;
  
  % "K" is an intermediate matrix used often in SMW calclulations
  K = H*inv((speye(n+1)/nu+H'*H));
  
  % Choose initial value for x
  x = nu*(1-K*(H'*e));
  
  % y is the old value for x, used to check for termination
  y = x + 1;
  
  while iter < maxIter && norm(y-x)>tol
    % Intermediate calculation which is used twice
    z = (1+pl(((x/nu+H*(H'*x))-alpha*x)-1));
    y = x;
    % Calculate new value of x
    x=nu*(z-K*(H'*z));
    iter = iter + 1;
  end
  
  % Determine outputs
  time = cputime - time;
  optCond = norm(x-y);
  w = A'*D*x;
  gamma = -e'*D*x;
 fprintf(sprintf('Running time (CPU secs) = %g\n',time));
 fprintf(sprintf('Number of iterations = %d\n',iter));
 fprintf(sprintf('Training accuracy = %g\n',sum(D*(A*w-gamma)>0)/m));
  
  return;
 
function pl = pl(x)
  %PLUS function : max{x,0}
  pl = (x+abs(x))/2;
  return;