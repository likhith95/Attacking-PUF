function xmin=cmaes2(N,in,lbl)  
  % (mu/mu_w, lambda)-CMA-ES 
  % CMA-ES: Evolution Strategy with Covariance Matrix Adaptation
  % for nonlinear function minimization. 
  %
  % This code is "an excerpt" from cmaes.m and implements the key 
  % parts of the algorithm. It is intendend to be used for READING
  % and UNDERSTANDING the basic flow and all details of the CMA-ES
  % *algorithm*. To run "serious" simulations better use the cmaes.m 
  % code: it is longer, but offers restarts, far better termination 
  % options, and, in particular, supposedly quite useful output.
  %
  % Author: Nikolaus Hansen, 2003-09. 
  % e-mail: hansen[at]lri.fr
  %
  % License: This code is released into the public domain (that is, 
  %   you may use and modify it however you like). 
  %
  % URL: http://www.lri.fr/~hansen/purecmaes.m
  % References: See end of file. Last change: April, 29, 2014

  % --------------------  Initialization --------------------------------  
  % User defined input parameters (need to be edited)
  strfitnessfct = 'frosenbrock';  % name of objective/fitness function
  % N = 11;               % number of objective variables/problem dimension
  xmean = rand(N,1);    % objective variables initial point
  sigma = 0.5;          % coordinate wise standard deviation (step size)
  stopfitness = 1e-10;  % stop if fitness < stopfitness (minimization)
  stopeval = 10000;%1e3*N^2;   % stop after stopeval number of function evaluations
  
  % Strategy parameter setting: Selection  
  lambda = 4+floor(3*log(N));  % population size, offspring number
  mu = lambda/2;               % number of parents/points for recombination
  weights = log(mu+1/2)-log(1:mu)'; % muXone array for weighted recombination
  mu = floor(mu);        
  weights = weights/sum(weights);     % normalize recombination weights array
  mueff=sum(weights)^2/sum(weights.^2); % variance-effectiveness of sum w_i x_i

  % Strategy parameter setting: Adaptation
  cc = (4 + mueff/N) / (N+4 + 2*mueff/N); % time constant for cumulation for C
  cs = (mueff+2) / (N+mueff+5);  % t-const for cumulation for sigma control
  c1 = 2 / ((N+1.3)^2+mueff);    % learning rate for rank-one update of C
  cmu = min(1-c1, 2 * (mueff-2+1/mueff) / ((N+2)^2+mueff));  % and for rank-mu update
  damps = 1 + 2*max(0, sqrt((mueff-1)/(N+1))-1) + cs; % damping for sigma 
                                                      % usually close to 1
  % Initialize dynamic (internal) strategy parameters and constants
  pc = zeros(N,1); ps = zeros(N,1);   % evolution paths for C and sigma
  B = eye(N,N);                       % B defines the coordinate system
  D = ones(N,1);                      % diagonal D defines the scaling
  C = B * diag(D.^2) * B';            % covariance matrix C
  invsqrtC = B * diag(D.^-1) * B';    % C^-1/2 
  eigeneval = 0;                      % track update of B and D
  chiN=N^0.5*(1-1/(4*N)+1/(21*N^2));  % expectation of 
                                      %   ||N(0,I)|| == norm(randn(N,1))
  out.dat = []; out.datx = [];  % for plotting output
%   j = 1;

  % -------------------- Generation Loop --------------------------------
  counteval = 0;  % the next 40 lines contain the 20 lines of interesting code 
  while counteval < stopeval
    % Generate and evaluate lambda offspring
    for k=1:lambda
      arx(:,k) = xmean + sigma * B * (D .* randn(N,1)); % m + sig * Normal(0,C) 
      arfitness(k) = feval(strfitnessfct, arx(:,k), in, lbl);%, in', lbl]); % objective function call
%       arfitness(k) = feval(strfitnessfct, [arx(:,k)]);%, in']);%; lbl(k)]);%, lbl);
%       e(k) = feval(feval(strfitnessfct,arx(:,k)),in);
      counteval = counteval+1;
%       j = j+1;
    end
    
    % Sort by fitness and compute weighted mean into xmean
    [arfitness, arindex] = sort(arfitness);  % minimization
    xold = xmean;
    xmean = arx(:,arindex(1:mu)) * weights;  % recombination, new mean value
    
    % Cumulation: Update evolution paths
    ps = (1-cs) * ps ... 
          + sqrt(cs*(2-cs)*mueff) * invsqrtC * (xmean-xold) / sigma; 
    hsig = sum(ps.^2)/(1-(1-cs)^(2*counteval/lambda))/N < 2 + 4/(N+1);
    pc = (1-cc) * pc ...
          + hsig * sqrt(cc*(2-cc)*mueff) * (xmean-xold) / sigma; 

    % Adapt covariance matrix C
    artmp = (1/sigma) * (arx(:,arindex(1:mu)) - repmat(xold,1,mu));  % mu difference vectors
    C = (1-c1-cmu) * C ...                   % regard old matrix  
         + c1 * (pc * pc' ...                % plus rank one update
                 + (1-hsig) * cc*(2-cc) * C) ... % minor correction if hsig==0
         + cmu * artmp * diag(weights) * artmp'; % plus rank mu update 

    % Adapt step size sigma
    sigma = sigma * exp((cs/damps)*(norm(ps)/chiN - 1)); 
    
    % Update B and D from C
    if counteval - eigeneval > lambda/(c1+cmu)/N/10  % to achieve O(N^2)
      eigeneval = counteval;
      C = triu(C) + triu(C,1)'; % enforce symmetry
      [B,D] = eig(C);           % eigen decomposition, B==normalized eigenvectors
      D = sqrt(diag(D));        % D contains standard deviations now
      invsqrtC = B * diag(D.^-1) * B';
    end
    
    % Break, if fitness is good enough or condition exceeds 1e14, better termination methods are advisable 
    if arfitness(1) <= stopfitness || max(D) > 1e7 * min(D)
      break;
    end

    % Output 
    more off;  % turn pagination off in Octave
    disp([num2str(counteval) ': ' num2str(arfitness(1)) ' ' ... 
          num2str(sigma*sqrt(max(diag(C)))) ' ' ...
          num2str(max(D) / min(D))]);
    % with long runs, the next line becomes time consuming
    out.dat = [out.dat; arfitness(1) sigma 1e5*D' ]; 
    out.datx = [out.datx; xmean'];
  end % while, end generation loop

  % ------------- Final Message and Plotting Figures --------------------
  disp([num2str(counteval) ': ' num2str(arfitness(1))]);
  xmin = arx(:, arindex(1)); % Return best point of last iteration.
                             % Notice that xmean is expected to be even
                             % better.
  figure(1); hold off; semilogy(abs(out.dat)); hold on;  % abs for negative fitness
  semilogy(out.dat(:,1) - min(out.dat(:,1)), 'k-');  % difference to best ever fitness, zero is not displayed
  title('fitness, sigma, sqrt(eigenvalues)'); grid on; xlabel('iteration');  
  figure(2); hold off; plot(out.datx);
  title('Distribution Mean'); grid on; xlabel('iteration')
  
% ---------------------------------------------------------------  
function f=frosenbrock(x,A_phi,D_rlb)

epsilon = abs(x(66));

t= sign(abs(A_phi*x(1:65)) - epsilon);
m=corrcoef(t,D_rlb);
f= 0.5-m(2,1)/2;
 
% % ---------------------------------------------------------------  
% %%% REFERENCES
% %
% % Hansen, N. and S. Kern (2004). Evaluating the CMA Evolution
% % Strategy on Multimodal Test Functions.  Eighth International
% % Conference on Parallel Problem Solving from Nature PPSN VIII,
% % Proceedings, pp. 282-291, Berlin: Springer. 
% % (http://www.bionik.tu-berlin.de/user/niko/ppsn2004hansenkern.pdf)
% % 
% % Further references:
% % Hansen, N. and A. Ostermeier (2001). Completely Derandomized
% % Self-Adaptation in Evolution Strategies. Evolutionary Computation,
% % 9(2), pp. 159-195.
% % (http://www.bionik.tu-berlin.de/user/niko/cmaartic.pdf).
% %
% % Hansen, N., S.D. Mueller and P. Koumoutsakos (2003). Reducing the
% % Time Complexity of the Derandomized Evolution Strategy with
% % Covariance Matrix Adaptation (CMA-ES). Evolutionary Computation,
% % 11(1).  (http://mitpress.mit.edu/journals/pdf/evco_11_1_1_0.pdf).
% %
