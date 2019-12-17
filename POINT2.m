clear all
close all
clc

load ('CMA_ES_CRPs.mat');

N = 10000/2;
tot = 10;
eps = rand(N,1);
n = 64;

response=1-2*response;

rtest = response(N+1:2*N);

reliability_labeled_data=1-2*reliability_labeled_data;



for ii = 1:10000
    for k = 1:n
        p(ii,k) = prod(challenge(ii,k:n));
    end
    p(ii,k+1) = 1;
end


xmin = cmaes2(n+2,p(1:N,:),reliability_labeled_data(1:N));
acc = size((find((rtest+sign(p(N+1:end,:)*xmin(1:n+1))~=0))),1)/N;

fprintf(sprintf('Accuracy = %g\n',acc));
