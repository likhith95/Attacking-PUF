clear all
close all
clc
load ('single_CRPs.mat')
response=Matrix_update(response)
[D,P]=Matrix_update2(challenge,10000)
lsvm(P,response(1:10000,:),-1,-1,-1,-1,-1,-1)