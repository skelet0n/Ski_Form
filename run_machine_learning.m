%wrapper for all methods
%
%

%load data
load images_stef_scans;
load images_scans;

Ystef_scans = ~Ystef_scans;
XX = [XXscans;XXstef_scans];
Y = [Yscans;Ystef_scans];

%rescale XX to have values between 0 and 1
XX = 1 - XX/255;
%XX(i,j) = 1 => pixel (i,j) is BLACK

%CODING:
%    Y(i) = 1 => image i has an X
%    Y(i) = 0 => image i DOES NOT have an X

%THREE METHODS:
%  1) RIDGE RIDGE REGRESSION, need CV to estimate \lambda
%
%  2) SEPERATING HYPERPLANES
%
%  3) FEED FORWARD NEURAL NET, need to train several times (local mins)


%% 1) RIDGE REGRESSION WITH CROSS VALIDATION (10-FOLD)
lambda_vals = 5:5;
global lambda

cv_err = zeros(numel(lambda_vals),1);
cv_std = zeros(numel(lambda_vals),1);

%loop through parameter values
for i =1:numel(lambda_vals)
    %set new parameter value
    lambda = lambda_vals(i);
    
    %run ridge regression with cv
    vals = crossval(@doRidge,XX,Y);
    
    %store mean and standard deviation
    cv_err(i) = mean(vals);
    cv_std(i) = std(vals);
end

figure(1);
errorbar(lambda_vals,cv_err,cv_std,'o');

%% 2) SEPERATING HYPERPLANES

%no parameters to estimate here,

%split data into 10 chunks and train and test on each
Ysvm = Y; Ysvm(Y==0) = -1;
vals = crossval(@doSVM,XX,Ysvm);

svm_err = mean(vals);


[err,svmstruct] = doSVM(XX,Ysvm,XX,Ysvm);
%calculate seperating hyperplane
alpha = svmstruct.Alpha;
sv = svmstruct.SupportVectors;
sv_ind = svmstruct.SupportVectorIndices;

w = zeros(1,256);
for i = 1:numel(alpha)
   w = w + alpha(i)*Ysvm(sv_ind(i))*XX(sv_ind(i),:); 
end



%% 3) FEED FORWARD NEURAL NETWORK

%CV on number of hidden layers?
hidden_unit_vals = 1:10;

%declare global variable for the # of hidden units
global hidden_units

cv_nn_err = zeros(numel(hidden_unit_vals),1);
cv_nn_std = zeros(numel(hidden_unit_vals),1);

%training Neural Net using nntrain starts at random weights
%=> solution is sensitive to initial guess and consequently
% performance results need to be average

%pick number of times to train neural network 
num_runs = 1;

%outer loop of hidden units
for i = 1:numel(hidden_unit_vals);
    hidden_units = hidden_unit_vals(i);
    
    %holder for average of NN performance
    vals_outer = zeros(10,1);
    
    %inner loop for averaging performance of neural net
    for j = 1:num_runs
        vals_inner = crossval(@doNN,XX,Y);
        
        %calculate average performance as you go
        vals_outer = vals_outer + vals_inner/num_runs;
    end
    
    %get mean and std of misclassifications
    cv_nn_err(i) = mean(vals_outer);
    cv_nn_std(i) = std(vals_outer);
    
    disp(i);
end

figure(2);
errorbar(hidden_unit_vals,cv_nn_err,cv_nn_std,'o');
