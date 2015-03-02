function [err, svmstruct] = doSVM(Xtrain,Ytrain,Xtest,Ytest)
%this function calculates the optimal seperating
%hyperplane on the training data (Xtrain, Ytrain)
%
%it then tests the data on Xtest, and returns 
%the misclassification error rate and svm structure

%add col of ones?
Xtrain = [ones(size(Xtrain,1),1) Xtrain];

%create options structure
opts = optimset('maxiter',50000);

%creat svm structure and train
svmstruct = svmtrain(Xtrain,Ytrain,'options',opts,'autoscale','false');

%classify test data
Yhat = svmclassify(svmstruct,[ones(size(Xtest),1) Xtest]);

%calculate misclassification rate
err = sum(Yhat~=Ytest)/numel(Ytest);