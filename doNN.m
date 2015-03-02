function [err, netff] = doNN(Xtrain,Ytrain,Xtest,Ytest)
%this function trains neural net on (Xtrain Ytrain)
% and tests on Xtest and calculates the missclassification error

%get number of hidden units to use
global hidden_units

%create neural network object with said number of hidden units
netff = feedforwardnet(hidden_units);

%Neural Network training tools require classification differently
T = zeros(2,size(Ytrain,1));
T(1,Ytrain==1) = 1;
T(2,Ytrain==0) = 1;
%XX is N x p => T is 2 x N, where:
%
%T(:,i) = [1 0]' <=> Y(i) = 1
%T(:,i) = [0 1]' <=> Y(i) = 0%


%don't pop up toolbox
netff.trainparam.showWindow = 0;

%train
netff = train(netff,Xtrain',T);

%get classification
classes_vec = netff(Xtest');

%change from [1 0], [0 1] format to just 0 or 1
classes = ~(vec2ind(classes_vec) - 1);

%get number of misclassifications
err = sum(classes'~=Ytest)/numel(Ytest);
