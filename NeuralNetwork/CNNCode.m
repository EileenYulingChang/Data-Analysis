%% MNIST Classifier with convolutional neural network
clear; close all; clc
load fashion_mnist.mat
X_train = im2double(X_train);
X_test = im2double(X_test);
X_train = reshape(X_train,[60000 28 28 1]);
X_train = permute(X_train,[2 3 4 1]);
X_test = reshape(X_test,[10000 28 28 1]);
X_test = permute(X_test,[2 3 4 1]);
X_valid = X_train(:,:,:,1:5000);
X_train = X_train(:,:,:,5001:end);
y_valid = categorical(y_train(1:5000))';
y_train = categorical(y_train(5001:end))';
y_test = categorical(y_test)';
layers = [    
    imageInputLayer([28 28 1],"Name","imageinput")    
    convolution2dLayer([5 5],32,"Name","conv_1","Padding","same")    
    reluLayer("Name","tanh_1")    
    averagePooling2dLayer([2 2],"Name","avgpool2d_1","Padding","same","Stride",2)    
    convolution2dLayer([5 5],64,"Name","conv_2")    
    reluLayer("Name","tanh_3")    
    averagePooling2dLayer([2 2],"Name","avgpool2d_2","Padding","same","Stride",2)    
    convolution2dLayer([5 5],128,"Name","conv_3")    
    reluLayer("Name","tanh_2")   
    averagePooling2dLayer([2 2],"Name","avgpool2d_3","Padding","same","Stride",2) 
    fullyConnectedLayer(128,"Name","fc_1")    
    reluLayer("Name","tanh_4")    
    fullyConnectedLayer(10,"Name","fc_2")    
    softmaxLayer("Name","softmax")    
    classificationLayer("Name","classoutput")];
options = trainingOptions('sgdm', ...    
    'MaxEpochs',15,...    
    'InitialLearnRate',9e-2, ...    
    'L2Regularization',5e-4, ...    
    'ValidationData',{X_valid,y_valid}, ...    
    'Verbose',false, ...
    'ValidationFrequency', 200, ...
    'Plots','training-progress');
net = trainNetwork(X_train,y_train,layers,options);
%% Confusion for training
figure(1)
y_pred = classify(net,X_train);
plotconfusion(y_train,y_pred)
%% Test classifier
figure(2)
y_pred = classify(net,X_test);
plotconfusion(y_test,y_pred)