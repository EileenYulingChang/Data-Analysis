clear all; clc; close all;

%loading videos
load('cam1_1.mat');
load('cam2_1.mat');
load('cam3_1.mat');

%find numbers of frame
numFrames1 = size(vidFrames1_1,4);
numFrames2 = size(vidFrames2_1,4);
numFrames3 = size(vidFrames3_1,4);
maxFrames = max([numFrames1 numFrames2 numFrames3]);

%find size of frame
[m1,n1] = size(vidFrames1_1(:,:,1,1));
[m2,n2] = size(vidFrames2_1(:,:,1,1));
[m3,n3] = size(vidFrames3_1(:,:,1,1));

for k = 1:maxFrames
    if k <= numFrames1
        mov1(k).cdata = vidFrames1_1(:,:,:,k);
        mov1(k).colormap = [];
    end
    if k <= numFrames2
        mov2(k).cdata = vidFrames2_1(:,:,:,k);
        mov2(k).colormap = [];
    end
    if k <= numFrames3
        mov3(k).cdata = vidFrames3_1(:,:,:,k);
        mov3(k).colormap = [];
    end    
end

X1=[];X2=[];X3=[];Y1=[];Y2=[];Y3=[];

for i = 1:maxFrames
    if i <= numFrames1
        abw = rgb2gray(frame2im(mov1(i)));
        abw(:,1:320) = 0; 
        abw(:,380:end) = 0;
        abw(1:200,:) = 0;
        [Max, Ind] = max(abw(:));
        [y1 x1] = ind2sub(size(abw), Ind);
        X1 = [X1 x1];
        Y1 = [Y1 y1];
    end
    if i <= numFrames2
        abw = rgb2gray(frame2im(mov2(i)));
        abw(:,1:260) = 0; 
        abw(:,330:end) = 0;
        [Max, Ind] = max(abw(:));
        [y2 x2] = ind2sub(size(abw), Ind);
        X2 = [X2 x2];
        Y2 = [Y2 y2];
    end
    if i <= numFrames3
        abw = rgb2gray(frame2im(mov3(i)));
        abw(1:250,:) = 0; 
        abw(310:end,:) = 0;
        abw(:, 1:260) = 0;
        [Max, Ind] = max(abw(:));
        [y3 x3] = ind2sub(size(abw), Ind);
        X3 = [X3 x3];
        Y3 = [Y3 y3];
    end
end

figure(1)
subplot(3,2,1)
plot(X1);xlabel('Frame');ylabel('Position in x');ylim([0 500]);title('Camera 1');

subplot(3,2,2)
plot(Y1);xlabel('Frame');ylabel('Position in y');ylim([0 500]);title('Camera 1');

subplot(3,2,3)
plot(X2);xlabel('Frame');ylabel('Position in x');ylim([0 500]);title('Camera 2');

subplot(3,2,4)
plot(Y2);xlabel('Frame');ylabel('Position in y');ylim([0 500]);title('Camera 2');

subplot(3,2,5)
plot(X3);xlabel('Frame');ylabel('Position in x');ylim([0 500]);title('Camera 3');

subplot(3,2,6)
plot(Y3);xlabel('Frame');ylabel('Position in y');ylim([0 500]);title('Camera 3');

[Min1 I1]=min(X1(1:50)); X1=X1(I1:I1+200); Y1=Y1(I1:I1+200);
[Min2 I2]=min(X2(1:50)); X2=X2(I2:I2+200); Y2=Y2(I2:I2+200);
[Min3 I3]=min(X3(1:50)); Y3=Y3(I3:I3+200); X3=X3(I3:I3+200);


A = [X1;Y1;X2;Y2;X3;Y3];
[m,n]=size(A); 
mn=mean(A,2); 
A=A-repmat(mn,1,n); 

[u,s,v] = svd(A);

Y = u'*A;

figure(2)
plot(1:6, (diag(s).^2)/sum(diag(s).^2), 'ko--', 'Linewidth', 2);
title("Case 1: Energy of each Diagonal Variance");
xlabel("Diagonal Variances"); ylabel("Energy Captured");



