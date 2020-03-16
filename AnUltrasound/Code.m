clear all; clc;close all;
load Testdata.mat;

L = 15; %spatial domain
n = 64; %Fourier modes

%Define grid vectors, create 3D cartesian grid
x2 = linspace(-L,L,n+1); % time discretization
x = x2(1:n); % only use the first n points (periodicity)
y = x;
z = x;
k = (2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; % frequency components
ks = fftshift(k);

[X,Y,Z] = meshgrid(x,y,z);
[Kx,Ky,Kz] = meshgrid(ks, ks, ks);

U = zeros(n,n,n); %the variable that going to calculate the average of frequency

for j = 1:20
    Un(:,:,:) = reshape(Undata(j,:),n,n,n);
    U = U + fftn(Un);
end

U = abs(fftshift(U))/20; %shift and take absoluate value of average frequency

%Finding the center frequency (max value in U)
maxi = -inf;
for a = 1:64
    for b = 1:64
        for c = 1:64
            if U(a,b,c) > maxi
                maxi = U(a,b,c);
                x1 = a;
                y1 = b;
                z1 = c;
            end
        end
    end
end

%Gaussian Filter Function 
tau = 0.2;
gfunc = exp(-tau*((Kx-ks(y1)).^2 + (Ky-ks(x1)).^2 +(Kz-ks(z1)).^2)/2);

figure(1)      
isosurface(Kx,Ky,Kz,abs(U)/max(U(:)),0.7)
set(gca,'FontSize', 18);
axis([ks(1) -ks(1) ks(1) -ks(1) ks(1) -ks(1)]), grid on;
xlabel('Kx');ylabel('Ky');zlabel('Kz');

%Gaussian Filter picture
figure(2)
slice(Kx,Ky,Kz,gfunc,2,-1,0);
set(gca,'FontSize', 18);
axis([ks(1) -ks(1) ks(1) -ks(1) ks(1) -ks(1)]), grid on;
xlabel('Kx');ylabel('Ky');zlabel('Kz');

marble_location = zeros(20,3);
for j = 1:20
    Un(:,:,:) = reshape(Undata(j,:),n,n,n);
    Ut = fftshift(fftn(Un));
    NewUt = Ut.*gfunc;
    Unf = abs(ifftn(fftshift(NewUt)));
    
    %find the indices that has max value
    maxi = -inf;
    for a = 1:64
        for b = 1:64
            for c = 1:64
                if Unf(a,b,c) > maxi
                    maxi = Unf(a,b,c);
                    x1 = X(1,a,1);
                    y1 = Y(b,1,1);
                    z1 = Z(1,1,c);
                end
            end
        end
    end
    
    marble_location(j,:) = [y1,x1,z1];
    
   %isosurface of trajectory of marble
    figure(3)
    isosurface(X,Y,Z,abs(Unf)/max(abs(Unf(:))),0.8);
    set(gca,'FontSize', 18);
    axis([-L L -L L -L L]), grid on;
    xlabel('X');ylabel('Y');zlabel('Z');
end

final_location = marble_location(end,:);

%isosurface of final marble position
figure(4)
isosurface(X,Y,Z,abs(Unf)/max(abs(Unf(:))),0.8);
set(gca,'FontSize', 18);
axis([-L L -L L -L L]), grid on;
xlabel('X');ylabel('Y');zlabel('Z');

%plot3 of the trajectory of the marble
figure(5)
plot3(marble_location(:,1), marble_location(:,2), marble_location(:,3),'k--o','LineWidth', 5);
set(gca,'FontSize', 18);
axis([-L L -L L -L L]), grid on;
xlabel('X');ylabel('Y');zlabel('Z');



