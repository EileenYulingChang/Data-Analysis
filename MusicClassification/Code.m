clc; close all; clear all;

Sgtvector1 = [];
a = 150;
for j = 1:15
    [y,Fs] = audioread(strcat('EdSheeran' , num2str(j) , '.wav'));
    %convert to single channel
    for i = 1:length(y)
        if (y(i,1) == 0 || y(i,2) == 0)
            monosong(i,1) = max(y(i,:));
        else 
            monosong(i,1) = (y(i,1) + y(i,2))/2;
        end
    end
    
    v1 = monosong'/2;
    tfinal = length(v1)/Fs; t = (1:length(v1))/Fs;
    t_step = 0:0.05:tfinal;
    for i = 1:length(t_step)
        g = exp(-a*(t-t_step(i)).^2);
        Sg1 = g.*v1; Sgt1 = fft(Sg1); Sgt1 = fftshift(Sgt1);
        Sgtvector1 = [Sgtvector1; abs(Sgt1)];
    end
    edsheeran(:,j) = sum(Sgtvector1);
    Sgtvector1 = [];
end

Sgtvector1 = [];
for j = 1:15
    [y,Fs] = audioread(strcat('Chain' , num2str(j) , '.wav'));
    %convert to single channel
    monosong = [];
    for i = 1:length(y)
        if (y(i,1) == 0 || y(i,2) == 0)
            monosong(i,1) = max(y(i,:));
        else 
            monosong(i,1) = (y(i,1) + y(i,2))/2;
        end
    end
    
    v1 = monosong'/2;
    tfinal = length(v1)/Fs; t = (1:length(v1))/Fs;
    t_step = 0:0.05:tfinal;
    for i = 1:length(t_step)
        g = exp(-a*(t-t_step(i)).^2);
        Sg1 = g.*v1; Sgt1 = fft(Sg1); Sgt1 = fftshift(Sgt1);
        Sgtvector1 = [Sgtvector1; abs(Sgt1)];
    end 
    chainsmoker(:,j) = sum(Sgtvector1);
    Sgtvector1 = [];
end


Sgtvector1 = [];
for j = 1:15
    [y,Fs] = audioread(strcat('Charlie' , num2str(j) , '.wav'));
    if length(y(1,:)) == 1
        monosong = y;
    else
        %convert to single channel
        monosong = [];
        for i = 1:length(y)
            if (y(i,1) == 0 || y(i,2) == 0)
                monosong(i,1) = max(y(i,:));
            else 
                monosong(i,1) = (y(i,1) + y(i,2))/2;
            end
        end
    end
    
    v1 = monosong'/2;
    tfinal = length(v1)/Fs; t = (1:length(v1))/Fs;
    t_step = 0:0.05:tfinal;
    for i = 1:length(t_step)
        g = exp(-a*(t-t_step(i)).^2);
        Sg1 = g.*v1; Sgt1 = fft(Sg1); Sgt1 = fftshift(Sgt1);
        Sgtvector1 = [Sgtvector1; abs(Sgt1)];
    end 
    puth(:,j) = sum(Sgtvector1);
    Sgtvector1 = [];
end


traindata = [edsheeran(:,1:10), chainsmoker(:,1:10), puth(:,1:10)];
testdata = [edsheeran(:,11:15), chainsmoker(:,11:15), puth(:,11:15)];

[U,S,V] = svd(traindata, 0);

ne=length(edsheeran(1,1:10)); nc=length(chainsmoker(1,1:10)); np=length(puth(1,1:10));
music = S*V';
suc = [];
succh = [];
succhar = [];
suced = [];
for feature = 1:9
    U1 = U(:, 1:feature);

    ed = music(1:feature, 1:ne);
    chain = music(1:feature, ne+1:ne+nc);
    char = music(1:feature, ne+nc+1:ne+nc+np);

    global_mean = mean(music(1:feature,:),2);

    med = mean(ed,2);
    mchain = mean(chain,2);
    mchar = mean(char,2);

    Sw = 0;
    for k = 1:ne
        Sw = Sw + (ed(:,k)-med)*(ed(:,k)-med)';
    end

    for k = 1:nc
        Sw = Sw + (chain(:,k)-mchain)*(chain(:,k)-mchain)';
    end

    for k = 1:np
        Sw = Sw + (char(:,k)-mchar)*(char(:,k)-mchar)';
    end

    Sb = 1/3*(med-global_mean)*(med-global_mean)';
    Sb = Sb + 1/3*(mchain-global_mean)*(mchain-global_mean)';
    Sb = Sb + 1/3*(mchar-global_mean)*(mchar-global_mean)';

    [V2, D] = eig(Sb,Sw);
    [~,ind] = max(abs(diag(D)));
    w = V2(:,ind); w = w/norm(w,2);

    ved = w'*ed;
    vchain = w'*chain;
    vchar = w'*char;

    if mean(ved) > mean(vchain)
        w = -w;
        ved = -ved;
        vchain = -vchain;
    end

    if mean(vchain) > mean(vchar)
        w = -w;
        vchain = -vchain;
        vchar = -vchar;
    end

    sorted = sort(ved);
    sortchain = sort(vchain);
    sortchar = sort(vchar);

    t1 = length(sortchain);
    t2 = 1;

    while sortchain(t1) > sortchar(t2)
        t1 = t1 - 1;
        t2 = t2 + 1;
        if (t1 == 0) || (t2 == 0)
            break
        end
    end
    
     if (t1 == 0) || (t2 == 0)
         continue
     end
    
    threshold1 = (sortchain(t1) + sortchar(t2))/2;


    t3 = length(sorted);
    t4 = 1;
    
    while sorted(t3) > sortchain(t4)
        t3 = t3 - 1;
        t4 = t4 + 1;
        if (t3 == 0) || (t4 == 0)
            break
        end
    end
    
    if (t3 == 0) || (t4 == 0)
        continue
    end
    
    threshold2 = (sorted(t3) + sortchain(t4))/2;
    
    if feature == 6
        figure(2)
        subplot(2,2,1)
        hist(sortchar); hold on; plot([threshold1 threshold1],[0 10],'r')
        set(gca,'Ylim',[0 10],'Fontsize',[14])
        title('Charlie Puth'); hold on;
        subplot(2,2,2)
        hist(sortchain); hold on; plot([threshold1 threshold1],[0 10],'r')
        set(gca,'Ylim',[0 10],'Fontsize',[14])
        title('Chainsmokers'); hold on;
        subplot(2,2,3)
        hist(sorted); hold on; plot([threshold2 threshold2],[0 10],'r')
        set(gca,'Ylim',[0 10],'Fontsize',[14])
        title('Ed Sheeran')
        subplot(2,2,4)
        hist(sortchain); hold on; plot([threshold2 threshold2],[0 10],'r')
        set(gca,'Ylim',[0 10],'Fontsize',[14])
        title('Chainsmokers')
    end

    testmat = U1'*testdata;
    pval = w'*testmat;

    hiddenlabels = [ones(5,1);ones(5,1)*2; ones(5,1)*3]';

    ResVec = [];
    for i = 1:length(pval)
        if pval(i) > threshold1
            ResVec(i) = 3;
        elseif pval(i) < threshold2
            ResVec(i) = 1;
        else
            ResVec(i) = 2;
        end
    end
    errNum = 0;
    for i = 1:length(ResVec)
        if ResVec(i) ~= hiddenlabels(i)
            errNum = errNum + 1;
        end
    end
    sucRate = 1-abs(errNum)/15
    suc(1,feature) = sucRate;

    errNum = 0;
    for i = 1:5
        if ResVec(i) ~= hiddenlabels(i)
            errNum = errNum + 1;
        end
    end
    sucRate = 1-abs(errNum)/5
    suced(1,feature) = sucRate;

    errNum = 0;
    for i = 6:10
        if ResVec(i) ~= hiddenlabels(i)
            errNum = errNum + 1;
        end
    end

    sucRate = 1-abs(errNum)/5
    succh(1,feature) = sucRate;

    errNum = 0;
    for i = 11:15
        if ResVec(i) ~= hiddenlabels(i)
            errNum = errNum + 1;
        end
    end
    
    sucRate = 1-abs(errNum)/5
    succhar(1,feature) = sucRate;
end

figure(3)
subplot(2,2,1)
plot(suc,'-o');
xlabel('Number of Feature'); ylabel('Successful Rate');
title('Overall');
set(gca, 'Ylim', [0 1]);
subplot(2,2,2)
plot(suced,'-o');
xlabel('Number of Feature'); ylabel('Successful Rate');
title('Ed Sheeran');
set(gca, 'Ylim', [0 1]);
subplot(2,2,3)
plot(succh,'-o');
xlabel('Number of Feature'); ylabel('Successful Rate');
title('Chainsmokers');
set(gca, 'Ylim', [0 1]);
subplot(2,2,4)
plot(succhar,'-o');
xlabel('Number of Feature'); ylabel('Successful Rate');
title('Charlie Puth');
set(gca, 'Ylim', [0 1]);


