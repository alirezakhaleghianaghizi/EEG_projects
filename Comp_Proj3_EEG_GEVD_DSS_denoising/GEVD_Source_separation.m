clc;
clear all;
cd('C:\Users\ASUS\OneDrive\Desktop\uni\simisterm7\eeg\Comp_HW3');
%%
Q1=load('Q1.mat');
X_org=Q1.X_org;
X1=Q1.X1;
X2=Q1.X2;
X3=Q1.X3;
X4=Q1.X4;
T1=Q1.T1;
T2=Q1.T2;

feq=100;
lab=[];
offset = max(abs(X_org(:))) ;
disp_eeg(X_org,offset,feq,lab);
title('X_org')

offset = max(abs(X1(:))) ;
disp_eeg(X1,offset,feq,lab);
title('X1')

offset = max(abs(X2(:))) ;
disp_eeg(X2,offset,feq,lab);
title('X2')
%%
offset = max((X3(:))) ;
disp_eeg((X3),offset,feq,lab);
title('X3')
%%
offset = max(abs(X4(:))) ;
disp_eeg(X4,offset,feq,lab);
title('X4')
%% alef GEVD triangolar
fs=100;
maxshift=4*fs;
L=length(X_org(1,:));
C_t=X_org(:,1:L-maxshift)*X_org(:,maxshift+1:L)';
C_t=(C_t+C_t')/2;
C_all=X_org*X_org';
[V1,D1] = eig(C_t,C_all);
[sVals,sIndex] = sort(diag(D1),'descend') ;
sV1 = V1(:,sIndex) ;
Z1=sV1'*X_org;
offset = max(abs(Z1(:))) ;
feq=100;
disp_eeg(Z1,offset,feq,[]);
title('source\_den2')
lenghtz=length(Z1(1,:));
Z1(2:8,:)=zeros(7,lenghtz);
X_den_GEVD=pinv(sV1')*Z1;
offset = max(abs(X_den_GEVD(:))) ;
% ElecName = Electrodes.labels ;
disp_eeg(X_den_GEVD,offset,feq,[]);
title('X\_den2')
PRMSE2=sqrt(sum(sum(pow2(X1-X_den_GEVD))))/sqrt(sum(sum(pow2(X1))));
%%  b not known the period
shift=0;
[c,lags] = xcorr(X_org(1,:));
[maxAutoCorr, ind1] = max(c(10300:10700));
shift=ind1;
for i=2:1:8
    [c,lags] = xcorr(X_org(i,:));
    [maxAutoCorr, ind1] = max(c(10300:10700));
    shift=shift+ind1;
end
shift=shift/8+300;
fs=100;
maxshift=ceil(shift);
L=length(X_org(1,:));
C_t=X_org(:,1:L-maxshift)*X_org(:,maxshift+1:L)';
C_t=(C_t+C_t')/2;
C_all=X_org*X_org';
[V1,D1] = eig(C_t,C_all);
[sVals,sIndex] = sort(diag(D1),'descend') ;
sV1 = V1(:,sIndex) ;
Z1=sV1'*X_org;
offset = max(abs(Z1(:))) ;
feq=100;
disp_eeg(Z1,offset,feq,[]);
title('source\_den2')
lenghtz=length(Z1(1,:));
Z1(2:8,:)=zeros(7,lenghtz);
X_den_GEVD=pinv(sV1')*Z1;
offset = max(abs(X_den_GEVD(:))) ;
% ElecName = Electrodes.labels ;
disp_eeg(X_den_GEVD,offset,feq,[]);
title('X\_den2')
PRMSE2=sqrt(sum(sum(pow2(X1-X_den_GEVD))))/sqrt(sum(sum(pow2(X1))));
%% j T1 on off
C_on=(X_org.*T1)*(X_org.*T1)';
C_all=X_org*X_org';
[V1,D1] = eig(C_on,C_all);
[sVals,sIndex] = sort(diag(D1),'descend') ;
sV1 = V1(:,sIndex) ;
Z1=sV1'*X_org;
offset = max(abs(Z1(:))) ;
feq=100;
disp_eeg(Z1,offset,feq,[]);
title('source\_den2')
lenghtz=length(Z1(1,:));
Z1(2:8,:)=zeros(7,lenghtz);
X_den_GEVD=pinv(sV1')*Z1;
offset = max(abs(X_den_GEVD(:))) ;
% ElecName = Electrodes.labels ;
disp_eeg(X_den_GEVD,offset,feq,[]);
title('X\_den2')
PRMSE2=sqrt(sum(sum(pow2(X2-X_den_GEVD))))/sqrt(sum(sum(pow2(X2))));
%% d T2 on off
C_on=(X_org.*T2)*(X_org.*T2)';
C_all=X_org*X_org';
[V1,D1] = eig(C_on,C_all);
[sVals,sIndex] = sort(diag(D1),'descend') ;
sV1 = V1(:,sIndex) ;
Z1=sV1'*X_org;
offset = max(abs(Z1(:))) ;
feq=100;
disp_eeg(Z1,offset,feq,[]);
title('source\_den2')
lenghtz=length(Z1(1,:));
Z1(2:8,:)=zeros(7,lenghtz);
X_den_GEVD=pinv(sV1')*Z1;
offset = max(abs(X_den_GEVD(:))) ;
% ElecName = Electrodes.labels ;
disp_eeg(X_den_GEVD,offset,feq,[]);
title('X\_den2')
PRMSE2=sqrt(sum(sum(pow2(X2-X_den_GEVD))))/sqrt(sum(sum(pow2(X2))));
%%  h freq 10-15
X_org_f=zeros(8,feq);
for i =1:1:8
    X_org_f(i,:)=abs(fft(X_org(i,:),feq));
end

mask=zeros(1,length(X_org_f(1,:)));
freqdom=10:1:15;
mask(freqdom)=ones(1,length(freqdom));
C_on=(X_org_f.*mask)*(X_org_f.*mask)';
C_all=X_org_f*X_org_f';
[V1,D1] = eig(C_on,C_all);
[sVals,sIndex] = sort(diag(D1),'descend') ;
sV1 = V1(:,sIndex) ;
Z1=sV1'*X_org;
Z1_f=zeros(8,feq);
for i =1:1:8
    Z1_f(i,:)=abs(fft(Z1(i,:),feq));
end
offset = max(abs(Z1_f(:))) ;
feq=100;
disp_eeg(Z1_f,offset,feq,[]);
title('source\_den2')

lenghtz=length(Z1(1,:));
Z1(2:8,:)=zeros(7,lenghtz);
X_den_GEVD=pinv(sV1')*Z1;
offset = max(abs(X_den_GEVD(:))) ;
% ElecName = Electrodes.labels ;
disp_eeg(X_den_GEVD,offset,feq,[]);
title('X\_den2')
PRMSE2=sqrt(sum(sum(pow2(X3-X_den_GEVD))))/sqrt(sum(sum(pow2(X3))));
%%  v freq 5-25

mask=zeros(1,length(X_org_f(1,:)));
freqdom=5:1:25;
mask(freqdom)=ones(1,length(freqdom));
C_on=(X_org_f.*mask)*(X_org_f.*mask)';
C_all=X_org_f*X_org_f';
[V1,D1] = eig(C_on,C_all);
[sVals,sIndex] = sort(diag(D1),'descend') ;
sV1 = V1(:,sIndex) ;
Z1=sV1'*X_org;
Z1_f=zeros(8,feq);
for i =1:1:8
    Z1_f(i,:)=abs(fft(Z1(i,:),feq));
end
offset = max(abs(Z1_f(:))) ;
feq=100;
disp_eeg(Z1_f,offset,feq,[]);
title('source\_den2')

lenghtz=length(Z1(1,:));
Z1(2:8,:)=zeros(7,lenghtz);
X_den_GEVD=pinv(sV1')*Z1;
offset = max(abs(X_den_GEVD(:))) ;
% ElecName = Electrodes.labels ;
disp_eeg(X_den_GEVD,offset,feq,[]);
title('X\_den2')
PRMSE2=sqrt(sum(sum(pow2(X3-X_den_GEVD))))/sqrt(sum(sum(pow2(X3))));
%% dss algorithms
cx=(X_org-mean(X_org')')*(X_org-mean(X_org')')'/length(X_org(1,:));
[V D]=eig(cx);
coeff = V;
D_d=diag(diag(D).^(-0.5))*V';
%EEG_y=zscore(D_d*EEG_sig);
X_org_w=(D_d*X_org);
%% alef dss
w1=rand(1,8);
i=0;
w2=w1+1;
while i<1000  && ~isequal(w1, w2)
    w1=w2;
    i=i+1;
    rp=w1*X_org_w;
    maxshift=4*fs;
    rp2=rp(1:maxshift);
    for j=1:1:ceil(length(rp)/maxshift)-2
        rp2=rp2+rp(j*maxshift+1:(j+1)*maxshift);
    end
    rp2=rp2/ceil(length(rp)/maxshift);
    rp(1:maxshift)=rp2;
    for j=1:1:ceil(length(rp)/maxshift)-2
        rp(j*maxshift+1:(j+1)*maxshift)=rp2;
    end
    w2=X_org_w*rp';
    w2=w2'/norm(w2);
    %disp(w2)
end

figure
plot(rp)
X_den_GEVD=pinv(D_d)*pinv(w1)*rp;
offset = max(abs(X_den_GEVD(:))) ;
% ElecName = Electrodes.labels ;
disp_eeg(X_den_GEVD,offset,feq,[]);
title('X\_den2')
PRMSE2=sqrt(sum(sum(pow2(X2-X_den_GEVD))))/sqrt(sum(sum(pow2(X2))));
%% b dss

w1=rand(1,8);
i=0;
w2=w1+1;
while i<1000  && ~isequal(w1, w2)
    w1=w2;
    i=i+1;
    rp=w1*X_org_w;
    maxshift=ceil(shift);
    rp2=rp(1:maxshift);
    for j=1:1:ceil(length(rp)/maxshift)-2
        rp2=rp2+rp(j*maxshift+1:(j+1)*maxshift);
    end
    rp2=rp2/ceil(length(rp)/maxshift);
    rp(1:maxshift)=rp2;
    for j=1:1:ceil(length(rp)/maxshift)-2
        rp(j*maxshift+1:(j+1)*maxshift)=rp2;
    end
    w2=X_org_w*rp';
    w2=w2'/norm(w2);
    %disp(w2)
end

figure
plot(rp)
X_den_GEVD=pinv(D_d)*pinv(w1)*rp;
offset = max(abs(X_den_GEVD(:))) ;
% ElecName = Electrodes.labels ;
disp_eeg(X_den_GEVD,offset,feq,[]);
title('X\_den2')
PRMSE2=sqrt(sum(sum(pow2(X2-X_den_GEVD))))/sqrt(sum(sum(pow2(X2))));
%% j dss
w1=rand(1,8);
i=0;
w2=w1+1;
while i<1000  && ~isequal(w1, w2)
    w1=w2;
    i=i+1;
    rp=w1*X_org_w;
    rp=rp.*T1;
    w2=X_org_w*rp';
    w2=w2'/norm(w2);
    %disp(w2)
end

figure
plot(rp)
X_den_GEVD=pinv(D_d)*pinv(w1)*rp;
offset = max(abs(X_den_GEVD(:))) ;
% ElecName = Electrodes.labels ;
disp_eeg(X_den_GEVD,offset,feq,[]);
title('X\_den2')
PRMSE2=sqrt(sum(sum(pow2(X2-X_den_GEVD))))/sqrt(sum(sum(pow2(X2))));
%% d dss
w1=rand(1,8);
i=0;
w2=w1+1;
while i<1000  && ~isequal(w1, w2)
    w1=w2;
    i=i+1;
    rp=w1*X_org_w;
    rp=rp.*T2;
    w2=X_org_w*rp';
    w2=w2'/norm(w2);
    %disp(w2)
end

figure
plot(rp)
X_den_GEVD=pinv(D_d)*pinv(w1)*rp;
offset = max(abs(X_den_GEVD(:))) ;
% ElecName = Electrodes.labels ;
disp_eeg(X_den_GEVD,offset,feq,[]);
title('X\_den2')
PRMSE2=sqrt(sum(sum(pow2(X2-X_den_GEVD))))/sqrt(sum(sum(pow2(X2))));
%% h dss

w1=rand(1,8);
i=0;
w2=w1+1;
while i<1000  && ~isequal(w1, w2)
    w1=w2;
    i=i+1;
    rp=w1*X_org_w;
    f_low = 10;
    f_high = 15;
    order = 4; 
    Wn = [f_low, f_high] / (feq/2);
    [b, a] = butter(order, Wn, 'bandpass');
    rp = filter(b, a, rp);
    w2=X_org_w*rp';
    w2=w2'/norm(w2);
    %disp(w2)
end

figure
plot(rp)
X_den_GEVD=pinv(D_d)*pinv(w1)*rp;
offset = max(abs(X_den_GEVD(:))) ;
% ElecName = Electrodes.labels ;
disp_eeg(X_den_GEVD,offset,feq,[]);
title('X\_den2')
PRMSE2=sqrt(sum(sum(pow2(X3-X_den_GEVD))))/sqrt(sum(sum(pow2(X3))));
%% v dss
w1=rand(1,8);
i=0;
w2=w1+1;
while i<1000  && ~isequal(w1, w2)
    w1=w2;
    i=i+1;
    rp=w1*X_org_w;
    f_low = 5;
    f_high = 25;
    order = 4; 
    Wn = [f_low, f_high] / (feq/2);
    [b, a] = butter(order, Wn, 'bandpass');
    rp = filter(b, a, rp);
    w2=X_org_w*rp';
    w2=w2'/norm(w2);
    %disp(w2)
end

figure
plot(rp)
X_den_GEVD=pinv(D_d)*pinv(w1)*rp;
offset = max(abs(X_den_GEVD(:))) ;
% ElecName = Electrodes.labels ;
disp_eeg(X_den_GEVD,offset,feq,[]);
title('X\_den2')
PRMSE2=sqrt(sum(sum(pow2(X3-X_den_GEVD))))/sqrt(sum(sum(pow2(X3))));
