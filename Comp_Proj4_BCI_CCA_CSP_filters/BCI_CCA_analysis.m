clc;
clear all;
cd('C:\Users\ASUS\OneDrive\Desktop\uni\simisterm7\eeg\Comp_HW4');
%%
SSVEP=load('SSVEP_EEG.mat');
SSVEP_Signal=SSVEP.SSVEP_Signal;
Event_samples=SSVEP.Event_samples;
Events=SSVEP.Events;
%%
fs=250;
fpass=[1,40];
order=30;
%% 
SSVEP_Signal_filterd=SSVEP_Signal;
for i=1:1:6
    SSVEP_Signal_filterd(i,:)=bandpass(SSVEP_Signal(i,:),[1 40],fs);
end
%%
evokes=zeros(15,6,5*fs);
for i =1:15
   evokes(i,:,:)= SSVEP_Signal_filterd(:,Event_samples(i):Event_samples(i)+5*fs-1);
end
%%
figure 
for i=1:15
    subplot(4,4,i)
    for j=1:6   
        [pxx,w] = pwelch(squeeze(evokes(i,j,:)),500,300,500,fs);
        plot(w/pi*fs/2/5000*128,10*log10(pxx))
        hold on
    end
    legend({'Pz','Oz','P7','P8',' O2','O1'})
        title(['frequency domain of experiment :' num2str(i)])

    xlabel("frequency");
    ylabel('magnitude');
end
%% cca
f_evok=unique(Events);
Y1=[];
Y2=[];
Y3=[];
Y4=[];
Y5=[];
t=length(evokes(1,1,:));
for  j=1:ceil(40/f_evok(1))
    Y1=[Y1;cos(2*pi*(1:t)*j*f_evok(1)/250);sin(2*pi*(1:t)*j*f_evok(1))/250];
end
for  j=1:ceil(40/f_evok(2))
    Y2=[Y2;cos(2*pi*(1:t)*j*f_evok(2)/250);sin(2*pi*(1:t)*j*f_evok(2))/250];
end
for  j=1:ceil(40/f_evok(3))
    Y3=[Y3;cos(2*pi*(1:t)*j*f_evok(3)/250);sin(2*pi*(1:t)*j*f_evok(3))/250];
end
for  j=1:ceil(40/f_evok(4))
    Y4=[Y4;cos(2*pi*(1:t)*j*f_evok(4)/250);sin(2*pi*(1:t)*j*f_evok(4))/250];
end
for  j=1:ceil(40/f_evok(5))
    Y5=[Y5;cos(2*pi*(1:t)*j*f_evok(5)/250);sin(2*pi*(1:t)*j*f_evok(5))/250];
end
%%
cors=zeros(15,5);

for i=1:15
    [A,B,r,U,V] = canoncorr(squeeze(evokes(i,:,:))',Y1');
    cors(i,1)=max(r);
    [A,B,r,U,V] = canoncorr(squeeze(evokes(i,:,:))',Y2');
    cors(i,2)=max(r);
    [A,B,r,U,V] = canoncorr(squeeze(evokes(i,:,:))',Y3');
    cors(i,3)=max(r);
    [A,B,r,U,V] = canoncorr(squeeze(evokes(i,:,:))',Y4');
    cors(i,4)=max(r);
    [A,B,r,U,V] = canoncorr(squeeze(evokes(i,:,:))',Y5');
    cors(i,5)=max(r);
    
end
%% accuracy
corss=cors';
[max1,ind]=max(corss);
predfreq=zeros(1,15);
for i=1:15
    
    predfreq(i)=f_evok(ind(i));
end
acc=sum(predfreq==Events)/15;
%% cca with one channel
evokes2=zeros(15,1,5*fs);
for i =1:15
   evokes2(i,:,:)= SSVEP_Signal_filterd(1,Event_samples(i):Event_samples(i)+5*fs-1);
end
cors2=zeros(15,5);

for i=1:15
    [A,B,r,U,V] = canoncorr(squeeze(evokes2(i,:,:)),Y1');
    cors2(i,1)=max(r);
    [A,B,r,U,V] = canoncorr(squeeze(evokes2(i,:,:)),Y2');
    cors2(i,2)=max(r);
    [A,B,r,U,V] = canoncorr(squeeze(evokes2(i,:,:)),Y3');
    cors2(i,3)=max(r);
    [A,B,r,U,V] = canoncorr(squeeze(evokes2(i,:,:)),Y4');
    cors2(i,4)=max(r);
    [A,B,r,U,V] = canoncorr(squeeze(evokes2(i,:,:)),Y5');
    cors2(i,5)=max(r);
    
end
%disp(corr(U,V))
%%
corss=cors2';
[max1,ind]=max(corss);
predfreq=zeros(1,15);
for i=1:15
    
    predfreq(i)=f_evok(ind(i));
end
acc2=sum(predfreq==Events)/15;
%% cca with 6 channel and but 2 sec
evokes3=zeros(15,6,2*fs);
for i =1:15
   evokes3(i,:,:)= SSVEP_Signal_filterd(:,Event_samples(i):Event_samples(i)+2*fs-1);
end
cors3=zeros(15,5);
Y1=[];
Y2=[];
Y3=[];
Y4=[];
Y5=[];
t=length(evokes3(1,1,:));
for  j=1:ceil(40/f_evok(1))
    Y1=[Y1;cos(2*pi*(1:t)*j*f_evok(1)/250);sin(2*pi*(1:t)*j*f_evok(1))/250];
end
for  j=1:ceil(40/f_evok(2))
    Y2=[Y2;cos(2*pi*(1:t)*j*f_evok(2)/250);sin(2*pi*(1:t)*j*f_evok(2))/250];
end
for  j=1:ceil(40/f_evok(3))
    Y3=[Y3;cos(2*pi*(1:t)*j*f_evok(3)/250);sin(2*pi*(1:t)*j*f_evok(3))/250];
end
for  j=1:ceil(40/f_evok(4))
    Y4=[Y4;cos(2*pi*(1:t)*j*f_evok(4)/250);sin(2*pi*(1:t)*j*f_evok(4))/250];
end
for  j=1:ceil(40/f_evok(5))
    Y5=[Y5;cos(2*pi*(1:t)*j*f_evok(5)/250);sin(2*pi*(1:t)*j*f_evok(5))/250];
end
for i=1:15
    [A,B,r,U,V] = canoncorr(squeeze(evokes3(i,:,:))',Y1');
    cors3(i,1)=max(r);
    [A,B,r,U,V] = canoncorr(squeeze(evokes3(i,:,:))',Y2');
    cors3(i,2)=max(r);
    [A,B,r,U,V] = canoncorr(squeeze(evokes3(i,:,:))',Y3');
    cors3(i,3)=max(r);
    [A,B,r,U,V] = canoncorr(squeeze(evokes3(i,:,:))',Y4');
    cors3(i,4)=max(r);
    [A,B,r,U,V] = canoncorr(squeeze(evokes3(i,:,:))',Y5');
    cors3(i,5)=max(r);
    
end
%disp(corr(U,V))
%%
corss=cors3';
[max1,ind]=max(corss);
predfreq=zeros(1,15);
for i=1:15
    
    predfreq(i)=f_evok(ind(i));
end
acc3=sum(predfreq==Events)/15;
%%
%% cca with 1 channel and but 2 sec
evokes3=zeros(15,1,2*fs);
for i =1:15
   evokes3(i,:,:)= SSVEP_Signal_filterd(1,Event_samples(i):Event_samples(i)+2*fs-1);
end
cors3=zeros(15,5);
Y1=[];
Y2=[];
Y3=[];
Y4=[];
Y5=[];
t=length(evokes3(1,1,:));
for  j=1:ceil(40/f_evok(1))
    Y1=[Y1;cos(2*pi*(1:t)*j*f_evok(1)/250);sin(2*pi*(1:t)*j*f_evok(1))/250];
end
for  j=1:ceil(40/f_evok(2))
    Y2=[Y2;cos(2*pi*(1:t)*j*f_evok(2)/250);sin(2*pi*(1:t)*j*f_evok(2))/250];
end
for  j=1:ceil(40/f_evok(3))
    Y3=[Y3;cos(2*pi*(1:t)*j*f_evok(3)/250);sin(2*pi*(1:t)*j*f_evok(3))/250];
end
for  j=1:ceil(40/f_evok(4))
    Y4=[Y4;cos(2*pi*(1:t)*j*f_evok(4)/250);sin(2*pi*(1:t)*j*f_evok(4))/250];
end
for  j=1:ceil(40/f_evok(5))
    Y5=[Y5;cos(2*pi*(1:t)*j*f_evok(5)/250);sin(2*pi*(1:t)*j*f_evok(5))/250];
end
for i=1:15
    [A,B,r,U,V] = canoncorr(squeeze(evokes3(i,:,:)),Y1');
    cors3(i,1)=max(r);
    [A,B,r,U,V] = canoncorr(squeeze(evokes3(i,:,:)),Y2');
    cors3(i,2)=max(r);
    [A,B,r,U,V] = canoncorr(squeeze(evokes3(i,:,:)),Y3');
    cors3(i,3)=max(r);
    [A,B,r,U,V] = canoncorr(squeeze(evokes3(i,:,:)),Y4');
    cors3(i,4)=max(r);
    [A,B,r,U,V] = canoncorr(squeeze(evokes3(i,:,:)),Y5');
    cors3(i,5)=max(r);
    
end
%disp(corr(U,V))
%%
corss=cors3';
[max1,ind]=max(corss);
predfreq=zeros(1,15);
for i=1:15
    
    predfreq(i)=f_evok(ind(i));
end
acc4=sum(predfreq==Events)/15;