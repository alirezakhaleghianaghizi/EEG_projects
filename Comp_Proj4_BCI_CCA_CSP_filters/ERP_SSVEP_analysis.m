clc;
clear all;
cd('C:\Users\ASUS\OneDrive\Desktop\uni\simisterm7\eeg\Comp_HW4');
%%
ERP=load('ERP_EEG.mat');
ERP_sig=ERP.ERP_EEG;
%%
figure
for N=100:100:2500
    mean_N=mean(ERP_sig(:,1:N)');
    plot((1:240)/240*1000,mean_N)
    hold on
    title('ERP')

    xlabel("time(ms)");
    ylabel('magnitude');
end
%%
figure
max_doms=zeros(2550,1);
for N=1:2550
    max_dom=max(mean(ERP_sig(:,1:N)'));
    max_doms(N)=max_dom;
end
plot(1:2550,max_doms)
    title('max ERP magnitud')

    xlabel("number of evokes that averaged");
    ylabel('magnitude');
%%
figure
rms=zeros(2549,1);
prev_algo=ERP_sig(:,1)';
for N=2:2550
    next_algo=mean(ERP_sig(:,1:N)');
    dif=sqrt(mean((next_algo-prev_algo).^2)/N);
    rms(N-1)=dif;
end
plot(1:2549,rms)
    title('dif ERP algorithm')

    xlabel("number of evokes used in algorithm");
    ylabel('rms dif');
%%
figure
N0=600;
N0_algo=mean(ERP_sig(:,1:N0)');
plot((1:240)/240*1000,N0_algo)
hold on 

all_algo=mean(ERP_sig(:,1:2550)');
plot((1:240)/240*1000,all_algo)

N03=N0/3;
N03_algo=mean(ERP_sig(:,1:N03)');
plot((1:240)/240*1000,N03_algo)

N0_rand_algo=mean(ERP_sig(:,randi(2550,[1,N0]))');
plot((1:240)/240*1000,N0_rand_algo)

    title('ERP')

    xlabel("time(ms)");
    ylabel('magnitude');
N03_rand_algo=mean(ERP_sig(:,randi(2550,[1,N03]))');
plot((1:240)/240*1000,N03_rand_algo)
legend((num2str(N0)),(num2str(2550)),num2str(N03),[num2str(N0) 'random'],[num2str(N03) 'random'])
%%
