clc;
clear all;
cd('C:\Users\ASUS\OneDrive\Desktop\uni\simisterm7\eeg\Comp_HW2\Ex2');
%%
EEG=load('Ex2.mat');
EEG_org=EEG.X_org;
elec=load('Electrodes.mat');
lab=elec.Electrodes.labels;
%%
offset = max(abs(EEG_org(:))) ;
feq = 256 ;
disp_eeg(EEG_org,offset,feq,lab);
title('X\_or')

EEG_nois1=EEG.X_noise_1;
offset = max(abs(EEG_nois1(:))) ;
disp_eeg(EEG_nois1,offset,feq,lab);
title('X\_noisy\_1')

EEG_nois2=EEG.X_noise_2;

offset = max(abs(EEG_nois2(:))) ;
disp_eeg(EEG_nois2,offset,feq,lab);
title('X\_noisy\_2')

EEG_nois3=EEG.X_noise_3;

offset = max(abs(EEG_nois3(:))) ;
disp_eeg(EEG_nois3,offset,feq,lab);
title('X\_noisy\_3')

EEG_nois4=EEG.X_noise_4;

offset = max(abs(EEG_nois4(:))) ;
disp_eeg(EEG_nois4,offset,feq,lab);
title('X\_noisy\_4')

EEG_nois5=EEG.X_noise_5;

offset = max(abs(EEG_nois5(:))) ;
disp_eeg(EEG_nois5,offset,feq,lab);
title('X\_noisy\_5')
%%
for snr=-20:2:-10
    EEG_noised=noisingEEG(EEG_org,EEG_nois4,snr);
    offset = max(abs(EEG_noised(:))) ;
    disp_eeg(EEG_noised,offset,feq,lab);
    title(['SNR' num2str(snr) 'db'])
end
%%
for snr=-20:2:-10
    EEG_noised=noisingEEG(EEG_org,EEG_nois5,snr);
    offset = max(abs(EEG_noised(:))) ;
    disp_eeg(EEG_noised,offset,feq,lab);
    title(['SNR' num2str(snr) 'db'])
end

%%
X_org=EEG_org;
rrmse=zeros(1,length(-20:2:-10));
i=1;
for snr=-20:2:-10
    EEG_noised=noisingEEG(EEG_org,EEG_nois4,snr);
    
    [F,W,K]=COM2R(EEG_noised,32);
    Z=pinv(F)*EEG_noised;
    offset = max(abs(Z(:))) ;
    disp_eeg(Z,offset,feq,lab);
    title(['SNR' num2str(snr) 'db'])
    if snr==-20
    lenghtz=length(Z(1,:));
    Z(1:24,:)=zeros(24,lenghtz);
    Z(26:32,:)=zeros(7,lenghtz);
    end
    if snr==-18
    lenghtz=length(Z(1,:));
    Z(1:21,:)=zeros(21,lenghtz);
    Z(23:32,:)=zeros(10,lenghtz);
    end
    if snr==-16
    lenghtz=length(Z(1,:));
    Z(1:19,:)=zeros(19,lenghtz);
    Z(21:32,:)=zeros(12,lenghtz);
    end
    if snr==-14
    lenghtz=length(Z(1,:));
    Z(1:13,:)=zeros(13,lenghtz);
    Z(15:32,:)=zeros(18,lenghtz);
    end
    if snr==-12
    lenghtz=length(Z(1,:));
    Z(1:6,:)=zeros(6,lenghtz);
    Z(8:32,:)=zeros(25,lenghtz);
    end
    if snr==-10
    lenghtz=length(Z(1,:));
    Z(1:2,:)=zeros(2,lenghtz);
    Z(4:32,:)=zeros(29,lenghtz);
    end
    X_den=F*Z;
    offset = max(abs(X_den(:))) ;
    disp_eeg(X_den,offset,feq,lab);
    title(['X_den' 'SNR' num2str(snr) 'db'])
    X=zeros(4,lenghtz);
    X(1,:)=X_den(13,:);
    X(2,:)=X_den(24,:);
    X(3,:)=X_org(13,:);
    X(4,:)=X_org(24,:);
    offset = max(abs(X(:))) ;
    disp_eeg(X,offset,feq,lab);
    title('chanel 13 & 14 den and the org of them sorted')
    PRMSE1=sqrt(sum(sum(pow2(X_org-X_den))))/sqrt(sum(sum(pow2(X_org))));
    rrmse(i)=PRMSE1;
    i=i+1;
end

%%
rrmse2=zeros(1,length(-20:2:-10));
i=1;
for snr=-20:2:-10
    EEG_noised=noisingEEG(EEG_org,EEG_nois5,snr);
    
    [F,W,K]=COM2R(EEG_noised,32);
    Z=pinv(F)*EEG_noised;
    offset = max(abs(Z(:))) ;
    disp_eeg(Z,offset,feq,lab);
    title(['SNR' num2str(snr) 'db'])
    if snr==-20
    lenghtz=length(Z(1,:));
    Z(1:18,:)=zeros(18,lenghtz);
    Z(20:32,:)=zeros(13,lenghtz);
    end
    if snr==-18
    lenghtz=length(Z(1,:));
    Z(1:18,:)=zeros(18,lenghtz);
    Z(20:32,:)=zeros(13,lenghtz);
    end
    if snr==-16
    lenghtz=length(Z(1,:));
    Z(1:12,:)=zeros(12,lenghtz);
    Z(14:32,:)=zeros(19,lenghtz);
    end
    if snr==-14
    lenghtz=length(Z(1,:));
    Z(1:10,:)=zeros(10,lenghtz);
    Z(12:32,:)=zeros(21,lenghtz);
    end
    if snr==-12
    lenghtz=length(Z(1,:));
    Z(1:6,:)=zeros(6,lenghtz);
    Z(8:32,:)=zeros(25,lenghtz);
    end
    if snr==-10
    lenghtz=length(Z(1,:));
    Z(1:3,:)=zeros(3,lenghtz);
    Z(5:32,:)=zeros(28,lenghtz);
    end
    X_den=F*Z;
    offset = max(abs(X_den(:))) ;
    disp_eeg(X_den,offset,feq,lab);
    title(['X_den' 'SNR' num2str(snr) 'db'])
    X=zeros(4,lenghtz);
    X(1,:)=X_den(13,:);
    X(2,:)=X_den(24,:);
    X(3,:)=X_org(13,:);
    X(4,:)=X_org(24,:);
    feq = 256 ;
    offset = max(abs(X(:))) ;
    disp_eeg(X,offset,feq,lab);
    title('chanel 13 & 14 den and the org of them sorted')
    PRMSE1=sqrt(sum(sum(pow2(X_org-X_den))))/sqrt(sum(sum(pow2(X_org))));
    rrmse2(i)=PRMSE1;
    i=i+1;
end

%%
rrmse3=zeros(1,length(-20:2:-10));
i=1;
for snr=-20:2:-10
    EEG_noised=noisingEEG(EEG_org,EEG_nois5,snr);
    
    [coeff,~,latent]=pca(EEG_noised');
    V=coeff;
    D=diag(latent);
    D_d=diag(diag(D).^(-0.5))*V';
    Z=(D_d*EEG_noised);
    offset = max(abs(Z(:))) ;
    disp_eeg(Z,offset,feq,lab);
    title(['SNR' num2str(snr) 'db'])
    if snr==-20
    lenghtz=length(Z(1,:));
    Z(1:23,:)=zeros(23,lenghtz);
    Z(25:28,:)=zeros(4,lenghtz);
    Z(30:32,:)=zeros(3,lenghtz);
    end
    if snr==-18
    lenghtz=length(Z(1,:));
    Z(1:22,:)=zeros(22,lenghtz);
    Z(24:28,:)=zeros(5,lenghtz);
    Z(30:32,:)=zeros(3,lenghtz);
    end
    if snr==-16
    lenghtz=length(Z(1,:));   
    Z(1:22,:)=zeros(22,lenghtz);
    Z(24:28,:)=zeros(5,lenghtz);
    Z(30:32,:)=zeros(3,lenghtz);
    end
    if snr==-14
    lenghtz=length(Z(1,:));
    Z(1:22,:)=zeros(22,lenghtz);
    Z(24:28,:)=zeros(5,lenghtz);
    Z(30:32,:)=zeros(3,lenghtz);
    end
    if snr==-12
    lenghtz=length(Z(1,:));
    Z(1:21,:)=zeros(21,lenghtz);
    Z(23:28,:)=zeros(6,lenghtz);
    Z(30:32,:)=zeros(3,lenghtz);
    end
    if snr==-10
    lenghtz=length(Z(1,:));
    Z(1:21,:)=zeros(21,lenghtz);
    Z(23:26,:)=zeros(4,lenghtz);
    Z(28:32,:)=zeros(5,lenghtz);
    end
    X_den=inv(D_d)*Z;
    offset = max(abs(X_den(:))) ;
    disp_eeg(X_den,offset,feq,lab);
    title(['X_den' 'SNR' num2str(snr) 'db'])
    X=zeros(4,lenghtz);
    X(1,:)=X_den(13,:);
    X(2,:)=X_den(24,:);
    X(3,:)=X_org(13,:);
    X(4,:)=X_org(24,:);
    offset = max(abs(X(:))) ;
    disp_eeg(X,offset,feq,lab);
    title('chanel 13 & 14 den and the org of them sorted')
    PRMSE1=sqrt(sum(sum(pow2(X_org-X_den))))/sqrt(sum(sum(pow2(X_org))));
    rrmse3(i)=PRMSE1;
    i=i+1;
end
%%
rrmse4=zeros(1,length(-20:2:-10));
i=1;
for snr=-20:2:-10
    EEG_noised=noisingEEG(EEG_org,EEG_nois4,snr);
    
    [coeff,~,latent]=pca(EEG_noised');
    V=coeff;
    D=diag(latent);
    D_d=diag(diag(D).^(-0.5))*V';
    Z=(D_d*EEG_noised);
    offset = max(abs(Z(:))) ;
    disp_eeg(Z,offset,feq,lab);
    title(['SNR' num2str(snr) 'db'])
    if snr==-20
    lenghtz=length(Z(1,:));
    Z(1:30,:)=zeros(30,lenghtz);
    Z(32,:)=zeros(1,lenghtz);
    end
    if snr==-18
    lenghtz=length(Z(1,:));
    Z(1:30,:)=zeros(30,lenghtz);
    Z(32,:)=zeros(1,lenghtz);
    end
    if snr==-16
    lenghtz=length(Z(1,:));   
    Z(1:30,:)=zeros(30,lenghtz);
    Z(32,:)=zeros(1,lenghtz);
    end
    if snr==-14
    lenghtz=length(Z(1,:));
    Z(1,:)=zeros(1,lenghtz);
    Z(3:30,:)=zeros(28,lenghtz);
    Z(32,:)=zeros(1,lenghtz);
    end
    if snr==-12
    lenghtz=length(Z(1,:));
    Z(1,:)=zeros(1,lenghtz);
    Z(3:30,:)=zeros(28,lenghtz);
    Z(32,:)=zeros(1,lenghtz);
    end
    if snr==-10
    lenghtz=length(Z(1,:));
    Z(1,:)=zeros(1,lenghtz);
    Z(3:30,:)=zeros(28,lenghtz);
    Z(32,:)=zeros(1,lenghtz);
    end
    X_den=inv(D_d)*Z;
    offset = max(abs(X_den(:))) ;
    disp_eeg(X_den,offset,feq,lab);
    title(['X_den' 'SNR' num2str(snr) 'db'])
    X=zeros(4,lenghtz);
    X(1,:)=X_den(13,:);
    X(2,:)=X_den(24,:);
    X(3,:)=X_org(13,:);
    X(4,:)=X_org(24,:);
    offset = max(abs(X(:))) ;
    disp_eeg(X,offset,feq,lab);
    title('chanel 13 & 14 den and the org of them sorted')
    PRMSE1=sqrt(sum(sum(pow2(X_org-X_den))))/sqrt(sum(sum(pow2(X_org))));
    rrmse4(i)=PRMSE1;
    i=i+1;
end
%%
figure;
subplot(2,2,1)
plot(-20:2:-10,rrmse)
title('ICA noise 4')
subplot(2,2,2)
plot(-20:2:-10,rrmse2)
title('ICA noise 5')
subplot(2,2,3)
plot(-20:2:-10,rrmse3)
title('pCA noise 5')

subplot(2,2,4)
plot(-20:2:-10,rrmse4)
title('pCA noise 4')
%%
function addN=noisingEEG(EEG_org,EEG_noise,snr)
Ps=sum(sum(EEG_org.*EEG_org));
Pn=sum(sum(EEG_noise.*EEG_noise));
k=sqrt(Ps/10^(snr/10)/Pn);
addN=EEG_org+k*EEG_noise;
end
