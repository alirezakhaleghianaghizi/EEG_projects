clc;
clear all;
cd('C:\Users\ASUS\OneDrive\Desktop\uni\simisterm7\eeg\Comp_HW2\Ex2');
%%
EEG=load('Ex2.mat');
EEG_org=EEG.X_org;
elec=load('Electrodes.mat');
lab=[];
offset = max(abs(EEG_org(:))) ;
feq = 256 ;
disp_eeg(EEG_org,offset,feq,lab);
title('X\_org')

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
title('noise')

EEG_nois5=EEG.X_noise_5;

offset = max(abs(EEG_nois5(:))) ;
disp_eeg(EEG_nois5,offset,feq,lab);
title('X\_noisy\_5')
%%
for snr=-20:2:-10
    EEG_noised=noisingEEG(EEG_org,EEG_nois4,snr);
    title(['SNR' num2str(snr) 'db'])
end

offset = max(abs(EEG_noised(:))) ;
disp_eeg(EEG_noised,offset,feq,lab);
title('X\_noised')
%%
X_org=EEG_org;
rrmse=zeros(1,length(-20:2:-10));
i=1;
%GEVD or generalized eigenvalue decomposition
sparki1=(8.9375:1/feq:10.65625)*feq;
sparki2=(11.25:1/feq:12.25)*feq;
sparki3=(15:1/feq:15.75)*feq;
sparki4=(17.65625:1/feq:18.75)*feq;
sparki5=(21:1/feq:22.9375)*feq;
sparki6=(24.5:1/feq:25.5)*feq;
sparki7=(31.75:1/feq:32.65625)*feq;
sparki8=(35.75:1/feq:36.9375)*feq;
sparki9=(39.5:1/feq:40)*feq;


for snr=-20:2:-10
    EEG_noised=noisingEEG(EEG_org,EEG_nois4,snr);
    sparks=EEG_noised(:,[sparki1,sparki2,sparki3,sparki4,sparki5,sparki6,sparki7,sparki8,sparki9]);
    cov_sparks=cov(sparks');
    cov_all_times=cov(EEG_noised');
    [V1,D1] = eig(cov_sparks,cov_all_times);
    [sVals,sIndex] = sort(diag(D1),'descend') ;
    sV1 = V1(:,sIndex) ;
    Z=sV1'*EEG_noised;
    offset = max(abs(Z(:))) ;
    disp_eeg(Z,offset,feq,[]);
    title('source\_den2')
    if snr==-20
    lenghtz=length(Z(1,:));
    Z(2:32,:)=zeros(31,lenghtz);
    end
    if snr==-18
    lenghtz=length(Z(1,:));
    Z(2:32,:)=zeros(31,lenghtz);
    end
    if snr==-16
    lenghtz=length(Z(1,:));
    Z(3:32,:)=zeros(30,lenghtz);
    end
    if snr==-14
    lenghtz=length(Z(1,:));
    Z(4:32,:)=zeros(29,lenghtz);
    end
    if snr==-12
    lenghtz=length(Z(1,:));
    Z(4:32,:)=zeros(29,lenghtz);
    end
    if snr==-10
    lenghtz=length(Z(1,:));
    Z(4:32,:)=zeros(29,lenghtz);
    end
    X_den=pinv(sV1')*Z;
    offset = max(abs(X_den(:))) ;
    disp_eeg(X_den,offset,feq,lab);
    title(['X_den' 'SNR' num2str(snr) 'db'])
    X=zeros(2,lenghtz);
    X(1,:)=X_den(13,:);
    %X(2,:)=X_den(24,:);
    X(2,:)=X_org(13,:);
    %X(4,:)=X_org(24,:);
    offset = max(abs(X(:))) ;
    disp_eeg(X,offset,feq,lab);
    title('')
    PRMSE1=sqrt(sum(sum(pow2(X_org-X_den))))/sqrt(sum(sum(pow2(X_org))));
    rrmse(i)=PRMSE1;
    i=i+1;
end

%%
rrmse2=zeros(1,length(-20:2:-10));
i=1;
for snr=-20:2:-10
    EEG_noised=noisingEEG(EEG_org,EEG_nois5,snr);
    
    sparks=EEG_noised(:,[sparki1,sparki2,sparki3,sparki4,sparki5,sparki6,sparki7,sparki8,sparki9]);
    cov_sparks=cov(sparks');
    cov_all_times=cov(EEG_noised');
    [V1,D1] = eig(cov_sparks,cov_all_times);
    [sVals,sIndex] = sort(diag(D1),'descend') ;
    sV1 = V1(:,sIndex) ;
    Z=sV1'*EEG_noised;
    offset = max(abs(Z(:))) ;
    disp_eeg(Z,offset,feq,[]);
    title('source\_den2')
    if snr==-20
    lenghtz=length(Z(1,:));
    Z(3:32,:)=zeros(30,lenghtz);
    end
    if snr==-18
    lenghtz=length(Z(1,:));
    Z(3:32,:)=zeros(30,lenghtz);
    end
    if snr==-16
    lenghtz=length(Z(1,:));
    Z(3:32,:)=zeros(30,lenghtz);
    end
    if snr==-14
    lenghtz=length(Z(1,:));
    Z(4:32,:)=zeros(29,lenghtz);
    end
    if snr==-12
    lenghtz=length(Z(1,:));
    Z(4:32,:)=zeros(29,lenghtz);
    end
    if snr==-10
    lenghtz=length(Z(1,:));
    Z(4:32,:)=zeros(29,lenghtz);
    end
    X_den=pinv(sV1')*Z;
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
    rrmse2(i)=PRMSE1;
    i=i+1;
end

%%
cx=(X_org-mean(X_org')')*(X_org-mean(X_org')')'/length(X_org(1,:));
[V D]=eig(cx);
coeff = V;
D_d=diag(diag(D).^(-0.5))*V';
%EEG_y=zscore(D_d*EEG_sig);
X_org_w=(D_d*X_org);

rrmse3=zeros(1,length(-20:2:-10));
m=1;
T1=zeros(1,length(EEG_org(1,:)));
T1([sparki1,sparki2,sparki3,sparki4,sparki5,sparki6,sparki7,sparki8,sparki9])=ones(1,length([sparki1,sparki2,sparki3,sparki4,sparki5,sparki6,sparki7,sparki8,sparki9]));
for snr=-20:2:-10
    EEG_noised=noisingEEG(EEG_org,EEG_nois5,snr); 
    cx=(EEG_noised-mean(EEG_noised')')*(EEG_noised-mean(EEG_noised')')'/length(EEG_noised(1,:));
    [V D]=eig(cx);
    coeff = V;
    D_d=diag(diag(D).^(-0.5))*V';
    %EEG_y=zscore(D_d*EEG_sig);
    X_org_w=(D_d*EEG_noised);
    W=[];
    for k=1:1:5
        w1=rand(1,32);
        i=0;
        w2=w1+1;
        while i<1000  && ~isequal(w1, w2)
            w1=w2;
            i=i+1;
            rp=w1*X_org_w;
            rp=rp.*T1;
            w2=X_org_w*rp';
            if k>1
               w2=w2-W*W'*w2; 
            end
            w2=w2'/norm(w2);
            %disp(w2)
        end
        W=[W,w2'];
    end
    
    Z=W'*X_org_w;
    offset = max(abs(Z(:))) ;
    disp_eeg(Z,offset,feq,[]);
    title('source\_den2')
    if snr==-20
    lenghtz=length(Z(1,:));
    Z(3:5,:)=zeros(3,lenghtz);
    end
    if snr==-18
    lenghtz=length(Z(1,:));
    Z(3:5,:)=zeros(3,lenghtz);
    end
    if snr==-16
    lenghtz=length(Z(1,:));
    Z(3:5,:)=zeros(3,lenghtz);
    end
    if snr==-14
    lenghtz=length(Z(1,:));
    Z(4:5,:)=zeros(2,lenghtz);
    end
    if snr==-12
    lenghtz=length(Z(1,:));
    Z(4:5,:)=zeros(2,lenghtz);
    end
    if snr==-10
    lenghtz=length(Z(1,:));
    Z(4:5,:)=zeros(2,lenghtz);
    end
    X_den=pinv(D_d)*pinv(W')*Z;
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
    rrmse3(m)=PRMSE1;
    m=m+1;
end
%%
rrmse4=zeros(1,length(-20:2:-10));
m=1;
for snr=-20:2:-10
    EEG_noised=noisingEEG(EEG_org,EEG_nois4,snr); 
    cx=(EEG_noised-mean(EEG_noised')')*(EEG_noised-mean(EEG_noised')')'/length(EEG_noised(1,:));
    [V D]=eig(cx);
    coeff = V;
    D_d=diag(diag(D).^(-0.5))*V';
    %EEG_y=zscore(D_d*EEG_sig);
    X_org_w=(D_d*EEG_noised);
    W=[];
    for k=1:1:5
        w1=rand(1,32);
        i=0;
        w2=w1+1;
        while i<1000  && ~isequal(w1, w2)
            w1=w2;
            i=i+1;
            rp=w1*X_org_w;
            rp=rp.*T1;
            w2=X_org_w*rp';
            if k>1
               w2=w2-W*W'*w2; 
            end
            w2=w2'/norm(w2);
            %disp(w2)
        end
        W=[W,w2'];
    end
    
    Z=W'*X_org_w;
    offset = max(abs(Z(:))) ;
    disp_eeg(Z,offset,feq,[]);
    title('source\_den2')
    if snr==-20
    lenghtz=length(Z(1,:));
    Z(2:5,:)=zeros(4,lenghtz);
    end
    if snr==-18
    lenghtz=length(Z(1,:));
    Z(2:5,:)=zeros(4,lenghtz);
    end
    if snr==-16
    lenghtz=length(Z(1,:));
    Z(3:5,:)=zeros(3,lenghtz);
    end
    if snr==-14
    lenghtz=length(Z(1,:));
    Z(4:5,:)=zeros(2,lenghtz);
    end
    if snr==-12
    lenghtz=length(Z(1,:));
    Z(4:5,:)=zeros(2,lenghtz);
    end
    if snr==-10
    lenghtz=length(Z(1,:));
    Z(4:5,:)=zeros(2,lenghtz);
    end
    X_den=pinv(D_d)*pinv(W')*Z;
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
    rrmse4(m)=PRMSE1;
    m=m+1;
end
%%
figure;
subplot(2,2,1)
plot(-20:2:-10,rrmse)
title('rrmse gevd noize 4')
subplot(2,2,2)
plot(-20:2:-10,rrmse2)
title('rrmse gevd noize 5')
subplot(2,2,3)
plot(-20:2:-10,rrmse3)
title('rrmse dss noize 5')
subplot(2,2,4)
plot(-20:2:-10,rrmse4)
title('rrmse dss noize 4')
%%
figure
subplot(2,1,1)
plot(-20:2:-10,rrmse)
title('BSS with GEVD ')
xlabel("SNR(dB)");
ylabel('RRMSE');
subplot(2,1,2)
plot(-20:2:-10,rrmse4)
title('BSS with DSS')

xlabel("SNR(dB)");
ylabel('RRMSE');
%%
function addN=noisingEEG(EEG_org,EEG_noise,snr)
Ps=sum(sum(EEG_org.*EEG_org));
Pn=sum(sum(EEG_noise.*EEG_noise));
k=sqrt(Ps/10^(snr/10)/Pn);
addN=EEG_org+k*EEG_noise;
end
