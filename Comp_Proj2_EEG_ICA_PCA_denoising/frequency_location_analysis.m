clc;
clear all;
cd('C:\Users\ASUS\OneDrive\Desktop\uni\simisterm7\eeg\Comp_HW2\Ex3');
%%

elec=load('Electrodes.mat');
lab=elec.Electrodes.labels;

EEG1=load('NewData1.mat');
EEG1_sig=EEG1.EEG_Sig;
offset = max(abs(EEG1_sig(:))) ;
feq = 250 ;
disp_eeg(EEG1_sig,offset,feq,lab);
title('XEEG\_1')

EEG2=load('NewData2.mat');
EEG2_sig=EEG2.EEG_Sig;
offset = max(abs(EEG2_sig(:))) ;
disp_eeg(EEG2_sig,offset,feq,lab);
title('XEEG\_2')

EEG3=load('NewData3.mat');
EEG3_sig=EEG3.EEG_Sig;
offset = max(abs(EEG3_sig(:))) ;
disp_eeg(EEG3_sig,offset,feq,lab);
title('XEEG\_3')

EEG4=load('NewData4.mat');
EEG4_sig=EEG4.EEG_Sig;
%EEG4_sig(17:18,:)=0*EEG4_sig(17:18,:);
offset = max(abs(EEG4_sig(:))) ;
disp_eeg(EEG4_sig,offset,feq,lab);
title('XEEG\_4')
%%

zz=zeros(5,length(EEG3_sig(3,:)));
zz(1:3,:)=EEG3_sig([1,2,14],:);

zz(4,:)=EEG2_sig(3,1:length(EEG3_sig(3,:)));
zz(5,:)=EEG1_sig(1,:);
offset = max(abs(zz(:))) ;
disp_eeg(zz,offset,feq,lab);
title('some artifacts and noises')
%%
[F1,W1,K1]=COM2R(EEG1_sig,32);
Z1=pinv(F1)*EEG1_sig;
offset = max(abs(Z1(:))) ;
disp_eeg(Z1,offset,feq,lab);
title('source1')


[F2,W2,K2]=COM2R(EEG2_sig,32);
Z2=pinv(F2)*EEG2_sig;
offset = max(abs(Z2(:))) ;
disp_eeg(Z2,offset,feq,lab);
title('source2')


[F3,W3,K3]=COM2R(EEG3_sig,32);
Z3=pinv(F3)*EEG3_sig;
offset = max(abs(Z3(:))) ;
disp_eeg(Z3,offset,feq,lab);
title('source3')


[F4,W4,K4]=COM2R(EEG4_sig,32);
Z4=pinv(F4)*EEG4_sig;
offset = max(abs(Z4(:))) ;
disp_eeg(Z4,offset,feq,lab);
title('source4')
%%
electr=load('Electrodes.mat');
electr=electr.Electrodes;
%%
elox_x=electr.X;
eloc_y=electr.Y;
%%
figure

freq1=Z1;
for i=1:1:length(Z1(:,1))
    subplot(3,7,i)
    pwelch(Z1(i,:));  
    title(['source' num2str(i)  'EEG1'])
end
figure

for i=1:1:length(Z2(:,1))
    subplot(3,7,i)
    pwelch(Z2(i,:));  
    title(['source' num2str(i)  'EEG2'])
end
figure


for i=1:1:length(Z3(:,1))
    subplot(3,7,i)
    pwelch(Z3(i,:));  
    title(['source' num2str(i)  'EEG3'])
end
figure
for i=1:1:length(Z4(:,1))
    subplot(3,7,i)
    pwelch(Z4(i,:));  
    title(['source' num2str(i)  'EEG4'])
end
%%
figure
for i=1:1:length(F1(1,:))
    subplot(3,7,i)
    plottopomap(elox_x,eloc_y,electr.labels,F1(:,i))
    title(['source' num2str(i)  'EEG1'])
end

figure
for i=1:1:length(F2(1,:))
    subplot(3,7,i)
    plottopomap(elox_x,eloc_y,electr.labels,F2(:,i))
    title(['source' num2str(i)  'EEG2'])
end
figure
for i=1:1:length(F3(1,:))
    subplot(3,7,i)
    plottopomap(elox_x,eloc_y,electr.labels,F3(:,i))
    title(['source' num2str(i)  'EEG3'])
end
figure
for i=1:1:length(F4(1,:))
    subplot(3,7,i)
    plottopomap(elox_x,eloc_y,electr.labels,F4(:,i))
    title(['source' num2str(i)  'EEG4'])
end
%%
indices1=[1,2,3,5,6,7,8,9,11,12,13,14,15,16,17,18,19,20,21];
indices3=[2,3,4,5,6,8,9,10,11,12,13,14,15,16,17,18,19,20,21];
indices2=[3,4,5,8,9,10,11,12,13,14,16,17,18,19,20,21];
indices4=[9,10,11,12,15,16,18,17,20];
X_den1=F1(:,indices1)*Z1(indices1,:);
X_den2=F2(:,indices2)*Z2(indices2,:);
X_den3=F3(:,indices3)*Z3(indices3,:);
X_den4=F4(:,indices4)*Z4(indices4,:);


offset = max(abs(X_den1(:))) ;
disp_eeg(X_den1,offset,feq,lab);
title('X_den1')

offset = max(abs(X_den2(:))) ;
disp_eeg(X_den2,offset,feq,lab);
title('X_den2')

offset = max(abs(X_den3(:))) ;
disp_eeg(X_den3,offset,feq,lab);
title('X_den3')

offset = max(abs(X_den4(:))) ;
disp_eeg(X_den4,offset,feq,lab);
title('X_den4')