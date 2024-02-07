clc;
clear all;
cd('C:\Users\ASUS\OneDrive\Desktop\uni\simisterm7\eeg\Comp_HW4');
%%
EX3=load('Ex3.mat');
TestData=EX3.TestData;
TrainData=EX3.TrainData;
TrainLabel=EX3.TrainLabel;

class1=TrainData(:,:,TrainLabel==0);
class2=TrainData(:,:,TrainLabel==1);
C1=zeros(30);

for i=1:length(class1(1,1,:))
    sq_c1=squeeze(class1(:,:,i));
    
   C1= C1+sq_c1*sq_c1';
end
C1=C1/length(class1(1,1,:));

C2=zeros(30);

for i=1:length(class2(1,1,:))
    sq_c2=squeeze(class2(:,:,i));
   C2= C2+sq_c2*sq_c2';
end
C2=C2/length(class2(1,1,:));
%%
[V1,D1] = eig(C1,C2);
[sVals,sIndex] = sort(diag(D1),'descend') ;
sV1 = V1(:,sIndex) ;
figure
for i=1:1:9
    subplot(3,6,2*i-1)  
    Z1=sV1(:,[1,length(sV1(1,:))])'*class1(:,:,i);
    
    plot((1:256)/256,Z1(1,:))
    title(['csp1 class 1' num2str(i)])

    xlabel("time(ms)");
    ylabel('magnitude');
    subplot(3,6,2*i) 
    plot((1:256)/256,Z1(2,:))
    title(['csp2 class 1' num2str(i)])

    xlabel("time(ms)");
    ylabel('magnitude');
    
end
figure
for i=1:1:9
    subplot(3,6,2*i-1)  
    Z2=sV1(:,[1,length(sV1(1,:))])'*class2(:,:,i);
    plot((1:256)/256,Z2(1,:))
        title(['csp1 class 2' num2str(i)])

    xlabel("time(ms)");
    ylabel('magnitude');
    subplot(3,6,2*i) 
    plot((1:256)/256,Z2(2,:))
        title(['csp2 class 2' num2str(i)])

    xlabel("time(ms)");
    ylabel('magnitude');
end
%%
X=sV1(:,[1,length(sV1(1,:))])';
offset = max(abs(X(:))) ;
feq = 1 ;
% ElecName = Electrodes.labels ;
disp_eeg(X,offset,feq,[]);
%%
electr=load('AllElectrodes.mat');
electr=electr.AllElectrodes;
electr=electr([37,7,5,38,40,42,10,47,45,15,13,48,50,52,18,32,55,23,22,21,20,31,57,58,59,60,26,63,27,64]);
%%

labels={'AFz';'F7';'F3';'Fz';'F4';'F8';'FC3';'FCz';'FC4';'T7';'C3';'Cz';'C4';'T8';'CP3';'CPz';'CP4';'P7';'P5';'P3';'P1';'Pz';'P2';'P4';'P6';'P8';'PO3';'PO4';'O1';'O2'};
X1=[79.0255388591416;49.8713779489202;57.5510633930990;60.7384809484625;57.5840261068105;49.9265268118817;30.9552849531915;32.9278836352560;30.9552849531915;5.17649253748256e-15;3.86812533613566e-15;5.20474889637625e-15;3.86788221025119e-15;5.17649253748256e-15;-30.9552849531915;-32.9278836352560;-30.9552849531915;-49.8713779489202;-54.0378881132511;-57.5510633930990;-59.9127302448179;-60.7384809484625;-59.8744127660117;-57.5840261068105;-54.0263340465386;-49.9265268118817;-76.1527667684845;-76.1527667684845;-80.7840137690914;-80.7840137690914];
Y=[0;68.4233350269540;48.2004273175388;0;-48.1425964684523;-68.3835902976096;59.2749781760892;0;-59.2749781760892;84.5385386396573;63.1712807125907;0;-63.1673101655785;-84.5385386396573;59.2749781760892;-4.03250272966127e-15;-59.2749781760892;68.4233350269539;63.0582218645482;48.2004273175389;26.0420933899754;-7.43831862786072e-15;-26.0254380421476;-48.1425964684523;-63.0447391225751;-68.3835902976096;31.4827967984807;-31.4827967984807;26.1330144040702;-26.1330144040702];
figure
subplot(1,2,1)
plottopomap(X1,Y,labels,X(1,:)')
title('CSP filter 1')
subplot(1,2,2)
plottopomap(X1,Y,labels,X(2,:)')
title('CSP filter 1')
%%
cv = cvpartition(length(TrainLabel), 'KFold', 4);
k_values = 1:15;
accuracy = zeros(length(k_values), 1);
% Perform cross-validation for different k values
for i = 1:length(k_values)
    current_k = k_values(i);
    
    % Initialize accuracy for the current k
    current_accuracy = 0;
    
    % Perform cross-validation
    for fold = 1:cv.NumTestSets
        train_idx = training(cv, fold);
        test_idx = test(cv, fold);
        
        % Train k-NN model
        TrainData1=TrainData(:,:,train_idx);
        TrainLabel1=TrainLabel(train_idx);
        class1=TrainData1(:,:,TrainLabel1==0);
        class2=TrainData1(:,:,TrainLabel1==1);
        C1=zeros(30);

        for kk=1:length(class1(1,1,:))
            sq_c1=squeeze(class1(:,:,kk));

           C1= C1+sq_c1*sq_c1';
        end
        C1=C1/length(class1(1,1,:));

        C2=zeros(30);

        for kk=1:length(class2(1,1,:))
            sq_c2=squeeze(class2(:,:,kk));
           C2= C2+sq_c2*sq_c2';
        end
        C2=C2/length(class2(1,1,:));
        [V1,D1] = eig(C1,C2);
        [sVals,sIndex] = sort(diag(D1),'descend') ;
        sV1 = V1(:,sIndex) ;
     
        features1=[];
        for nn=1:length(class1(1,1,:))
            Z1=sV1(:,[1:current_k,length(sV1(1,:))-current_k+1:length(sV1(1,:))])'*class1(:,:,nn);
            
            features1=[features1;diag(cov(Z1'))'];
        end
        features2=[];
        for nn=1:length(class2(1,1,:))
            Z1=sV1(:,[1:current_k,length(sV1(1,:))-current_k+1:length(sV1(1,:))])'*class2(:,:,nn);
            
            features2=[features2;diag(cov(Z1'))'];
        end
        features=[features1;features2];
        labels=[zeros(length(features1),1);ones(length(features2),1)];
        knn_model = fitcknn(features, labels, 'NumNeighbors', 3);
        TestData1=TrainData(:,:,test_idx);
        features3=[];
        for nn=1:length(TestData1(1,1,:))
            Z1=sV1(:,[1:current_k,length(sV1(1,:))-current_k+1:length(sV1(1,:))])'*TestData1(:,:,nn);
            %disp(size(Z1))
            features3=[features3;diag(cov(Z1'))'];
        end
        TestLabel1=TrainLabel(test_idx);
        % Predict on the test set
        y_pred = predict(knn_model, features3)';
        
        % Calculate accuracy for the current fold
        current_accuracy = current_accuracy + sum(y_pred == TestLabel1) / numel(TestLabel1);
    end
    
    % Average accuracy across folds for the current k
    accuracy(i) = current_accuracy / cv.NumTestSets;
end

% Find the best k (the one with the highest accuracy)
[best_accuracy, best_k_idx] = max(accuracy);
best_k = k_values(best_k_idx);

%%
class1=TrainData(:,:,TrainLabel==0);
class2=TrainData(:,:,TrainLabel==1);
features1=[];
for nn=1:length(class1(1,1,:))
    Z1=sV1(:,[1:best_k,length(sV1(1,:))-best_k+1:length(sV1(1,:))])'*class1(:,:,nn);

    features1=[features1;diag(cov(Z1'))'];
end
features2=[];
for nn=1:length(class2(1,1,:))
    Z1=sV1(:,[1:best_k,length(sV1(1,:))-best_k+1:length(sV1(1,:))])'*class2(:,:,nn);

    features2=[features2;diag(cov(Z1'))'];
end
features=[features1;features2];
labels=[zeros(length(features1),1);ones(length(features2),1)];
knn_model = fitcknn(features, labels, 'NumNeighbors', 3);
features3=[];
for nn=1:length(TestData(1,1,:))
    Z1=sV1(:,[1:best_k,length(sV1(1,:))-best_k+1:length(sV1(1,:))])'*TestData(:,:,nn);
    %disp(size(Z1))
    features3=[features3;diag(cov(Z1'))'];
end

TestLabel = predict(knn_model, features3)';
