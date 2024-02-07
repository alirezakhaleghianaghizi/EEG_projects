clc;
clear all;
cd('C:\Users\ASUS\OneDrive\Desktop\uni\simisterm7\eeg\Comp_HW2');
EEG=load('EX1.mat');
EEG_sig=EEG.EEG_Sig;
%%
f=200;
[N K] = size(EEG_sig);
t=(1:K)/f;

for i = 1:size(EEG_sig, 1)
    subplot(3,1,i)
    plot(t,EEG_sig(i, :));
    xlabel('t');
    ylabel('Values');
    title(['chanel ' num2str(i) ' of EEG' ]);
    legend(['ch ' num2str(i)]);
end
%%
figure
scatter3(EEG_sig(1,:),EEG_sig(2,:),EEG_sig(3,:),".",'b');
cx=(EEG_sig-mean(EEG_sig')')*(EEG_sig-mean(EEG_sig')')'/length(EEG_sig(1,:));
[V D]=eig(cx);
coeff = V;

hold on;
quiver3(0, 0, 0, D(1, 1) * V(1, 1),  D(1, 1) * V(2, 1),  D(1, 1) * V(3, 1), 'y', 'LineWidth', 2);
quiver3(0, 0, 0, D(2, 2) * V(1, 2),  D(2, 2) * V(2, 2), D(2, 2) * V(3, 2), 'g', 'LineWidth', 2);
quiver3(0, 0, 0, D(3, 3) * V(1, 3), D(3, 3) * V(2, 3), D(3, 3)* V(3, 3), 'r', 'LineWidth', 2);

hold off;
%%

D_d=diag(diag(D).^(-0.5))*V';
%EEG_y=zscore(D_d*EEG_sig);
EEG_y=(D_d*EEG_sig);
figure
for i = 1:size(EEG_y, 1)
    subplot(3,1,i)
    plot(t,EEG_y(i, :));
    xlabel('t');
    ylabel('Values');
    title(['chanel ' num2str(i) ' of EEG' ]);
    legend(['ch ' num2str(i)]);
end
figure
scatter3(EEG_y(1,:),EEG_y(2,:),EEG_y(3,:),'b');
%%
cy=(EEG_y-mean(EEG_y')')*(EEG_y-mean(EEG_y')')'/length(EEG_y(1,:));
%%
figure
scatter3(EEG_sig(1,:),EEG_sig(2,:),EEG_sig(3,:),".",'b');
[coeff,~,latent]=pca(EEG_sig');
V=coeff;
D=diag(latent);
hold on;
quiver3(0, 0, 0, D(1, 1) * V(1, 1),  D(1, 1) * V(2, 1),  D(1, 1) * V(3, 1), 'y', 'LineWidth', 2);
quiver3(0, 0, 0, D(2, 2) * V(1, 2),  D(2, 2) * V(2, 2), D(2, 2) * V(3, 2), 'g', 'LineWidth', 2);
quiver3(0, 0, 0, D(3, 3) * V(1, 3), D(3, 3) * V(2, 3), D(3, 3)* V(3, 3), 'r', 'LineWidth', 2);

hold off;
%%

D_d=diag(diag(D).^(-0.5))*V';
%EEG_y=zscore(D_d*EEG_sig);
EEG_y=(D_d*EEG_sig);
figure
for i = 1:size(EEG_y, 1)
    subplot(3,1,i)
    plot(t,EEG_y(i, :));
    xlabel('t');
    ylabel('Values');
    title(['chanel ' num2str(i) ' of EEG' ]);
    legend(['ch ' num2str(i)]);
end
figure
scatter3(EEG_y(1,:),EEG_y(2,:),EEG_y(3,:),'b');
%%
cy=(EEG_y-mean(EEG_y')')*(EEG_y-mean(EEG_y')')'/length(EEG_y(1,:));
%%

[U, S, V] = svd(cx);0

% Display matrices and related concepts
disp('Left Singular Vectors (U):');
%disp(U);

disp('Singular Values (S):');
%disp(S);

disp('Right Singular Vectors (V):');
%disp(V);

% Compute explained variance
explained_variance = diag(S).^2 / sum(diag(S).^2);
disp('Explained Variance:');
disp(explained_variance);
%%
svd_normalized=(diag(diag(S(:,1:3)).^(-0.5))*U'*EEG_sig);
figure
scatter3(svd_normalized(1,:),svd_normalized(2,:),svd_normalized(3,:),".",'b');

%%
cy=(svd_normalized-mean(svd_normalized')')*(svd_normalized-mean(svd_normalized')')'/length(svd_normalized(1,:));

%%
[coeff,~,latent]=pca(svd_normalized');
V=coeff;
D=diag(latent);
hold on;
quiver3(0, 0, 0, D(1, 1) * V(1, 1),  D(1, 1) * V(2, 1),  D(1, 1) * V(3, 1), 'y', 'LineWidth', 2);
quiver3(0, 0, 0, D(2, 2) * V(1, 2),  D(2, 2) * V(2, 2), D(2, 2) * V(3, 2), 'g', 'LineWidth', 2);
quiver3(0, 0, 0, D(3, 3) * V(1, 3), D(3, 3) * V(2, 3), D(3, 3)* V(3, 3), 'r', 'LineWidth', 2);

hold off;