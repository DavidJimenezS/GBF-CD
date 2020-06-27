%This code was made from the article 
% "Rayleigh-Rice mixture parameter estimation via EM algorithm for change 
% detection in multispectral images "

% if you used remember to cited bibtex here:
% @article{zanetti2015rayleigh,
%   title={Rayleigh-Rice mixture parameter estimation via EM algorithm for change detection in multispectral images},
%   author={Zanetti, Massimo and Bovolo, Francesca and Bruzzone, Lorenzo},
%   journal={IEEE Transactions on Image Processing},
%   volume={24},
%   number={12},
%   pages={5004--5016},
%   year={2015},
%   publisher={IEEE}
% }
%%
clear all
close all
clc

warning off

addpath(genpath(strcat(pwd,'\Data\')))
addpath(genpath(strcat(pwd,'\Functions\')))

%% Initialization
addpath(genpath('C:\Users\lenoco\Desktop\CD databases'));
%________________ lake data ___________

%load('mulargia_95_96.mat')
%load('fire_2013.mat')
%load('alaska_dataset')
%load('Madeirinha_dataset')
%load('katios_dataset')
load('dique_dataset')
%load('SF_dataset')
%load('Wenchuan_dataset')
%load('canada_dataset')
%load('california_flood')
%load('contest_dataset')
% 
% Aspot = (cropsaf);
% 
% I = (cropsbf);

I = double(before);
Aspot = double(after);
clear lake lakef before after;

%________________ contest data ___________
%load('data_contest.mat')


%_________________Magnitude of difference _______________


ro = sqrt((I(:) - Aspot(:)).^2);

ro(ro == 0) = eps;




%% Population of W1 and W2 for ML initial parameters

T = (max(ro) - min(ro))/2;

idx_w1 = ro <= T;
idx_w2 = ro > T;
W1 = ro(idx_w1);
W2 = ro(idx_w2);

change_map_rR = ro;
change_map_rR(idx_w1) = 0;
change_map_rR(idx_w2) = 1;


[m_c n_c] = size(Aspot);

figure, imshow(reshape(ro,m_c,n_c)), colorbar
figure, imshow(reshape(change_map_rR,m_c,n_c)), colorbar

N = length(ro);


alpha = length(W1)/N;
b_n = raylfit(W1);
pd = fitdist(W2,'Rician');
nu = pd.s;
sigma_c = pd.sigma;

%% Compute log-likelihood
p_x_w1 = rayleigh(ro,b_n);
p_x_w2 = rician(ro,sigma_c,nu);

%Solve problem with low probabilities round to zero in order to avoid
%Nan or Inf due to de fisr order modified Bessel function and
%log-likelihood
p_x_w1(p_x_w1 == 0) = eps;
p_x_w2(p_x_w2 == 0) = eps;
p_x_w2(isnan(p_x_w2)) = eps;
p_x_w2(isinf(p_x_w2)) = eps;

%joint distribution
p_ro_Psi = (alpha*p_x_w1) + ...
    ((1 - alpha)*p_x_w2);


clear idx_w1 idx_w2

%imagesc(reshape(p_ro_Psi,m_c,n_c))

log_likelihood(1) = sum(log(p_ro_Psi));

log_likeliold = 100;

figure;
hold on

Q_psi_psihat = zeros(N,2);
k = 1;

while abs(log_likelihood(k) - log_likeliold)/log_likeliold > 10e-8
    
    
    
    %%Expectation step
    p_x_w1 = rayleigh(ro,b_n);
    p_x_w2 = rician(ro,sigma_c,nu);
    
    %Solve problem with low probabilities round to zero
    p_x_w1(p_x_w1 == 0) = eps;
    p_x_w2(p_x_w2 == 0) = eps;
    p_x_w2(isnan(p_x_w2)) = eps;
    p_x_w2(isinf(p_x_w2)) = eps;
    
    
    p_w1_x_psi = (alpha*p_x_w1)./p_ro_Psi;
    p_w2_x_psi = ((1 - alpha)*p_x_w2)./p_ro_Psi;
    
    Q_psi_psihat(:,1) = p_w1_x_psi;%.*log(alpha*p_x_w1);
    Q_psi_psihat(:,2) = p_w2_x_psi;%.*log((1 - alpha)*p_x_w2);
    
    
    %% Maximitation step
    alpha_old = alpha;
    b_n_old = b_n;
    nu_old = nu;
    sigma_c_old = sigma_c;
    
    alpha = (1/N)*sum(p_w1_x_psi);
    
    b_n = sqrt((sum(p_w1_x_psi.*(ro.^2)))/(2*(sum(p_w1_x_psi))));
    
    I_0 = besseli(0,(ro.*nu_old)./(sigma_c_old^2));
    I_1 = besseli(1,(ro.*nu_old)./(sigma_c_old^2));
    nu = ((sum(p_w2_x_psi.*(I_1./I_0).*ro))/(sum(p_w2_x_psi)));
    
    aux = (ro.^2) + (nu_old.^2) - (2*ro.*nu_old.*(I_1./I_0));
    sigma_c = sqrt((sum(p_w2_x_psi.*aux))/(2*(sum(p_w2_x_psi))));
    
    
    %% Log-likelihood
    p_x_w1 = rayleigh(ro,b_n);
    p_x_w2 = rician(ro,sigma_c,nu);
    
    %Solve problem with low probabilities round to zero in order to avoid
    %Nan or Inf due to de fisr order modified Bessel function and
    %log-likelihood
    p_x_w1(p_x_w1 == 0) = eps;
    p_x_w2(p_x_w2 == 0) = eps;
    p_x_w2(isnan(p_x_w2)) = eps;
    p_x_w2(isinf(p_x_w2)) = eps;
    
    %joint distribution
    p_ro_Psi = (alpha*p_x_w1) + ...
        ((1 - alpha)*p_x_w2);
    
    log_likeliold = log_likelihood(k);
    
    k = k + 1
    
    log_likelihood(k) = sum(log(p_ro_Psi));
    
    plot(log_likelihood,'r','linewidth',3)
    title('Log-likelihood','FontSize',13)
    xlabel('Iterations','FontSize',13)
    ylabel('$\mathcal{L}(\mathbf{ \Psi }) = \sum log(p(\mathbf{X}|\mathbf{\Psi}))$','FontSize',13,'interpreter','latex')
    set(gca,'FontSize',12)
    grid on
    
    drawnow
    
end

%% BDR (Bayesian Decision Rule)
[val index] = max(Q_psi_psihat,[],2);
clear val;

change_map_rR = ro;
change_map_rR(index == 1) = 0; %No change
change_map_rR(index == 2) = 1; %Change
change_map_rR = reshape(change_map_rR,m_c,n_c);

figure, imshow(change_map_rR), colorbar

[counts ,centers] = hist(ro,100);
counts_n = counts./max(counts);
figure, bar(centers,counts_n)
hold on
plot(ro,p_ro_Psi./max(p_ro_Psi),'rx')
title('Joint distribution and Histogram of $\rho$','FontSize',13,'interpreter','latex')
xlabel('$\rho$','FontSize',13,'interpreter','latex')
ylabel('Normalized probabilities','FontSize',13,'interpreter','latex')
set(gca,'FontSize',12)
grid on
