function Change_map = gbf_cd(dataset,n)
% graph based data fusion for change detection
%
% FORMAT change_map = gbf_cd(dataset,n)
%
% dataset     - Contains an availabla dataset in string format or
%               a cell array with size 2 with before and after images
% n           - Number of sample nodes (pixels).
% Change_map  - Binary image of the change zone detected
%__________________________________________________________________________
% Copyright (C) 2019 Graph based data fusion
% David Alejandro Jimenez Sierra

%--------------------Available datasets-----------------------
%   sardania_dataset
%   omodeo_dataset
%   alaska_dataset
%   Madeirinha_dataset
%   katios_dataset
%   dique_dataset
%   SF_dataset
%   Wenchuan_dataset
%   canada_dataset
%   california_flood
%   contest_dataset
%-------------------------------------------------------------
%% Validation of inputs

if isa(dataset,'cell')
    before = dataset{1};
    after = dataset{2};
elseif isa(dataset,'string')
    load(dataset)
else
    msgbox('The dataset must be a cell by 1x2 with before and after images or a string with the name of an available dataset',...
        'Error','Error');
    return;
end

if isa(n,'double')
    if floor(n) ~= n
        msgbox('The number os samples n must be an integer',...
        'Error','Error');
        return;
    end
end
warning off

%% Data read and parameters

% sampling
addpath(genpath([pwd , '/RR_ImageSampling/']));
% MI
addpath(genpath([pwd , '/minf/']));
%Data
addpath(genpath([pwd , '/Data/']));
%Functions
addpath(genpath([pwd , '/Functions/']));


e_ditsAB = 3;

if isa(before,'uint8')
    before = double(before);
    after = double(after);
end

disp('Normalization of data')
% ------------ for lake & fire data----------
Aspot = after; %(lakef');%logt2_clipped(:,:,2) + abs(min(min(logt2_clipped(:,:,2))));%(lakef');
maxA = max(max(after));
Aspot = (Aspot./maxA);

I = before; %(lake'); % t1_L8_clipped(:,:,3) + abs(min(min(t1_L8_clipped(:,:,3))));%
maxA = max(max(before));
I = (I./maxA);

%----- Check for Nan values in the data
idx = find(isnan(I));
if ~isempty(idx)
    I(idx) = 0;
end

idx = find(isnan(Aspot));
if ~isempty(idx)
    Aspot(idx) = 0;
end

clear before after idx

clear logt2_clipped t1_L8_clipped
%%

[m_c n_c] = size(Aspot);

modalities = 2;

%Numeber of samples

Xl_AA = cell(modalities,1);
Xl = cell(modalities,1);

%% Extract random pixels from Images
%-----------if you prefer to get the locations manually--------------------

%imshow(Aspot(:,:,1),[]);

%[xi,yi] = getpts;
%locations = round([xi yi]);

% Manual
%locations = sub2ind([m_c  n_c],locations(:,2),locations(:,1)); %manual

% --------------Get n random locations (Automatically)--------------------
[ locations, v_j ] = Sampling_Grid(I(:,:,1), n, false );
%[ locations, v_j ] = Sampling_Jittered(I(:,:,1), n, true);
%[locations, v_j] = Sampling_Uniform(I(:,:,1), n, true );
%[ locations, v_bc ] = Sampling_BestCandidate(I(:,:,1), n, ceil(n/4), true );

%Autom
locations = sub2ind([m_c  n_c],locations(2,:),locations(1,:));


n = length(locations);

% Get values at those locations:

Xl{1} = (I(:,:,1)); %Ir
Xl{2} = (Aspot(:,:,1)); %Ig

Xl_AA{1} = Xl{1}(locations); %Ir
Xl_AA{2} = Xl{2}(locations); %Ig

% local graph
wl = cell(modalities,1);

complement = setdiff(1:(m_c*n_c), locations);


%------another metric taking into account the index of the pixel

% indexdisa = pdist2(locations',locations');
% max_distance = max(nonzeros(indexdisa));
% indexdisa = indexdisa./max_distance;
% indexistda = std(indexdisa(:));
% % %
% indexdisb = pdist2((complement)',locations');
% max_distance = max(nonzeros(indexdisb));
% indexdisb = indexdisb./max_distance;
% indexistdb = std(indexdisb(:));
%
% metrica = @(x,kernelstda) exp(-x ./ (2*(kernelstda ^ 2))).*...
%  exp(-indexdisa ./ (2*(indexistda ^ 2)));
% %
% metricb = @(x,kernelstda) exp(-x ./ (2*(kernelstda ^ 2))).*...
% exp(-indexdisb ./ (2*(indexistdb ^ 2)));

%-------------Kernel heat-------------------------
metric = @(x,kernelstda) exp(-(x.^2) ./ ((kernelstda ^ 2)));

disp('Computing graphs...')
for i = 1 : modalities
    
    %% L2 norm
    
    distlAA = pdist2(Xl_AA{i}',Xl_AA{i}','euclidean');
    distlAB = pdist2(Xl{i}(complement)',Xl_AA{i}','euclidean').^e_ditsAB;
    
    
    %% L1 norm
    %                 x=repmat(Xl_AA{i}',[size(Xl_AA{i}',1),1]);
    %                 y=repmat(Xl_AA{i},[size(Xl_AA{i},2),1]);
    %                 y = y(:);
    %                 distlAA = reshape(sum(abs(x-y),2).^2,n,n);
    %                 x=repmat(Xl_AA{i},[size(Xl{i}(complement),2),1]);
    %                 x = x(:);
    %                 y=repmat(Xl{i}(complement)',[size(Xl_AA{i}',1),1]);
    %                 distlAB = reshape(sum(abs(x-y),2).^2,length(complement),n);
    %
    %% Normalization and Degree
    
    D1 = distlAA*ones(n,1) + distlAB'*ones(length(complement),1);%diag(sum(distlAA, 2));
    distlAA = distlAA./repmat((D1),1,n);
    D2 = distlAB*ones(n,1) + (distlAB*pinv(distlAA))*(distlAB'*ones(length(complement),1));%sum(distlAB, 2);
    distlAB = distlAB./repmat((D2),1,n);
    
    
    clear x y D1 D2
    
    %% Metric computation
    
    %mean data
    kernelstd = mean(distlAB(:));
       
    distlAA = metric(distlAA,kernelstd);
    distlAB = metric(distlAB,kernelstd);
    
        wl{i} = [distlAA ;distlAB];
    
    clear distlAA distlAB distlAAsmooth
    
    disp(['graph ',num2str(i),' complete...'])
end
clear Xl Xl_AA
%% Multimodal Weights

n = length(locations);
WL = min(cat(3,wl{1} , wl{2}),[],3);

figure,
imagesc(WL), colorbar
title('Multimodal Weights')
set(gca,'FontSize',12)
drawnow

disp('Fusing graphs done...')


% One shot method of Nystrom
W_AA = WL(1:n,1:end);
W_BA = WL(n+1:end,1:end);


%% Nystrom aproximation

disp('Computing eigenvectors and eigenvalues...')

W_AA_sqrtinv = pinv(sqrtm(W_AA));%W_AA^-0.5;
%
S = W_AA + (W_AA_sqrtinv*(W_BA'*W_BA)*W_AA_sqrtinv);

[U_s,D_s] = eig(S);

Uhat_W = [ W_AA ; (W_AA_sqrtinv*W_BA')']*(U_s*pinv(sqrtm(D_s)));%%%[U_s; W_BA*U_s*(pinv(D_s))];%[W_AA;W_BA*W_AA_sqrtinv]*(U_s*pinv(sqrtm(D_s)));%%%

clear S W_AA_sqrtinv W_AA W_BA

figure; imshow(I(:,:,1),[])%, colorbar
title('Before')
set(gca,'FontSize',12)
%export_fig(strcat(pwd,'\Figures\bf_fire'),'-pdf','-transparent',h)
drawnow

figure; imshow((Aspot(:,:,1)),[])%, colorbar
title('After')
set(gca,'FontSize',12)
%export_fig(strcat(pwd,'\Figures\af_fire'),'-pdf','-transparent',h)
drawnow

%E = zeros(1,n);
MI = zeros(1,n);
prior1  = (Aspot(:,:,1) - I)./((Aspot(:,:,1) + I));
prior2  = (I - Aspot(:,:,1))./((Aspot(:,:,1) + I));
prior = imbinarize(prior1) + imbinarize(prior2);


%prior = imbinarize(prior);
%figure, imshow((prior(:,:,1)),[]), colorbar
%title('Prior')
%set(gca,'FontSize',12)

disp('Selecting the final change map by maximizing MI....')

for i = 2 : n
    Iaux = Uhat_W(:,i)*sqrt(D_s(i,i));
    
    A = Iaux(1:n);
    AB = Iaux(n+1:end);
    Iaux(locations) = A;
    Iaux(complement) = AB;
    Iaux = ((reshape(Iaux,m_c,n_c)));
    
    Iaux = imbinarize(abs(Iaux));
    
    %     figure, imshow(Iaux,[]),colorbar%, colormap hot
    %     title(['Eigenvector sample ' , num2str(i) ])
    %     set(gca,'FontSize',12)
    
    MI(i) =  mi(prior,Iaux);
    clear Iaux;
end

figure, plot(MI,'linewidth',2)
title('MI of eigenvectors')
xlabel('$u_{i} \sqrt{d_{i}}$','FontSize',13,'interpreter','latex')
ylabel('$MI(I_{u_{i}},I_{Prior})$','FontSize',13,'interpreter','latex')
set(gca,'FontSize',12)

i = find(MI == max(MI),1,'first');

hold on

plot(i,MI(i),'ro','MarkerSize',12)
drawnow

Iaux = Uhat_W(:,i)*sqrt(D_s(i,i));

A = Iaux(1:n);
AB = Iaux(n+1:end);
Iaux(locations) = A;
Iaux(complement) = AB;
Iaux = ((reshape(Iaux,m_c,n_c)));

Change_map = imbinarize(abs(Iaux));

figure, imshow(Change_map,[]);
title(['Eigenvector sample ' , num2str(i) ])
set(gca,'FontSize',12)
drawnow

%      saveas(h,strcat(pwd,'/Runs_images/omodeo_data_n_',num2str(traials_n),'.png'));
%      clearvars -except modalities m_c n_c Aspot I e_ditsAB traials_n
%      close all
end