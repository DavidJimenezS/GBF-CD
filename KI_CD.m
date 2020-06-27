%%
% This code is based on the Kittler– Illingsworth method
% please refer to the next bibtex:
%@article{kittler1986minimum,
%  title={Minimum error thresholding},
%  author={Kittler, Josef and Illingworth, John},
%  journal={Pattern recognition},
%  volume={19},
%  number={1},
%  pages={41--47},
%  year={1986},
%  publisher={Elsevier}
%}
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
%load('dique_dataset')
%load('SF_dataset')
%load('Wenchuan_dataset')
%load('canada_dataset')
load('california_flood')
%load('contest_dataset')



I = double(before);
Aspot = double(after);
clear lake lakef before after;

%________________ contest data ___________
%load('data_contest.mat')


%_________________Magnitude of difference _______________


ro = sqrt((I(:) - Aspot(:)).^2);

T = kittler(ro);

idx_w1 = ro <= T;
idx_w2 = ro > T;
W1 = ro(idx_w1);
W2 = ro(idx_w2);

change_map = ro;
change_map(idx_w1) = 0;
change_map(idx_w2) = 1;
[m_c n_c] = size(Aspot);
figure, imshow(reshape(ro,m_c,n_c)), colorbar
change_map_KI = reshape(change_map,m_c,n_c);
figure, imshow(change_map_KI), colorbar
