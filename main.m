clear all
close all
clc

warning off

addpath(genpath(strcat(pwd,'\Data\')))
addpath(genpath(strcat(pwd,'\Functions\')))

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

dataset = "sardania_dataset";

%if you want to try your own dataset it must be like:
%   dataset{1} = before; image with one chanel
%   dataset{2} = after;  image with one chanel

n = 93;

Change_map = gbf_cd(dataset,n);

%% Metrics and error map

%all available datasets has the gt
load('sardania_dataset','gt');

[MA, FA, Precision, recall, kappa, OE] = cohensKappa(gt,Change_map);
