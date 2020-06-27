function [MA, FA, P, R, kappa, OE, label_map] = cohensKappa(y, yhat)
% Metrics and map error
%
% FORMAT [MA, FA, P, R, kappa, OE, label_map] = cohensKappa(y, yhat)
% 
% y             - Graound truth
% yhat          - Change_map predicted 
% MA            - Missed Alarms
% FA            - False Alarms
% P             - Precision
% R             - Recall
% kappa         - Cohen's Kappa coefficient
% OE            - Over all error
% label_map     - Error map (green values = correct
%                            Blue values  = MA
%                            Red values   = FA )    
%
%_______________________________________________________________________
% Copyright (C) 2019 
% David Alejandro Jimenez Sierra

%% Validation of inputs
    if ~isa(y,'logical') || ~isa(yhat,'logical')
        msgbox('The input data must be logical',...
        'Error','Error');
        return;
    end
%% 

    y_v = y(:);
    yhat_v = yhat(:);
    C = confusionmat(y_v, yhat_v); % compute confusion matrix
    
    %% computation of MA FA precision, recall and OE
    total_cd_pixels = length(find(y_v));
    total_ncd_pixels = length(y_v) - total_cd_pixels;
    
    MA = (C(2,1)/total_cd_pixels)*100;
    FA = (C(1,2)/total_ncd_pixels)*100;
    
    P = C(2,2)/(C(2,2) + C(2,1));%1 - FA/(total_cd_pixels - MA + FA);
    R = C(2,2)/(C(2,2) + C(1,2));%1 - (MA/total_cd_pixels);
    
    OE = ((C(1,2) + C(2,1))/sum(C(:)))*100;
    
    %% computation of kappa
    n = sum(C(:)); % get total N
    C = C./n; % Convert confusion matrix counts to proportion of n
    r = sum(C,2); % row sum
    s = sum(C); % column sum
    expected = r*s; % expected proportion for random agree
    po = sum(diag(C)); % Observed proportion correct
    pe = sum(diag(expected)); % Proportion correct expected
    kappa = (po-pe)/(1-pe); % Cohen's kappa
    
    table(MA, FA, P, R, kappa, OE)
    %% Plotting map error
    Err = y - yhat;
    Err(Err == -1) = 3; %FA
    Err(Err == 1) = 2; %MA
    aux_1 = y + yhat;
    Err(aux_1 == 2)= 1;

    label_map = label2rgb(Err,[0 1 0;0 0 1; 1 0 0],[1 1 1]);
    figure; imshow(label_map)
    % title('Error map using rR-EM','FontSize',13,'interpreter','latex')
    set(gca,'FontSize',12)
end