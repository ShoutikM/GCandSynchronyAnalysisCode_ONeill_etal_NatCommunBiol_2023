%%
clear all
clc
addpath('~/Code/');    

%%
load('MEAdata_BDNFdoseResp.mat');

% exp##_MEAinfo(##).div##.spikeMatrix specifies a particular recording
% Modify as necessary to select a different recording
% Spike Times are in msec
spkMat = (1e-3)*exp7_MEAinfo(5).div07.spikeMatrix;
bin_size = 0.300; % burstlets

preprocessing % run the preprocessing to obtain time series data from exp
clearvars -except spks;

fn_save = 'Synchrony_exp7_5_div07.mat';

valid_ele = [2:4,6:7,9:56,58:63];

spks = double(spks>=4);                    % burstlets
tmp = sum(spks,2);
neur_num = unique(find(sum(spks,2)>=30)); % Selection of most active units for tractability
n = spks(neur_num,:);

L = size(n,1);
L_star = 2^L - 1;
T = size(n,2);
alpha = 1e-3;
sigma = 0.05;

%% Data partition for cross-validation
tmp1 = downsample(n',2)'; N1 = size(tmp1,2);
tmp2 = downsample(n',2,1)'; N2 = size(tmp2,2); % Even-Odd

incl_idx = [];
X_star = [];
id_list = 1:L_star;
for idx = id_list
    bi_idx = de2bi(idx,L);
    idx0 = find(bi_idx==0); idx1 = find(bi_idx==1);
    
    tmp1 = prod(n(idx1,:),1);
    if sum(tmp1) == 0
        id_gr = id_list(id_list>=idx);
        for ii=1:numel(idx1)
            id_gr_bits(ii,1:numel(id_gr)) = bitget(id_gr, idx1(ii));
        end
        elim_id = find(sum(id_gr_bits)==numel(idx1));
        id_list( ismember(id_list, elim_id + idx-1) ) = [];
    else
        tmp = max(0, prod(n(idx1,:),1)-double(sum(n(idx0,:),1)>0))';
        if sum(tmp)>1
            X_star = [X_star tmp];
            incl_idx = [incl_idx idx];
        end        
    end
end

ord_idx = sum( de2bi(incl_idx, L) ,2);
pw2idx = find(ord_idx == 1);

X1_star = downsample(X_star,2);
X2_star = downsample(X_star,2,1);

if ~exist('n_star','var')
    n_star = X_star';
end
if ~exist('ng','var')
    ng = sum(n_star, 1);
end

%% Static Analysis
M = length(incl_idx);

% Gradient Descent
theta_star = zeros(M,1);
for gditer=1:2000
    grad = sum(X_star)' - T*( ( exp(theta_star) )/( 1+sum(exp(theta_star)) ) );
    theta_star = theta_star + alpha*grad;
end

lambdastar_est = exp(theta_star)/(1+sum(exp(theta_star)));
lambdagstar_est = sum(lambdastar_est);

%% Statistical Inference (Static)
% Set null nalues for out-of-support parameters
lambda_est = zeros(L,1);
bi_incl_idx = de2bi(incl_idx, L);
for l=1:L
    tmp = find( bi_incl_idx(:,l) );
    lambda_est(l) = sum( lambdastar_est(tmp) );
end

bi_incl_idx = de2bi(incl_idx, L);
odds = zeros(size(incl_idx))';
for l=1:length(incl_idx)
    k_tmp = find( bi_incl_idx(l,:) );
    odds(l) = prod( lambda_est(k_tmp) ) / prod( 1 - lambda_est(k_tmp) );
end

% Iterated Deviance-based Hypothesis Tests to find highest order of Synch
ord = 2:L;
j_stat_ord = zeros(size(ord));
Ms = zeros(size(ord));
Devs = zeros(size(ord));
nus = zeros(size(ord));

ordKidx = [];
for k = ord
    tmp = sort(find(ord_idx == k), 'ascend');
    ordKidx = 1:M; ordKidx(tmp) = [];
    if isempty(ordKidx)
        continue;
    end
    
    [h, pval, Jstat, Dev, Md, nu] = SynchTest_static(theta_star, X_star, odds, ordKidx, sigma, alpha);
    j_stat_ord(k) = h*(1-pval); % or Jstat
    Ms(k) = Md;
    Devs(k) = Dev;
    nus(k) = nu;
end

figure;
stem(ord, double(j_stat_ord(2:end)>0), 'r');
hold on; bar(ord, j_stat_ord(2:end), 'k'); ylim([0 1]); hold off

%%
% save(fn_save);
