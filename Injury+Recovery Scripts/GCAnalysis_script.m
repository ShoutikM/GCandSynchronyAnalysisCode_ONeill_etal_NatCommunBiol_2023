clear all
clc
addpath('~/Code/');

%% Load Data
load('MEAdata_glu+BDNF.mat');

% exp##_MEAinfo(##).div##.spikeMatrix specifies a particular recording
% Modify as necessary to select a different recording
spkMat = (1e-3)*exp09_MEAinfo(10).div14.spikeMatrix;
bin_size = 0.300; % burstlets

preprocessing % run the preprocessing to obtain time series data from exp
clearvars -except spks;

fn_save = 'GC_exp09_10_div14.mat';

if ispc
    fn_save = strrep(fn_save, '/', '\');
end

valid_ele = [2:4,6:7,9:56,58:63];

spks = double(spks>=4);                    % burstlets
neur_num = unique(find(sum(spks,2)>=10));
n = spks(neur_num,:);

L = size(n,1);
T = size(n,2);

win_n = [1 2 4 5]; % burstlets - 300ms
p = length(win_n);
M = L*p+1;
alpha = 1e-3;

%% Form Covariate Matrices
tmp1 = downsample(n',2)'; N1 = size(tmp1,2);
tmp2 = downsample(n',2,1)'; N2 = size(tmp2,2);

theta_tmp = zeros(M,1);

X1 = getDesMat(L, tmp1, p, win_n);
X2 = getDesMat(L, tmp2, p, win_n);

L = size(n,1);
T = size(n,2);
X = getDesMat(L, n, p, win_n);

%% GC Analysis with OMP
GClink = zeros(1,L,L);
Dev = zeros(1,L,L);
Md = zeros(1,L,L);

theta_full = zeros(L,M);
for l=1:L
    % Full Model
    %     Cross-validation to determine s_full
	reg = floor(100*(sum(n(l,:))/T)*(1-(sum(n(l,:))/T))); reg=reg+double(reg==0);
    [~, LL1_tmp] = omp_cv(tmp1(l,:), X1, tmp2(l,:), X2, ceil(M/2), alpha/reg);
    [~, LL2_tmp] = omp_cv(tmp2(l,:), X2, tmp1(l,:), X1, ceil(M/2), alpha/reg);
    
    LL = LL1_tmp+LL2_tmp;
    [~, s_full] = max(LL');
    
    % Full model
    theta_full(l,:) = (omp(n(l,:), X, s_full, alpha/reg))';
    S_full = find(theta_full(l,:));
    
    % Reduced Models and GC computation
    cbar = 1:L; cbar(l) = [];
    for c = cbar
        tht = theta_full(l,:); tht(1+(c-1)*p+1:1+c*p) = [];
        s_red_max = length(find(tht));
        
        X1_red = X1; X1_red(:,1+(c-1)*p+1:1+c*p)=[];
        X2_red = X2; X2_red(:,1+(c-1)*p+1:1+c*p)=[];
        X_red = X; X_red(:,1+(c-1)*p+1:1+c*p)=[];
        
        % Cross-validation to determine s_red
        reg = floor(100*(sum(n(l,:))/T)*(1-(sum(n(l,:))/T))); reg=reg+double(reg==0);
        [~, LL1_tmp] = omp_cv(tmp1(l,:), X1_red, tmp2(l,:), X2_red, s_red_max, alpha/reg);
        [~, LL2_tmp] = omp_cv(tmp2(l,:), X2_red, tmp1(l,:), X1_red, s_red_max, alpha/reg);
        
        LL = LL1_tmp+LL2_tmp;
        [~, s_red] = max(LL');
        
        % Reduced Model
        theta_red = (omp(n(l,:), X_red, s_red, alpha/reg))';
        S_red = find(theta_red);
        
        % GC links and related computations
        ll_diff = ( n(l,:)*X*theta_full(l,:)' - sum(log(1+exp(X*theta_full(l,:)'))) )...
                - ( n(l,:)*X_red*theta_red' - sum(log(1+exp(X_red*theta_red'))) );
        
        Dev(1,l,c) = 2*ll_diff;
        
        Md(1,l,c) = length(S_full) - length(S_red);
    end
end

%%% Non-centrality param estimation
nu = max(Dev-Md, 0);

%%% Statistical Inference Test: FDR Control
% Formatting considerations for compatibility
wknft = cell(L,1);
for l=1:L
    tmp = theta_full(l,:);
    tmp(2+(l-1)*p:1+l*p)=[];
    tht=[tmp(1), theta_full(l,2+(l-1)*p:1+l*p), tmp(2:end)];
    wknft{l} = tht;
end
Whc = win_n;
% Benjamini-Hochburg-Yuketieli Test
sigma = 0.1; % FDR rate
[GCs,~] = FDRcontrolBY(Dev,nu,sigma,Md,wknft,Whc); % J-statistics
GCs = GCs .* (Md>0);

GCs = squeeze(GCs);
thresh = 0.98;
GCLinkWeights = GCs.*(abs(GCs)>=thresh); % Retain the most reliably detected links

%%
% save(fn_save);
