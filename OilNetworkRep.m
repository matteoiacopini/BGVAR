% This code implements the analysis in
% Casarin, R., Iacopini, M., Molina, G. Ter Horst, E., Espinasa, R., Sucre, C., and Rigobon, R. (2019).
% Multilayer Network Analysis of Oil Linkages
% Econometrics Journal, 2019.
%========================================================================

close all; clc; clear all;
% Add path of functions, data and figures
addpath('functions');
addpath('data');
addpath('results');
%-------------------------LOAD VARIABLES----------------------------------
load('data/ListVar.mat');
ShortID = ID(2:end);
clear ID Short
%-----------------------------LOAD DATA-----------------------------------
% dat1: oil prices (pt)
% dat2: oil production (yt)
% dat3: number of rigs (rt)
% dat4: economic world activity index build by Lutz Killian (dt).
load('data/AllData.mat')
pt = dat1(:,1);
yt = dat2;
rt = dat3;
dt = dat4(:,1);
differ = 1;
if (differ == 1)
   Data = [((log(pt(2:end,:)))-(log(pt(1:end-1,:))))*100,...
            dt(2:end,1),...
            ((yt(2:end,:))-(yt(1:end-1,:)))*100,...
            ((rt(2:end,:))-(rt(1:end-1,:)))*100];
else
   Data = [pt dt yt rt];
end


%% Main algorithm
%----------------------------PRELIMINARIES--------------------------------
ny  = size(Data,2);        % number of response variables
ny0 = 2;                   % common factors
ny1 = 12;                  % nodes in layer 1
ny2 = 12;                  % nodes in layer 2
nsimu     = 100;         % Gibbs sampler, number of Simulations
nsimuInit = 100;           % Initialization step, number of Simulations
sdz = 0;                   % 1. Standardize dataset  0. Demean dataset
lag = 1;                   % lag order

%-------------------SAMPLING MIN & MAR STRUCTURE -------------------------
% Initialization
disp('---------- Initialization ------------')
MAR0 = SAMPLE_BGMAR_DAG(Data, nsimuInit, sdz, lag, ny);
disp('MAR initialized')
MIN0 = SAMPLE_BGMIN_DAG(Data, nsimuInit, sdz, ny);
disp('MIN initialized')

% Estimation
disp('---------- Estimation ------------')
MAR = SAMPLE_BGMAR_DAGblocks(Data, nsimu, sdz, lag, ny0,ny1,ny2,MAR0.DAG);
disp('MAR estimated');
MIN = SAMPLE_BGMIN_DAGblocks(Data, nsimu, sdz, ny0,ny1,ny2,MIN0.DAG);
disp('MIN estimated');



%% Plots of the results
disp('---------- Plot and Save ------------');
figure(1)
heatmap((MAR.DAG-diag(diag(MAR.DAG))));
colormap([1 1 1;0.2 0.4 0.9]);
title('Lagged dependence (All layers): Adjacency matrix');

figure(2)
heatmap(max(MIN.DAG,MIN.DAG'));
colormap([1 1 1;0.2 0.4 0.9]);
title('Contemporaneous dependence (All layers): Adjacency matrix');

%%
figure(3)
Gmar = digraph((MAR.DAG-diag(diag(MAR.DAG))));
plot(Gmar)
title('Lagged dependence (All layers): Graph');

figure(4)
Gmin = graph(max(MIN.DAG,MIN.DAG'));
plot(Gmin)
title('Contemporaneous dependence (All layers): Graph');


%%
figure(5)
lb    = ShortID;
lb    = lb(ny0+1:ny0+ny1);
Gmar1 = digraph((MAR.DAG(ny0+1:ny0+ny1,ny0+1:ny0+ny1)-diag(diag(MAR.DAG(ny0+1:ny0+ny1,ny0+1:ny0+ny1)))),lb);
plot(Gmar1)
title('Lagged dependence (Production Layer): Graph');

figure(6)
lb    = ShortID;
lb    = lb(ny0+ny1+1:ny0+ny1+ny2);
Gmar2 = digraph((MAR.DAG(ny0+ny1+1:ny0+ny1+ny2,ny0+ny1+1:ny0+ny1+ny2)-...
                 diag(diag(MAR.DAG(ny0+ny1+1:ny0+ny1+ny2,ny0+ny1+1:ny0+ny1+ny2)))),lb);
plot(Gmar2)
title('Lagged dependence (Riggs Layer): Graph');

%%
figure(7)
lb = ShortID(ny0+1:ny0+ny1+ny2);  % ShortID;
for jjj=1:ny1
   lb{jjj} = ['Prod',lb{jjj}];
end
for jjj=1:ny2
   lb{ny1+jjj} = ['Rig',lb{ny1+jjj}];
end
A12 = [zeros(ny1,ny1),(MAR.DAG(ny0+1:ny0+ny1,ny0+ny1+1:ny0+ny1+ny2)-...
       diag(diag(MAR.DAG(ny0+1:ny0+ny1,ny0+ny1+1:ny0+ny1+ny2))));...
       zeros(ny2,ny1+ny2)];
Gmar12 = digraph(A12,lb);
plot(Gmar12)
title('Lagged dependence (Riggs to Oil Inter-Layer): Graph');


figure(8)
lb = ShortID(ny0+1:ny0+ny1+ny2);   % ShortID;
for jjj=1:ny1
   lb{jjj} = ['Prod',lb{jjj}];
end
for jjj=1:ny2
   lb{ny1+jjj} = ['Rig',lb{ny1+jjj}];
end
A21 = [zeros(ny1,ny1+ny2);...
       (MAR.DAG(ny0+ny1+1:ny0+ny1+ny2,ny0+1:ny0+ny1)-...
       diag(diag(MAR.DAG(ny0+ny1+1:ny0+ny1+ny2,ny0+1:ny0+ny1)))),zeros(ny2,ny2)];
Gmar21 = digraph(A21,lb);
plot(Gmar21)
title('Lagged dependence (Oil to Riggs Inter-Layer): Graph');



%% Save the results
if (differ == 1)
   suff = '_diff';
else
   suff = '';
end
a = clock;
fname = sprintf('Oil_%d_%d_%d_%d', a(2),a(3),a(4),a(5));
save(['results/' fname '_' suff 'New.mat']);


%%
w     = zeros(ny,ny);
elMAR = adj2edgeL(MAR.DAG,w);
elMIN = adj2edgeL(MIN.DAG,w);
% n = ny/2;

% n rigs and n production level
filename = ['results/GlabelNew',suff,'.csv'];
fileID   = fopen(filename,'w');
frm      = '%s, %s, %s\n';
fprintf(fileID,frm,'id','label','partition');
for i=1:size(ShortID,1)
   fprintf(fileID,frm,num2str(i),num2str(cell2mat(ShortID(i))),num2str((i-2<=0)*0+((i-2>0)&(i-2<=n))*1+(i-2>n)*2));
end
fclose(fileID);

%%
filename = ['results/GmarNew',suff,'.csv'];
fileID   = fopen(filename,'w');
frm      = '%s, %s, %s, %s, %s\n';
fprintf(fileID,frm,'Source','Target','Type','Id','Weight');
for i=1:size(elMAR,1)
   il = (elMAR(i,1)-2<=n)*(elMAR(i,2)-2<=n)*1+(elMAR(i,1)-2<=n)*(elMAR(i,2)-2>n)*2+...
        (elMAR(i,1)-2>n)*(elMAR(i,2)-2<=n)*3+(elMAR(i,1)-2>n)*(elMAR(i,2)-2>n)*4;
   fprintf(fileID,frm,num2str(elMAR(i,1)),num2str(elMAR(i,2)),'Directed',num2str(i),num2str(il));
end
fclose(fileID);

%%
filename = ['results/GminNew',suff,'.csv'];
fileID   = fopen(filename,'w');
frm      = '%s, %s, %s, %s, %s\n';
fprintf(fileID,frm,'Source','Target','Type','Id','Weight');
for i=1:size(elMIN,1)
   il = (elMIN(i,1)-2<=n)*(elMIN(i,2)-2<=n)*1+(elMIN(i,1)-2<=n)*(elMIN(i,2)-2>n)*2+...
        (elMIN(i,1)-2>n)*(elMIN(i,2)-2<=n)*3+(elMIN(i,1)-2>n)*(elMIN(i,2)-2>n)*4;
   fprintf(fileID,frm,num2str(elMIN(i,1)),num2str(elMIN(i,2)),'Undirected',num2str(i),num2str(il));
end
fclose(fileID);
