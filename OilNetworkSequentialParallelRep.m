%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is a parallel implementation of the sequential estimation analysis in
% Casarin, R., Iacopini, M., Molina, G. Ter Horst, E., Espinasa, R., Sucre, C., and Rigobon, R. (2019).
% Multilayer Network Analysis of Oil Linkages
% Econometrics Journal, 2019.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clc; clear;
% Add path of functions, data and figures
addpath('functions');
addpath('data');
addpath('results');
%-------------------------LOAD VARIABLES----------------------------------
load('data/ListVar.mat');
ShortID=ID(2:end);
clear ID Short
%-----------------------------LOAD DATA-----------------------------------
%dat1: oil prices (pt)
%dat2: oil production (yt)
%dat3: number of rigs (rt)
%dat4: economic world activity index build by Lutz Killian (dt).
load('data/AllData.mat')
pt=dat1(:,1);
yt=dat2;
rt=dat3;
dt=dat4(:,1);
differ=1;
if differ==1
    Data=[((log(pt(2:end,:)))-(log(pt(1:end-1,:))))*100,...
        dt(2:end,1),...
        ((yt(2:end,:))-(yt(1:end-1,:)))*100,...
        ((rt(2:end,:))-(rt(1:end-1,:)))*100];
else
    Data=[pt dt yt rt];
end
%%
%----------------------------PRELIMINARIES--------------------------------
ny = size(Data,2);        % number of response variables
ny0=2;                        % common factors
ny1=12;                      % nodes in layer 1
ny2=12;                      % nodes in layer 2
nsimu  = 10000;           % Gibbs sampler, number of Simulation
nsimuInit  = 100;          % Initialization step, number of Simulation
sdz = 0;                       % 1. Standardize dataset  0. Demean dataset
lag=1;
%%%%%%%%%%
wnd=60;
T=size(Data,1);
MARall=cell(T-wnd,1);
MINall=cell(T-wnd,1);
DataSub=NaN(T-wnd,wnd,ny);
for t=wnd+1:T
    DataSub(t-wnd,:,:)=Data(t-wnd+1:t,:);
end
obj=parpool(6);
parfor t=wnd+1:T
    disp(['########## Iteration: ',num2str(t-wnd),' out of ',num2str(T-wnd)]);
    DataS=squeeze(DataSub(t-wnd,:,:));
    [MAR,MIN]=Gibbs(DataS, nsimuInit, nsimu, sdz, lag, ny,ny0,ny1,ny2);
    MARall{t-wnd}=MAR.DAG;
    MINall{t-wnd}=MIN.DAG;
end
obj.delete;

%-------------------SAVE RESULTS ----------------------------------------
%%
fname = sprintf('Oil_Sequential');
save(['results/' fname 'New.mat']);
%%
MARInDg=zeros(T-wnd,ny-2);
MAROuDg=zeros(T-wnd,ny-2);
MINDg=zeros(T-wnd,ny-2);

MARInDgY=zeros(T-wnd,n);
MAROuDgY=zeros(T-wnd,n);
MINDgY=zeros(T-wnd,n);

MARInDgR=zeros(T-wnd,n);
MAROuDgR=zeros(T-wnd,n);
MINDgR=zeros(T-wnd,n);

MARInDgYR=zeros(T-wnd,n);
MAROuDgYR=zeros(T-wnd,n);
MINDgYR=zeros(T-wnd,n);

MARInDgP=zeros(T-wnd,1);
MAROuDgP=zeros(T-wnd,1);
MINDgP=zeros(T-wnd,1);

MARInDgYP=zeros(T-wnd,1);
MAROuDgYP=zeros(T-wnd,1);
MINDgYP=zeros(T-wnd,1);

MARInDgRP=zeros(T-wnd,1);
MAROuDgRP=zeros(T-wnd,1);
MINDgRP=zeros(T-wnd,1);

MARInDgD=zeros(T-wnd,1);
MAROuDgD=zeros(T-wnd,1);
MINDgD=zeros(T-wnd,1);

MARInDgYD=zeros(T-wnd,1);
MAROuDgYD=zeros(T-wnd,1);
MINDgYD=zeros(T-wnd,1);

MARInDgRD=zeros(T-wnd,1);
MAROuDgRD=zeros(T-wnd,1);
MINDgRD=zeros(T-wnd,1);


for t=wnd+1:T
    % overall in and out excluding price and ec activity
    MARInDg(t-wnd,:)=sum(MARall{t-wnd}(3:26,3:26));
    MAROuDg(t-wnd,:)=sum(MARall{t-wnd}(3:26,3:26)');
    MINDg(t-wnd,:)=sum(MINall{t-wnd}(3:26,3:26));
    % in and out within OIL excluding price and ec activity
    MARInDgY(t-wnd,:)=sum(MARall{t-wnd}(3:14,3:14));
    MAROuDgY(t-wnd,:)=sum(MARall{t-wnd}(3:14,3:14)');
    MINDgY(t-wnd,:)=sum(MINall{t-wnd}(3:14,3:14));
    % in and out within RIG excluding price and ec activity
    MARInDgR(t-wnd,:)=sum(MARall{t-wnd}(15:26,15:26));
    MAROuDgR(t-wnd,:)=sum(MARall{t-wnd}(15:26,15:26)');
    MINDgR(t-wnd,:)=sum(MINall{t-wnd}(15:26,15:26));
    % in and out for RIG and OIL
    MARInDgYR(t-wnd,:)=sum(MARall{t-wnd}(3:14,15:26));
    MAROuDgYR(t-wnd,:)=sum(MARall{t-wnd}(15:26,3:14)');
    MINDgYR(t-wnd,:)=sum(MINall{t-wnd}(3:14,15:26));
    % in and out for PRICE
    MARInDgP(t-wnd,:)=sum(MARall{t-wnd}(1,3:26));
    MAROuDgP(t-wnd,:)=sum(MARall{t-wnd}(3:26,1)');
    MINDgP(t-wnd,:)=sum(MINall{t-wnd}(1,3:26));
    % in and out for PRICE and OIL
    MARInDgYP(t-wnd,:)=sum(MARall{t-wnd}(3:14,1));
    MAROuDgYP(t-wnd,:)=sum(MARall{t-wnd}(1,3:14)');
    MINDgYP(t-wnd,:)=sum(MINall{t-wnd}(3:14,1));
    % in and out for PRICE and RIG
    MARInDgRP(t-wnd,:)=sum(MARall{t-wnd}(15:26,1));
    MAROuDgRP(t-wnd,:)=sum(MARall{t-wnd}(1,15:26)');
    MINDgRP(t-wnd,:)=sum(MINall{t-wnd}(15:26,1));
    % in and out for EC ACTIVITY
    MARInDgD(t-wnd,:)=sum(MARall{t-wnd}(2,3:26));
    MAROuDgD(t-wnd,:)=sum(MARall{t-wnd}(3:26,2)');
    MINDgD(t-wnd,:)=sum(MINall{t-wnd}(2,3:26));
    % in and out for EC ACTIVITY and OIL
    MARInDgYD(t-wnd,:)=sum(MARall{t-wnd}(3:14,2));
    MAROuDgYD(t-wnd,:)=sum(MARall{t-wnd}(2,3:14)');
    MINDgYD(t-wnd,:)=sum(MINall{t-wnd}(3:14,2));
    % in and out for EC ACTIVITY and RIG
    MARInDgRD(t-wnd,:)=sum(MARall{t-wnd}(15:26,2));
    MAROuDgRD(t-wnd,:)=sum(MARall{t-wnd}(2,15:26)');
    MINDgRD(t-wnd,:)=sum(MINall{t-wnd}(15:26,2));
end
MARToDg=MARInDg+MAROuDg;
MARToDgY=MARInDgY+MAROuDgY;
MARToDgR=MARInDgR+MAROuDgR;
MARToDgP=MARInDgP+MAROuDgP;
MARToDgD=MARInDgD+MAROuDgD;
MARToDgYR=MARInDgYR+MAROuDgYR;
MARToDgYP=MARInDgYP+MAROuDgYP;
MARToDgRP=MARInDgRP+MAROuDgRP;
MARToDgYD=MARInDgYD+MAROuDgYD;
MARToDgRD=MARInDgRD+MAROuDgRD;
%%%
MINToDg=MINDg;
MINToDgY=MINDgY;
MINToDgR=MINDgR;
MINToDgP=MINDgP;
MINToDgD=MINDgD;
MINToDgYR=MINDgYR;
MINToDgYP=MINDgYP;
MINToDgRP=MINDgRP;
MINToDgYD=MINDgYD;
MINToDgRD=MINDgRD;
%%
sm=1;
p=3;
lab={'2004:12 ','2006:08 ','2008:04 ','2009:12 ','2011:08 ','2013:04 ','2014:12 ','2016:08 '};
n2=2*n;
fs=10;

%%%
figure(1)
plot((wnd+1:T),smth(sum(MARToDg,2)/(n2*(n2-1)),sm,p),'k-');
hold on;
plot((wnd+1:T),smth(sum(MINToDg,2)/(n2*(n2-1)),sm,p),'k--');
hold off;
leg=legend('$(Y_t,R_t)\leftarrow (Y_{t-1},R_{t-1})$:  Lagged Total Degree','$(Y_t,R_t)\leftrightarrow (Y_{t},R_{t})$: Contemporaneous Total Degree');
set(leg,'Interpreter','Latex');
xlim([wnd T]);
ylim([0 1]);
set(gca,'XTick',[60 80 100 120 140 160 180 200]);
set(gca,'XTickLabel',lab);
set(gca,'FontSize',fs);
%%%

figure(2)
%plot((wnd+1:T),smth(sum(MARToDg,2)/(n2*(n2-1)),sm,p),'k-');
plot((wnd+1:T),smth(sum(MARToDgY,2)/(n*(n-1)),sm,p),'k-','Marker','o');
hold on;
plot((wnd+1:T),smth(sum(MARToDgR,2)/(n*(n-1)),sm,p),'k-','Marker','o','MarkerFaceColor','k');
%plot((wnd+1:T),smth(sum(MARToDgYR,2)/(n*n),sm,p),'k-','Marker','s');
hold off;
leg=legend('$Y_t\leftarrow Y_{t-1}$: Oil Lagged Total Degree','$R_t\leftarrow R_{t-1}$: Rig Lagged Total Degree','location','SouthEast');
set(leg,'Interpreter','Latex');
xlim([wnd T]);
ylim([0 1]);
set(gca,'XTick',[60 80 100 120 140 160 180 200]);
set(gca,'XTickLabel',lab);
set(gca,'FontSize',fs);
%%%

figure(22)
plot((wnd+1:T),smth(sum(MARToDgYR,2)/(n*n*2),sm,p),'k-','Marker','s');
hold on;
plot((wnd+1:T),smth(sum(MARInDgYR,2)/(n*n*2),sm,p),'k-');
plot((wnd+1:T),smth(sum(MAROuDgYR,2)/(n*n*2),sm,p),'k--');
hold off;
leg=legend('$Y_t \leftarrow R_{t-1}$, $R_t \leftarrow Y_{t-1}$: Oil-Rig Lagged Total Degree','$Y_t \leftarrow R_{t-1}$: Oil-Rig Lagged In-Degree','$R_t \leftarrow Y_{t-1}$: Oil-Rig Lagged Out-Degree','location','SouthEast');
set(leg,'Interpreter','Latex');
xlim([wnd T]);
ylim([0 0.40]);
set(gca,'XTick',[60 80 100 120 140 160 180 200]);
set(gca,'XTickLabel',lab);
set(gca,'FontSize',fs);
%%%

figure(3)
%plot((wnd+1:T),smth(sum(MINToDg,2)/(n2*(n2-1)),sm,p),'k-');
plot((wnd+1:T),smth(sum(MINToDgY,2)/(n*(n-1)),sm,p),'k-','Marker','o');
hold on;
plot((wnd+1:T),smth(sum(MINToDgR,2)/(n*(n-1)),sm,p),'k-','Marker','o','MarkerFaceColor','k');
plot((wnd+1:T),smth(sum(MINToDgYR,2)/(n*n*2),sm,p),'k-','Marker','s');
hold off;
leg=legend('$Y_t \leftrightarrow Y_t$ Oil Contemporaneous Total Degree','$R_t \leftrightarrow R_{t}$: Rig Contemporaneous Total Degree','$R_t \leftrightarrow Y_t$: Oil-Rig Contemporaneous Total Degree');
set(leg,'Interpreter','Latex');
xlim([wnd T]);
ylim([0 0.40]);
set(gca,'XTick',[60 80 100 120 140 160 180 200]);
set(gca,'XTickLabel',lab);
set(gca,'FontSize',fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Price %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(100)
plot((wnd+1:T),smth(sum(MARToDgP,2)/(2*n2),sm,p),'k-');
hold on;
plot((wnd+1:T),smth(sum(MARToDgYP,2)/(2*n),sm,p),'k-','Marker','o');
plot((wnd+1:T),smth(sum(MARToDgRP,2)/(2*n),sm,p),'k-','Marker','o','MarkerFaceColor','k');
hold off;
leg=legend('$P_t \leftarrow (Y_{t-1},R_{t-1})$, $(Y_{t},R_{t}) \leftarrow P_{t-1}$: Price Lagged Total Degree','$P_t \leftarrow Y_{t-1}$, $Y_t \leftarrow P_{t-1}$: Price-Oil Lagged Total Degree','$P_t \leftarrow R_{t-1}$, $R_t \leftarrow P_{t-1}$: Price-Rig Lagged Total Degree');
set(leg,'Interpreter','Latex');
xlim([wnd T]);
ylim([0 0.8]);
set(gca,'XTick',[60 80 100 120 140 160 180 200]);
set(gca,'XTickLabel',lab);
set(gca,'FontSize',fs);

%%%
figure(101)
plot((wnd+1:T),smth(sum(MARToDgYP,2)/(2*n),sm,p),'k-','Marker','o');
hold on;
plot((wnd+1:T),smth(sum(MAROuDgYP,2)/(2*n),sm,p),'k-');
plot((wnd+1:T),smth(sum(MARInDgRP,2)/(2*n),sm,p),'k--');
hold off;
leg=legend('$P_t \leftarrow Y_{t-1}$, $Y_t \leftarrow P_{t-1}$: Price-Oil Lagged Total Degree','$Y_t \leftarrow P_{t-1}$: Price-Oil Lagged Out-Degree','$P_t \leftarrow Y_{t-1}$: Price-Oil Lagged In-Degree');
set(leg,'Interpreter','Latex');
xlim([wnd T]);
ylim([0 0.8]);
set(gca,'XTick',[60 80 100 120 140 160 180 200]);
set(gca,'XTickLabel',lab);
set(gca,'FontSize',fs);

%%%
figure(102)
plot((wnd+1:T),smth(sum(MARToDgRP,2)/(2*n),sm,p),'k-','Marker','o','MarkerFaceColor','k');
hold on;
plot((wnd+1:T),smth(sum(MAROuDgRP,2)/(2*n),sm,p),'k-');
plot((wnd+1:T),smth(sum(MARInDgRP,2)/(2*n),sm,p),'k--');
hold off;
leg=legend('$P_t \leftarrow R_{t-1}$, $R_t \leftarrow P_{t-1}$: Price-Rig Lagged Total Degree','$R_t \leftarrow P_{t-1}$: Price-Rig Lagged Out-Degree','$P_t \leftarrow R_{t-1}$: Price-Rig Lagged In-Degree');
set(leg,'Interpreter','Latex');
xlim([wnd T]);
ylim([0 0.8]);
set(gca,'XTick',[60 80 100 120 140 160 180 200]);
set(gca,'XTickLabel',lab);
set(gca,'FontSize',fs);
%%%
lab={'2004:12','2006:08','2008:04','2009:12','2011:08','2013:04','2014:12','2016:08'};
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EcoAct %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(103)
plot((wnd+1:T),smth(sum(MARToDgD,2)/(2*n2),sm,p),'k-');
hold on;
plot((wnd+1:T),smth(sum(MARToDgYD,2)/(2*n),sm,p),'k-','Marker','o');
plot((wnd+1:T),smth(sum(MARToDgRD,2)/(2*n),sm,p),'k-','Marker','o','MarkerFaceColor','k');
hold off;
leg=legend('$D_t \leftarrow (Y_{t-1},R_{t-1})$, $(Y_{t},R_{t}) \leftarrow D_{t-1}$: EcoAct Lagged Total Degree','$D_t \leftarrow Y_{t-1}$, $Y_t \leftarrow D_{t-1}$: EcoAct-Oil Lagged Total Degree','$D_t \leftarrow R_{t-1}$, $R_t \leftarrow D_{t-1}$: EcoAct-Rig Lagged Total Degree');
set(leg,'Interpreter','Latex');
xlim([wnd T]);
ylim([0 0.8]);
set(gca,'XTick',[60 80 100 120 140 160 180 200]);
set(gca,'XTickLabel',lab);
set(gca,'FontSize',fs);

%%%
figure(104)
plot((wnd+1:T),smth(sum(MARToDgYD,2)/(2*n),sm,p),'k-','Marker','o');
hold on;
plot((wnd+1:T),smth(sum(MAROuDgYD,2)/(2*n),sm,p),'k-');
plot((wnd+1:T),smth(sum(MARInDgRD,2)/(2*n),sm,p),'k--');
hold off;
leg=legend('$D_t \leftarrow Y_{t-1}$, $Y_t \leftarrow D_{t-1}$: EcoAct-Oil Lagged Total Degree','$Y_t \leftarrow D_{t-1}$: EcoAct-Oil Lagged Out-Degree','$D_t \leftarrow Y_{t-1}$: EcoAct-Oil Lagged In-Degree');
set(leg,'Interpreter','Latex');
xlim([wnd T]);
ylim([0 0.8]);
set(gca,'XTick',[60 80 100 120 140 160 180 200]);
set(gca,'XTickLabel',lab);
set(gca,'FontSize',fs);

%%%
figure(105)
plot((wnd+1:T),smth(sum(MARToDgRD,2)/(2*n),sm,p),'k-','Marker','o','MarkerFaceColor','k');
hold on;
plot((wnd+1:T),smth(sum(MAROuDgRD,2)/(2*n),sm,p),'k-');
plot((wnd+1:T),smth(sum(MARInDgRD,2)/(2*n),sm,p),'k--');
hold off;
leg=legend('$D_t \leftarrow R_{t-1}$, $R_t \leftarrow D_{t-1}$: EcoAct-Rig Lagged Total Degree','$R_t \leftarrow D_{t-1}$: EcoAct-Rig Lagged Out-Degree','$D_t \leftarrow R_{t-1}$: EcoAct-Rig Lagged In-Degree');
set(leg,'Interpreter','Latex');
xlim([wnd T]);
ylim([0 0.8]);
set(gca,'XTick',[60 80 100 120 140 160 180 200]);
set(gca,'XTickLabel',lab);
set(gca,'FontSize',fs);

%%%
lab={'2004:12','2006:08','2008:04','2009:12','2011:08','2013:04','2014:12','2016:08'};
%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% data price and ec activity
figure(400)
plot((wnd+1:T),Data((wnd+1:T),1),'k-');
hold on;
plot((wnd+1:T),Data((wnd+1:T),2),'k--');
hold off;
leg=legend('$P_t$','$D_t$');
set(leg,'Interpreter','Latex');
xlim([wnd T]);
set(gca,'XTick',[60 80 100 120 140 160 180 200]);
set(gca,'XTickLabel',lab);
set(gca,'FontSize',fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Save file sequence
%%
%%
for j=1:size(MARall,1)
    A=MARall{j}-diag(diag(MARall{j}));
    adj=(A>0);
    edges=find(adj>0);
    m=size(A,1);
    [ind,jnd]=ind2sub([m,m],edges);
    elMAR=[ind,jnd];
    %%%
    filename = ['results/GmarSeqCom',num2str(j),'.csv'];
    fileID = fopen(filename,'w');
    frm='%s, %s, %s, %s, %s\n';
    fprintf(fileID,frm,'Source','Target','Type','Id','Weight');
    for i=1:size(ind)
        il=(elMAR(i,1)-2<=n)*(elMAR(i,2)-2<=n)*1+(elMAR(i,1)-2<=n)*(elMAR(i,2)-2>n)*2+...
            (elMAR(i,1)-2>n)*(elMAR(i,2)-2<=n)*3+(elMAR(i,1)-2>n)*(elMAR(i,2)-2>n)*4;
        fprintf(fileID,frm,num2str(elMAR(i,1)),num2str(elMAR(i,2)),'Directed',num2str(i),num2str(il));
    end
    fclose(fileID);
    %%%
end
%%
cent=zeros(size(MARall,1),5);
mes=zeros(size(MARall,1),4);
for j=1:size(MARall,1)
    A=MARall{j}-diag(diag(MARall{j}));
    adj=(A>0);
    C=digraph(adj);
    mes(j,1)=max(centrality(C,'hubs'));
    mes(j,2)=max(centrality(C,'authorities'));
    %
    B=graph((adj+adj')/2);
    cent(j,1)=quantile(centrality(B,'degree'),0.95);
    cent(j,2)=quantile(centrality(B,'degree'),0.05);
    cent(j,3)=mean(centrality(B,'degree'));
    cent(j,4)=max(centrality(B,'eigenvector'));
    cent(j,5)=min(centrality(B,'eigenvector'));
end

figure(400)
plot((wnd+1:T),smth(cent(:,3),sm,p),'k-');
hold on;
plot((wnd+1:T),smth(cent(:,1),sm,p),'k--');
plot((wnd+1:T),smth(cent(:,2),sm,p),'k--');
hold off;
leg=legend('$d_t$','$q_{0.975,t}$','$q_{0.05,t}$');
set(leg,'Interpreter','Latex');
xlim([wnd T]);
set(gca,'XTick',[60 80 100 120 140 160 180 200]);
set(gca,'XTickLabel',lab);
set(gca,'FontSize',fs);

figure(401)
plot((wnd+1:T),smth(mes(:,1),sm,p),'k-');
hold on;
plot((wnd+1:T),smth(mes(:,2),sm,p),'k--');
hold off;
leg=legend('$h_t$','$a_t$');
set(leg,'Interpreter','Latex');
xlim([wnd T]);
ylim([0.04,0.1]);
set(gca,'XTick',[60 80 100 120 140 160 180 200]);
set(gca,'XTickLabel',lab);
set(gca,'FontSize',fs);

%%
figure(402)
plot((wnd+1:T),smth(cent(:,4),sm,p),'k-');
hold on;
plot((wnd+1:T),smth(cent(:,5),sm,p),'k--');
hold off;
leg=legend('$\bar{e}_t$','$\underbar{e}_t$');
set(leg,'Interpreter','Latex');
xlim([wnd T]);
%ylim([0.04,0.1]);
set(gca,'XTick',[60 80 100 120 140 160 180 200]);
set(gca,'XTickLabel',lab);
set(gca,'FontSize',fs);


