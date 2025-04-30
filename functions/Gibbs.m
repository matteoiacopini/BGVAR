function [MAR,MIN]=Gibbs(DataSub, nsimuInit, nsimu, sdz, lag, ny,ny0,ny1,ny2)
% Initialization
    MAR0 = SAMPLE_BGMAR_DAG(DataSub, nsimuInit, sdz, lag, ny);
    MIN0 = SAMPLE_BGMIN_DAG(DataSub, nsimuInit, sdz, ny);
    % Estimation
    %-------------------SAMPLING MIN & MAR STRUCTURE -------------------------
    MAR = SAMPLE_BGMAR_DAGblocks(DataSub, nsimu, sdz, lag, ny0,ny1,ny2,MAR0.DAG);
    MIN = SAMPLE_BGMIN_DAGblocks(DataSub, nsimu, sdz, ny0,ny1,ny2,MIN0.DAG);
