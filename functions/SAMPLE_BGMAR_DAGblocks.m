function M = SAMPLE_BGMAR_DAGblocks(data, nsimu, sdz, lag, ny0,ny1,ny2,MAR0)
% Input
% data:  time series, one for each node
% nsimu: number of Gibbs draws from the posterior
% sdz:   standardize data (sdz=1), demean data (sdz=0)
% lag:   takes values 1,2,... and indicates the VAR order
% ny0:   commont factors (e.g. oil and economic activity)
% ny1:   nodes in the layer 1
% ny2:   nodes in the layer 2
% MAR0:  initial value of the network
% Output
% M:     estimated laged Granger multilayer network
%====== PRELIMINARIES;
percent_burn = 0.5;
nburn = fix(percent_burn*nsimu);        
nsave = nsimu - nburn;

D = data;
[~,nx] = size(D); 
ny=ny0+ny1+ny2;

%====== PROCESS DATA;
[Sigma_pr, Sigma_pst, nu_pr, nu_pst, Kn, nvar] = PROC_DATA(D, ny,sdz,lag);

% ========= INITIALIZATION;
LogL  = zeros(nsave+1,ny);

logL  = zeros(1,ny);
for yi = 1:ny
   logL(1,yi) = LOG_SCORE(yi, [], Sigma_pr, Sigma_pst, nu_pr, nu_pst,...
       Kn, nvar, lag);
end
    
LogL(1,:)  = logL;
Nx = nx*lag;
Dag = MAR0;
DG  = Dag(:);
scn = 5e2;

%====================== Start Sampling ============================
tic;  
ct = 1;
for t = 1:nsimu
    % Local Move: Blocks G11 and G12 
    for yi = ny0+1:ny0+ny1
        j = randperm(Nx,1);
        parents = Dag(yi,:);
        new_parents = parents;
        new_parents(1,j) = 1 - new_parents(1,j);
        
        xj = find(new_parents);    
        nlogL = LOG_SCORE(yi, xj, Sigma_pr, Sigma_pst, nu_pr, nu_pst,...
            Kn, nvar, lag); 
        
        R = exp(nlogL - logL(yi));
        u = rand; 
        if u < min(1,R) % accept the move:
            Dag(yi,:) = new_parents;
            logL(yi) = nlogL;  
        end
    end
    % Local Move: Block G22 and G21
    for yi = ny0+ny1+1:ny
        j = randperm(Nx,1);
        parents = Dag(yi,:);
        new_parents = parents;
        new_parents(1,j) = 1 - new_parents(1,j);
        
        xj = find(new_parents);    
        nlogL = LOG_SCORE(yi, xj, Sigma_pr, Sigma_pst, nu_pr, nu_pst,...
            Kn, nvar, lag); 
        
        R = exp(nlogL - logL(yi));
        u = rand; 
        if u < min(1,R) % accept the move:
            Dag(yi,:) = new_parents;
            logL(yi) = nlogL;  
        end
    end
    % Global Move: Blocks G11, G22, G12 and G21
    for yi = 1:ny
        j = randperm(Nx,1);
        parents = Dag(yi,:);
        new_parents = parents;
        new_parents(1,j) = 1 - new_parents(1,j);
        
        xj = find(new_parents);    
        nlogL = LOG_SCORE(yi, xj, Sigma_pr, Sigma_pst, nu_pr, nu_pst,...
            Kn, nvar, lag); 
        
        R = exp(nlogL - logL(yi));
        u = rand; 
        if u < min(1,R) % accept the move:
            Dag(yi,:) = new_parents;
            logL(yi) = nlogL;  
        end
    end
    if t > nburn
        ct = ct + 1;        
        DG = DG + Dag(:);       LogL(ct,:) = logL;
    end
end
Time = toc;

[PostG, DAG, R1, R2] = CONVERGENCE(DG, LogL, ny, Nx, lag);

M.DAG   = DAG;
M.PostG = PostG;
M.PSRF   = R1;
M.MPSRF = R2;
M.Time    = Time;