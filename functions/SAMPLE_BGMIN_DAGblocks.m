function M = SAMPLE_BGMIN_DAGblocks(data, nsimu, sdz, ny0,ny1,ny2,MIN0)
% Input
% data:  time series, one for each node
% nsimu: number of Gibbs draws from the posterior
% sdz:   standardize data (sdz=1), demean data (sdz=0)
% ny0:   commont factors (e.g. oil and economic activity)
% ny1:   nodes in the layer 1
% ny2:   nodes in the layer 2
% MIN0:  initial value of the network
%
% Output
% M: estimated contemporaneous multilayer network
%====== PRELIMINARIES;
percent_burn = 0.5;
nburn = fix(percent_burn*nsimu);
nsave = nsimu - nburn;

D = data;
[~,nx] = size(D);
ny=ny0+ny1+ny2;
lag = 0;
D = detrend(D);
%====== PROCESS DATA;
[Sigma_pr, Sigma_pst, nu_pr, nu_pst, Kn, nvar] = PROC_DATA(D, ny, sdz);

% ========= INITIALIZATION;
VxAll = cell(1,1);
for j = 1:ny
    nbrs = 1:ny;     
    nbrs(j)=[];
    VxAll{j,1} =  nbrs;
end

for j = 1:ny1
    nbrs = ny0+1:ny0+ny1;     
    nbrs(j)=[];
    Vx1{j,1} =  nbrs;
end

for j = 1:ny2
    nbrs = ny0+ny1+1:ny0+ny1+ny2;     
    nbrs(j)=[];
    Vx2{j,1} =  nbrs;
end


LogL  = zeros(nsave+1,ny);

logL  = zeros(1,ny);
for yi = 1:ny
    logL(1,yi) = LOG_SCORE(yi, [], Sigma_pr, Sigma_pst, nu_pr, nu_pst,...
        Kn, nvar, lag);
end

LogL(1,:)  = logL;
Dag = MIN0;
DG  = Dag(:);
scn = 5e2;

%====================== Start Sampling ============================
tic;
ct = 1;
for t = 1:nsimu
    %%%% Block G11
    py = randperm(ny1);
    Nx = ny1-1;
    Vx=Vx1;
    for i = 1:ny1
        yi = py(i);
        cand_pa = Vx{yi,1};   % candidate parents
        idj = randperm(Nx,1);
        xi = cand_pa(idj);
        
        Dag_n = Dag;
        e_yx = Dag_n(yi,xi);
        e_xy = Dag_n(xi,yi);
        % first step verification
        if e_yx == 0 && e_xy == 1;
            Dag_n(yi,xi) = 1;  Dag_n(xi,yi) = 0;     % Reverse
            rev = 1;
        else
            Dag_n(yi,xi) = 1 - Dag_n(yi,xi);         % Add or delete
            rev = 0;
        end
        
        % second step verification
        if graphisdag(sparse(Dag_n)) == 1;
            nlogL = logL;
            xj = find(Dag_n(yi,:));
            nlogL(yi) = LOG_SCORE(yi, xj, Sigma_pr, Sigma_pst, nu_pr,...
                nu_pst, Kn, nvar, lag);
            if rev == 1;
                yj = find(Dag_n(xi,:));
                nlogL(xi) = LOG_SCORE(xi, yj, Sigma_pr, Sigma_pst, ...
                    nu_pr, nu_pst, Kn, nvar, lag);
            end
            
            R = exp(sum(nlogL) - sum(logL));
            u = rand;
            if u < min(1,R) % accept the move:
                Dag = Dag_n;
                logL = nlogL;
            end
        end
    end
    %%%% Block G22
    py = randperm(ny2);
    Nx = ny2-1;
    Vx=Vx2;
    for i = 1:ny2
        yi = py(i);
        cand_pa = Vx{yi,1};   % candidate parents
        idj = randperm(Nx,1);
        xi = cand_pa(idj);
        Dag_n = Dag;
        e_yx = Dag_n(yi,xi);
        e_xy = Dag_n(xi,yi);
        % first step verification
        if e_yx == 0 && e_xy == 1;
            Dag_n(yi,xi) = 1;  Dag_n(xi,yi) = 0;     % Reverse
            rev = 1;
        else
            Dag_n(yi,xi) = 1 - Dag_n(yi,xi);         % Add or delete
            rev = 0;
        end
        % second step verification
        if graphisdag(sparse(Dag_n)) == 1
            nlogL = logL;
            xj = find(Dag_n(yi,:));
            nlogL(yi) = LOG_SCORE(yi, xj, Sigma_pr, Sigma_pst, nu_pr,...
                nu_pst, Kn, nvar, lag);
            if rev == 1
                yj = find(Dag_n(xi,:));
                nlogL(xi) = LOG_SCORE(xi, yj, Sigma_pr, Sigma_pst, ...
                    nu_pr, nu_pst, Kn, nvar, lag);
            end
            
            R = exp(sum(nlogL) - sum(logL));
            u = rand;
            if u < min(1,R) % accept the move:
                Dag = Dag_n;
                logL = nlogL;
            end
        end
    end
    %%%%%%%%%
    %%%% Blocks G12 and G21
    py = randperm(ny);
    Nx = ny-1;
    Vx=VxAll;
    for i = 1:ny
        yi = py(i);
        cand_pa = Vx{yi,1};   % candidate parents
        idj = randperm(Nx,1);
        xi = cand_pa(idj);
        Dag_n = Dag;
        e_yx = Dag_n(yi,xi);
        e_xy = Dag_n(xi,yi);
        % first step verification
        if e_yx == 0 && e_xy == 1;
            Dag_n(yi,xi) = 1;  Dag_n(xi,yi) = 0;     % Reverse
            rev = 1;
        else
            Dag_n(yi,xi) = 1 - Dag_n(yi,xi);         % Add or delete
            rev = 0;
        end
        % second step verification
        if graphisdag(sparse(Dag_n)) == 1;
            nlogL = logL;
            xj = find(Dag_n(yi,:));
            nlogL(yi) = LOG_SCORE(yi, xj, Sigma_pr, Sigma_pst, nu_pr,...
                nu_pst, Kn, nvar, lag);
            if rev == 1;
                yj = find(Dag_n(xi,:));
                nlogL(xi) = LOG_SCORE(xi, yj, Sigma_pr, Sigma_pst, ...
                    nu_pr, nu_pst, Kn, nvar, lag);
            end
            
            R = exp(sum(nlogL) - sum(logL));
            u = rand;
            if u < min(1,R) % accept the move:
                Dag = Dag_n;
                logL = nlogL;
            end
        end
    end
    if t > nburn
        ct = ct + 1;
        DG = DG + Dag(:);       LogL(ct,:) = logL;
    end
end
Time = toc;



[PostG, DAG, R1] = CONVERGENCE(DG, LogL, ny, ny, lag);

M.DAG   = DAG;
M.PostG = PostG;
M.PSRF   = R1;
M.Time    = Time;