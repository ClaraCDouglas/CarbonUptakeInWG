function [ev_index_1,tda_1,ev_index_2,tda_2,pev,r] = calc_pigup_EOF_2D(index3d_1,index3d_2,time)
%function [ev_index_1,tda_1,ev_index_2,tda_2,pev] = calc_pigup_EOF_2D(index3d_1,index3d_2,time)
% calculation as described (somewhat opaquely) by Bjornsson & Venegas " A
% Manual for EOF & SVD analysis of Climate Data"
%
% first things first convert 3D grid into 2D
% gridcells on rows time on cols

[r3d1,c3d1,z3d1] = size(index3d_1);
[r3d2,c3d2,z3d2] = size(index3d_2);

% remember how you do this so you can put it back!
% this is column major
for i = 1:1:z3d1
    index2d_1(:,i) = reshape(index3d_1(:,:,i),r3d1*c3d1,1);
    % puts it back
    % cdd_new(:,:,i) = reshape(cdd2d(:,i),r3d,c3d);
end

for i = 1:1:z3d2
    index2d_2(:,i) = reshape(index3d_2(:,:,i),r3d2*c3d2,1);
    % puts it back
    % cdd_new(:,:,i) = reshape(cdd2d(:,i),r3d,c3d);
end

% now we need to remove linear trend and mean at each location
[r2d1,~] = size(index2d_1);
for i = 1:1:r2d1
    tmp_index = index2d_1(i,:);
    tmp_time  = time;
    
    % remove any NaN
    iix            = find(isnan(tmp_index));
    tmp_index(iix) = [];
    tmp_time(iix)  = [];
    
    if length(tmp_index) > 1
        % calculate mean
        m_idx          = nanmean(tmp_index);
    
        % calulate linear trend of de-meaned data
        tmp_index      = tmp_index - m_idx;
        m_pp           = polyfit(tmp_time',tmp_index,1);
        m_trend        = polyval(m_pp,time');
        slopes_1(i)    = m_pp(1);
    
        % remove mean and linear trend
        eof_index_1(i,:) = (index2d_1(i,:) - m_idx) - m_trend;
    
        % calulate std
        m_std          = nanstd(eof_index_1(i,:));
    
        % nomalize
        %if abs(m_std) > 0
        %    eof_index_1(i,:) = eof_index_1(i,:) ./ m_std;
        %end
    
        % this is very quick and very dirty
        % replace any NaN with zero
    end
    
    eof_index_1(i,iix) = 0;
    %disp(['No of NaN : ',num2str(length(iix))]);
    
end

% now we need to remove linear trend and mean at each location
[r2d2,~] = size(index2d_2);
for i = 1:1:r2d2
    tmp_index = index2d_2(i,:);
    tmp_time  = time;
    
    % remove any NaN
    iix            = find(isnan(tmp_index));
    tmp_index(iix) = [];
    tmp_time(iix)  = [];
    
    if length(tmp_index) > 1
        % calculate mean
        m_idx          = nanmean(tmp_index);
    
        % calulate linear trend of de-meaned data
        tmp_index      = tmp_index - m_idx;
        m_pp           = polyfit(tmp_time',tmp_index,1);
        m_trend        = polyval(m_pp,time');
        slopes_2(i)      = m_pp(1);
    
        % remove mean and linear trend
        eof_index_2(i,:) = (index2d_2(i,:) - m_idx) - m_trend;
    
        % calulate std
        m_std          = nanstd(eof_index_2(i,:));
    
        % nomalize
        %if abs(m_std) > 0
        %     eof_index_2(i,:) = eof_index_2(i,:) ./ m_std;
        %end
    
        % this is very quick and very dirty
        % replace any NaN with zero
    end
    
    eof_index_2(i,iix) = 0;
    %disp(['No of NaN : ',num2str(length(iix))]);
    
end

eof_index = eof_index_1' * eof_index_2;
% right now the difficult bit
% first a singular value decompsition of our data (this is as described in
% E&T P333 (the zero is for compatibility with E&T's description
[U,Lam,V] = svd(eof_index);

% now the values in U are the eigenvectors (so it's the same dims as the
% input) so we can out this back
% this is column major

%NB Z3 (time axis must be same for both 1 & 2
for i = 1:1:z3d1
    PC1               = eof_index_1 * U(:,i);
    PC2               = eof_index_2 * V(:,i);
    
    r(i)              = calc_corr_val(PC1,PC2);
    
    ev_index_1(:,:,i) = reshape(PC1,r3d1,c3d1);
    ev_index_2(:,:,i) = reshape(PC2,r3d2,c3d2);
    
    tda_1(:,i)        =  U(:,i);
    tda_2(:,i)        =  V(:,i);
end

% the precentage variance explained by each mode
expVar   = diag(Lam).^2;
pev      =(expVar/sum(expVar))*100;

end
    
    
    
