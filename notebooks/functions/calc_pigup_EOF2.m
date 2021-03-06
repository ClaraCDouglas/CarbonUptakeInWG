function [ev_index,tda,pev,trends] = calc_pigup_EOF2(index3d,time)
%function [ev_index,tda,pev,trends] = calc_pigup_EOF2(index3d,time)
% calculation as described (somewhat opaquely) by Bjornsson & Venegas " A
% Manual for EOF & SVD analysis of Climate Data"

% first things first convert 3D grid into 2D
% gridcells on rows time on cols

[r3d,c3d,z3d] = size(index3d);

% remember how you do this so you can put it back!
% this is column major
for i = 1:1:z3d
    index2d(:,i) = reshape(index3d(:,:,i),r3d*c3d,1);
    % puts it back
    % cdd_new(:,:,i) = reshape(cdd2d(:,i),r3d,c3d);
end

% now we need to remove linear trend and mean at each location
[r2d,c2d] = size(index2d);
for i = 1:1:r2d
    tmp_index = index2d(i,:);
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
        slopes(i)      = m_pp(1);
    
        % remove mean and linear trend
        eof_index(i,:) = (index2d(i,:) - m_idx) - m_trend;
    
        % calulate std
        %m_std          = nanstd(eof_index(i,:));
    
        % nomalize
        %eof_index(i,:) = eof_index(i,:) ./ m_std;
    
        % this is very quick and very dirty
        % replace any NaN with zero
    end
    
    eof_index(i,iix) = 0;
    %disp(['No of NaN : ',num2str(length(iix))]);
    
end

% see if there's any trends in this data
trends = []%reshape(slopes,r3d,c3d); % doesn't work if some cells are removed because of NaNs - changes the size of the array?

% right now the difficult bit
% first a singular value decompsition of our data (this is as described in
% E&T P333 (the zero is for compatibility with E&T's description
[C,Lam,CC] = svd(eof_index);

% now the values in U are the eigenvectors (so it's the same dims as the
% input) so we can out this back
% this is column major

for i = 1:1:z3d
    
    PCi             = eof_index * CC(:,i);
    ev_index(:,:,i) = reshape(PCi,r3d,c3d);
    tda(:,i)        =  CC(:,i);
    
    % this should now be of same dimension as eof_index
    %ccmi            = PCi * tda(:,i)';
   
    %for r = 1:1:r2d
    %    rrmi(r) = calc_corr_val(eof_index(r,:),ccmi(r,:));
    %end
    
    %r_map(:,:,i)   = reshape(rrmi,r3d,c3d);
end

% the precentage variance explained by each mode
expVar   = diag(Lam).^2;
pev      =(expVar/sum(expVar))*100;

end
    
    
    
