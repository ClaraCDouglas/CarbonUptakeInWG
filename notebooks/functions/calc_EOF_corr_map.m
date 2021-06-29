function r_map = calc_EOF_corr_map(index3d,eof3d,tda,no_mod)

[r3d,c3d,z3d] = size(index3d);

% remember how you do this so you can put it back!
% this is column major
for i = 1:1:z3d
    index2d(:,i)   = reshape(index3d(:,:,i),r3d*c3d,1);
    eof_index(:,i) = reshape(eof3d(:,:,i),r3d*c3d,1);
    % puts it back
    % cdd_new(:,:,i) = reshape(cdd2d(:,i),r3d,c3d);
end

[r2d,c2d] = size(index2d);
if no_mod <= c2d
    for i = 1:1:no_mod
    
        % this should now be of same dimension as index2d
        ccmi            = eof_index(:,i) * tda(:,i)';
   
        for r = 1:1:r2d
            rrmi(r) = calc_corr_val(index2d(r,:),ccmi(r,:));
        end
    
        r_map(:,:,i)   = reshape(rrmi,r3d,c3d);
    end
else
    r_map = NaN;
    error('Not that many modes')
end

end