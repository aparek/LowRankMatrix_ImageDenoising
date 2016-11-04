function [ estCube, wtCube, outPatch, outWt] = Est_Patch( inCube, inWtCube, inpPatch, srchWindow, blcksize, overlap, threshold, is2d,lam)

% threshold - number of similar blocks
estCube = inCube;
wtCube = inWtCube;
if is2d
    [len,wid] = size(inCube);
    I = 1:overlap:len;
    J = 1:overlap:wid;
    K = 1;
else
    [m,n,p] = size(inCube);
    I = 1:overlap:m;
    J = 1:overlap:n;
    K = 1:1:p;
end

Q = cell(length(I)*length(J)*length(K),1);
metricVal = zeros(length(I)*length(J)*length(K),1);

% Assemble all blocks in the given search Window
idx = 1;
for i = I
    for j = J
        for k = K
            if is2d
                
                inLoc_st = [i-blcksize(1)/2, j-blcksize(2)/2];
                inLoc_end = [i+blcksize(1)/2-1, j+blcksize(2)/2-1];
                
                if inLoc_st(1) < 1
                    inLoc_st(1) = i;
                    inLoc_end(1) = i + blcksize(1) - 1;
                end
                
                if inLoc_st(2) < 1
                    inLoc_st(2) = j;
                    inLoc_end(2) = j + blcksize(2) - 1;
                end
                
                if inLoc_end(1) > len
                    inLoc_st(1) = len - blcksize(1) + 1;
                    inLoc_end(1) = len;
                end
                
                if inLoc_end(2)> wid
                    inLoc_st(2) = wid - blcksize(2) + 1;
                    inLoc_end(2) = wid;
                end
                
                Q{idx}{1} = srchWindow(inLoc_st(1):inLoc_end(1), inLoc_st(2):inLoc_end(2));
                metricVal(idx) = (norm(Q{idx}{1}(:)-inpPatch(:))^2)/(blcksize(1) * blcksize(2));
                Q{idx}{2} = inLoc_st;
                Q{idx}{3} = inLoc_end;
                
            else
                inLoc_st = [i-blcksize(1)/2, j-blcksize(2)/2, k-blcksize(3)/2];
                inLoc_end = [i+blcksize(1)/2-1, j+blcksize(2)/2-1, k+blcksize(3)/2-1];
                
                if inLoc_st(1)< 1
                    inLoc_st(1) = i;
                    inLoc_end(1) = i + blcksize(1) - 1;
                end
                
                if inLoc_st(2) < 1
                    inLoc_st(2) = j;
                    inLoc_end(2) = j + blcksize(2) - 1;
                end
                
                if inLoc_st(3) < 1
                    inLoc_st(3) = k;
                    inLoc_end(3) = k + blcksize(3) - 1;
                end
                
                if inLoc_end(1) >= m
                    inLoc_st(1) = m - blcksize(1) + 1;
                    inLoc_end(1) = m;
                end
                
                if inLoc_end(2) >= n
                    inLoc_st(2) = n - blcksize(2) + 1;
                    inLoc_end(2) = n;
                end
                
                if inLoc_end(3) >= p
                    inLoc_st(3) = p - blcksize(3) + 1;
                    inLoc_end(3) = p;
                end
                
                Q{idx}{1} = srchWindow(inLoc_st(1):inLoc_end(1), inLoc_st(2):inLoc_end(2), inLoc_st(3):inLoc_end(3));
                
                metricVal(idx) = norm(Q{idx}{1}(:)-inpPatch(:))^2;
                Q{idx}{2} = inLoc_st;
                Q{idx}{3} = inLoc_end;
            end
            
            idx = idx + 1;
        end
    end
end


if is2d
    StackedMatrix = zeros(blcksize(1)*blcksize(2), threshold);
    StackedMatrix(:,1) = inpPatch(:);
else
    StackedMatrix = zeros(blcksize(1)*blcksize(2)*blcksize(3), threshold);
    StackedMatrix(:,1) = inpPatch(:);
end

% Find indices of the (threshold) similar patches
Indx = zeros(threshold,1);
metricVal(metricVal == 0) = Inf;
for i = 2:threshold
    Indx(i) = find(metricVal == min(metricVal),1);
    StackedMatrix(:,i) = Q{Indx(i)}{1}(:);
    metricVal(Indx(i)) = Inf;
end

% denoise the stacked matrix
[U,S,V] = svd(StackedMatrix,'econ');
Denoised = U * diag(thresh(diag(S), lam, 0.45/lam, 'firm')) * V';

% Place the estimated patches back in the cube or image

if is2d
    outPatch = reshape(Denoised(:,1), [blcksize(1), blcksize(2)]);
    outWt = ones(blcksize(1), blcksize(2));
    % Place all the other columns
    for i = 2:length(Indx)
        loc = Q{Indx(i)}{2};
        loc_end = Q{Indx(i)}{3};
        estCube(loc(1):loc_end(1),...
            loc(2):loc_end(2)) = (inCube(loc(1):loc_end(1),...
            loc(2):loc_end(2))+ ...
            reshape(Denoised(:,i), [blcksize(1), blcksize(2)]));
        
        wtCube(loc(1):loc_end(1),...
            loc(2):loc_end(2)) = (inWtCube(loc(1):loc_end(1),...
            loc(2):loc_end(2))+ ...
            ones(blcksize(1), blcksize(2)));
        
    end
else
    % First column of patch matrix gets placed in the cube
    outPatch = reshape(Denoised(:,1), [blcksize(1), blcksize(2), blcksize(3)]);
    outWt = ones(blcksize(1), blcksize(2), blcksize(3));
    % Place all the other columns
    for i = 2:length(Indx)
        loc = Q{Indx(i)}{2};
        loc_end = Q{Indx(i)}{3};
        estCube(loc(1):loc_end(1),...
            loc(2):loc_end(2),...
            loc(3):loc_end(3)) = (inCube(loc(1):loc_end(1),...
            loc(2):loc_end(2),...
            loc(3):loc_end(3)) + ...
            reshape(Denoised(:,i), [blcksize(1), blcksize(2), blcksize(3)]));
        
        wtCube(loc(1):loc_end(1),...
            loc(2):loc_end(2),...
            loc(3):loc_end(3)) = (inWtCube(loc(1):loc_end(1),...
            loc(2):loc_end(2),...
            loc(3):loc_end(3)) + ones(blcksize(1), blcksize(2), blcksize(3)));
    end
end




