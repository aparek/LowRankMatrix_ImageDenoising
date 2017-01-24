function [ estCube ] = lowRank3D(Y, blcksize, overlap, threshold, searchSize, is2d,lam )

if is2d
    [len,wid] = size(Y);
    estCube = zeros(len,wid);
    wtCube =  zeros(len,wid);
else
    [len, wid, ht] = size(Y);
    estCube = zeros(len,wid,ht);
    wtCube =  zeros(len,wid,ht);
end

I = 1:overlap:len;
J = 1:overlap:wid;
if is2d
    K = 1;
else
    K = 1:1:ht;
end


% Compute all the patches with index
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
                
                
                inpPatch = Y(inLoc_st(1):inLoc_end(1), inLoc_st(2):inLoc_end(2));
                strt_loc = [i-searchSize(1)/2, j-searchSize(2)/2]';
                end_loc = [i+searchSize(1)/2-1, j+searchSize(2)/2-1]';
                
                
                if strt_loc(1) < 1
                    strt_loc(1) = i;
                    end_loc(1) = i+searchSize(1)-1;
                end
                
                if strt_loc(2) < 1
                    strt_loc(2) = j;
                    end_loc(2) = j + searchSize(2)-1;
                end
                
                
                if end_loc(1) > len
                    strt_loc(1) = len - searchSize(1) + 1;
                    end_loc(1) = len;
                end
                
                if end_loc(2) > wid
                    strt_loc(2) = wid - searchSize(2) + 1;
                    end_loc(2) = wid;
                end
                
                srchWindow = Y(strt_loc(1):end_loc(1), strt_loc(2):end_loc(2));
                inCube = estCube(strt_loc(1):end_loc(1), strt_loc(2):end_loc(2));
                inWtCube = wtCube(strt_loc(1):end_loc(1), strt_loc(2):end_loc(2));
                
                [estCube(strt_loc(1):end_loc(1), ...
                    strt_loc(2):end_loc(2)), ...
                    wtCube(strt_loc(1):end_loc(1), ...
                    strt_loc(2):end_loc(2)),...
                    outPatch, outWt] = Est_Patch(inCube, inWtCube, inpPatch, ...
                    srchWindow, blcksize, overlap, threshold, is2d,lam);
                
                estCube(inLoc_st(1):inLoc_end(1),...
                    inLoc_st(2):inLoc_end(2)) = (estCube(inLoc_st(1):inLoc_end(1),...
                                                 inLoc_st(2):inLoc_end(2)) + outPatch);
                
                wtCube(inLoc_st(1):inLoc_end(1),...
                    inLoc_st(2):inLoc_end(2)) = (wtCube(inLoc_st(1):inLoc_end(1),...
                                                    inLoc_st(2):inLoc_end(2)) + outWt);
%                 fprintf('Denoised all cubes at [%d, %d] \n',i,j)
                
            else
                inLoc_st = [i-blcksize(1)/2, j-blcksize(2)/2, k-blcksize(3)/2];
                inLoc_end = [i+blcksize(1)/2-1, j+blcksize(2)/2-1, k+blcksize(3)/2-1];
                
                if inLoc_st(1) < 1
                    inLoc_st(1) = i;
                    inLoc_end(1) = i + blcksize(1) - 1;
                end
                
                if inLoc_st(2)< 1
                    inLoc_st(2) = j;
                    inLoc_end(2) = j + blcksize(2) - 1;
                end
                
                if inLoc_st(3) < 1
                    inLoc_st(3) = k;
                    inLoc_end(3) = k + blcksize(3) - 1;
                end
                
                if inLoc_end(1) >= len
                    inLoc_st(1) = len - blcksize(1) + 1;
                    inLoc_end(1) = len;
                end
                
                if inLoc_end(2) >= wid
                     inLoc_st(2) = wid - blcksize(2) + 1;
                     inLoc_end(2) = wid;
                end
                
                if inLoc_end(3) >= ht
                     inLoc_st(3) = ht - blcksize(3) + 1;
                    inLoc_end(3) = ht;
                end
                
                inpPatch = Y(inLoc_st(1):inLoc_end(1), inLoc_st(2):inLoc_end(2), inLoc_st(3):inLoc_end(3));
                strt_loc = [i-searchSize(1)/2,j-searchSize(2)/2,k-searchSize(3)/2]';
                end_loc = [i+searchSize(1)/2-1, j+searchSize(2)/2-1,  k+searchSize(3)/2-1]';
                
                if strt_loc(1) < 1
                    strt_loc(1) = i;
                    end_loc(1) = i+searchSize(1)-1;
                end
                
                if strt_loc(2) < 1
                    strt_loc(2) = j;
                    end_loc(2) = j + searchSize(2)-1;
                end
                
                if strt_loc(3) < 1
                    strt_loc(3) = k;
                    end_loc(3) = k + searchSize(3)-1;
                end
                
                if end_loc(1) >= len
                    strt_loc(1) = len-searchSize(1)+1;
                    end_loc(1) = len;
                end
                
                if end_loc(2) >= wid
                    strt_loc(2) = wid-searchSize(2)+1;
                    end_loc(2) = wid;
                end
                
                if end_loc(3) >= ht
                    strt_loc(3) = ht-searchSize(3)+1;
                    end_loc(3) = ht;
                end
                
                srchWindow = Y(strt_loc(1):end_loc(1), strt_loc(2):end_loc(2), strt_loc(3):end_loc(3));
                inCube = estCube(strt_loc(1):end_loc(1), strt_loc(2):end_loc(2), strt_loc(3):end_loc(3));
                inWtCube = wtCube(strt_loc(1):end_loc(1), strt_loc(2):end_loc(2), strt_loc(3):end_loc(3));
                
                [estCube(strt_loc(1):end_loc(1), ...
                    strt_loc(2):end_loc(2), ...
                    strt_loc(3):end_loc(3)), ...
                    wtCube(strt_loc(1):end_loc(1), ...
                    strt_loc(2):end_loc(2), ...
                    strt_loc(3):end_loc(3)), ...
                    outPatch, outWt] = Est_Patch(inCube, inWtCube, inpPatch,srchWindow, blcksize, overlap, threshold, is2d,lam);
                
                estCube(inLoc_st(1):inLoc_end(1),...
                    inLoc_st(2):inLoc_end(2),...
                    inLoc_st(3):inLoc_end(3)) = (estCube(inLoc_st(1):inLoc_end(1),...
                                                    inLoc_st(2):inLoc_end(2),...
                                                    inLoc_st(3):inLoc_end(3)) + outPatch);
                
                
                wtCube(inLoc_st(1):inLoc_end(1),...
                    inLoc_st(2):inLoc_end(2),...
                    inLoc_st(3):inLoc_end(3)) = (wtCube(inLoc_st(1):inLoc_end(1),...
                                                inLoc_st(2):inLoc_end(2),...
                                                inLoc_st(3):inLoc_end(3)) + outWt);
                
                fprintf('Denoised all cubes at [%d, %d, %d] \n',i,j,k)
            end
        end
    end
end

estCube = estCube./(wtCube+eps);
end

