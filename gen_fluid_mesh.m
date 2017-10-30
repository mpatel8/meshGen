clc; clear all ;
dbg = 0;

% generate fluid mesh
% [ id (1),  x,y,z (2:4)  p,t,r (5:7)  x,y (8) ]

% read in lumen mesh [ node-set NL and S4R-element-set EL ]
szDir = 'mesh';
%szFlLumen = 'seg-1-36.inp';
szFlLumen = 'seg-5.inp';
szFlLumen = 'seg-36-end.inp';
%szFlLumen = 'seg-10-end.inp';
%szFlLumen = 'seg-all.inp';
%szFlLumen = 'lumen-S_A0039R_W9A1_FSI_Y2017M10D23H16M41.inp';
szFlLumen = 'S_A0039R_W9A1_FEA_Y2017M10D05H12M13.inp';

[ mNL, mEL ] = read_input(strcat(szDir, '\', szFlLumen ));

% determine how many nodes per segment 'N' using first set of nodes that have
% same X value [IDEALLY MAKE AVAILABLE VIA MESH INPUT FILE - OR VIA MATLAB CODE ]
i = 1; 
x_org = mNL(i, 2:2);
while abs(mNL(i+1, 2:2)-x_org) < 1e-4
    i = i + 1;

end
cNLS = i
cNLS = 36;

% look through each N nodes and get plane definition
% if plane samilar to prior then its rigid-wide-end (check PointOnPlane function
    % %% TBA %% note which N-batch this starts, RS, ends, RE, number-segments RN, and plane, RN
    % %% TBA %%collate these nodes into list R
% store each N-batch in a separate list, NS
% store each plane info, PS
% note total number of actual segments T
cLS = size(mNL, 1)/cNLS; cLSRS = 0; cLSRE = 0; mmNLS = []; mPLS = [];cLSENDWIDE=10000;
for iNLS = 1:cLS
    mmNLS(:, :, iNLS) = [ mNL((iNLS-1)*cNLS+1:iNLS*cNLS, 1:4) zeros(cNLS, 5)] ;
    p = fitPlane(mmNLS(:, 2:4, iNLS));
    if iNLS > 1
        if abs(distancePoints3d(mean(mmNLS(:, 2:4, iNLS-1)), mean(mmNLS(:, 2:4, iNLS)))) < 0.001
            if cLSRS == 0
                cLSENDWIDE = iNLS-1;
                cLSRS = iNLS;
            end
        elseif cLSRS ~= 0 && cLSRE == 0
            if abs(distancePoints3d(mean(mmNLS(:, 2:4, iNLS-1)), mean(mmNLS(:, 2:4, iNLS)))) > 0.001
                cLSRE = iNLS-1;
            end
        end
    end
                
    mPLS = [ mPLS; p ];
    mmNLS(:, 8:9, iNLS) = planePosition(mmNLS(:, 2:4, iNLS), mPLS(iNLS, :));
    mmNLS(:, 5:7, iNLS) = cart2sph2d(horzcat(mmNLS(:, 8:9, iNLS), zeros(cNLS, 1)));
end
if cLSENDWIDE == 0
   cLSENDWIDE = cLS; 
else
    cLSTAPER = cLSENDWIDE - (cLSRE-cLSRS);
end
nNIDmax = mNL(end,1);
mNF = []; mEF = []; 

bBLadd = 1;
if bBLadd == 1
% define boundary layer as defined by first internal segment on rigid wide end or 10%
% Splits the boundary into 5 layers, each 20% wider than before
% Boundary layer is stored in BL and NLS is switched to inner layer
mmNBL1 = mmNLS(:,:,1:cLSRS-1); % outer layer
% mNF = reshape(permute(mmNBL1(:, 1:4,:),[1 3 2]),[],4,1);

cBLgwth = 1.2; % growth factor of boundary layer
cBLnum = 5;
nBLtotal = 0;
for i = 1:cBLnum
    nBLtotal = nBLtotal+cBLgwth^(i-1);
end
cBLpct=0.1; % total pct of 
if cLSRE - cLSRS > 0
    mBLpct = [ mmNLS(:, 7, cLSRS-1) - mmNLS(:, 7, cLSRS)]./mmNLS(:, 7, cLSRS-1);
else
    mBLpct = repmat(cBLpct, cNLS, 1);
end
mBLLpct=[]; mBLLpctPrev=0;
for iLAYER = 1:cBLnum
    mBLLpctPrev = [ mBLLpctPrev + (cBLgwth^(iLAYER-1))*mBLpct/nBLtotal ];
    %mmBLLpct(:, iLAYER) = 1-prevBLLpct;
    mBLLpct = [ mBLLpct, 1-mBLLpctPrev ];
end
 
mmNBL=[]; mmEBL = []; 
nEFmax = 0;
for iNLS = 1:cLSENDWIDE
    
    % Node boundary layer
    for iLAYER = 1:cBLnum+1
        if iLAYER == 1 % outer layer 
            mmNBL(:,:, iLAYER, iNLS) = mmNLS(:, :, iNLS);
        elseif iNLS == cLSENDWIDE && iLAYER == cBLnum+1 % inside layer at wide end
            mmNBL(:,:, iLAYER, iNLS) = mmNLS(:, :, cLSRS);
        else
            mSPH = [ mmNLS(:, 5:6, iNLS), mmNLS(:, 7, iNLS).*mBLLpct(:,iLAYER-1) ];
            mCYL = sph2cart2d(mSPH);
            m3D = planePoint(mPLS(iNLS, :), mCYL(:, 1:2));
            mmNBL(:,:, iLAYER, iNLS) = [ mmNLS(:, 1, 1)+nNIDmax, m3D, mSPH, mCYL(:, 1:2) ];
        end
        mNF = [ mNF; mmNBL(:,1:4, iLAYER, iNLS) ];

        % Add the elemebts
        if iNLS > 1 && iLAYER > 1
            mmEBL(:, :, iLAYER, iNLS-1) = gen_outer_elem(nEFmax, ...
                                            mmNBL(:,1, iLAYER, iNLS-1), mmNBL(:,1, iLAYER-1, iNLS-1), ...
                                            mmNBL(:,1, iLAYER, iNLS), mmNBL(:,1, iLAYER-1, iNLS));
            mEF = [ mEF; mmEBL(:, :, iLAYER, iNLS-1) ];
            nEFmax = max(mEF(:, 1));
        end
        if nNIDmax < mNF(end,1), nNIDmax = mNF(end,1); end;
    end
end


mPLSadj=[];
for iNLS = 1:cLSENDWIDE
    mmNLSadj(:,:,iNLS) = mmNBL(:,:,cBLnum+1, iNLS);
    mPLSadj = [ mPLSadj; mPLS(iNLS,:) ];
end
for iNLS = 1:cLS-cLSRS
    mmNLSadj(:,:,iNLS+cLSENDWIDE) = mmNLS(:,:,iNLS+cLSRS);
    mPLSadj = [ mPLSadj; mPLS(iNLS+cLSRS,:) ];
end
cLSENDWIDEadj = cLSENDWIDE;
cLSRSadj = cLSRS;
cLSREadj = cLSRE-1;
cLSTAPERadj = cLSTAPER+1;
cLSadj = cLS-1;


else
    mPLSadj=[];
    for iNLS = 1:cLS
        mmNLSadj(:,:,iNLS) = mmNLS(:,:,iNLS);
    end
    cLSENDWIDEadj = cLSENDWIDE;
    cLSRSadj = cLSRS;
    cLSREadj = cLSRE;
    cLSTAPERadj = cLSTAPER;
    cLSadj = cLS;
end

nNLIS = 10^(round(log10(nNIDmax)+0.5));
nNLIB = 10*nNLIS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1 == 1    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% read in associalted idealised mesh (2 sets of nodes and 1 set of elements
% store 1st N-batch of ideal nodes in temmplate set NI1
% store 2nd N-batch of ideal nodes in associated set NI2
% store element definitions, EI
% determing how many 
    % interal nodes 1st segment NII
    % numbering start of 1st set of internal nodes NIIS
    % numbering start of 2nd set of internal nodes NIIE
szFlIdeal = strcat(szDir, '\', 'ideal-', int2str(cNLS), '-one-layer.inp');
[ mNI, mEI ] = read_input (szFlIdeal);
mNI1 = mNI(1:size(mNI, 1)/2, :);
mNI2 = mNI(1+size(mNI, 1)/2:size(mNI, 1), :);
cNI = size(mNI1, 1);
cNII = cNI - cNLS;
cNI1IS = mNI1(cNLS+1, 1);
cNI1IE = mNI1(size(mNI1, 1), 1);
cNI2IS = mNI( cNI+cNLS+1, 1);
cEI = size(mEI, 1);

% list circumferential nodes, NIC
% list internal internal nodes, NIN
mNI1C = mNI1(1:cNLS, :);
mNI1N = mNI1(cNLS+1:size(mNI1, 1), :);

% for NIT, for each node, calc radius NIR and acngle NIA
mI_p = fitPlane(mNI1(:, 2:4));
mNI1_2d = planePosition(mNI1(:,2:4), mI_p);
mNI1_sp = cart2sph2d(horzcat(mNI1_2d, zeros(cNI, 1)));
mNIC_sp = mNI1_sp(1:cNLS, :);
mNIA_sp = mNI1_sp(cNLS+1:cNI, :);

nRotate = mNIC_sp(1, 2);
for i = 1:cNLS
    nAngle = mNIC_sp(i, 2)-nRotate;
    if nAngle < 0
        nAngle = 360.0 + nAngle;
    elseif nAngle > 360.0
        nAngle = nAngle - 360.0;
    end
    mNIC_sp(i, 2) = nAngle;
end

for i = 1:cNII
    nAngle = mNIA_sp(i, 2)-nRotate;
    if nAngle < 0
        nAngle = 360.0 + nAngle;
    elseif nAngle > 360.0
        nAngle = nAngle - 360.0;
    end
    mNIA_sp(i, 2) = nAngle;
end

% for internal nodes, get two NIC whose NIA <= & > its NIA, store NINC
mNIA_sp = horzcat(mNIA_sp, zeros(size(mNIA_sp, 1), 3));
for i = 1:cNII
    bFound = 0;
    for j = 1:cNLS-1
        if sign(mNIA_sp(i,2)-mNIC_sp(j,2)) ~= sign(mNIA_sp(i,2)-mNIC_sp(j+1,2))
            bFound = 1;
            break;
        end
    end
    if bFound == 0
        mNIA_sp(i, 4:6) = [ mNI1N(i, 1)-(cNI1IS-1), cNLS, 1 ];
    else
        mNIA_sp(i, 4:6) = [ mNI1N(i, 1)-(cNI1IS-1), j, j+1 ];
    end
end

% for N-batch in S, NLS
    % if SN is between RS-(RE-RS) and RE then skip
    % add nodes to fluid segment node list, FSN
    % get plane definition P
    % Transform SN to 2d plane SNT
    % 
    % for each node in SNT, calculate radius SNTR and angle SNTG
    % for each nodes in NIN
        % calculate equivalent node locaton NIT_i on P based on the two equivalent nodes on SNT
        % rotate node back to global coordinates NITC_i
        % append NITC_i node id, coordinates in FSN
mTrans = [];
%mNF = []; mEF = [];
mmNFS = []; mmEFS = [];
nLSadj=cLSRSadj-(cLSREadj-cLSRSadj)-1; cFS=0; cFSO=0; nNEW=0;
for iNLS = 1:cLSadj
    % collate full segment info 
    mNSFP_add=[]; % plane to transform the nodes into
    mmNSFinf = []; % ids of outer nodes and their 2D coords reference
    nS = 0; nR = 0;
    if iNLS < cLSTAPERadj || iNLS > cLSREadj % Normal main mesh or constriction section
        nS = nS+1;
        mmNSFinf(:,:,nS) = [ mmNLSadj(1:cNLS, 1, iNLS), planePosition(mmNLSadj(:, 2:4, iNLS), mPLSadj(iNLS, :)) ];
        mNSFP_add = [ mNSFP_add; mPLSadj(iNLS, :) ];

    elseif iNLS <= cLSENDWIDEadj % between taperring section near end
        % Add taper points based on end points (inner first then on segment)
        % collate full segment info for each
        nS = nS+1; nNEW = nNEW + 1; nENDS = cLSENDWIDEadj+(iNLS-cLSTAPERadj)+1;

        mmNSFinf(:,:,nS) = [ mmNLSadj(1:cNLS, 1, cLSadj)+cNLS*nNEW+nNIDmax, planePosition(mmNLSadj(:, 2:4, nENDS), mPLSadj(cLSENDWIDEadj, :)) ];
        mNSFP_add = [ mNSFP_add; mPLSadj(iNLS-1, 1:3)+0.5*(mPLSadj(iNLS, 1:3)-mPLSadj(iNLS-1, 1:3)) mPLSadj(iNLS, 4:9) ];

        nS = nS+1;
        if nENDS < cLSREadj
            nNEW = nNEW + 1; 
            mmNSFinf(:,:,nS) = [ mmNLSadj(1:cNLS, 1, cLSadj)+cNLS*nNEW+nNIDmax, planePosition(mmNLSadj(:, 2:4, nENDS), mPLSadj(cLSENDWIDEadj, :)) ];
        else
            mmNSFinf(:,:,nS) = [ mmNLSadj(1:cNLS, 1, nENDS), planePosition(mmNLSadj(:, 2:4, nENDS), mPLSadj(cLSENDWIDEadj, :)) ];
        end
        mNSFP_add = [ mNSFP_add; mPLSadj(iNLS, :) ];

        % if required add additional mesh points
        for iPT = 1:iNLS-cLSTAPERadj
            nR = nR+1;
            if iNLS < cLSENDWIDEadj
                nNEW = nNEW + 1;
                mmNSFinf(:,:,nR+nS) = [ mmNLSadj(1:cNLS, 1, cLSadj)+cNLS*nNEW+nNIDmax, planePosition(mmNLSadj(:, 2:4, nENDS-iPT), mPLSadj(cLSENDWIDEadj, :)) ];
            else
                mmNSFinf(:,:,nR+nS) = [ mmNLSadj(1:cNLS, 1, nENDS-iPT), planePosition(mmNLSadj(:, 2:4, nENDS-iPT), mPLSadj(cLSENDWIDEadj, :)) ];
            end
            mNSFP_add = [ mNSFP_add; mPLSadj(iNLS, :) ];
        end
        nR = nR+1;
        mmNSFinf(:,:,nR+nS) = [ mmNLSadj(1:cNLS, 1, iNLS), planePosition(mmNLSadj(:, 2:4, iNLS), mPLSadj(iNLS, :)) ];
        mNSFP_add = [ mNSFP_add; mPLSadj(iNLS, :) ];
        % if first additional point log first for straddle mesh ring
        % collate remaining additional points to add ring of elements
        % if no additional points log outer segment point for straddle mesh ring
        % otherwise collate outer segment points to add ring of elements
    end
    
    for iS = 1:nS
        [ mNFSN_2d ] = gen_internal_nodes(mmNSFinf(:,2:3,iS), mNIC_sp, mNIA_sp, mTrans, iNLS, dbg);
        % add nodes to respective collections
        cFS = cFS+1;
        mmNFS(1:cNLS, :, cFS) = [ mmNSFinf(:,1,iS), planePoint(mNSFP_add(iS,:), mmNSFinf(:,2:3,iS)) ];
        mmNFS(cNLS+1:cNLS+cNII, :, cFS) = [ mNIA_sp(1:cNII, 4)+nNLIB+nNLIS*(cFS-1), planePoint(mNSFP_add(iS,:), mNFSN_2d) ];
        mNF = [ mNF; mmNFS(:, :, cFS) ];
        
        % for N-Batch between RS-(RE-RS) and RE
        if cFS > 1
            for iEFS = cFS-1:cFS-1
                mEFNS = [];
                for iEI = 1:size(mEI,1)
                    mEFN = [mEI(iEI,1)+nNLIS*(iEFS)]; % Element ID
                    for iE = 2:size(mEI,2)
                        if mEI(iEI, iE) < cNI1IS % circumfrential node
                            if mEI(iEI, iE) < 200 % circumfrential nodes of 1st or 2nd segment
                                mEFN = horzcat(mEFN, mmNFS(mEI(iEI, iE)-100, 1, iEFS));
                            else
                                mEFN = horzcat(mEFN, mmNFS(mEI(iEI, iE)-200, 1, iEFS+1));
                            end
                        else % internal node
                            if mEI(iEI, iE) < 20000
                                mEFN = horzcat(mEFN, mmNFS(mEI(iEI, iE)-10000+cNLS, 1, iEFS));
                            else
                                mEFN = horzcat(mEFN, mmNFS(mEI(iEI, iE)-20000+cNLS, 1, iEFS+1));
                            end
                        end
                    end
                    mEFNS = [mEFNS; mEFN];
                end
                mmEFS(:, :, iEFS) = mEFNS;
                mEF = [ mEF; mEFNS ];
            end
        end
    end    
    
    for iR = nS+1:nS+nR
        cFSO = cFSO + 1;
        mmNFSO(1:cNLS, :, cFSO) = [ mmNSFinf(:,1,iR), planePoint(mNSFP_add(iR,:), mmNSFinf(:,2:3,iR)) ];
        mNF = [ mNF; mmNFSO(:, :, cFSO) ];
        
        if iR-nS == 1
            mmEFSO(:, :, cFSO) = gen_outer_elem(mEF(end,1), mmNFS(1:cNLS,1,cFS-1), mmNFS(1:cNLS,1,cFS-2), mmNFS(1:cNLS,1,cFS), mmNFSO(1:cNLS, :, cFSO));
        elseif iR-nS == 2
            mmEFSO(:, :, cFSO) = gen_outer_elem(mEF(end,1), mmNFS(1:cNLS,1,cFS-2), mmNFSO(1:cNLS, :, cFSO-nR), mmNFSO(1:cNLS, :, cFSO-1), mmNFSO(1:cNLS, :, cFSO));
        else
            mmEFSO(:, :, cFSO) = gen_outer_elem(mEF(end,1), mmNFSO(1:cNLS, :, cFSO-nR-1), mmNFSO(1:cNLS, :, cFSO-nR), mmNFSO(1:cNLS, :, cFSO-1), mmNFSO(1:cNLS, :, cFSO));
        end
        mEF = [ mEF; mmEFSO(:, :, cFSO) ];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Write out nodes and elements
szFlOut = strcat('fluid_', szFlLumen);
fid = fopen(strcat(szDir, '\', szFlOut),'wt');
fprintf(fid, '*Heading\n**\n');
fprintf(fid, '** PARTS\n**\n');
fprintf(fid, '*Part, name=blood\n');
fprintf(fid, '*Node\n');
for iNF = 1:size(mNF,1)
    fprintf(fid, '%u, %f, %f, %f\n', mNF(iNF, 1), mNF(iNF, 2), mNF(iNF, 3), mNF(iNF, 4));
end

fprintf(fid, '*Element, type=FC3D8\n');
for iE = 1:size(mEF,1)
    fprintf(fid, '%u, %u, %u, %u, %u, %u, %u, %u, %u\n', mEF(iE, 1), mEF(iE, 2), mEF(iE, 3), mEF(iE, 4), mEF(iE, 5), mEF(iE, 6), mEF(iE, 7), mEF(iE, 8), mEF(iE, 9));
end

if 1 == 0
    % write out key sets
    for i = 1:size(mEF,1)
        if i == 1
            fprintf(fid,'*Elset, elset=all\n');
        end

        if (i/16) == round(i/16) || i == size(mEF,1)
            fprintf(fid,'%8u\n', mEF(i,1));
        else
            fprintf(fid,'%8u,', mEF(i,1));
        end    
    end

    % inlet set - node and element
    for i = 1:size(mmNFS,1)
        if i == 1
            fprintf(fid,'*Nset, nset=set_n_inlet\n');
        end

        if (i/16) == round(i/16) || i == size(mmNFS,1)
            fprintf(fid,'%8u\n', mmNFS(i,1, 1));
        else
            fprintf(fid,'%8u,', mmNFS(i,1, 1));
        end    
    end
    for i = 1:size(mmEFS,1)
        if i == 1
            fprintf(fid,'*Elset, elset=set_e_inlet\n');
        end

        if (i/16) == round(i/16) || i == size(mmEFS,1)
            fprintf(fid,'%8u\n', mmEFS(i,1, 1));
        else
            fprintf(fid,'%8u,', mmEFS(i,1, 1));
        end    
    end

    % outlet set - node and element
    for i = 1:size(mmNFS,1)
        if i == 1
            fprintf(fid,'*Nset, nset=set_n_outlet\n');
        end

        if (i/16) == round(i/16) || i == size(mmNFS,1)
            fprintf(fid,'%8u\n', mmNFS(i,1, end));
        else
            fprintf(fid,'%8u,', mmNFS(i,1, end));
        end    
    end
    for i = 1:size(mmEFS,1)
        if i == 1
            fprintf(fid,'*Elset, elset=set_e_outlet\n');
        end

        if (i/16) == round(i/16) || i == size(mmEFS,1)
            fprintf(fid,'%8u\n', mmEFS(i,1, end));
        else
            fprintf(fid,'%8u,', mmEFS(i,1, end));
        end    
    end

    for i = 1:size(mNL,1)
        if i == 1
            fprintf(fid,'*Nset, nset=wall\n');
        end

        if (i/16) == round(i/16) || i == size(mNL,1)
            fprintf(fid,'%8u\n', mNL(i,1));
        else
            fprintf(fid,'%8u,', mNL(i,1));
        end    
    end

end

fprintf(fid, '*End Part\n**\n**\n');
fprintf(fid, '** ASSEMBLY\n**\n');
fprintf(fid, '*Assembly, name=Assembly\n**\n');
fprintf(fid, '*Instance, name=blood-1, part=blood\n');
fprintf(fid, '*End Instance\n**\n');
fprintf(fid, '*End Assembly\n');
fclose(fid)
