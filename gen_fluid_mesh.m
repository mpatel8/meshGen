clc; clear all ;
dbg = 0;

% generate fluid mesh
% [ id (1),  x,y,z (2:4)  r,t,p (5:7)  r,t (8:9) ]

% read in lumen mesh [ node-set NL and S4R-element-set EL ]
szDir = 'mesh';
szFlLumen = 'seg-1-36.inp';
szFlLumen = 'seg-5.inp';
szFlLumen = 'seg-36-end.inp';
szFlLumen = 'seg-10-end.inp';
szFlLumen = 'seg-all.inp';

[ mNL, mEL ] = read_input(strcat(szDir, '\', szFlLumen ));

% determine how many nodes per segment 'N' using first set of nodes that have
% same X value [IDEALLY MAKE AVAILABLE VIA MESH INPUT FILE - OR VIA MATLAB CODE ]
% i = 1; 
% x_org = mNL(i, 2:2);
% while abs(mNL(i+1, 2:2)-x_org) < 1e-4
%     i = i + 1;

% end
% %cNLS = i;
cNLS = 36;

% look through each N nodes and get plane definition
% if plane samilar to prior then its rigid-wide-end (check PointOnPlane function
    % %% TBA %% note which N-batch this starts, RS, ends, RE, number-segments RN, and plane, RN
    % %% TBA %%collate these nodes into list R
% store each N-batch in a separate list, NS
% store each plane info, PS
% note total number of actual segments T
cLS = size(mNL, 1)/cNLS; cLSRS = 0; cLSRE = 0; mmNLS = []; mPLS = [];cLSENDWIDE=10000
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
end
cLSTAPER = cLSENDWIDE - (cLSRE-cLSRS);
nNLIS = 10^(round(log10(max(mNL(:, 1)))+0.5));
nNLIB = 10*nNLIS;

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
mNIA_sp = horzcat(mNIA_sp, zeros(99, 3));
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
mNF = []; mmNFS = [];
mEF = []; mmEFS = [];
nLSadj=cLSRS-(cLSRE-cLSRS)-1; cFS=0; cFSO=0; nNEW=0;
for iNLS = 1:cLS
    % collate full segment info 
    mNSFP_add=[]; % plane to transform the nodes into
    mmNSFinf = []; % ids of outer nodes and their 2D coords reference
    nS = 0; nR = 0;
    if iNLS < cLSTAPER || iNLS > cLSRE % Normal main mesh or constriction section
        nS = nS+1;
        mmNSFinf(:,:,nS) = [ mmNLS(1:cNLS, 1, iNLS), planePosition(mmNLS(:, 2:4, iNLS), mPLS(iNLS, :)) ];
        mNSFP_add = [ mNSFP_add; mPLS(iNLS, :) ];
        
    elseif iNLS <= cLSENDWIDE % between taperring section near end
        % Add taper points based on end points (inner first then on segment)
        % collate full segment info for each
        nS = nS+1; nNEW = nNEW + 1; nENDS = cLSENDWIDE+(iNLS-cLSTAPER)+1;
        mmNSFinf(:,:,nS) = [ mmNLS(1:cNLS, 1, cLS)+cNLS*nNEW, planePosition(mmNLS(:, 2:4, nENDS), mPLS(cLSENDWIDE, :)) ];
        mNSFP_add = [ mNSFP_add; mPLS(iNLS-1, 1:3)+0.5*(mPLS(iNLS, 1:3)-mPLS(iNLS-1, 1:3)) mPLS(iNLS, 4:9) ];

        nS = nS+1;
        if nENDS < cLSRE
            nNEW = nNEW + 1;
            mmNSFinf(:,:,nS) = [ mmNLS(1:cNLS, 1, cLS)+cNLS*nNEW, planePosition(mmNLS(:, 2:4, nENDS), mPLS(cLSENDWIDE, :)) ];
        else
            mmNSFinf(:,:,nS) = [ mmNLS(1:cNLS, 1, nENDS), planePosition(mmNLS(:, 2:4, nENDS), mPLS(cLSENDWIDE, :)) ];
        end
        mNSFP_add = [ mNSFP_add; mPLS(iNLS, :) ];

        % if required add additional mesh points
        for iPT = 1:iNLS-cLSTAPER
            nR = nR+1;
            if iNLS < cLSENDWIDE
                nNEW = nNEW + 1;
                mmNSFinf(:,:,nR+nS) = [ mmNLS(1:cNLS, 1, cLS)+cNLS*nNEW, planePosition(mmNLS(:, 2:4, nENDS-iPT), mPLS(cLSENDWIDE, :)) ];
            else
                mmNSFinf(:,:,nR+nS) = [ mmNLS(1:cNLS, 1, nENDS-iPT), planePosition(mmNLS(:, 2:4, nENDS-iPT), mPLS(cLSENDWIDE, :)) ];
            end
            mNSFP_add = [ mNSFP_add; mPLS(iNLS, :) ];
        end
        nR = nR+1;
        mmNSFinf(:,:,nR+nS) = [ mmNLS(1:cNLS, 1, iNLS), planePosition(mmNLS(:, 2:4, iNLS), mPLS(iNLS, :)) ];
        mNSFP_add = [ mNSFP_add; mPLS(iNLS, :) ];
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

% Write out nodes and elements
szFlOut = 'test.inp';
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

fprintf(fid, '*End Part\n**\n**\n');
fprintf(fid, '** ASSEMBLY\n**\n');
fprintf(fid, '*Assembly, name=Assembly\n**\n');
fprintf(fid, '*Instance, name=blood-1, part=blood\n');
fprintf(fid, '*End Instance\n**\n');
fprintf(fid, '*End Assembly\n');
fclose(fid)
