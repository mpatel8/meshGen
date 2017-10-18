clc; clear all ;
dbg = 0;

% generate fluid mesh
% [ id (1),  x,y,z (2:4)  r,t,p (5:7)  r,t (8:9) ]

% read in lumen mesh [ node-set NL and S4R-element-set EL ]
szDir = 'mesh';
szFlLumen = 'seg-1-36.inp';
szFlLumen = 'seg-5.inp';
szFlLumen = 'seg-36-end.inp';
%szFlLumen = 'seg-1.inp';
%szFlLumen = 'seg-10-end.inp';

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
cLS = size(mNL, 1)/cNLS;
cLSRS = 0; cLSRE = 0;
mNLS = []; mPLS = [];
for iNLS = 1:cLS
    mNLS(:, :, iNLS) = [ mNL((iNLS-1)*cNLS+1:iNLS*cNLS, 1:4) zeros(cNLS, 5)] ;
    p = fitPlane(mNLS(:, 2:4, iNLS));
    if iNLS > 1
        if abs(distancePoints3d(mean(mNLS(:, 2:4, iNLS-1)), mean(mNLS(:, 2:4, iNLS)))) < 0.001
            if cLSRS == 0
                cLSRS = iNLS;
            end
        elseif cLSRS ~= 0 && cLSRE == 0
            if abs(distancePoints3d(mean(mNLS(:, 2:4, iNLS-1)), mean(mNLS(:, 2:4, iNLS)))) > 0.001
                cLSRE = iNLS-1;
            end
        end
    end
                
  mPLS = [ mPLS; p ];
  mNLS(:, 8:9, iNLS) = planePosition(mNLS(:, 2:4, iNLS), mPLS(iNLS, :));
end
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
mNIC_sp = mNI1_sp(1:cNLS, :, :);
mNIA_sp = mNI1_sp(cNLS+1:cNI, :, :);

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
%        if sign(a360(mNIA_sp(i, 2))-a360(mNIC_sp(j, 2))) ~= sign(a360(mNIA_sp(i, 2))-a360(mNIC_sp(j+1, 2)))
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
mNFS = []; mNF = [];
nLSadj = cLSRS - (cLSRE - cLSRS) - 1; cFS = 0; nMid = 0;
for iNLS = 1:cLS
    if iNLS < cLSRS || iNLS > cLSRE
        cFS = cFS + 1;
        if iNLS >= nLSadj && iNLS < cLSRS % in the adjusted zone
            iNLS_adj = iNLS+nLSadj; nMid = nMid + 1;
            mNLS_2d_adj = planePosition(mNLS(:, 2:4, iNLS_adj), mPLS(cLSRS-1, :));
            p_mid = [ mPLS(iNLS-1, 1:3)+0.5*(mPLS(iNLS, 1:3)-mPLS(iNLS-1, 1:3)) mPLS(iNLS, 4:9) ];
        else
            iNLS_adj = iNLS;
            mNLS_2d_adj = mNLS(:, 8:9, iNLS_adj);
            p_mid = [];
        end

        [ mNFSN_2d ] = gen_internal_nodes(mNLS_2d_adj, mNIC_sp, mNIA_sp, mTrans, iNLS, dbg);

        nNLI_base = nNLIB+nNLIS*(iNLS-1);

        % IF in replase section
        % - add a plane based on XYZ = XYZ-1 + (XYZ-2 - XYZ-1)/2, dXdYdZ, XnYnZn
        % - add mNFS based on this plane NOTE need new IDs for outer nodes
        % increment cFS by one
        if ~isempty(p_mid)
            mNFS(1:cNLS, :, cFS) = [ mNLS(1:cNLS, 1, cLS)+cNLS*nMid, planePoint(p_mid, mNLS_2d_adj) ];
            mNFS(cNLS+1:cNLS+cNII, :, cFS) = [ mNIA_sp(1:cNII, 4)+nNLIB+nNLIS*(cFS-1), planePoint(p_mid, mNFSN_2d) ];
            mNF = [ mNF; mNFS(:, :, cFS) ];
            cFS = cFS+1; nMid = nMid + 1;
            mNFS(1:cNLS, :, cFS) = [ mNLS(1:cNLS, 1, cLS)+cNLS*nMid, planePoint(mPLS(iNLS, :), mNLS_2d_adj) ];
            mNFS(cNLS+1:cNLS+cNII, :, cFS) = [ mNIA_sp(1:cNII, 4)+nNLIB+nNLIS*(cFS-1), planePoint(mPLS(iNLS,:), mNFSN_2d) ];
            mNF = [ mNF; mNFS(:, :, cFS) ];
        else
            mNFS(1:cNLS, :, cFS) = [ mNLS(1:cNLS, 1, iNLS_adj), planePoint(mPLS(iNLS, :), mNLS_2d_adj) ];
            mNFS(cNLS+1:cNLS+cNII, :, cFS) = [ mNIA_sp(1:cNII, 4)+nNLIB+nNLIS*(cFS-1), planePoint(mPLS(iNLS,:), mNFSN_2d) ];
            mNF = [ mNF; mNFS(:, :, cFS) ];
        end
    end
end

% for N-Batch between RS-(RE-RS) and RE
mEF = []; mEFS = [];
for iEFS = 1:cFS-1
    mEFNS = [];
    for iEI = 1:size(mEI,1)
        mEFN = [mEI(iEI,1)+nNLIS*(iEFS)]; % Element ID
        for iE = 2:size(mEI,2)
            if mEI(iEI, iE) < cNI1IS % circumfrential node
                if mEI(iEI, iE) < 200 % circumfrential nodes of 1st or 2nd segment
                    mEFN = horzcat(mEFN, mNFS(mEI(iEI, iE)-100, 1, iEFS));
                else
                    mEFN = horzcat(mEFN, mNFS(mEI(iEI, iE)-200, 1, iEFS+1));
                end
            else % internal node
                if mEI(iEI, iE) < 20000
                    mEFN = horzcat(mEFN, mNFS(mEI(iEI, iE)-10000+cNLS, 1, iEFS));
                else
                    mEFN = horzcat(mEFN, mNFS(mEI(iEI, iE)-20000+cNLS, 1, iEFS+1));
                end
            end
        end
        mEFNS = [mEFNS; mEFN];
    end
    mEFS(:, :, iEFS) = mEFNS;
    mEF = [ mEF; mEFNS ];
end

% Fill in the outer sections of the end
cFSO = 0; nLAYER = cLSRE-cLSRS;
for iNFSO = 1:nLAYER
    nPTSadd = 0;
    mNMID_2d = planePosition(mNLS(1:cNLS, 2:4, cLSRE-iNFSO), mPLS(cLSRS-1, :));
    for iMID = 1:iNFSO
        cFSO = cFSO + 1;  nMid = nMid + 1; nPTSadd = nPTSadd+1;

        mNMID = planePoint(mPLS(cLSRS-iNFSO+(iMID-1), :), mNMID_2d);
        mNFSO(1:cNLS, :, cFSO) = [ mNLS(1:cNLS, 1, cLS)+cNLS*nMid mNMID ];
        mNF = [ mNF; mNFSO(:, :, cFSO) ];

        if iMID == 1
            mEFSO(:, :, cFSO) = gen_outer_elem(mEF(end,1), mNFS(1:cNLS,1,cLSRE-(2*iNFSO-1)), mNFS(1:cNLS,1,cLSRE-(2*iNFSO)), mNFS(1:cNLS,1,cLSRE-(iNFSO-1)*2), mNFSO(1:cNLS, :, cFSO));
            mEF = [ mEF; mEFSO(:, :, cFSO) ];
            m34 = [ mNFS(1:cNLS,1,cLSRE-(iNFSO-1)*2), mNFSO(1:cNLS, :, cFSO) ];
         else
             mEFSO(:, :, cFSO) = gen_outer_elem(mEF(end,1), m34(:,1), m34(:,2), mNFSO(1:cNLS, 1, cFSO-iNFSO), mNFSO(1:cNLS, 1, cFSO));
             mEF = [ mEF; mEFSO(:, :, cFSO) ];
             m34 = [ mNFSO(1:cNLS, 1, cFSO-iNFSO), mNFSO(1:cNLS, 1, cFSO) ];
        end
    end
end

% Add remining outer end segments
m34 = [mNFS(1:cNLS,1,cLSRS-(cLSRE-cLSRS)-1), mNFS(1:cNLS,1,cLSRS-(cLSRE-cLSRS)-2)]; nPTSadd = 0;
for iNFSO = cLSRS-(cLSRE-cLSRS)-1:cLSRS-1       
    cFSO = cFSO + 1; nPTSadd = nPTSadd + 1;
    mNFSO(1:cNLS, :, cFSO) = [ mNLS(1:cNLS, 1:4, iNFSO) ];
    mNF = [ mNF; mNFSO(:, :, cFSO) ];

    if nPTSadd == 1
        mEFSO(:, :, cFSO) = gen_outer_elem(mEF(end,1), mNFS(1:cNLS,1,cLSRS-(cLSRE-cLSRS)-1), mNFS(1:cNLS,1,cLSRS-(cLSRE-cLSRS)-2), mNFS(1:cNLS,1,cLSRS-(cLSRE-cLSRS)), mNFSO(1:cNLS, :, cFSO));
        mEF = [ mEF; mEFSO(:, :, cFSO) ];
        m34 = [ mNFS(1:cNLS,1,cLSRE-(iNFSO-1)*2), mNFSO(1:cNLS, :, cFSO) ];
     else
         mEFSO(:, :, cFSO) = gen_outer_elem(mEF(end,1), m34(:,1), m34(:,2), mNFSO(1:cNLS, 1, iNFSO-1), mNFSO(1:cNLS, 1, cFSO));
         mEF = [ mEF; mEFSO(:, :, cFSO) ];
         m34 = [ mNFSO(1:cNLS, 1, iNFSO-1), mNFSO(1:cNLS, 1, cFSO) ];
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
