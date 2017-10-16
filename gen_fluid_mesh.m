clc; clear all ;
debug = 0;

% generate fluid mesh

% read in lumen mesh [ node-set NL and S4R-element-set EL ]
szDir = 'mesh';
szFlLumen = 'seg-1-36.inp';
szFlLumen = 'seg-5.inp';
%szFlLumen = 'seg-all-noEndFace.inp';
%szFlLumen = 'seg-36-end.inp';

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
    mNLS(:, :, iNLS) = mNL((iNLS-1)*cNLS+1:iNLS*cNLS, 1:4);
    p = fitPlane(mNLS(:, 2:4, iNLS));
    if iNLS > 1
        if abs(distancePoints3d(mean(mNLS(:, 2:4, iNLS-1)), mean(mNLS(:, 2:4, iNLS)))) < 0.001
            if cLSRS == 0
                cLSRS = iNLS-1;
            end
        elseif cLSRS ~= 0 && cLSRE == 0
            if abs(distancePoints3d(mean(mNLS(:, 2:4, iNLS-1)), mean(mNLS(:, 2:4, iNLS)))) > 0.001
                cLSRE = iNLS-1;
            end
        end
    end
                
  mPLS(:, :, iNLS) = p;
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
mNIA_sp_assoc = horzcat(mNIA_sp, zeros(99, 2));
mNINC = [];
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
        mNINC(:,:,i) = [ cNLS, 1 ]; 
        mNIA_sp_assoc(i, 4:5) = [ cNLS, 1 ];
    else
        mNINC(:,:,i) = [ j, j+1 ];
        mNIA_sp_assoc(i, 4:5) = [ j, j+1 ];
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
for iNLS = 1:cLS
    %if iNLS >= cLSRS-(cLSRE-cLSRS) || iNLS >= cLSRE
        
        
    %mL_p = fitPlane(mNLS(:, 2:4, iNLS));
    mNLS_2d = planePosition(mNLS(:, 2:4, iNLS), mPLS(:,:,iNLS));
    mNLS_sp = cart2sph2d(horzcat(mNLS_2d, zeros(cNLS, 1)));
    if debug, figure; hold on, end;
    if debug, drawPoint3d(sph2cart2d(mNLS_sp)), end

    nRotate = mNLS_sp(1, 2);
    for i = 1:cNLS
        nAngle = mNLS_sp(i, 2)-nRotate;
        if nAngle < 0
            nAngle = 360.0 + nAngle;
        elseif nAngle > 360.0
            nAngle = nAngle - 360.0;
        end
        mNLS_sp(i, 2) = nAngle;
    end

    mNLI_sp = []; mNLI_2d = []; mNLI_id = [];
    for j = 1:cNII
        nID = mNINC(:, :, j);
        dT = (mNIA_sp(j, 2)-mNIC_sp(nID(1), 2));      if abs(dT) > 2.0*360.0/cNLS, dT = dT - sign(dT)*360; end;
        dI = (mNIC_sp(nID(2), 2)-mNIC_sp(nID(1), 2)); if abs(dI) > 2.0*360.0/cNLS, dI = dI - sign(dI)*360; end;
        dA = (mNLS_sp(nID(2), 2)-mNLS_sp(nID(1), 2)); if abs(dA) > 2.0*360.0/cNLS, dA = dA - sign(dA)*360; end;
        nT = mNLS_sp(nID(1), 2)+dA*(dT/dI);
        if j < 2
            if debug, [ j nT mNINC(:, :, j) mNIA_sp(j, 2) mNIC_sp(nID(1), 2) mNIC_sp(nID(2), 2) mNLS_sp(nID(1), 2) mNLS_sp(nID(2), 2)], end
        end
        if debug, mTrans = [ mTrans; [ iNLS, j nT dI dT dA mNINC(:, :, j) mNIA_sp(j, 2) mNIC_sp(nID(1), 2) mNIC_sp(nID(2), 2) mNLS_sp(nID(1), 2) mNLS_sp(nID(2), 2)] dI dA], end

        nT = nT + nRotate;
        if nT > 180.0
            nT = nT - 360.0;
        elseif nT < -180.0
            nT = 360.0 + nT;
        end

        
        % dR1 = (mNIA_sp(j, 2)-mNIC_sp(nID(2), 2))/dI;
        % dR2 = (mNIC_sp(nID(1), 2)-mNIA_sp(j, 2))/dI;
        % nR = mNIA_sp(j, 3)*(dR1*mNLS_sp(nID(1), 3) + dR2*mNLS_sp(nID(2), 3));

        nR = mNIA_sp(j,3)*(mNLS_sp(nID(1),3)*dT/dI + mNLS_sp(nID(2),3)*(1-dT/dI));

        mNLI_sp = [ mNLI_sp; [ mNIA_sp(j, 1) nT nR ] ];
        mNLI_id = [ mNLI_id; nNLIB+nNLIS*(iNLS-1)+mNI1N(j, 1)-(cNI1IS-1) ];
        mNLI_2d = [ mNLI_2d; sph2cart2d(mNLI_sp(j,:)) ];
        if debug, drawPoint3d(mNLI_2d(j, :)), end
    end
    %mNLI(:, :, iNLS) = horzcat(mNLI_id, planePoint(mL_p, sph2cart2d(mNLI_sp(:,:))));
    mNLI(:, :, iNLS) = horzcat(mNLI_id, planePoint(mPLS(:,:,iNLS), sph2cart2d(mNLI_sp(:,:))));
    mNFS(:, :, iNLS) = [ mNLS(:, :, iNLS); mNLI(:, :, iNLS) ];
end

% for N-Batch between RS-(RE-RS) and RE
mEF = [];
for iELS = 1:cLS-1
    mEFS=[];
    for iEI = 1:size(mEI,1)
        mEFN = [mEI(iEI,1)+nNLIS*(iELS)];
        for iE = 2:size(mEI,2)
            if mEI(iEI, iE) < cNI1IS
                if mEI(iEI, iE) < 200
                    mEFN = horzcat(mEFN, mNLS(mEI(iEI, iE)-100, 1, iELS));
                else
                    mEFN = horzcat(mEFN, mNLS(mEI(iEI, iE)-200, 1, iELS+1));
                end
            elseif mEI(iEI, iE) < 20000
                mEFN = horzcat(mEFN, mEI(iEI, iE)-10000+(nNLIB+nNLIS*(iELS-1)));
            else
                mEFN = horzcat(mEFN, mEI(iEI, iE)-20000+(nNLIB+nNLIS*(iELS)));
            end
        end
                            
    end
    mEF(:, :, iELS) = mEFS;
end

% for N-batch in S-1
    % get FSN and FSN+1
    % calculate and store element references in FSE
szFlOut = 'test.inp';
fid = fopen(strcat(szDir, '\', szFlOut),'wt');
fprintf(fid, '*Heading\n');
fprintf(fid, '**\n');
fprintf(fid, '** PARTS\n');
fprintf(fid, '**\n');
fprintf(fid, '*Part, name=blood\n');
fprintf(fid, '*Node\n');
for iLS = 1:cLS
    for iNLS = 1:cNI
        fprintf(fid, '%u, %f, %f, %f\n', mNFS(iNLS, 1, iLS), mNFS(iNLS, 2, iLS), mNFS(iNLS, 3, iLS), mNFS(iNLS, 4, iLS));
    end
end
% fclose(fid)
% for iLS = 1:cLS
%     dlmwrite(szFlOut,mNFS(:, :, iLS),'precision','%.6f'), '-append');
% end
% fid = fopen(szFlOut,'a+');
fprintf(fid, '*Element, type=FC3D8\n');
% for iLS = 1:cLS-1
%     dlmwrite(szFlOut,mEF(:, :, iLS), '-append');
% end
for iLS = 1:cLS-1
    for iE = 1:cEI
        fprintf(fid, '%u, %u, %u, %u, %u, %u, %u, %u, %u\n', mEF(iE, 1, iLS), mEF(iE, 2, iLS), mEF(iE, 3, iLS), mEF(iE, 4, iLS), mEF(iE, 5, iLS), mEF(iE, 6, iLS), mEF(iE, 7, iLS), mEF(iE, 8, iLS), mEF(iE, 9, iLS));
    end
end
fprintf(fid, '*End Part\n');
fprintf(fid, '**\n');
fprintf(fid, '**\n');
fprintf(fid, '** ASSEMBLY\n');
fprintf(fid, '**\n');
fprintf(fid, '*Assembly, name=Assembly\n');
fprintf(fid, '**\n');
fprintf(fid, '*Instance, name=blood-1, part=blood\n');
fprintf(fid, '*End Instance\n');
fprintf(fid, '**\n');
fprintf(fid, '*End Assembly\n');
fclose(fid)

