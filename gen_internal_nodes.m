% clc; clear all ;

% [ id,  x,y,x  r,t,p,  r,t ]

function [ mNFSN_2d ] = gen_internal_nodes(mNLS_2d, mNIC_sp, mNIA_sp, mTrans, iNLS, dbg)
    cNLS = size(mNLS_2d, 1); cNII = size(mNIA_sp, 1);
    mFSN_sp = cart2sph2d(horzcat(mNLS_2d, zeros(cNLS, 1)));
    if dbg, figure; hold on, end;
    if dbg, drawPoint3d(sph2cart2d(mFSN_sp)), end;

    % Rotate the points to aling 1st to 0 degrees
    nRotate = mFSN_sp(1, 2);
    for i = 1:cNLS
        nAngle = mFSN_sp(i, 2)-nRotate;
        if nAngle < 0
            nAngle = 360.0 + nAngle;
        elseif nAngle > 360.0
            nAngle = nAngle - 360.0;
        end
        mFSN_sp(i, 2) = nAngle;
    end

    mNFSN_2d = [];
    for j = 1:cNII
        nID = mNIA_sp(j, 5:6);
        dT = (mNIA_sp(j, 2)-mNIC_sp(nID(1), 2));      if abs(dT) > 2.0*360.0/cNLS, dT = dT - sign(dT)*360; end;
        dI = (mNIC_sp(nID(2), 2)-mNIC_sp(nID(1), 2)); if abs(dI) > 2.0*360.0/cNLS, dI = dI - sign(dI)*360; end;
        dA = (mFSN_sp(nID(2), 2)-mFSN_sp(nID(1), 2)); if abs(dA) > 2.0*360.0/cNLS, dA = dA - sign(dA)*360; end;
        nT = mFSN_sp(nID(1), 2)+dA*(dT/dI);
        if j < 2
            if dbg, [ j nT mNIA_sp(j, 5:6) mNIA_sp(j, 2) mNIC_sp(nID(1), 2) mNIC_sp(nID(2), 2) mNLS_2d(nID(1), 2) mNLS_2d(nID(2), 2)], end
        end
        if dbg, mTrans = [ mTrans; [ iNLS, j nT dI dT dA mNIA_sp(j, 5:6) mNIA_sp(j, 2) mNIC_sp(nID(1), 2) mNIC_sp(nID(2), 2) mNLS_2d(nID(1), 2) mNLS_2d(nID(2), 2)] dI dA], end

        % correct the rotation applied to circumfrential nodes to get actual angle
        nT = nT + nRotate;
        if nT > 180.0
            nT = nT - 360.0;
        elseif nT < -180.0
            nT = 360.0 + nT;
        end
        nR = mNIA_sp(j,3)*(mFSN_sp(nID(1),3)*dT/dI + mFSN_sp(nID(2),3)*(1-dT/dI));

        mN_sp = sph2cart2d([ mNIA_sp(j, 1) nT nR ]);
        mNFSN_2d = [ mNFSN_2d; mN_sp(1:2) ];
        if dbg, drawPoint3d(sph2cart2d([ mNIA_sp(j, 1) nT nR ])), end;
    end
end
