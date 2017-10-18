% clc; clear all ;

function [mEFSO] = gen_outer_elem (nEmax, m1, m2, m3, m4)
    mEFSO = [];
    for iE = 1:size(m1,1)
        if iE < size(m1, 1)
            mEFO = [ nEmax+iE, ...
                    m1(iE,1), m1(iE+1,1), m2(iE+1,1), m2(iE,1),  ...
                    m3(iE,1), m3(iE+1,1), m4(iE+1,1), m4(iE,1) ];
        else
            mEFO = [ nEmax+iE, ...
                    m1(iE,1), m1(1,1), m2(1,1), m2(iE,1), ...
                    m3(iE,1), m3(1,1), m4(1,1), m4(iE,1) ];
        end
        mEFSO = [ mEFSO; mEFO ];
    end
