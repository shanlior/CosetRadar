function [M1, A1, A2, B , B2] = BuildRadarMatrices2(sSimParams)
% Build matrices used to solve CS.
% Output:
%   M1      - size (nKappaSingleTrans x nKappaSingleTrans)

% Input:
%   sSimParams
% Remarks:


% M1:
fftIdx      = sSimParams.kappaSingleTrans ;
fftIdx      = mod(fftIdx, sSimParams.sigLen) + 1;
M1Diag      = sSimParams.H0(fftIdx);
M1          = diag(M1Diag) / sSimParams.T;

tauVec      = 0 : sSimParams.sResolution.timeBinSuperResolution : sSimParams.T - sSimParams.sResolution.timeBinSuperResolution ;
% each column corresponds to ragne grid point
kVec        = sSimParams.kappaSingleTrans; %- ceil(sSimParams.NTag/2);
kVec        = kVec(:);
A1          = exp(-1j * 2 * pi / sSimParams.T * kVec * tauVec);

% each column corresponds to ragne grid point
fVec        = repmat(sSimParams.freqDivisionLocation,1,sSimParams.N);
fVec        =fVec(:);
A2          = exp(-1j * 2 * pi *fVec  * tauVec);

%A(m)=A1*diag(A2(m,:)))

sinThetaVec = -1 : sSimParams.sResolution.sinThetaBin : 1 - sSimParams.sResolution.sinThetaBin;
sinThetaVec = sinThetaVec(:);
% each raw corresponds to azimuth grid point
B          = exp(2j * pi * sinThetaVec * sSimParams.sVec);

% each raw corresponds to azimuth grid point
B2          = exp(2j * pi * sinThetaVec *(fVec.'.*sSimParams.sVec)/sSimParams.fc);

%b(m)=B(:,m)

