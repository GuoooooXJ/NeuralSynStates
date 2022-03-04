function out=spANFilter(CoefRs,msFiltP,msPSDP,data,OutType)
%msFiltP: parameters for bandstop filter
%msPSDP: parameters for power spectra analysis

Rs=CoefRs(1);
Rp=Rs/10;

Fs=msFiltP(1);
Fc=msFiltP(2);
Ds=msFiltP(3);
Dp=msFiltP(4);

Wp=[Fc-Dp Fc+Dp]/Fs*2;
Ws=[Fc-Ds Fc+Ds]/Fs*2;

[N, Wn] = buttord(Wp, Ws, Rp, Rs);
[B,A] = butter(N,Wn,'stop');
X=filtfilt(B,A,data);%zero phase filtering. 

%[1] Oppenheim, A.V., and R.W. Schafer, Discrete-Time Signal Processing, Prentice-Hall, 1989, pp.311-312.
%[2] Mitra, S.K., Digital Signal Processing, 2nd ed., McGraw-Hill, 2001, Sections 4.4.2 and 8.2.5. 
%[3] Gustafsson, F., Determining the initial states in forward-backward filtering, IEEE Transactions on Signal Processing,
%    April 1996, Volume 44, Issue 4, pp.988¨C992.

NFFT=msPSDP(1);
NOVERLAP=msPSDP(2);
WinLen=msPSDP(3);
WinType=msPSDP(4);
STF=msPSDP(5);
EDF=msPSDP(6);

switch WinType
case 1
    w=window(@hann,WinLen);
case 2
    w=window(@hamming,WinLen);
case 3
    w=window(@bartlett,WinLen);
case 4
    w=window(@rectwin,WinLen);
end

X1=detrend(X);
[P1,F] = pwelch(X1,w,NOVERLAP,NFFT,Fs,'onesided');

stF0=find(F>(Fc-Ds), 1 ); 
edF0=find(F<(Fc+Ds), 1, 'last' ); 

stF1=find(F>STF, 1 ); 
edF1=find(F<(Fc-Dp*2), 1, 'last' ); 

stF2=find(F>(Fc+Dp*2), 1 ); 
edF2=find(F<EDF, 1, 'last' ); 

switch OutType
    case 0
        v1=mean(P1(stF0:edF0));
        v2=(mean(P1(stF1:edF1))+mean(P1(stF2:edF2)))/2;
        out=(v1-v2);
    case 1
        out=X;
    case 2
        [P0,F] = pwelch(data,w,NOVERLAP,NFFT,Fs,'onesided');
        out=[F P0 P1];
end