%20171009, adaptive notch filter. call function spANFilter
% 
% fs=500;%nfft=2048;winPSD=0.6;winPl=round(winPSD*fs);overlapPSD=0.4;overlapPl=round(overlapPSD*fs);
% stp=1;%round(handles.RawData.sttime*handles.RawData.fs);
% edp=25*fs;%round(handles.RawData.endtime*handles.RawData.fs);
% ftype='.txt';
% filepathway='D:\MATLAB_work\1707_TS\ts20160101cyx\';%
% txtfile=[filepathway,'*',ftype];
% file=dir(txtfile);
% filename=file(1).name;
% filename(end-3:end)=[];%when read in txt files. the file contains the filetype .txt at the end of "file.name"
% impfilename=[filepathway,filename,ftype];
% sori=importdata(impfilename,' ',1);
% s=sori.data(stp:edp,:);s(end,:)=[];

function data = spAdaptiveNotchFilter(s, fs, notchband, referband, centrefreq)
%function input: 
%nargin==6, last variable is filename, signal after filtering will be saved into .txt
% s; fs; filename(optional)
% notchband=3;%bandwidth
% referband=6;%frequency band
% centrefreq=50;%frequency need to be notch filtered


%% function
[nl,nchn]=size(s);
f_stop = notchband/2;
f_pass = f_stop+1;
filtpara=[fs, centrefreq, f_stop, f_pass];

% Rp is the adjusted parameter
% Rs is Rp/10;

psd_nfft=16384;
psd_nwin=round(1*fs);
psd_noverlap=round(0.4*fs);
psd_wintype=1;
psd_stfreq=centrefreq-referband;
psd_edfreq=centrefreq+referband;
psdpara=[psd_nfft,psd_noverlap,psd_nwin,psd_wintype,psd_stfreq,psd_edfreq];

CoefRs(1)=0.0001;
LB=0.00005;
UB=10;
% parameter settings %

rescoef=zeros(nl,nchn);
for ichn=1:nchn
    ls=s(:,ichn);
    [params,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(@spANFilter,CoefRs,LB,UB,[],filtpara,psdpara,ls,0);    
    rescoef(:,ichn)=spANFilter(params,filtpara,psdpara,ls,1);
end
% handles.ResultData.coef=rescoef;

data=rescoef;

end %end function
