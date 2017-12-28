%script to generate scaling for Duke and Ella Virtual Family
% models to make samples from a normal distribution based 
% on FAA data from 1993 published here: 
% http://www.dtic.mil/dtic/tr/fulltext/u2/a268661.pdf

%% start with original dims of Duke and Ella:
%Sellion-Menton Lengths:
DukeSM=12/2.54;
EllaSM=10.5/2.54;
%Head breadth:
DukeHB=(71-41)*.5/2.54;
EllaHB=(70-41)*.5/2.54;
%Head Length:
DukeHL=(55-14)*.5/2.54;
EllaHL=(55-16)*.5/2.54;
%% parameters of normal distribution of head measurements
%SM lengths(in):
femSMavg=4.37;
femSMstd=0.25;
malSMavg=4.76;
malSMstd=0.29;
%Head Length:
femHLavg=7.4;
femHLstd=0.29;
malHLavg=7.86;
malHLstd=0.30;
%Head Breadth:
femHBavg=5.74;
femHBstd=0.22;
malHBavg=6.00;
malHBstd=0.21;
%% generate a bunch of head sizes for M/F in 3 non-correlated dims
Nsub=100;
NF=Nsub/2;
NM=Nsub/2;
femHBs = femHBavg + femHBstd.*randn(NF,1);
malHBs = malHBavg + malHBstd.*randn(NF,1);

femHLs = femHLavg + femHLstd.*randn(NF,1);
malHLs = malHLavg + malHLstd.*randn(NF,1);

femSMs = femSMavg + femSMstd.*randn(NF,1);
malSMs = malSMavg + malSMstd.*randn(NF,1);
%% figure out factor to convert current models to the average:
DukeHBfac=malHBavg/DukeHB;
EllaHBfac=femHBavg/EllaHB;
DukeHLfac=malHLavg/DukeHL;
EllaHLfac=femHLavg/EllaHL;
DukeSMfac=malSMavg/DukeSM;
EllaSMfac=femSMavg/EllaSM;
%% then convert the head normal distrib. into factors to mult. by
femHBfacs=femHBs./femHBavg;
malHBfacs=malHBs./malHBavg;
femHLfacs=femHLs./femHLavg;
malHLfacs=malHLs./malHLavg;
femSMfacs=femSMs./femSMavg;
malSMfacs=malSMs./malSMavg;
%% add in the Duke-Ella avg factors:
EllaLRfacs=femHBfacs.*EllaHBfac;
DukeLRfacs=malHBfacs.*DukeHBfac;
EllaAPfacs=femHLfacs.*EllaHLfac;
DukeAPfacs=malHLfacs.*DukeHLfac;
EllaHFfacs=femSMfacs.*EllaSMfac;
DukeHFfacs=malSMfacs.*DukeSMfac;

keyboard

%% save out
save('normal_head_dist','Duke*facs','Ella*facs');