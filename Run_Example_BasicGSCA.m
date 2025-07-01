%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Illustration for GSCA.Basic_Prime package                                %
%   Author: Gyeongcheol Cho                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:                                                            %
%   - This code aims to illustrate how to use GSCA.Basic_Prime package.    %
%   - The dataset is a replica of the ACSI data used in Cho & Hwang       %
%     (2024).                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References                                                              %
%     * Hwang, H. & Takane, Y. (2004). Generalized structured component   %
%         analysis. Psychometrika, 69(1), 81-99.                          %
%     * Cho, G., Hwang, H. Generalized Structured Component Analysis      %
%         Accommodating Convex Components: A Knowledge-Based Multivariate %
%         Method with Interpretable Composite Indexes. Psychometrika 89,  %
%         241–266 (2024). https://doi.org/10.1007/s11336-023-09944-3      %   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

help BasicGSCA()

Data=readtable('ACSI_774_Replica.csv');

%% If missing values are included 
Z0=Data{:,:};
W0=[1 1 1 0 0 0 0 0 0 0 0 0 0 0 ; ...
   0 0 0 1 1 1 0 0 0 0 0 0 0 0 ; ...
   0 0 0 0 0 0 1 1 0 0 0 0 0 0 ; ...
   0 0 0 0 0 0 0 0 1 1 1 0 0 0 ; ...
   0 0 0 0 0 0 0 0 0 0 0 1 0 0 ; ...
   0 0 0 0 0 0 0 0 0 0 0 0 1 1 ]';
C0=zeros(6,14);%W0';
B0=[0 1 1 1 0 0;...
    0 0 1 1 0 0;...
    0 0 0 1 0 0;...
    0 0 0 0 1 1;...
    0 0 0 0 0 1;...
    0 0 0 0 0 0];
ind_sign=[1,4,7,9,12,13];
N_Boot=1000;
Max_iter = 1000;
Min_limit = 10^(-8);
Flag_C_Forced = true;
Flag_Parallel = false;
Opt_Missing = 0;
Results=BasicGSCA(Z0,W0,C0,B0,ind_sign,N_Boot,Max_iter,Min_limit,Flag_C_Forced,Flag_Parallel,Opt_Missing);
INI=Results.INI;
TABLE=Results.TABLE;
ETC=Results.ETC;

INI
INI.GoF
INI.Converge
INI.iter
INI.W
INI.C
INI.B
INI.R2_m
INI.R2_s
TABLE
TABLE.W
TABLE.C
TABLE.B
ETC

%% If missing values are included in data
Percentage_deleted=.10;
[N,J]=size(Z0);
ind_missing=rand(N,J)<=Percentage_deleted;
Z0_miss=Z0;
Z0_miss(ind_missing)=NaN;
Opt_Missing = 4; 
Results=BasicGSCA(Z0_miss,W0,C0,B0,ind_sign,N_Boot,Max_iter,Min_limit,Flag_C_Forced,Flag_Parallel,Opt_Missing);
INI=Results.INI;
TABLE=Results.TABLE;
ETC=Results.ETC;

INI
INI.GoF
INI.Converge
INI.iter
INI.W
INI.C
INI.B
INI.R2_m
INI.R2_s
TABLE
TABLE.W
TABLE.C
TABLE.B
ETC