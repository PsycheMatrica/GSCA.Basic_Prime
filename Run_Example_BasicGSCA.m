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
[INI,TABLE,ETC]=BasicGSCA(Data{:,:},W0,C0,B0,ind_sign,N_Boot,Max_iter,Min_limit,Flag_C_Forced,Flag_Parallel);
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