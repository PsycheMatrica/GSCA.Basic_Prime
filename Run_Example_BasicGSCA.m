%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Illustration for BasicGSCA Prime package                                %
%   Author: Gyeongcheol Cho                                               %
%   Last Revision Date: September 24, 2024                                %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:                                                            %
%   - This code aims to illustrate how to use BasicGSCA Prime package.    %
%   - The dataset is a replica of the ACSI data used in Cho (submitted).  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References                                                              %
%     * Hwang, H. & Takane, Y. (2004). Generalized structured component   %
%         analysis. Psychometrika, 69(1), 81-99.                          %
%     * Cho, G. (submitted). Predictor exclusion threshold: A criterion   %
%         for determining predictor relevance in regularized generalized  % 
%         structured component analysis.                                  %   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

help BasicGSCA()

Data=readtable('ACSI_Comp_Replica.csv');
W0=[1 1 1 0 0 0 0 0 0 0 0 0 0 0 ; ...
   0 0 0 1 1 1 0 0 0 0 0 0 0 0 ; ...
   0 0 0 0 0 0 1 1 0 0 0 0 0 0 ; ...
   0 0 0 0 0 0 0 0 1 1 1 0 0 0 ; ...
   0 0 0 0 0 0 0 0 0 0 0 1 0 0 ; ...
   0 0 0 0 0 0 0 0 0 0 0 0 1 1 ]';
C0=W0';
B0=[0 1 1 1 0 0;...
    0 0 1 1 0 0;...
    0 0 0 1 0 0;...
    0 0 0 0 1 1;...
    0 0 0 0 0 1;...
    0 0 0 0 0 0];
N_Boot=1000;
Max_iter = 1000;
Min_limit = 10^(-8);
Flag_Parallel = true;
[INI,TABLE,ETC]=BasicGSCA(Data{:,:},W0,C0,B0,N_Boot,Max_iter,Min_limit);
INI
INI.GoF
INI.W
INI.C
INI.B
INI.R2_m
Ini.R2_s
TABLE
TABLE.W
TABLE.C
TABLE.B
ETC