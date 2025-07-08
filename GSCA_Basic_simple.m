function [est_W,est_C,est_B,Flag_Converge]=BasicGSCA_simple(Data,W,C,B,ind_sign,Max_iter,Min_limit,Flag_C_Forced)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BasicGSCA_simple() - Simplified Version of BasicGSCA() for Simulation   %                  
% Author: Gyeongcheol Cho                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:                                                        %
%   Data = an N by J matrix of scores for N individuals on J indicators   %
%   W = a J by P matrix of weight parameters                              %
%   C = a P by J matrix of loading parameters                             %
%   B = a P by P matrix of path coefficients                              %
%   ind_sign = a P by 1 vector whose p-th element represents the number   %
%               of the sign-fixing indicator for the p-th component       % 
%   N_Boot = Integer representing the number of bootstrap samples for     %
%            calculating standard errors (SE) and 95% confidence          %
%            intervals (CI)                                               %
%   Max_iter = Maximum number of iterations for the Alternating Least     % 
%              Squares (ALS) algorithm                                    %
%   Min_limit = Tolerance level for ALS algorithm                         %
%   Flag_C_Forced = Logical value to determine whether to force loadings  %
%                   for canonical components to be estimated. It won't    %
%                   change overall results but additionally allows for    %
%                   estimating loadings for canonical components if some  %
%                   of the components are canonical. If every component   %
%                   is nomological, the results remain unchanged          %
%                   regardless of its value.                              %                        %
% Output arguments:                                                      
%     est_W: a J by P matrix of weight estimates                          %
%     est_C: a P by J matrix of loading estimates                         %
%     est_B: a P by P matrix of path coefficient estimates                %
%     Flag_Converge = Logical value indicating whether the ALS algorithm  %
%                 converges within the maximum number of iterations       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Z=Data;  
    [N,J]=size(Z);
    P=size(B,1);
    T=J+P;
    W0=W~=0;
    C0=C~=0; 
    B0=B~=0;
    ind_Cdep=sum(C0,1)>0; Jy = sum(ind_Cdep,2); loc_Cdep=find(ind_Cdep); ind_Cdep_post=sum(C0_post,1)>0;
    ind_Bdep=sum(B0,1)>0; Py = sum(ind_Bdep,2); loc_Bdep=find(ind_Bdep);
    ind_Adep=[ind_Cdep, ind_Bdep]; ind_Adep_post=[ind_Cdep_post, ind_Bdep];
    C0_post= C0;
    if Flag_C_Forced; C0_post=W0'; end
    loc_w_t=cell(1,P);
    loc_b_t=cell(1,P);
    for p=1:P
        loc_w_t{1,p}=find(W(:,p))';
        loc_b_t{1,p}=find(B(:,p))';
    end
    loc_c_t=cell(1,J);
    for j=1:J
        loc_c_t{1,j}=find(C(:,j))';
    end
    W(W0)=1;
    [est_W,est_C,est_B,~,Flag_Converge]=ALS_Basic(Z,W,W0,C0,B0,ind_sign,ind_Adep,ind_Adep_post,Min_limit,Max_iter,Flag_C_Forced,C0_post,N,J,P,T,Jy,Py,loc_Cdep,loc_Bdep);      
end