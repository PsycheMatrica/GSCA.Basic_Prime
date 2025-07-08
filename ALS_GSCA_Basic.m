function [W,C,B,vec_err,Flag_Converge,iter,Z,CVr]=ALS_Basic(Z,Flag_LS,W,W0,C0,B0,ind_sign,ind_Adep,ind_Adep_post,Min_limit,Max_iter,Flag_C_Forced,C0_post,N,J,P,T,Jy,Py,loc_Cdep,loc_Bdep)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALS_Basic() - MATLAB function to implement the basic ALS algorithm for  %
%               Generalized Structured Component Analysis (GSCA).         %
% Author: Gyeongcheol Cho                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set initial value    
    iter=0;
    improve=100000000;
    dif_c=100000000;   
    eye_J=eye(J);
    flag=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~Flag_LS
% Normalize & Sample Covariance    
    Z=zscore(Z,1)/sqrt(N);    
    cov_Z=Z'*Z;
    for p=1:P
        w_p=W(W0(:,p),p);
        cov_zp=cov_Z(W0(:,p),W0(:,p));
        norm_cv=w_p'*cov_zp*w_p;
        W(W0(:,p),p)=W(W0(:,p),p)/sqrt(norm_cv);
    end
    V=[eye_J, W]; 
% Reduce the computational cost
    Zr=Z;
    if N>J
        [R,flag]=chol(cov_Z);
        if flag==0; Zr=R; end
    end
    CVr = Zr*W;
    C=double(C0); B=double(B0);
    while improve > Min_limit && iter < Max_iter
        iter=iter+1;
        dif_p=dif_c;
    % (2-1) to estimate A given W ... SS(ZV-ZWA)=SS(Psi-GamA)
        [C,B]=Estimation_A(Zr,CVr,C,B,C0,B0,loc_Cdep,loc_Bdep,Jy,Py);
        A=[C,B];
    % (2-2) to estimate W given A 
        [W,V]=Estimation_W(Zr,cov_Z,V,A,W,W0,P,J,T);
        CVr=Zr*W;
        
        ERROR=[Zr CVr] - CVr*A;
        vec_err=sum(ERROR(:,ind_Adep).^2,1);
        dif_c=sum(vec_err,2);   
        improve=dif_p-dif_c; % sum(error_past^2) - sum(error_current^2) 
    end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
% Ini-Step1: Z
    ind_Miss_Z0=isnan(Z);
    loc_Miss_J=find(sum(ind_Miss_Z0,1));
    Jmiss=size(loc_Miss_J,2);
    for jm = 1:Jmiss
        j = loc_Miss_J(jm);
        Z(ind_Miss_Z0(:,j),j)=mean(Z(~ind_Miss_Z0(:,j),j),1);
    end
% Ini-Step2: W and CV
    Z=zscore(Z,1)/sqrt(N);    
    cov_Z=Z'*Z;
    for p=1:P
        w_p=W(W0(:,p),p);
        cov_zp=cov_Z(W0(:,p),W0(:,p));
        norm_cv=w_p'*cov_zp*w_p;
        W(W0(:,p),p)=W(W0(:,p),p)/sqrt(norm_cv);
    end
    V=[eye_J, W]; 
    Zr=Z;
    if N>J
        [R,flag]=chol(cov_Z);
        if flag==0; Zr=R; end
    end
    CVr = Zr*W;
% Ini-Step3: [C,B]
    C=double(C0); B=double(B0);
  %  [C,B]=Estimation_A(Zr,CVr,C,B,C0,B0,loc_Cdep,loc_Bdep,Jy,Py);
    eye_J = eye(J);
    while improve > Min_limit && iter < Max_iter
        iter=iter+1;
        dif_p=dif_c;
    % (2-1) to estimate missing values 
        Q = V - W*[C,B];  
        for jm = 1:Jmiss
            j = loc_Miss_J(jm);
            zj = Z(:,j);
            qj = Q(j,:);
            %{
            Qj0=Q; Qj0(j,:)=0;
            Zj0=Z; Zj0(:,j)=0;            
            PI=-Zj0*Qj0;
            THETA=kron(qj',eye(N));
            THETA_hat=THETA(:,~ind_Miss_Z0(:,j));            
            THETA_mj=THETA(:,ind_Miss_Z0(:,j));
            zj(ind_Miss_Z0(:,j),1)=(THETA_mj'*THETA_mj)\(THETA_mj'*(PI(:)-THETA_hat*zj(~ind_Miss_Z0(:,j),1)));
            Z(:,j) = zj;
            %}
            eye_J_j0=eye_J;eye_J_j0(j,j) = 0;  
            zj(ind_Miss_Z0(:,j),1) = -Z(ind_Miss_Z0(:,j),:)*(eye_J_j0*Q)*(qj'/(qj*qj'));
            zj=zj-mean(zj,1);
            Z(:,j) = (zj/sqrt(zj'*zj))';
         end
         cov_Z=Z'*Z;
         Zr=Z;        
         if N>J
             [R,flag]=chol(cov_Z);
             if flag==0; Zr=R; end
         end
         for p=1:P
             w_p=W(W0(:,p),p);
             cov_zp=cov_Z(W0(:,p),W0(:,p));
             norm_cv=w_p'*cov_zp*w_p;
             W(W0(:,p),p)=W(W0(:,p),p)/sqrt(norm_cv);
         end
         V(:,(J+1):end)=W;
         CVr = Zr*W;
    % (2-2) to estimate A given W ... SS(ZV-ZWA)=SS(Psi-GamA)
        [C,B]=Estimation_A(Zr,CVr,C,B,C0,B0,loc_Cdep,loc_Bdep,Jy,Py);
    % (2-3) to estimate W given A    
        [W,V]=Estimation_W(Zr,cov_Z,V,[C,B],W,W0,P,J,T);
        CVr=Zr*W;

        ERROR=[Zr CVr] - CVr*[C,B];
        vec_err=sum(ERROR(:,ind_Adep).^2,1);
        dif_c=sum(vec_err,2);   
        improve=dif_p-dif_c;  %[dif_p,dif_c,improve] 
     end  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 Flag_Converge=true;
 if (iter == Max_iter) && (improve > Min_limit); Flag_Converge=false; end
 if Flag_C_Forced
    loc_Cdep=1:J;
    C0=C0_post;
    Jy=J;
 end
 [C,B]=Estimation_A(Zr,CVr,C,B,C0,B0,loc_Cdep,loc_Bdep,Jy,Py);
 
 for p=1:P
    if ind_sign(1,p)>0
        if Zr(:,ind_sign(1,p))'*CVr(:,p)<0
            CVr(:,p)=-CVr(:,p);
            W(:,p)=-W(:,p);
            C(p,:)=-C(p,:);
            B(:,p)=-B(:,p);
            B(p,:)=-B(p,:);
        end
    end
 end
 ERROR=[Zr CVr] - CVr*[C,B];
 vec_err=sum(ERROR(:,ind_Adep_post).^2,1);   
 if flag==0; CVr=Z*W*sqrt(N);end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [W,V]=Estimation_W(Zr,cov_Z,V,A,W,W0,P,J,T)
    for p=1:P
        t=J+p;

        V_m_t=V; V_m_t(:,t)=0; % _m_ : minus
        A_m_p=A; A_m_p(p,:)=0; % _p_ : existent
        a_p=A(p,:);

        W_m_p=W; W_m_p(:,p)=0;

        Lam_m_p=W_m_p*A_m_p;

        i_t=zeros(1,T);
        i_t(1,t)=1;

        beta=i_t - a_p;
        delta=Lam_m_p - V_m_t;

        Ksi=kron(beta',Zr);
        Ksi=Ksi(:,W0(:,p));
        vec_y=Zr*delta;
        vec_y=vec_y(:);
        theta=((Ksi'*Ksi)\Ksi')*vec_y;

        norm_cv=theta'*cov_Z(W0(:,p),W0(:,p))*theta;
        theta=theta/sqrt(norm_cv);
        W(W0(:,p),p)=theta;
        V(W0(:,p),J+p)=theta;
    end
end
function [C,B]=Estimation_A(Zr,CV,C,B,C0,B0,loc_Cdep,loc_Bdep,Jy,Py)
    %{ 
        oneshot
        % Psi = Z * V;
        % CV = Z * W; % N by P 
        % Phi=kron(eye(T),CV);
        % Phi=Phi(:,ind_A);
        % a=(Phi'*Phi)\Phi'*Psi(:);
        % A(ind_A)=a;
    %}
    if Jy>0
        for j=loc_Cdep
            CVx = CV(:,C0(:,j));
            C(C0(:,j),j)=(CVx'*CVx)\(CVx'*Zr(:,j));
        end
    end
    if Py>0
        for q=loc_Bdep
            CVx = CV(:,B0(:,q));             
            bq_est=(CVx'*CVx)\(CVx'*CV(:,q));
            B(B0(:,q),q)=bq_est;
        end
    end
end