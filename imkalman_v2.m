%u1,u2初始概率
% x0_ikf_cv1,x0_ikf_cv2初始状态
% P1,P2初始协方差
% R1,R2噪声方差
% u_a_1，u_a_2下一时刻的概率
% P_a_1，P_a_2下一时刻的协方差
% X_a_1，X_a_2预测值
function [u_a_1,u_a_2,P_a_1,P_a_2,X_a_1,X_a_2]=imkalman_v2(F,T,H,Q,R1,R2,Z,x0_ikf_cv1,x0_ikf_cv2,P0_ikf_cv1,P0_ikf_cv2,u1,u2,Pi)
c1=Pi(1,1)*u1+Pi(2,1)*u2;
c2=Pi(1,2)*u1+Pi(2,2)*u2;
u11=Pi(1,1)*u1/c1;u12=Pi(1,2)*u1/c2;
u21=Pi(2,1)*u2/c1;u22=Pi(2,2)*u2/c2;
x1_m = x0_ikf_cv1*u11+x0_ikf_cv2*u21;
x2_m = x0_ikf_cv1*u12+x0_ikf_cv2*u22;
p1_k_1=(P0_ikf_cv1+(x0_ikf_cv1-x1_m)*(x0_ikf_cv1-x1_m)')*u11+(P0_ikf_cv2+(x0_ikf_cv2-x1_m)*(x0_ikf_cv2-x1_m)')*u21;
p2_k_1=(P0_ikf_cv1+(x0_ikf_cv1-x2_m)*(x0_ikf_cv1-x2_m)')*u12+(P0_ikf_cv2+(x0_ikf_cv2-x2_m)*(x0_ikf_cv2-x2_m)')*u22;
R_x=T*Q*T';
xke1=F*x1_m;
xke2=F*x2_m;
pke1=F*p1_k_1*F'+R_x;
pke2=F*p2_k_1*F'+R_x;

Zpre1=H*xke1;
Zpre2=H*xke2;
dzz1=Z-Zpre1;
dzz2=Z-Zpre2;
Sv1=H*pke1*H'+R1;
Sv2=H*pke2*H'+R2;
sigma1=c1*1/(2*pi)^(length(Zpre1)/2)/(abs(det(Sv1)))^0.5*exp(-0.5*dzz1'*inv(Sv1)*dzz1);
sigma2=c2*1/(2*pi)^(length(Zpre2)/2)/(abs(det(Sv2)))^0.5*exp(-0.5*dzz2'*inv(Sv2)*dzz2);

sum_sigma=sigma1+sigma2;
if(sum_sigma==0)
%     sigma1=c1;
%     sigma2=c2;
%% 6月23日
    sigma1=1;
    sigma2=0;
    R1=dzz1*dzz1';
else
    sigma1=sigma1/(sum_sigma);
    sigma2=sigma2/(sum_sigma);
end


Gk1=pke1*H'*inv(H*pke1*H'+R1);
Gk2=pke2*H'*inv(H*pke2*H'+R2);

X_a_1=xke1+Gk1*(Z-H*xke1);
X_a_2=xke2+Gk2*(Z-H*xke2);

P_a_1=pke1-Gk1*H*pke1;
P_a_2=pke2-Gk2*H*pke2;

u_a_1=sigma1;u_a_2=sigma2;
end
