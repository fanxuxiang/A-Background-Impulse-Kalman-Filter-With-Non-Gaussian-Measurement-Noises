function [P_k,P_k_k_1,X_k]=kalman(F,T,H,Q,R,Z,X0,P0)
X_k_k_1 = F * X0;  %k-1ʱ�̶�kʱ��xֵ��Ԥ��
P_k_k_1 = F*P0*F' + T*Q*T'; %k-1ʱ�̶�kʱ��pֵ��Ԥ��
K_k = P_k_k_1 * H' * inv(H*P_k_k_1*H' + R);%kʱ��kalman�˲�����
X_k = X_k_k_1+K_k*(Z - H*X_k_k_1);
P_k = P_k_k_1 - P_k_k_1* H' * inv(H*P_k_k_1*H' + R) * H * P_k_k_1;
end
