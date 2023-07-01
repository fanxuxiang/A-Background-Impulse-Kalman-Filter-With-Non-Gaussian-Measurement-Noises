close all;
clear all;
q=1;
Len=1000;
%   X = stblrnd(alpha,beta,gamma,delta,M,N,..) %Generates an M by N by.. array of S(alpha,beta,gamma,delta) random variables .
MC=100;
for mn=1:MC
% v = stblrnd(1.2,0,3,0,q,Len);%alpha ·Ö²¼ÔëÉù
        v1 = randn(q,Len) * .1;
        v2 = randn(q,Len) * 50;
        for ii=1:q
            numrand = rand(1,Len);
            v(ii,:) = (numrand>=0.8).*v2(ii,:) + (numrand<0.8).*v1(ii,:);
        end % »ìºÏ¸ßË¹³å»÷ÔëÉù
[pai_1,R1_1,R2_1]=EM_GMM2(v);
R_1=var(v);
pai(:,mn)=pai_1;R1(mn)=R1_1;R2(mn)=R2_1;R(mn)=R_1;
end
pai_mn=mean(pai,2);
R1_mn=mean(R1);
R2_mn=mean(R2);
R_mn=mean(R);
plot(v);
xlabel('Samples'),ylabel('Amplitude');