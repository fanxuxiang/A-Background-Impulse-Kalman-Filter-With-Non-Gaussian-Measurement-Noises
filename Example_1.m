%对于文章Example1
%注意蒙特卡洛实验次数、噪声参数、MCKF与MEEKF核宽的选择
clear all;clc;
% tic
H = [1 0 0 0;0 0 1 0];
q=2;Len=1000;
N = Len;T = 1;
MC=100;
p=4;
ex1=zeros(MC,N);ex2=zeros(MC,N);
ey1=zeros(MC,N);ey2=zeros(MC,N);% 滤波误差初始化
%蒙特卡罗仿真
for mn=1:MC
    %% 初始x数据产生
    x0 = [1000,10,1000,-10]';
    p0=[ 100 0 0 0;0 1 0 0;0 0 100 0;0 0 0 1;];
    xA = [];
    zA = [x0(1);x0(3)];
    %model-1,匀速运动
    A1 = [1,T,0,0;
        0,1,0,0;
        0,0,1,T;
        0,0,0,1];
    G1=[T^2/2,    0;
        T,      0;
        0,      T^2/2;
        0,      T] ;
    Q1=[0.1^2 0;
        0 0.1^2];
    %model-2,匀速转弯模型
    omgia1=pi/270;
    %  omgia1=pi/90;
    A2=[1 sin(omgia1)/omgia1 0 (cos(omgia1)-1)/omgia1;
        0 cos(omgia1) 0 -sin(omgia1);
        0 -(cos(omgia1)-1)/omgia1 1 sin(omgia1)/omgia1;
        0 sin(omgia1) 0 cos(omgia1);];
    G2=[1/2 0;
        1 0;
        0 1/2;
        0 1;];
    Q2=[0.0144^2 0;
        0 0.0144^2];
    % 产生真实数据
    x = x0;
    for k = 1:350%匀速直线
        x = A1*x + G1*sqrt(Q1)*[randn,randn]';
        xA =[xA x];
    end
%             for k = 1:220%匀速圆周转弯
%                 x = A2*x + G2*sqrt(Q2)*[randn,randn]';
%                 xA =[xA x];
%             end
    for k = 1:650%匀速直线
        x = A1*x + G1*sqrt(Q1)*[randn,randn]';
        xA =[xA x];
    end
    
    %% 量测噪声产生
        v1 = randn(q,Len) * .1;
        v2 = randn(q,Len) * 100;
        for ii=1:q
            numrand = rand(1,Len);
            v(ii,:) = (numrand>=0.9).*v2(ii,:) + (numrand<0.9).*v1(ii,:);
        end % 混合高斯冲击噪声
        v1_noise = randn(q,Len) * .1;
        v2_noise = randn(q,Len) * 100;
        for ii=1:q
            numrand = rand(1,Len);
            v_noise(ii,:) = (numrand>=0.9).*v2_noise(ii,:) + (numrand<0.9).*v1_noise(ii,:);
        end % 接收混合高斯冲击噪声
    
%     X = stblrnd(alpha,beta,gamma,delta,M,N,..) %Generates an M by N by.. array of S(alpha,beta,gamma,delta) random variables .
%         v = stblrnd(1.2,0,2,0,q,Len);%alpha 分布噪声
%         v_noise = stblrnd(1.2,0,2,0,q,Len);%接收alpha 分布噪声
    
%     v=raylrnd(30,q,Len);%瑞利分布噪声
%     v=v-mean(v')';

%     original_x=rand(q,Len);
%     v=tan((original_x-1/2)*pi);%柯西分布噪声    
%       v = randn(q,Len) * 20;% 高斯噪声
%       v_noise = randn(q,Len) * 20;% 高斯噪声
      
    [pai,R1,R2]=EM_GMM2(v_noise);  
    %     R1=v1*v1'/Len;
    %     R2=v2*v2'/Len;
    
    R=v*v'/Len;
%     R=v_noise*v_noise'/Len;  
    %     R1=R;
    %     R2=R;%高斯噪声
    %% 量测产生
    for i=1:Len
        zA(:,i) = H*xA(:,i) +v(:,i);
    end
    %% 初始状态给定
%       x0_hj=[700,10,700,-10]';
    x0_hj=[1000,10,1000,-10]';
    X_ikv(:,1)=x0_hj;
    X_kv(:,1)=x0_hj;
    X_MCKF_CV(:,1)=x0_hj;
    X_MEEKF_CV(:,1)=x0_hj;
    %% CV模型卡尔曼滤波
    tic %KF时间
    x0_kf_cv=x0_hj;
    P0_kf_cv=p0;%初始协方差
    for i=2:Len
        [P_kv,P_k_k_1v,X1_kv]=kalman(A1,G1,H,Q1,R,zA(:,i),x0_kf_cv,P0_kf_cv);  %%
        X_kv(:,i)=X1_kv;
        x0_kf_cv=X1_kv;
        P0_kf_cv=P_kv;
        Err_kf(i)=norm(X_kv(:,i)-xA(:,i));
    end
    toc
    tim_KF=toc;
    tim_KF_save(mn,:)=tim_KF;
    %% MCKFCV 单模型滤波
    tic
    x0_mckf_cv = x0_hj;
    P0_mckf_cv = p0;
    for i=2:Len
        [Pkk_MCKF3_CV,X_MCKF3_CV]=MCKF3(A1,G1,H,Q1,R,zA(:,i),x0_mckf_cv,P0_mckf_cv,p);
        X_MCKF_CV(:,i)=X_MCKF3_CV;
        P0_mckf_cv=Pkk_MCKF3_CV;
        x0_mckf_cv=X_MCKF3_CV;
        Err_mckf(i)=norm(X_MCKF_CV(:,i)-xA(:,i));
    end
    toc
    tim_MCKF=toc;
    tim_MCKF_save(mn,:)=tim_MCKF;
    %% MEEKF CV单模型滤波
    tic
    x0_meekf_cv = x0_hj;
    P0_meekf_cv = p0;
    for i=2:Len
        [Pkk_MEEKF3_CV,X_MEEKF3_CV]=MEEKF(A1,G1,H,Q1,R,zA(:,i),x0_meekf_cv,P0_meekf_cv,p);
        X_MEEKF_CV(:,i)=X_MEEKF3_CV;
        P0_meekf_cv=Pkk_MEEKF3_CV;
        x0_meekf_cv=X_MEEKF3_CV;
        Err_meekf(i)=norm(X_MEEKF_CV(:,i)-xA(:,i));
    end
    toc
    tim_MEEKF=toc;
    tim_MEEKF_save(mn,:)=tim_MEEKF;
    %% 方法3，利用后验概率加权求和，加迭代
    tic
%     Pi=[0.9,0.1;0.9,0.1];
    Pi=[pai(1),pai(2);pai(1),pai(2);];
    u1=pai(1);
    u2=pai(2);
    x0_ikf_cv=x0_hj;
    x0_ikf_cv1=x0_hj;
    x0_ikf_cv2=x0_hj;
    P0_ikf_cv1=p0;%初始协方差
    P0_ikf_cv2=p0;
    M4=[0;0];
    M2=0;
    for i=2:Len
        [u_a_1,u_a_2,P_a_1,P_a_2,X_a_1,X_a_2]=imkalman_v2(A1,G1,H,Q1,R1,R2,zA(:,i),x0_ikf_cv1,x0_ikf_cv2,P0_ikf_cv1,P0_ikf_cv2,u1,u2,Pi);
        X_ikv(:,i)=u_a_1*X_a_1+u_a_2*X_a_2;
        x0_ikf_cv1=X_a_1;
        x0_ikf_cv2=X_a_2;
        P0_ikf_cv1=P_a_1;
        P0_ikf_cv2=P_a_2;
        u1=u_a_1;
        u2=u_a_2;
        Err_ikf(i)=norm(X_ikv(:,i)-xA(:,i));
%         %% 峭度
%         eee(:,i)=zA(:,i)-H*X_ikv(:,i);
%         M4=0.9*M4+0.1*eee(:,i).^4;
%         M2=0.9*M2+0.1*eee(:,i).^2;
%         k4(:,i)=M4./M2.^2-3;
    end
    toc
    tim_BIKF=toc;
    tim_BIKF_save(mn,:)=tim_BIKF;
    %% 误差计算
    ex_kf_cv(mn,:)=X_kv(1,:)-xA(1,:);
    ey_kf_cv(mn,:)=X_kv(3,:)-xA(3,:);%kf
    ex_ikf_cv(mn,:)=X_ikv(1,:)-xA(1,:);
    ey_ikf_cv(mn,:)=X_ikv(3,:)-xA(3,:);%ikf
    ex_mckf_cv(mn,:)=X_MCKF_CV(1,:)-xA(1,:);
    ey_mckf_cv(mn,:)=X_MCKF_CV(3,:)-xA(3,:);%MCKF_CV
    ex_meekf_cv(mn,:)=X_MEEKF_CV(1,:)-xA(1,:);
    ey_meekf_cv(mn,:)=X_MEEKF_CV(3,:)-xA(3,:);%MEEKF_CV
    Err_ikf_mn(mn,:)=Err_ikf;
    Err_kf_mn(mn,:)=Err_kf;
    Err_mckf_mn(mn,:)=Err_mckf;
    Err_meekf_mn(mn,:)=Err_meekf;
end

% [pai_e,R_e1,R_e2]=EM_GMM2(eee);

mex_kf_cv=mean(ex_kf_cv,1);mey_kf_cv=mean(ey_kf_cv,1);%KF-CV
mex_ikf_cv=mean(ex_ikf_cv,1);mey_ikf_cv=mean(ey_ikf_cv,1);%IKF-CV
mex_mckf_cv=mean(ex_mckf_cv,1);mey_mckf_cv=mean(ey_mckf_cv,1);%MCKF_CV
mex_meekf_cv=mean(ex_meekf_cv,1);mey_meekf_cv=mean(ey_meekf_cv,1);%MEEKF_CV
me_Err_ikf_mn=mean(Err_ikf_mn,1);
me_Err_kf_mn=mean(Err_kf_mn,1);
me_Err_mckf_mn=mean(Err_mckf_mn,1);
me_Err_meekf_mn=mean(Err_meekf_mn,1);

me_tim_KF=mean(tim_KF_save);
me_tim_MCKF=mean(tim_MCKF_save);
me_tim_MEEKF=mean(tim_MEEKF_save);
me_tim_BIKF=mean(tim_BIKF_save);
%% 绘图
t=1:Len;
%% 绘制各算法跟踪效果
figure(1)
plot(xA(1,t),xA(3,t),'k' ,zA(1,t),zA(2,t),'m.',xA(1,t)+mex_kf_cv(t),xA(3,t)+mey_kf_cv(t), 'g--'...
    ,xA(1,t)+mex_mckf_cv(t),xA(3,t)+mey_mckf_cv(t),'b--',xA(1,t)+mex_meekf_cv(t),xA(3,t)+mey_meekf_cv(t),'y--'...
    ,xA(1,t)+mex_ikf_cv(t),xA(3,t)+mey_ikf_cv(t),'r--','linewidth',1);
legend('true', 'measurment','KF','MCKF','MEEKF','BIKF');
title('Tracking results')
xlabel('x/(m)'),ylabel('y/(m)');

%% 误差图
figure(2)
subplot(2,1,1)
plot(t,mex_kf_cv(t),'g.-',t,mex_ikf_cv(t),'r.-',t,mex_mckf_cv(t),'b.-',t,mex_meekf_cv(t),'y.-');
title('X方向滤波误差')
legend('CV模型常规Kalman滤波','CV模型混合高斯Kalman滤波','CV模型MCKF滤波','CV模型MEEKF滤波');
subplot(2,1,2)
plot(t,mey_kf_cv(t),'g.-',t,mey_ikf_cv(t),'r.-',t,mey_mckf_cv(t),'b.-',t,mey_meekf_cv(t),'y.-');
legend('CV模型常规Kalman滤波','CV模型混合高斯Kalman滤波','CV模型MCKF滤波','CV模型MEEKF滤波');
title('Y方向滤波误差')
xlabel('t(s)'),ylabel('位置误差(m)');

%%  计算RMSE-均方根误差
EX_kf_cv=sqrt(sum(ex_kf_cv.^2,1)/MC);
EY_kf_cv=sqrt(sum(ey_kf_cv.^2,1)/MC);
EX_ikf_cv=sqrt(sum(ex_ikf_cv.^2,1)/MC);
EY_ikf_cv=sqrt(sum(ey_ikf_cv.^2,1)/MC);
EX_mckf_cv=sqrt(sum(ex_mckf_cv.^2,1)/MC);
EY_mckf_cv=sqrt(sum(ey_mckf_cv.^2,1)/MC);
EX_meekf_cv=sqrt(sum(ex_meekf_cv.^2,1)/MC);
EY_meekf_cv=sqrt(sum(ey_meekf_cv.^2,1)/MC);
figure(3)
subplot(1,2,1)
plot(t,20*log10(EX_kf_cv(t)),'g.-',t,20*log10(EX_mckf_cv(t)),'b.-',t,20*log10(EX_meekf_cv(t)),'y.-',t,20*log10(EX_ikf_cv(t)),'r.-');
% plot(t,EX_kf_cv(t),'g.-',t,EX_mckf_cv(t),'b.-',t,EX_meekf_cv(t),'y.-',t,EX_ikf_cv(t),'r.-');
% plot(t,20*log10(EX_ikf_cv(t)),'r.-');
title('X-direction location RMSE')
xlabel('time/(s)'),ylabel('RMSE/(db)');
legend('KF','MCKF','MEEKF','BIKF');
% legend('BIKF without dealing with abnormal impulse');
mean_x_cv=mean(EX_kf_cv);
mean_x_icv=mean(EX_ikf_cv);
mean_x_mckfcv=mean(EX_mckf_cv);
mean_x_meekfcv=mean(EX_meekf_cv);
subplot(1,2,2)
plot(t,20*log10(EY_kf_cv(t)),'g.-',t,20*log10(EY_mckf_cv(t)),'b.-',t,20*log10(EY_meekf_cv(t)),'y.-',t,20*log10(EY_ikf_cv(t)),'r.-');
% plot(t,EY_kf_cv(t),'g.-',t,EY_mckf_cv(t),'b.-',t,EY_meekf_cv(t),'y.-',t,EY_ikf_cv(t),'r.-');
% plot(t,20*log10(EY_ikf_cv(t)),'r.-');
legend('KF','MCKF','MEEKF','BIKF');
% legend('BIKF without dealing with abnormal impulse');
title('Y-direction location RMSE')
xlabel('time/(s)'),ylabel('RMSE/(db)');
mean_y_cv=mean(EY_kf_cv);
mean_y_icv=mean(EY_ikf_cv);
mean_y_mckfcv=mean(EY_mckf_cv);
mean_y_meekfcv=mean(EY_meekf_cv);
%%
figure(4)
plot(t,Err_kf(t),'g',t,Err_ikf(t),'r',t,Err_mckf(t),'b',t,Err_meekf(t),'y');
legend('CV模型常规Kalman滤波','CV模型混合高斯Kalman滤波','CV模型MCKF滤波','CV模型MEEKF滤波');
%%
figure(5)
plot(t,me_Err_kf_mn(t),'g',t,me_Err_ikf_mn(t),'r',t,me_Err_mckf_mn(t),'b',t,me_Err_meekf_mn(t),'y');
legend('CV模型常规Kalman滤波','CV模型混合高斯Kalman滤波','CV模型MCKF滤波','CV模型MEEKF滤波');
% %% 
% figure(6)
% plot(t,k4(2,t),t,v(2,t));
% toc