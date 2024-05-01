function [A, Q, C, R, pi, v]=train_kalman(trial, p)



%% training phase
A_sum1=0;
A_sum2=0;
%try identity for A
for j=1:p
    x=trial(p).spikes;
    z=trial(p).handPos;
    T= length(z);
    for i=2:T
        A1=z(:,i)*z(:,i-1)';
        A2=z(:,i-1)*z(:,i-1)';
        A_sum1=A_sum1+A1;
        A_sum2=A_sum2+A2;
    end
end
A=A_sum1*inv(A_sum2);

Q_sum=0;
for j=1:p
    x=trial(p).spikes;
    z=trial(p).handPos;
    T= length(z);
    for i=2:T
        Q1=(z(:,i)-A*z(:,i-1))*(z(:,i)-A*z(:,i-1))';
        Q_sum=Q_sum+Q1;
    end
end
Q=(1/(T-1))*Q_sum;


C1_sum=0;
C2_sum=0;
for j=1:p
    x=trial(p).spikes;
    z=trial(p).handPos;
    T= length(z);
    for i=1:T
        C1=x(:,i)*z(:,i)';
        C2=z(:,i)*z(:,i)';
        C1_sum=C1_sum+C1;
        C2_sum=C2_sum+C2;
    end
end
C=C1_sum*inv(C2_sum);

start_point=zeros(3,p);
start_var=zeros(3,3,p);
R_sum=0;
for j=1:p
    x=trial(p).spikes;
    z=trial(p).handPos;
    T= length(z);
    for i=1:T
        R=(x(:,i)-C*z(:,i))*(x(:,i)-C*z(:,i))';
        R_sum=R_sum+R;
    end
    start_point(:,j)=z(:,1);
    start_var(:,:,j)=cov(z');
end
R=(1/T)*R_sum;

pi=mean(start_point,2);
v=mean(start_var,3);


end