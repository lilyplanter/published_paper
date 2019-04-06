function [B] = PLS( X,Y )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
r=size(X,2);%��¼X�ж�����
q=size(Y,2);%��¼Y�ж�����

%���ڼ���BPLS���ĸ���������
W=zeros(r,q);%��w1һ���жѵ�����
P=zeros(r,q);
Q=zeros(q,q);
b=zeros(q);
time=0;%ѭ������
yuzhi=norm(Y);
while(time<q*0.9&&norm(Y)<=yuzhi),%ѭ������������߲в�����Сʱ����ѭ��
time=time+1;  %ÿ��ѭ����һ
X0=X;
Y0=Y;
%[~,i] = max(var(Y0));%i��ʾ���������еı��
u1=Y0(:,1);
u1=u1/norm(u1);

stop=1;
while(stop>0.001),%ѭ��ֱ��u1���䣬��Ӧ�����е�step2
w1=X0'*u1;
w1=w1/norm(w1);
t1=X0*w1;
t1=t1/norm(t1);

c1=Y0'*t1;
c1=c1/norm(c1);
u1new=Y0*c1;
diff=u1-u1new;
stop=diff'*diff;
u1=u1new;
end;
%��Ӧ�����е�step3

W(:,time)=w1;
P(:,time)=X0'*t1;
C(:,time)=c1;
dx=t1*(P(:,time)');
X=X0-dx;

b(time)=u1'*t1/(t1'*t1);
Y=Y0-b(time)*t1*c1';

end;
B=W(:,1:time)/(P(:,1:time)'*W(:,1:time))*diag(b(1:time))*(C(:,1:time)');%����BPLS

end


