clear
%�������
L=10E-12;
PHI_0=2.067833757E-15;
C=50E-15;
I0=100E-6;
beta_L=2*I0*L/PHI_0;
R=7.256511161;
beta_C=(2*pi/PHI_0)*I0*R*R*C;

%��Ե���ƫ��
ib=2;
%��-0.5��0.5֮��ֳ�N��������ܹ�N+1��phi_bֵ��phi_bΪ�Ա�����ÿһ��phi_b��ֵ���ǽ���һ�η�����
%�ܹ���΢�ַ���N��
N=5000;
phi_b=zeros(N+1,1);
phi_b_max=5;
for ite=1:N+1
    phi_b(ite)=-0.5+(ite-1)*phi_b_max/N;
end

%����ʱ�䲽���ͽ�洢�ռ�
t=5;
n=5000;
h=t/n;
delta1d=zeros(N+1,n+1);
delta2d=zeros(N+1,n+1);
delta1=zeros(N+1,n+1);
delta2=zeros(N+1,n+1);
for ite0=1:N+1 %ÿѭ��һ�ξͼ���һ�η�����
    for ite=2:n+1
        k1a=f1(delta1(ite0,ite-1),delta2(ite0,ite-1),delta1d(ite0,ite-1),beta_L,beta_C,phi_b(ite0));
        k2a=f1(delta1(ite0,ite-1),delta2(ite0,ite-1),delta1d(ite0,ite-1)+h/2*k1a,beta_L,beta_C,phi_b(ite0));
        k3a=f1(delta1(ite0,ite-1),delta2(ite0,ite-1),delta1d(ite0,ite-1)+h/2*k2a,beta_L,beta_C,phi_b(ite0));
        k4a=f1(delta1(ite0,ite-1),delta2(ite0,ite-1),delta1d(ite0,ite-1)+h*k3a,beta_L,beta_C,phi_b(ite0));
        
        k1b=f2(delta1(ite0,ite-1),delta2(ite0,ite-1),delta2d(ite0,ite-1),beta_L,beta_C,phi_b(ite0),ib);
        k2b=f2(delta1(ite0,ite-1),delta2(ite0,ite-1),delta2d(ite0,ite-1)+h/2*k1b,beta_L,beta_C,phi_b(ite0),ib);
        k3b=f2(delta1(ite0,ite-1),delta2(ite0,ite-1),delta2d(ite0,ite-1)+h/2*k2b,beta_L,beta_C,phi_b(ite0),ib);
        k4b=f2(delta1(ite0,ite-1),delta2(ite0,ite-1),delta2d(ite0,ite-1)+h*k3b,beta_L,beta_C,phi_b(ite0),ib);
        
        delta1d(ite0,ite)=delta1d(ite0,ite-1)+h*(1/6)*(k1a+2*k2a+2*k3a+k4a);
        delta2d(ite0,ite)=delta2d(ite0,ite-1)+h*(1/6)*(k1b+2*k2b+2*k3b+k4b);
        delta1(ite0,ite)=delta1(ite0,ite-1)+h/2*(delta1d(ite0,ite)+delta1d(ite0,ite-1));
        delta2(ite0,ite)=delta2(ite0,ite-1)+h/2*(delta2d(ite0,ite)+delta2d(ite0,ite-1));
    end
end
plot(phi_b,delta2d(:,5001))
