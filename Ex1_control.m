clear; clc; clf;
t=[0];
dt=1e-4;
tau=0.5;
T=15;
h=0.0001;
t_total=[0:h:T];
x_total=zeros(2,length(t_total));
y_total=zeros(1,length(t_total));
y_total(1,1)=1.21;
u=zeros(1,length(t_total));
u(1,1)=-3.0755;

a=0.1;
 A=[-0.7,0;0,-0.8];B=[1.8,0;0,2];C=[0.1;0.1];D=[1.2,0.1;0.1,1.2];
         P=[5.35,-0.0018;-0.0018,5.35];
         H=[-1.864,-1.864];

a1=0.5;a2=1;b1=0.85;b2=2;b3=0.85;b4=0.3;
x_total(1,1)=1.05;
x_total(2,1)=0.6;
d=3.5;
m=5;
for k=1:1:7
     for i=2:length(t_total)
         t_now=(i-1)*h;
         M=round((t(k)+d)/h);
       if t_now>t(k)+d && [x_total(1,i-1),x_total(2,i-1)]*P*[x_total(1,i-1),x_total(2,i-1)]'>=exp(0.1)*[x_total(1,M-1),x_total(2,M-1)]*P*[x_total(1,M-1),x_total(2,M-1)]'
           t(k+1)=t_now;
       elseif t_now<t(k)+d
           t(k+1)=t(k)+m;
       end
       if t(k+1)>=t(k)+m
           t(k+1)=t(k)+m;
       end
     if t(k)<=t_now&&t_now<t(k)+d
         taut=t_now-tau;
           if taut<=0
               M1=0;
           elseif taut>0
           M1=round(taut/h);
           end
         x_total(1,i)=x_total(1,i-1)+((A(1,1)+a*sin(0.5*t_now))*x_total(1,i-1)+A(1,2)*x_total(2,i-1)+(B(1,1)+a*sin(0.5*t_now))*(b4*(abs(b1*x_total(1,i-1)+b2)-abs(b3*x_total(1,i-1)-b2)))+B(1,2)*(b4*(abs(b1*x_total(2,M1+1)+b2)-abs(b3*x_total(2,M1+1)-b2)))+(C(1)+a*sin(0.5*t_now))*H(1)*x_total(1,i-1)+(C(1)+a*sin(0.5*t_now))*H(2)*x_total(2,i-1))*h;
         x_total(2,i)=x_total(2,i-1)+(A(2,1)*x_total(1,i-1)+(A(2,2)+a*sin(0.5*t_now))*x_total(2,i-1)+B(2,1)*(b4*(abs(b1*x_total(1,i-1)+b2)-abs(b3*x_total(1,i-1)-b2)))+(B(2,2)+a*sin(0.5*t_now))*(b4*(abs(b1*x_total(2,M1+1)+b2)-abs(b3*x_total(2,M1+1)-b2)))+C(2)*H(1)*x_total(1,i-1)+(C(2))*H(2)*x_total(2,i-1))*h;
         y_total(1,i)=sqrt(x_total(1,i)^2+x_total(2,i)^2);
         u(1,i)=H*[x_total(1,i);x_total(2,i)];
             end
  if abs(t_now-t(k))<=1e-6
  x_total(1,i)=(D(1,1)+a*sin(0.5*t_now))*x_total(1,i-1)+D(1,2)*x_total(2,i-1);
  x_total(2,i)=(D(2,1))*x_total(1,i-1)+(D(2,2)+a*sin(0.5*t_now))*x_total(2,i-1);
  y_total(1,i)=sqrt(x_total(1,i)^2+x_total(2,i)^2);
  u(1,i)=H*[x_total(1,i);x_total(2,i)];
   end
       
       if t(k)+d<=t_now&&t_now<t(k+1)
           taut=t_now-tau;
           if taut<=0
               M1=0;
           elseif taut>0
           M1=round(taut/h);
           end
           x_total(1,i)=x_total(1,i-1)+((A(1,1)+a*sin(0.5*t_now))*x_total(1,i-1)+A(1,2)*x_total(2,i-1)+(B(1,1)+a*sin(0.5*t_now))*(b4*(abs(b1*x_total(1,i-1)+b2)-abs(b3*x_total(1,i-1)-b2)))+B(1,2)*(b4*(abs(b1*x_total(2,M1+1)+b2)-abs(b3*x_total(2,M1+1)-b2))))*h;
         x_total(2,i)=x_total(2,i-1)+(A(2,1)*x_total(1,i-1)+(A(2,2)+a*sin(0.5*t_now))*x_total(2,i-1)+B(2,1)*(b4*(abs(b1*x_total(1,i-1)+b2)-abs(b3*x_total(1,i-1)-b2)))+(B(2,2)+a*sin(0.5*t_now))*(b4*(abs(b1*x_total(2,M1+1)+b2)-abs(b3*x_total(2,M1+1)-b2))))*h;
         y_total(1,i)=sqrt(x_total(1,i)^2+x_total(2,i)^2);
   
        end
  if abs(t_now-t(k))<=1e-6
  x_total(1,i)=(D(1,1)+a*sin(0.5*t_now))*x_total(1,i-1)+D(1,2)*x_total(2,i-1);
  x_total(2,i)=(D(2,1))*x_total(1,i-1)+(D(2,2)+a*sin(0.5*t_now))*x_total(2,i-1);
  y_total(1,i)=sqrt(x_total(1,i)^2+x_total(2,i)^2);
  end  
     y1(1,i)=0.45*t_now/t_now;
     
end
end

for i=1:6
    t1(i)=t(i+1);
y(1,i)=0.3*(t1(i)/t1(i));
end

   m=0:0.01:15; l=0:0.01:5;
   h=1.15*m;   j=7.5.*l./l;
   n=h./m;     
a=0:0.01:2;
b=9.5.*a./a;
  plot(t_total,y_total(1,:),'g',t1,y(1,i),'x k',m,n,'--',j,l,'--',b,a,'--');%,t_total,y1(1,i),'--'
legend('||\iota(t)||','t_k');
xlabel('t'),ylabel('||\iota(t)||')
title('(b)')
 xlim([0,12]);ylim([0.2,1.4]);
 annotation('arrow',[0.28 0.28],[0.85,0.76])
 text(1,1.3,'\rho=1.15')
 
 annotation('arrow',[0.5 0.61],[0.25,0.25])
 str='$$\mathcal{T}-\pi=7.5$$';
 text(5.4,0.46,str,'interpreter','latex')
 
 annotation('arrow',[0.85 0.75],[0.3,0.3])
 str='$$\mathcal{T}=9.5$$';
 text(10,0.53,str,'interpreter','latex')
 
 figure
 plot(t_total,x_total(1,:),'r',t_total,x_total(2,:),'b',t1,y(1,i),'x k',b,a,'--');%,m,n,'--',j,l,'--',t_total,y1(1,i),'--'
 xlim([0,10]);ylim([0,1.2]);
 legend('\iota_1(t)','\iota_2(t)','t_k');
xlabel('t'),ylabel('\iota_1(t),\iota_2(t)')
title('(a)')
annotation('arrow',[0.75 0.85],[0.6,0.6])
 str='$$\mathcal{T}=9.5$$';
 text(7.9,0.75,str,'interpreter','latex')
 figure
  plot(t_total,u(1,:),'m',t,y(1,i),'x');
  xlim([0,10]);ylim([-4,0.1]);
  legend('u(t)');
xlabel('t'),ylabel('u(t)')
title('(c)')
 
  
  