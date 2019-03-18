% % % % ---------------ȫ��ѭ��--------------
clc;
clear all;


rR=1:1:30;%rmax
rr=0.05:0.05:0.5;%rmin
% % % % % % % % %--------------˫�������̬�ֲ�,�����--------------------------
rmin=0.01;
rmax=10;


rm1=0.18;lnof1=0.46;vf1=1.5;rm2=3.2;lnof2=0.66;vf2=1;%˫����й�ҵ7/1.45-0.010  ��ָ��
%  rm1=0.14;lnof1=0.44;vf1=1.8;rm2=3.4;lnof2=0.74;vf2=1;%˫������ȼ��2/1.50-0.010  ��ָ��
% rm1=0.15;lnof1=0.46;vf1=1.4;rm2=3.6;lnof2=0.73;vf2=1;%˫������ȼ��3/1.50-0.010  ��ָ��
% rm1=0.15;lnof1=0.48;vf1=1.6;rm2=3.5;lnof2=0.75;vf2=1;%˫������ȼ��4/1.50-0.010  ��ָ��
% rm1=0.16;lnof1=0.45;vf1=1.8;rm2=3.6;lnof2=0.76;vf2=1;%˫������ȼ��5/1.50-0.010  ��ָ��
% rm1=0.16;lnof1=0.46;vf1=2;rm2=3.7;lnof2=0.78;vf2=1;%˫������ȼ��6/1.50-0.010  ��ָ��




f_r_LL=rmin:0.01:rmax;
n1_1=-((log(f_r_LL)-log(rm1)).^2)/(2*lnof1^2);
n2_1=exp(n1_1); 
n3_1=lnof1*sqrt(2*pi);
n1=(vf1*n2_1)/n3_1;
n1_2=-((log(f_r_LL)-log(rm2)).^2)/(2*lnof2^2);
n2_2=exp(n1_2); 
n3_2=lnof2*sqrt(2*pi);
n2=(vf2*n2_2)/n3_2;
nn1=n1+n2;
nn=nn1./f_r_LL;


% f_r_LL=rmin:0.01:rmax;
% a=89.972;b=0.019443;c=0.33154;d=0.84;%����1%0.5-25
% n=a.*((2.*f_r_LL).^b).*exp(-c.*(2.*f_r_LL).^d);%����1

% a=17797;b=5.21;c=8.21;d=0.4;%����2 0.5-33
% n=a.*(f_r_LL.^b).*exp(-c.*(f_r_LL.^d));%����2
% 
% nn=4/3*pi.*(f_r_LL.^3).*n;%�����dvdr
% nn1=nn;

figure(1);
semilogx(f_r_LL,nn1,'b');%
% %  loglog(f_r,nn,'-');
% %  plot(r,nn,'-');
hold on;

% FQext_i_Z=xlsread('1.42-0.004.xlsx','A2:K1001');

% FQext_i_Z=load('O:\����52-7+5(1.50-0.003i).txt');
FQext_i_Z=xlsread('C:\Users\Julian\Desktop\�����ͼ����\�����ͼ����\1.45-0.010i.xlsx');%��
% FQext_i=FQext_i_Z(1:1000,[1,3,5,8,9,10]);%340\400\870\1020\1627
FQext_i=FQext_i_Z(1:3300,1:end);%340\500\1020\1627\2200

FQext1=FQext_i'; %%%������Ϊ��������Ϊ�뾶
rmin1=round(rmin*100);
rmax1=round(rmax*100);
% FQext_LL=FQext1([3,4,5,6,7,8,9],rmin1:rmax1);
FQext_LL=FQext1(:,rmin1:rmax1);%12*1000
[line_LL,row_LL]=size(FQext_LL);

% % % % % % % % % %----------------------��ʼ��kp�������---------------------
KP_LL=zeros(line_LL,row_LL);
r_2=3/4;
for k=1:12
 KP_LL(k,1:end)=r_2.*FQext_LL(k,1:end)./f_r_LL; 
end
FK1=0.01.*KP_LL*nn';%%%%%%%%%%%%%%%%%����*0.001��%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% e_x=[1.05,1.05,1.05,1.05,1.05,1.05,1.05];
e_x=[1,1,1,1,0.95,1,1,1,1,1,1,1];
FK=FK1.*e_x';%����ϵ��ok
FQext_LL=FQext1;

% % % % % % % % %----------------�������ʲ�ѭ��--------------------------
data=zeros(3300,12,630);
FQext_XH1=zeros(12,3300,630);
pathname=uigetdir('C:\Users\Julian\Desktop\�����ͼ����\�����ͼ����\*.xlsx','��ѡ���ļ���');
if pathname==0
    msgbox=('û����ȷѡ���ļ���');
    return; 
end;
filexlsx=ls(strcat(pathname,'/*.xlsx'));
files=cellstr(filexlsx);
len=length(files);
for i=1:len
filesname(i)=strcat(pathname,'/',files(i));
data(:,:,i)=xlsread(cell2mat(filesname(i)),'A1:L3300');
FQext_XH1_10(:,:,i)=(data(:,:,i))';
end
FQext_XH1=FQext_XH1_10(:,:,:);
% FQext_XH1=FQext_XH1_10;

i_i=1; 
error=[];
r_max_c=[];
r_min_c=[0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5];
for i=1:630
% i=4;

for r_min=r_min_c
    for r_max=1:1:30
      r_max_c=[r_max_c,r_max];  
% r_min=0.05;r_max=5;
f_r=r_min:0.01:r_max;       
r_min1=round(r_min*100);
r_max1=round(r_max*100);
FQext=FQext_XH1(:,r_min1:r_max1,i);
% FQext=FQext_LL(:,r_min1:r_max1);
[line,row]=size(FQext);%12*r

% % % % % % % % % % % ----------------------kp�������---------------------ok
KP=zeros(line,row);
r_2=3/4;
for k=1:line
 KP(k,1:end)=r_2.*FQext(k,1:end)./f_r; 
end

% % % % %  ----------------------B��������(��������)-----------------------ok
r_linelog=log10(r_min);
r_rowlog=log10(r_max);
r=logspace(r_linelog,r_rowlog,14);
B=zeros(line,row);
a1=[0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4];
for j=1:1:12
%     a1=0.41;
%     B(j,1:end)=(1./(a1.*sqrt(2.*pi))).*exp(-(log(f_r)-log(r(j))).^2./(2.*a1.^2));
B(j,1:end)=(1./(a1(j).*sqrt(2.*pi))).*exp(-(log(f_r)-log(r(j+1))).^2./(2.*a1(j).^2));
end

% semilogx(f_r,B(1,:));
% hold on 
% semilogx(f_r,B(2,:));
% hold on 
% semilogx(f_r,B(3,:));
% hold on 
% semilogx(f_r,B(4,:));
% hold on 
% semilogx(f_r,B(5,:));
% hold on 
% semilogx(f_r,B(6,:));
% hold on 
% semilogx(f_r,B(7,:));
% hold on 
% semilogx(f_r,B(8,:));
% hold on 
% semilogx(f_r,B(9,:));
% hold on 
% semilogx(f_r,B(10,:));
% hold on 
% semilogx(f_r,B(11,:));
% hold on 
% semilogx(f_r,B(12,:));%��������ok

A_pn=0.01.*KP*B';%ok
% A_pn=KP*B';

% % % % % % % % % % % % ------------------�Ľ����������ճ���---------------------------
H=[1 -2 1 0 0 0 0 0 0 0 0 0;
   -2 5 -4 1 0 0 0 0 0 0 0 0;
   1 -4 6 -4 1 0 0 0 0 0 0 0;
   0 1 -4 6 -4 1 0 0 0 0 0 0;
   0 0 1 -4 6 -4 1 0 0 0 0 0
   0 0 0 1 -4 6 -4 1 0 0 0 0
   0 0 0 0 1 -4 6 -4 1 0 0 0
   0 0 0 0 0 1 -4 6 -4 1 0 0
   0 0 0 0 0 0 1 -4 6 -4 1 0
   0 0 0 0 0 0 0 1 -4 6 -4 1
   0 0 0 0 0 0 0 0 1 -4 5 -2
   0 0 0 0 0 0 0 0 0 1 -2 1];
A_A=A_pn'*A_pn;

error_kb=[];
for b=20:1:28 
for K=1:1:25
    larg=(2^K)*10^(-b);
    C1=pinv(A_A+larg*H);
    C=C1*A_pn'*FK;
% % % HH=C'*B; 
% % % % HH1=HH.*f_r;
v = 0;
for i_U=1:12;
	if C(i_U)<0
		 C(i_U)=0;
	end     
	v=v+C(i_U).*B(i_U,:);
end;
HH=v; 
ZK_abs=0.01.*KP*HH';
error1_kb=((sum(abs(FK-ZK_abs)))/12)*100;
error_kb=[error_kb,error1_kb];
end
end
error_fz=reshape(error_kb,25,9); 
[er_line,er_row]=find(error_fz==min(error_fz(:)));
K1=er_line(1);
b1=er_row(1)+19;
larg_min=(2^K1)*10^(-b1);
C1_min=pinv(A_A+larg_min*H);
C_min=C1_min*A_pn'*FK;
% % % % HH_min=C_min'*B;
v_min = 0;
for i_P=1:12;
	if C_min(i_P)<0
		 C_min(i_P)=0;
	end     
	v_min=v_min+C_min(i_P).*B(i_P,:);
end;
HH_min=v_min;
ZK=0.01.*KP*HH_min';
HH_quanbu(r_min1:r_max1,i_i)=HH_min';
[l_min,l]=min(error_kb);
error=[error,l_min];
[dis_min,dis_l]=min(error);%%%%%�ҵ��о�����С�ģ��оݵ�ֵ�Լ���Ӧ�ı�ţ�������֤��ΪʲôҪƽ����
i_i=i_i+1; 
     end  
end 
end

% % % refrac_real=[];
% % % refrac_imagin=[];
rmaxrange_qj=[];
rminrange_qj=[];
% [a,b]=min(error);

% [error_paixu,error_I]=sort(error);
[a,b]=min(error);
% cd=length(error_paixu);
% N=round(cd*0.01);%ȡ�ٷֱȵ�����
% y_liu=[];
% for i_value=1:N 
%         i_p=error_I(i_value);
%          y_liu=[ y_liu,i_p];
% H_cun(:,i_value)=HH_quanbu(:,i_p); 
% end
H_zuizhong1=HH_quanbu(:,b);
% H_zuizhong1=HH_quanbu(:,b); 
% % % % % % % % % % % % ---------------�������ʺͰ뾶����------------------------- 
% xx=ceil(b/30);%rmin
% x=rr(xx);
% y=rem(b,30);
% if y==0
%     y=30;%rmax
% end
xuhaok=b;
refrac_index1=ceil(xuhaok/300);%%%%�������������С����뾶�����ĳ˻�   246
% xhk=length(y_liu);
% for i_k=1:xhk  
file_zuizhong=files(refrac_index1);
% end;
% % % % % -------------------------������-----------------
refrac=(file_zuizhong)
% refrac_real12(:,1)=str2num(refrac(:,1:4));
% refrac_real1=roundn(sum(refrac_real12)/xhk,-2);
% % % % refrac_real=[refrac_real,refrac_real1];
% refrac_imagin12(:,1)=str2num(refrac(:,6:10));
% refrac_imagin1=roundn(sum(refrac_imagin12)/xhk,-3);
% refrac_imagin=[refrac_imagin,refrac_imagin1]; 
% % % % ----------�뾶����----------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xuhaok1=xuhaok-floor(xuhaok/100)*100;%%%ע��

%rmin
xx=rem(xuhaok,300);
if xx==0
    rxiao=rr(end);%rmin
else
    x=ceil(xx/30);  
    rxiao=rr(x);%rmin
end
%rmax
yy=rem(xuhaok,30);
if yy==0
    rda=rR(end);%rmax
else
    rda=rR(yy);%rmax
end


% rminrange_qj1=max(rminrange_qj12);
% rmaxrange_qj1=min(rmaxrange_qj12);
% rmaxrange_qj=[rmaxrange_qj,rmaxrange_qj1];
% rminrange_qj=[rminrange_qj,rminrange_qj1];
% % % % % % % % % % % % % % ------------------------------��ͼ---------------------
% r_zuixiao1=r_min_c(rminrange_qj1);%%%%������ƽ������С�뾶��λ��
% r_zuixiao=round(100*r_zuixiao1);
% r_zuida1=r_max_c(rmaxrange_qj1);  %%%%%������ƽ�������뾶��λ��
% r_zuida=round(r_zuida1*100);
% f_rzuizhong=rxiao:0.01:rda;
if rxiao <= rmin
    rxiao=rmin;
end
f_rzuizhong=rxiao:0.01:rda;
% H_zuizhong=H_zuizhong1(rxiao:rda);
H_zuizhong=H_zuizhong1(rxiao*100:rda*100);
% ggg=H_zuizhong*KP'
H_zuizhonglnr=H_zuizhong.*f_rzuizhong';
% H_zuizhonglnr=H_zuizhong;
semilogx(f_rzuizhong,H_zuizhonglnr,'-r');
hold on;
xlabel('Radius(um)');
ylabel('dv(r)/dlnr(um^3cm^-3)');
legend('Given','Retrived');
nn1=nn1';
figure(2)
errorvv=(H_zuizhonglnr-nn1(rxiao*100:rda*100))./nn1(rxiao*100:rda*100)*100;
semilogx(f_rzuizhong,errorvv);

errorvpj=sum(abs(errorvv))./length(f_rzuizhong);

nnqujian=nn(rxiao*100:rda*100);
fenzi1=sum(H_zuizhong.*0.01);%��õ��Ƿ��ݳ�������Ч�뾶��
fenmu1=sum(H_zuizhong./f_rzuizhong'.*0.01);
reff1=fenzi1/fenmu1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%������ǳ�ʼ״̬�£��ٶ�����ֵ������Ч�뾶%%%%%
fenzi=sum(nnqujian.*0.01);
fenmu=sum(nnqujian./f_rzuizhong.*0.01);
reff=fenzi/fenmu;
errorr=((reff1-reff)/reff)*100;%�������Ч�뾶��������
%%%%%%%%%%%%%%%%%%%�����Ũ�����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s=3*sum(H_zuizhong./f_rzuizhong'.*0.01);%����
s1=3*sum(nnqujian./f_rzuizhong.*0.01);%��ʼ
errors=((s-s1)/s1)*100;
%%%%%%%%%%%%%%%%%%%��Ũ�����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n=VVV.*(3./(4*pi.*power(r,3)));%������Ũ����
% N=sum(n.*0.01);%������Ũ��
% n1=v0'.*(3./(4*pi.*power(r,3)));%��ʼ
% N1=sum(n1*0.01);%��ʼ��Ũ��
% errorn=((N-N1)/N1)*100;
n_confz=sum(0.01*((3*H_zuizhong)./(4*pi*f_rzuizhong'.^3)));
n_conys=sum(0.01*((3*nnqujian)./(4*pi*f_rzuizhong.^3)));
errorn=((n_confz-n_conys)/n_conys)*100;
%%%%%%%%%%%%%%%%%%%%%%%%���Ũ��%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v1=sum(H_zuizhong.*0.01);%����
v2=sum(nnqujian.*0.01);%��ʼ
errorv=((v1-v2)/v2)*100;
% result=[delt1,delt2,error,errorv,errors,pingjunjueduiwucha,pingjunxiangduiwucha];
result=[errorvpj,errorr,errorv,errors,errorn]
