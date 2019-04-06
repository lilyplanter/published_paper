
%��ȡԭʼ����
 miRNAdisease = xlsread('1.miRNA-disease��������\miRNA-disease�������ֱ��.xlsx');
 mirnaname=xlsread('1.miRNA-disease��������\miRNA���ֱ��.xlsx');
 diseasename=xlsread('1.miRNA-disease��������\disease���ֱ��.xlsx');
 diseaseweight = textread('2.disease semantic similarity 1\�������������Լ�Ȩ����1.txt');
 diseasesimilarity1 = textread('2.disease semantic similarity 1\�������������Ծ���1.txt');
 diseasesimilarity2 = textread('3.disease semantic similarity 2\�������������Ծ���2.txt');
 mirnaweight = textread('4.miNA���������Ծ���\miRNA���������Լ�Ȩ����.txt');
 mirnasimilarity = textread('4.miNA���������Ծ���\miRNA���������Ծ���.txt');
 r=180; 
 p=220;
 q=170;
%����MiRNA-������ϵ����A
mirnanum=size(mirnaname,1);
diseasenum=size(diseasename,1);
associationnum=size(miRNAdisease,1);
A=zeros(mirnanum,diseasenum);

for i=1:1:5430
A(miRNAdisease(i,1),miRNAdisease(i,2))=1;
end;


localresult=zeros(5430,1);
globalresult=zeros(5430,1);
% ���н�����֤(5430��)
for i=1:1:5430
    %�޸ľ���A����������֤���������㣩
    mirna=miRNAdisease(i,1);
    disease=miRNAdisease(i,2);
    Al=A;
    Al(mirna,disease)=0;
    %����MiRNA�ͼ��������ƶȾ���
    A2=Al.^2;
    rd=1/(sum(sum(A2))/diseasenum);
    rm=1/(sum(sum(A2,2))/mirnanum);
    KD=zeros(diseasenum);
    KM=zeros(mirnanum);
    for j=1:1:diseasenum
    KD(j,:)=exp(-rd.*sum((repmat(Al(:,j),1,diseasenum)-Al).^2));
    end;
    for j=1:1:mirnanum
    KM(j,:)=exp(-rm.*sum((repmat(Al(j,:),mirnanum,1)-Al).^2,2));
    end;
    Sd=(diseasesimilarity1+diseasesimilarity2)./2.*diseaseweight+KD.*(diseaseweight-1).*(-1);
    Sm=mirnasimilarity.*mirnaweight+KM.*(mirnaweight-1).*(-1);
    
    %���
    [U,S,V] = svd (Al) ;
    S2=sqrtm(S(1:r,1:r));%��S���п���
    Ar=U(:,1:r)*S2;%mxr
    Ad=V(:,1:r)*S2;%nxr
    
    [U1,S1,~] = svd (Sm) ;
    S2=sqrtm(S1(1:p,1:p));%��S���п���
    Fr=U1(:,1:p)*S2;%mxp
    
    [U3,S3,~] = svd (Sd) ;
    S2=sqrtm(S3(1:q,1:q));%��S���п���
    Fd=U3(:,1:q)*S2;%nxq
    
    Br=PLS1(Fr,Ar);
    Bd=PLS1(Fd,Ad);
    final=Fr*Br*Bd'*Fd';
 %�õ�local��global����
    compare=final-Al.*10;
    localresult(i,1)=sum(compare(:,disease)>compare(mirna,disease))+sum(compare(:,disease)==compare(mirna,disease))/2+0.5;
    globalresult(i,1)=sum(sum(compare>compare(mirna,disease)))+sum(sum(compare==compare(mirna,disease)))/2+0.5;%global����
end;
xlswrite('newlocal.xlsx',localresult');
xlswrite('newglobal.xlsx',globalresult');



