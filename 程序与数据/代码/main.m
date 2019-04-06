
%读取原始数据
 miRNAdisease = xlsread('1.miRNA-disease关联数据\miRNA-disease关联数字编号.xlsx');
 mirnaname=xlsread('1.miRNA-disease关联数据\miRNA数字编号.xlsx');
 diseasename=xlsread('1.miRNA-disease关联数据\disease数字编号.xlsx');
 diseaseweight = textread('2.disease semantic similarity 1\疾病语义类似性加权矩阵1.txt');
 diseasesimilarity1 = textread('2.disease semantic similarity 1\疾病语义类似性矩阵1.txt');
 diseasesimilarity2 = textread('3.disease semantic similarity 2\疾病语义类似性矩阵2.txt');
 mirnaweight = textread('4.miNA功能类似性矩阵\miRNA功能类似性加权矩阵.txt');
 mirnasimilarity = textread('4.miNA功能类似性矩阵\miRNA功能类似性矩阵.txt');
 r=180; 
 p=220;
 q=170;
%构建MiRNA-疾病关系矩阵A
mirnanum=size(mirnaname,1);
diseasenum=size(diseasename,1);
associationnum=size(miRNAdisease,1);
A=zeros(mirnanum,diseasenum);

for i=1:1:5430
A(miRNAdisease(i,1),miRNAdisease(i,2))=1;
end;


localresult=zeros(5430,1);
globalresult=zeros(5430,1);
% 进行交叉验证(5430次)
for i=1:1:5430
    %修改矩阵A（把用于验证的数据置零）
    mirna=miRNAdisease(i,1);
    disease=miRNAdisease(i,2);
    Al=A;
    Al(mirna,disease)=0;
    %构造MiRNA和疾病的相似度矩阵
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
    
    %打分
    [U,S,V] = svd (Al) ;
    S2=sqrtm(S(1:r,1:r));%对S进行开方
    Ar=U(:,1:r)*S2;%mxr
    Ad=V(:,1:r)*S2;%nxr
    
    [U1,S1,~] = svd (Sm) ;
    S2=sqrtm(S1(1:p,1:p));%对S进行开方
    Fr=U1(:,1:p)*S2;%mxp
    
    [U3,S3,~] = svd (Sd) ;
    S2=sqrtm(S3(1:q,1:q));%对S进行开方
    Fd=U3(:,1:q)*S2;%nxq
    
    Br=PLS1(Fr,Ar);
    Bd=PLS1(Fd,Ad);
    final=Fr*Br*Bd'*Fd';
 %得到local和global排序
    compare=final-Al.*10;
    localresult(i,1)=sum(compare(:,disease)>compare(mirna,disease))+sum(compare(:,disease)==compare(mirna,disease))/2+0.5;
    globalresult(i,1)=sum(sum(compare>compare(mirna,disease)))+sum(sum(compare==compare(mirna,disease)))/2+0.5;%global排序
end;
xlswrite('newlocal.xlsx',localresult');
xlswrite('newglobal.xlsx',globalresult');



