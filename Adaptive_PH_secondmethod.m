close all;
clear all;
clc;format long
tic
%% material information
l=0.0125; k=1e-19;
gc=2.7e-3; %KN/mm（能量释放率）
lambda=121.15; mu=80.7692; %GPa（弹性体模量和剪切模量）
%Emod=mu*(2*mu+3*lambda)/(mu+lambda);nu=lambda/(2*(mu+lambda));%拉梅常数与弹性模量泊松比之间的转换关系
%lambda=(nu*Emod)/[(1+nu)*(1-2nu)];mu=Eomd/[2*(1+nu)];
%elasticity matrix (plane stress)
%Cmat=Emod/(1-nu^2)*[1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];
%elasticity matrix (plane strain)
stressState = 'PLANE_STRAIN';
Cmat=[lambda+2*mu,lambda,0;lambda lambda+2*mu 0;0 0 mu];
disp([num2str(toc),'  MATERIAL information'])
%% NURBS geometry and meshing
L = 1;         % Length
H = 1;         % Hight

numPatches = 4;

target_rel_error = 1e-3;
targetScale = 0.5;
fai=0.25;
% degree
p = 3;
q = 3; 

dimBasis = zeros(1, numPatches);

GIFTmesh = init2DGeometryGIFTMP(numPatches,L,H);
%[ GIFTmesh] = init2DGeometryGIFTMPrectangle(L,H,numPatches);

quadList = cell(numPatches,1);
PHTelem = cell(numPatches, 1);
for i=1:numPatches
    [PHTelem{i}, dimBasis(i)] = initPHTmesh(p,q);
    quadList{i} = 2:5;    
end

patchBoundaries = {1, 2, 2, 4; 2, 3, 3, 1; 3, 4, 4, 2};

tic %启动秒表计时器\
keep_refining = 1;
num_steps = 0;
while keep_refining
num_steps = num_steps + 1;
toc
disp(['Step ', num2str(num_steps)])
figure
plotPHTMeshMP( PHTelem, GIFTmesh )
[ PHTelem, dimBasis, quadList ] = checkConforming( PHTelem, dimBasis, patchBoundaries, p, q, quadList );
[ PHTelem, sizeBasis ] = zipConforming( PHTelem, dimBasis, patchBoundaries, p, q);
    sizeBasis
    
    %assemble the linear system
%disp('Assembling the linear system...')
%[ stiff ] = assembleGalerkinSysGIFTMP( PHTelem, GIFTmesh, sizeBasis, p, q, Cmat);
[iK,jK,iKd,jKd,iFd,jFd,dedofMat,edofMat] = edofMatconnect(PHTelem);
[Nd,Nu,Ndp,Bd,Bu,Bdp] = shape_function(PHTelem,GIFTmesh,p,q,edofMat,dedofMat);
nel=size(edofMat,1);
[ detJacobi ] = fejacobi( PHTelem, GIFTmesh, nel, p, q );

    %impose boundary conditions
toc
disp('Imposing boundary conditions...')
%[loaddofs,fixeddofs,freedofs,d,dfree,dfixed] =BoundryC2(PHTelem, p, q,sizeBasis);
[loaddofs,fixeddofs,freedofs] =BoundryC(PHTelem, p, q,sizeBasis);
%size(stiff)
%%
%开始进行位移和d的交错更新
% crack phase field solution（裂纹相场解）
u_increment = 1e-5; step_total = 1009; 
u = zeros(2*sizeBasis,1); u(loaddofs) = u_increment;
u = precomutation(u,nel,Cmat,loaddofs,freedofs,Bu,iK,jK,detJacobi);%u0 计算初始位置的计算
Fload = zeros(step_total+1,1); Uload = zeros(step_total+1,1); %存储加载的应力和位移数据点
Hn = zeros(16,nel);
niter=0;
numd=10; a=step_total/numd;
numD=zeros(sizeBasis,numd);
for step = 1:step_total 
    niter=niter+1;
    disp(['       Iter.: ' num2str(niter)])
    %for w=1:3高斯积分
    Hn = Hn_evaluation(u,Hn,Bu,edofMat,lambda,mu,nel);% H计算
    %d计算
    Bterm = kron(gc*l,ones(16,nel));
    Nterm = kron(gc/l,ones(16,nel))+2*Hn;
    sKKd_B=zeros(256,nel); sKKd_D=zeros(256,nel); sFD=zeros(16,nel);
    parfor iel=1:nel
        sKKd_B(:,iel)=Bdp(:,16*(iel-1)+1:16*iel).*kron(detJacobi(:,iel)',ones(256,1))*Bterm(:,iel);
        sKKd_D(:,iel)=Ndp(:,16*(iel-1)+1:16*iel).*kron(detJacobi(:,iel)',ones(256,1))*Nterm(:,iel);
        sFD(:,iel)=2*Nd(:,16*(iel-1)+1:16*iel)'.*kron(detJacobi(:,iel)',ones(16,1))*Hn(:,iel);
    end
    sKd_B=reshape(sKKd_B,256*nel,1);
    sKd_D=reshape(sKKd_D,256*nel,1);
    sKd=sKd_B+sKd_D;
    sFd=reshape(sFD,16*nel,1);
    Kd = sparse(iKd,jKd,sKd);
    Fd = sparse(iFd,jFd,sFd);
    %d(dfree) = Kd(dfree,dfree)\(Fd(dfree)-Kd(dfree,dfixed)*d(dfixed));
    d = Kd\Fd;
    d(d>1)=1;d(d<0)=0;
    %u计算
    epsilon_gauss=zeros(3,16*nel); d_gauss=zeros(1,16*nel); U=u(edofMat); DD=d(dedofMat);
    for iel=1:nel
       epsilon_gauss(:,16*(iel-1)+1:16*iel)=reshape(Bu(:,32*(iel-1)+1:32*iel)*U(iel,:)',3,16); %高斯点上的应变值
       d_gauss(:,16*(iel-1)+1:16*iel)=reshape(Nd(:,16*(iel-1)+1:16*iel)*DD(iel,:)',1,16);%高斯点上的d值
    end
    [R_plus,R_minus] = operator_R(epsilon_gauss);%R+ R-
    P_plus = zeros(9,16*nel); P_minus = P_plus;%  P+ P-
    parfor i=1:16*nel
       [P_plus(:,i),P_minus(:,i)] = operator_P(epsilon_gauss(:,i)); %  P+ P-
    end
    I = [1;1;0]; Ip = reshape(I*I',9,1);
    D = kron(((1-d_gauss).^2+k),ones(9,1)).*(kron(lambda*R_plus,Ip)+2*mu*P_plus)+kron(lambda*R_minus,Ip)+2*mu*P_minus;
    sK = zeros(32*32,nel);
    for iel=1:nel
       sK(:,iel)=reshape(kron(detJacobi(:,iel)',ones(32,3)).*Bu(:,32*(iel-1)+1:32*iel)'*blkdiag(reshape(D(:,16*(iel-1)+1),3,3),reshape(D(:,16*(iel-1)+2),3,3),...
       reshape(D(:,16*(iel-1)+3),3,3),reshape(D(:,16*(iel-1)+4),3,3),reshape(D(:,16*(iel-1)+5),3,3),reshape(D(:,16*(iel-1)+6),3,3),reshape(D(:,16*(iel-1)+7),3,3),...
       reshape(D(:,16*(iel-1)+8),3,3),reshape(D(:,16*(iel-1)+9),3,3),reshape(D(:,16*(iel-1)+10),3,3),reshape(D(:,16*(iel-1)+11),3,3),reshape(D(:,16*(iel-1)+12),3,3),...
       reshape(D(:,16*(iel-1)+13),3,3),reshape(D(:,16*(iel-1)+14),3,3),reshape(D(:,16*(iel-1)+15),3,3),reshape(D(:,16*(iel-1)+16),3,3))*Bu(:,32*(iel-1)+1:32*iel),32*32,1);%,reshape(D(:,25*(iel-1)+17),3,3),...
       %reshape(D(:,25*(iel-1)+18),3,3),reshape(D(:,25*(iel-1)+19),3,3),reshape(D(:,25*(iel-1)+20),3,3),reshape(D(:,25*(iel-1)+21),3,3),reshape(D(:,25*(iel-1)+22),3,3),...
       %reshape(D(:,25*(iel-1)+23),3,3),reshape(D(:,25*(iel-1)+24),3,3),reshape(D(:,25*(iel-1)+25),3,3))*Bu(:,32*(iel-1)+1:32*iel),32*32,1);reshape(kron(detJacobi(:,iel)',ones(32,3)).*
    end
    K = sparse(iK,jK,sK(:));
    u(loaddofs) = u_increment*step;
    u(freedofs) =  -K(freedofs,freedofs)\(K(freedofs,loaddofs)*u(loaddofs)); 
    %end
    Fload(step+1) = sum(K(loaddofs,:)*u);
    Uload(step+1) = u_increment*step;
    %%
    %存储d，进行后处理
    if mod(step,a) == 0 
     n=step/a;
     numD(:,n)=full(d);
    end
end
     toc
   
    disp('Plotting the solution...')
    plotSolPHTElasticMPVM( d, PHTelem, GIFTmesh, p, q )

    disp('find Phase Mesh...')
    [MP_quadRef] = Phase_mesh(numPatches,PHTelem,quadList,d,fai);
    MP_meshend = ones(1,numPatches);
    MP_indexQuad = cell(1, numPatches);
    for patchIndex = 1:numPatches
        MP_indexQuad{patchIndex} = find(MP_quadRef{patchIndex} > 0);
        if ~isempty(MP_indexQuad{patchIndex}) 
            numNewPatches = length(MP_indexQuad{patchIndex});
           disp(['In geometric patch ', num2str(patchIndex), 'Phase_mesh refining ',num2str(numNewPatches), ' quadruplets out of ', num2str(length(quadList{patchIndex}))])
            [quadList{patchIndex}, PHTelem{patchIndex}, dimBasis(patchIndex)] = refinePH_Mesh(MP_quadRef{patchIndex}, quadList{patchIndex}, PHTelem{patchIndex}, p, q, dimBasis(patchIndex));
        else
          MP_meshend(patchIndex) = 0;
        end
    end
keep_refining=sum(MP_meshend);
end
   %vtuFile = 'EdgeCrackMultiPatch.vtu';
   %PlotStressDispMP 
     
   %calculate the error norms 计算误差规范
    %disp('Calculating error norms...')
    %[J,l2relerr, h1relerr] = calcErrorNormsECMP( sol0, PHTelem, GIFTmesh, p, q, Cmat, force, Emod, nu, stressState);
    %l2relerr
    %h1relerr
    %disp(['Effectivity index: ',num2str(estErrorGlobTotal/h1relerr)])
    
    %adaptive refinement 自适应细分
  
%figure
%PlotStressDispMPSol( u, PHTelem, GIFTmesh, p, q, Cmat )
%%
%存储文件，进行后处理（Tecplot打开）
%post_treatment(coor,numD,connect)
%% 
%画u 
plot(Uload,Fload);% axis tight;
xlabel('displacement [mm]'); ylabel('load [KN]');
%title(['tension-y: load-deflection with nel=',num2str(nel),', du=',num2str(u_increment)]); pause(1e-6);
%saveas(gcf,['tension-y_load-deflection_nel_',num2str(nel),'_du_',num2str(u_increment),'.png']); 
%close all; 
% filename = 'test.mat';
% save(filename)