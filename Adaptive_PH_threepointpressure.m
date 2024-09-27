close all;
clear all;
clc;format long
tic
%% material information
k=1e-19;
gc=5e-3; %N/mm（能量释放率）
%lambda=121.15; mu=80.7692; %GPa（弹性体模量和剪切模量）
%Emod = 210000;nu = 0.3;
%Emod=mu*(2*mu+3*lambda)/(mu+lambda);nu=lambda/(2*(mu+lambda));%拉梅常数与弹性模量泊松比之间的转换关系
%lambda=(nu*Emod)/((1+nu)*(1-2*nu));mu=Emod/(2*(1+nu));
lambda = 12;mu=8;
%elasticity matrix (plane stress)
%Cmat=Emod/(1-nu^2)*[1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];
%elasticity matrix (plane strain)
stressState = 'PLANE_STRAIN';
Cmat=[lambda+2*mu lambda 0;lambda lambda+2*mu 0;0 0 mu];
disp([num2str(toc),'  MATERIAL information'])
%% NURBS geometry and meshing
L = 8;         % Length
H = 2;         % Hight

numPatches = 100;

target_rel_error = 1e-2;
targetScale = 0.5;
fai=0.8;
% degree
p = 3;
q = 3; 

dimBasis = zeros(1, numPatches);

GIFTmesh = init2DGeometry_threepointGIFTMP(numPatches,L,H);


quadList = cell(numPatches,1);
PHTelem = cell(numPatches, 1);
for i=1:numPatches
    [PHTelem{i}, dimBasis(i)] = initPHTmesh(p,q);
    quadList{i} = 2:5;    
end
load('threepoint_patchBoundary');
tic %启动秒表计时器\
keep_refining = 1;
%keep_E_refining = 1;
num_steps = 0;
level=1;
while keep_refining
num_steps = num_steps + 1;
toc
disp(['Step ', num2str(num_steps)])
figure
plotPHTMeshMP( PHTelem, GIFTmesh )

[ PHTelem, dimBasis, quadList ] = checkConforming( PHTelem, dimBasis, patchBoundaries, p, q, quadList );
[ PHTelem, sizeBasis ] = zipConforming( PHTelem, dimBasis, patchBoundaries, p, q);
    %sizeBasis
%% Staggered solution procedure
%开始进行位移和d的交错更新
%step=1500;
l=(4*0.4)/(2.^level);
%l=0.03;

u_inc = [-1e-3*ones(40,1);-5*1e-4*ones(1200,1)];%[1e-4*ones(70,1);1e-5*ones(800,1)];
u = zeros(2*sizeBasis,1);
%u_inc = 1e-5;u = zeros(2*sizeBasis,1);
Fload = zeros(length(u_inc)+1,1); Uload = zeros(length(u_inc)+1,1);
%Fload = zeros(step+1,1); Uload = zeros(step+1,1);
niter=0;
nel = numbermesh(PHTelem);
Hn1 = zeros(16,nel);

%ij = 1:step
for ij = 1:length(u_inc)
    niter=niter+1;
    disp(['       Iter.: ' num2str(niter)])
%%  assemble the linear system
 % crack phase field solution（裂纹相场解）
    disp('Assembling the linear system...')
    [dstiff,drhs,Hn1] = assembleKd_phase( PHTelem, GIFTmesh, sizeBasis, p, q,lambda,mu,gc,l,u,Hn1);
    toc
    disp('Imposing d_mesh boundary conditions...')
    [ d ] = imposeDirichletd_threepoint(dstiff, drhs, PHTelem, p, q,sizeBasis);
    %d = dstiff\drhs;
    %d(d>1)=1;d(d<0)=0;
    %u计算
    disp('Assembling the linear system...')
    [ ustiff] = assembleKu_phase( PHTelem, GIFTmesh, sizeBasis, p, q,lambda,mu,u,d,k);
    %impose  displacement boundary conditions
    toc
    disp('Imposing displacement boundary conditions...')
    [u,Fload,Uload ] = imposeDirichletu_threepoint(ustiff, PHTelem, p, q,sizeBasis,ij,u_inc,Fload,Uload,u);
    
        toc
    disp('Plotting the solution...')
    if ij==40||ij==240||ij==640||ij==1240||ij==1
         plotSolPHTElasticMPVM( d, PHTelem, GIFTmesh, p, q )
    end
end

 disp('find Phase Mesh...')
   [MP_quadRef,level] = Phase_mesh(numPatches,PHTelem,quadList,d,fai);
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

%% 
%画u 
plot(Uload,Fload);% axis tight;
xlabel('displacement [mm]'); ylabel('load [KN]');
title(['tension-y: load-deflection with nel=',num2str(nel),', du=',num2str(u_increment)]); pause(1e-6);
saveas(gcf,['tension-y_load-deflection_nel_',num2str(nel),'_du_',num2str(u_increment),'.png']); 
%close all; 
% filename = 'test.mat';
% save(filename)