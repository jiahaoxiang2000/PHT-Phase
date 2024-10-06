close all;
clear all;
clc;format long
tic

%% material information
k=1e-19;
gc=2.7; % KN/mm power relase rate
%lambda=121.15; mu=80.7692;
Emod = 210000;nu = 0.3;
%Emod=mu*(2*mu+3*lambda)/(mu+lambda);nu=lambda/(2*(mu+lambda)); %The transformation relationship between Lamay's constant and Poisson's ratio of elastic modulus
lambda=(nu*Emod)/((1+nu)*(1-2*nu)); mu=Emod/(2*(1+nu)); % GPa Elastomer modulus and shear modulus
% elasticity matrix (plane stress)
% Cmat=Emod/(1-nu^2)*[1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];
% elasticity matrix (plane strain)
stressState = 'PLANE_STRAIN';
Cmat=[lambda+2*mu lambda 0;lambda lambda+2*mu 0;0 0 mu];
disp([num2str(toc),'  MATERIAL information'])

%% NURBS geometry and meshing
L = 1;         % Length
H = 1;         % Hight

numPatches = 4;

target_rel_error = 1e-2;
targetScale = 0.5;
fai=0.8;
% degree
p = 3;
q = 3;

dimBasis = zeros(1, numPatches);

GIFTmesh = init2DGeometryGIFTMP(numPatches);
%[ GIFTmesh] = init2DGeometryGIFTMPrectangle(L,H,numPatches);

quadList = cell(numPatches,1);
PHTelem = cell(numPatches, 1);
for i=1:numPatches
    [PHTelem{i}, dimBasis(i)] = initPHTmesh(p,q);
    quadList{i} = 2:5;
end

% the format of patchBoundaries: patchA, edgeA, patchB, edgeB; the edge format: 1-down, 2-right, 3-up, 4-left.
patchBoundaries = {1, 2, 2, 4; 2, 3, 3, 1; 3, 4, 4, 2};

tic % start timing
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
    %???????¦Ë???d?????????
    %step=1500;
    l=(2*0.5)/(2.^level)-0.0003;
    %l=0.125;
    
    u_inc = 1e-5*ones(1500,1);
    u = zeros(2*sizeBasis,1);
    Fload = zeros(length(u_inc)+1,1); Uload = zeros(length(u_inc)+1,1);
    niter=0;
    nel = numbermesh(PHTelem);
    Hn1 = zeros(25,nel);
    for ij = 1:length(u_inc)
        niter=niter+1;
        disp(['       Iter.: ' num2str(niter)])
        %%  assemble the linear system
        % crack phase field solution??????????
        disp('Assembling the linear system...')
        [dstiff,drhs,Hn1] = assembleKd_phase( PHTelem, GIFTmesh, sizeBasis, p, q,lambda,mu,gc,l,u,Hn1);
        toc
        disp('Imposing d_mesh boundary conditions...')
        [ d ] = imposeDirichletd_mesh(dstiff, drhs, PHTelem, p, q,sizeBasis);
        %d(d>1)=1;d(d<0)=0;
        %u????
        disp('Assembling the linear system...')
        [ ustiff] = assembleKu_phase( PHTelem, GIFTmesh, sizeBasis, p, q,lambda,mu,u,d,k);
        %impose  displacement boundary conditions
        toc
        disp('Imposing displacement boundary conditions...')
        [ u,Fload,Uload ] = imposeDirichletu_mesh(ustiff, PHTelem, p, q,sizeBasis,ij,u_inc,Fload,Uload,u);
        %[ ustiff, urhs, bcdof, bcval ] = imposeDirichletESNeuMP_u(ustiff, urhs, PHTelem, GIFTmesh, p, q, ij,u_inc);
        %[ ustiff, urhs, bcdof, bcval,bcdof_up_x  ] = imposeDirichletESMP_u(ustiff, urhs, PHTelem, p, q, ij,u_inc);
        % u = ustiff\urhs;
        % Fload(ij+1) = sum(ustiff(bcdof_up_x,:)*u);
        %Uload(ij+1) = u_inc*ij;
        
        %%
        %?›¥d?????§Ü???
        %if mod(step,a) == 0
        %n=step/a;
        %numD(:,n)=full(d);
        %end
        disp('Plotting the solution...')
        if ij==200||ij==500||ij==1000||ij==1500||ij==1
            plotSolPHTElasticMPVM( d, PHTelem, GIFTmesh, p, q )
        end
    end
    toc
    
    %disp('Plotting the solution...')
    %plotSolPHTElasticMPVM( d, PHTelem, GIFTmesh, p, q )
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
%?›¥????????§Ü?????Tecplot???
%post_treatment(coor,numD,connect)
%%
%??u
plot(Uload,Fload);% axis tight;
xlabel('displacement [mm]'); ylabel('load [N]');
title(['tension-y: load-deflection with nel=',num2str(nel),', du=',num2str(u_increment)]); pause(1e-6);
%close all;
% filename = 'test.mat';
% save(filename)