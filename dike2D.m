clc
clear
close all
%% geometry input
free = [
    0.0 10.0;
    10.0 10.0;
    14.0 14.0;
    18.0 14.0;
    30.0 8.0;
    50.0 8.0;
    ];%[m]
wl = [
    0.0 12.0;
    14.0 12.0;
    30.0 7.5;
    50.0 7.5;
    ];%[m]
disc = [1.0 0.5];%[m]
%% soil input
clx = 15.0;%[m]
cly = 1.0;%[m]
E_soil = 1e5;%[kPa]
nu_soil = 0.30;%[-]
gamma_soil = 20.0;%[kN/m3]
phi_soil = [20.0 25.0];%[deg]
coh_soil = [15.0 2.2];%[kPa];
ks_soil = [0.01 0.005];%[m/d]
dilation_factor = 0.1;%[psi/phi[-]]
%% pile input
activePile = false;%[true|false]
pile = [3 6.0];%[node|depth[m]]
I_pile = 5e-5;%[m4]
E_pile = 200e6;%[kPa]
nu_pile = 0.20;%[-]
gamma_pile = 78.5;%[kN/m3]
phi_pile = 0.0;%[deg]
psi_pile = 0.0;%[deg]
coh_pile = 125e3;%[kPa]
%% anchor input
activeAnchor = false;%[true|false]
anchor = [12.0 5.0 8.0];%[height y-axis[m]|base anchor x-axis[m]|base anchor y-axis[m]]
E_anchor = 200e6;%[kPa]
d_anchor = 0.15;%[m]
HoH_anchor = 2.5;%[m]
%% interface input
activeInter = false;%[true|false]
inter = [3 4 3.0];%[nr1|nr2|depth[m]]
E_inter = 0.5e5;%[kPa]
nu_inter = 0.28;%[-]
gamma_inter = 20.0;%[kN/m3]
phi_inter = 18.0;%[deg]
psi_inter = 0.0;%[deg]
coh_inter = 8.5;%[kPa]
%% solver input
lim = 500;%[-]
tol = 1e-4;%[-]
def_scaling = 0.05;%[-]
mc = 1;%[-]
[equivalentPileThickness] = equivalentThicknessRectangle(I_pile);
[intersectionPointList] = findIntersectionPoints(free,wl);
[X,Y] = domainPartition(free,disc,pile,equivalentPileThickness,intersectionPointList,activePile);
[no,nc,el,ec] = meshGenerator(X,Y);
[no,freeSurfaceNode] = elevateMesh(free,X,no,el);
[anchorId,anchorSpringX,anchorSpringY] = anchorSpring(free,pile,anchor,E_anchor,d_anchor,HoH_anchor,no,activeAnchor);
[df] = degreeFreedom(el,ec);
[ad,fh,fs] = boundaryDirichlet(wl,X,no,nc,el,freeSurfaceNode);
[nodalWaterPressure] = boundaryNeumann(wl,no,nc,freeSurfaceNode);
[lip] = integrationPointLocation(no,el,ec);
[L] = randomField(lip(:,1),lip(:,2),clx,cly);
[E,nu] = elasticityParametersIntegrationPoint(free,pile,inter,E_soil,E_pile,E_inter,nu_soil,nu_pile,nu_inter,equivalentPileThickness,ec,lip,activePile,activeInter);
[Km_d] = displacementStiffness(no,nc,el,ec,df,ad,E,nu,anchorId,anchorSpringX,anchorSpringY);
[gamma] = specificWeightIntegrationPoint(free,pile,inter,gamma_soil,gamma_pile,gamma_inter,equivalentPileThickness,ec,lip,activePile,activeInter);
[nodalSelfWeight] = bodyWeight(no,nc,el,ec,gamma);
[H] = initialHead(wl,X,no,el,freeSurfaceNode);
normalFlowRate = zeros(mc,1);
FoS = normalFlowRate;
for i = 1:mc
    [ks,phi,coh] = randomFieldInstance(ks_soil,phi_soil,coh_soil,ec,L);
    [ks] = hydraulicConductivityIntegrationPoint(free,pile,ks,equivalentPileThickness,lip,activePile);
    [phi,psi,coh] = plasticityParametersIntegrationPoint(free,pile,inter,phi,phi_pile,phi_inter,dilation_factor,psi_pile,psi_inter,coh,coh_pile,coh_inter,equivalentPileThickness,lip,activePile,activeInter);
    [Kc] = conductivityStiffness(no,el,ec,ks);
    [H] = hydraulicSolver(H,Kc,no,el,fh,fs,lim,tol);
    [vee,normalFlowRate(i,1)] = flowDarcy(H,no,el,ec,ks);
    [nodalPorePressure] = hydroMechanicalCoupling(H,no,nc,el,ec,df);
    Fext = nodalSelfWeight + nodalWaterPressure + nodalPorePressure;
    [U,evpt1,evpt2,evpt3,evpt4,sigma,FoS(i,1)] = findSafetyFactor(no,nc,el,ec,df,phi,psi,coh,E,nu,Fext,Km_d,ad,lim,tol);
    fprintf('Monte Carlo Iteration %u\n',i)
end
displayResults(free,wl,pile,def_scaling,equivalentPileThickness,anchorId,no,nc,el,ec,nodalWaterPressure,lip,normalFlowRate,FoS,ks,H,U,evpt1,evpt2,evpt3,evpt4,sigma,activePile);
%% functions
function [X,Y] = domainPartition(free,disc,pile,equivalentPileThickness,intersectionPointList,activePile)
ey = max(free(:,2))/disc(2);
if activePile == true
    Y = unique([0:1.0/ey:1.0 (free(pile(1),2) - pile(2))/free(pile(1),2)]);
else
    Y = 0:1.0/ey:1.0;
end
X = free(:,1)';
for i = 1:size(free,1) - 1
    L = free(i + 1,1) - free(i,1);
    ex = ceil(L/disc(1));
    dx = L/ex;
    for j = 1:ex
        X = cat(2,X,free(i,1) + (j - 1)*dx);
    end
end
if activePile == true
    X = cat(2,X,free(pile(1),1) + equivalentPileThickness);
end
if isempty(intersectionPointList) == 0
    X = cat(2,X,intersectionPointList(:,1)');
end
X = unique(X);
end
function [no,nc,el,ec] = meshGenerator(X,Y)
ex = numel(X) - 1;
ey = numel(Y) - 1;
ec = ex*ey;
nx = ex + 1;
ny = ey + 1;
nc = (2*ex + 1)*(2*ey + 1) - ec;
el = zeros(ec,8);
no = zeros(nc,2);
for i = 1:ey
    for j = 1:ex
        el((i - 1)*ex + j,:) = [
            (i - 1)*nx + j,...
            i*nx + j,...
            i*nx + j + 1,...
            (i - 1)*nx + j + 1,...           
            nx*ny + ex + (i - 1)*(nx + ex) + j,...
            nx*ny + i*(nx + ex) + j,...
            nx*ny + ex + (i - 1)*(nx + ex) + j + 1,...
            nx*ny + (i - 1)*(nx + ex) + j,...
            ];
        no(el((i - 1)*ex + j,:),1) = [X(j) X(j) X(j + 1) X(j + 1) X(j) (X(j) + X(j + 1))/2 X(j + 1) (X(j) + X(j + 1))/2];
        no(el((i - 1)*ex + j,:),2) = [Y(i) Y(i + 1) Y(i + 1) Y(i) (Y(i) + Y(i + 1))/2 Y(i + 1) (Y(i) + Y(i + 1))/2 Y(i)];
    end
end
end
function [no,freeSurfaceNode] = elevateMesh(free,X,no,el)
nc = max(max(el(:,1:4)));
nx = numel(X);
ex = nx - 1;
freeSurface = interp1(free(:,1),free(:,2),X);
freeSurface = cat(2,freeSurface,(freeSurface(1:end - 1) + freeSurface(2:end))/2);
for i = 1:nx
    no(i:nx:nc,2) = freeSurface(i)*no(i:nx:nc,2);
end
for i = 1:nx
    no(nc + ex + i:ex + nx:size(no,1),2) = freeSurface(i)*no(nc + ex + i:ex + nx:size(no,1),2);
end
for i = 1:ex
    no(nc + i:ex + nx:size(no,1),2) = freeSurface(nx + i)*no(nc + i:ex + nx:size(no,1),2);
end
freeSurfaceNode = [nc - ex:nc size(no,1) - ex + 1:size(no,1)]';
[~,idSorted] = sort(no(freeSurfaceNode,1));
freeSurfaceNode = freeSurfaceNode(idSorted);
end
function [p,q,t] = stressInvariants(s1,s2,s3,s4)
    q = zeros(size(s1));
    t = zeros(size(s1));
    I1 = s1 + s2 + s4;
    I2 = s1.*s2 + s2.*s4 + s4.*s1 - s3.^2;
    I3 = s1.*s2.*s4 - s3.^2.*s4;
    J2 = 1/3*I1.^2 - I2;
    J3 = 2/27*I1.^3 - 1/3*I1.*I2 + I3;
    p = I1/3;
    id = J2 > eps;
    J2_id = J2(id);
    sqrt_J2_id = sqrt(3*J2_id);
    q(id) = sqrt_J2_id;
    tp = -3*sqrt(3)/2*J3(id)./J2_id.^1.5;
    tp = min(max(tp, -1.0), 1.0);
    t(id) = asin(tp)/3;
end
function [lip] = integrationPointLocation(no,el,ec)
[xi,yi] = integrationPoint();
lip = zeros(4*ec,2);
for i = 1:ec
    co = no(el(i,:),:);
    for j = 1:4
        [N] = shapeTensor8(xi(j),yi(j));
        lip((i - 1)*4 + j,:) = N*co;
    end
end
end
function [var_mean] = meanElement(var)
size_var = numel(var);
var_mean = (var(1:4:size_var) + var(2:4:size_var) + var(3:4:size_var) + var(4:4:size_var))/4;
end
function [ks] = hydraulicConductivityIntegrationPoint(free,pile,ks,equivalentPileThickness,lip,activePile)
if activePile == true
    id = selectPileIntegrationPoints(free,pile,equivalentPileThickness,lip);
    ks(id) = 0.0;
end
end
function [E,nu] = elasticityParametersIntegrationPoint(free,pile,inter,E_soil,E_pile,E_inter,nu_soil,nu_pile,nu_inter,equivalentPileThickness,ec,lip,activePile,activeInter)
E = E_soil*ones(4*ec,1);
nu = nu_soil*ones(4*ec,1);
if activeInter == true
    id = selectInterIntegrationPoints(free,inter,lip);
    E(id) = E_inter;
    nu(id) = nu_inter;
end
if activePile == true
    id = selectPileIntegrationPoints(free,pile,equivalentPileThickness,lip);
    E(id) = E_pile;
    nu(id) = nu_pile;
end
end
function [gamma] = specificWeightIntegrationPoint(free,pile,inter,gamma_soil,gamma_pile,gamma_inter,equivalentPileThickness,ec,lip,activePile,activeInter)
gamma = gamma_soil*ones(4*ec,1);
if activeInter == true
    id = selectInterIntegrationPoints(free,inter,lip);
    gamma(id) = gamma_inter;
end
if activePile == true
    id = selectPileIntegrationPoints(free,pile,equivalentPileThickness,lip);
    gamma(id) = gamma_pile;
end
end
function [phi,psi,coh] = plasticityParametersIntegrationPoint(free,pile,inter,phi,phi_pile,phi_inter,dilation_factor,psi_pile,psi_inter,coh,coh_pile,coh_inter,equivalentPileThickness,lip,activePile,activeInter)
psi = dilation_factor*phi;
if activeInter == true
    id = selectInterIntegrationPoints(free,inter,lip);
    phi(id) = phi_inter;
    psi(id) = psi_inter;
    coh(id) = coh_inter;
end
if activePile == true
    id = selectPileIntegrationPoints(free,pile,equivalentPileThickness,lip);
    phi(id) = phi_pile;
    psi(id) = psi_pile;
    coh(id) = coh_pile;
end
end
function [Lambda,Mu] = elasticityParameters(E,nu)
Lambda = E.*nu./((1 - 2*nu).*(1 + nu));
Mu = E./(2*(1 + nu));
end
function [D,dA] = globalDerivative(dN,co)
J = dN*co;
dA = det(J);
D = J\dN;
end
function [D1_1,D1_2,D1_3,D1_4,D1_5,D1_6,D1_7,D1_8,D2_1,D2_2,D2_3,D2_4,D2_5,D2_6,D2_7,D2_8,dA] = globalDerivativeVectorized(dN1_1,dN1_2,dN1_3,dN1_4,dN1_5,dN1_6,dN1_7,dN1_8,dN2_1,dN2_2,dN2_3,dN2_4,dN2_5,dN2_6,dN2_7,dN2_8,c1_1,c1_2,c2_1,c2_2,c3_1,c3_2,c4_1,c4_2,c5_1,c5_2,c6_1,c6_2,c7_1,c7_2,c8_1,c8_2)
J1_1 = c1_1*dN1_1 + c2_1*dN1_2 + c3_1*dN1_3 + c4_1*dN1_4 + c5_1*dN1_5 + c6_1*dN1_6 + c7_1*dN1_7 + c8_1*dN1_8;
J1_2 = c1_2*dN1_1 + c2_2*dN1_2 + c3_2*dN1_3 + c4_2*dN1_4 + c5_2*dN1_5 + c6_2*dN1_6 + c7_2*dN1_7 + c8_2*dN1_8;
J2_1 = c1_1*dN2_1 + c2_1*dN2_2 + c3_1*dN2_3 + c4_1*dN2_4 + c5_1*dN2_5 + c6_1*dN2_6 + c7_1*dN2_7 + c8_1*dN2_8;
J2_2 = c1_2*dN2_1 + c2_2*dN2_2 + c3_2*dN2_3 + c4_2*dN2_4 + c5_2*dN2_5 + c6_2*dN2_6 + c7_2*dN2_7 + c8_2*dN2_8;
dA = J1_1.*J2_2 - J1_2.*J2_1;
iJ1_1 = J2_2./dA;
iJ1_2 = -J1_2./dA;
iJ2_1 = -J2_1./dA;
iJ2_2 = J1_1./dA;
D1_1 = dN1_1*iJ1_1 + dN2_1*iJ1_2;
D1_2 = dN1_2*iJ1_1 + dN2_2*iJ1_2;
D1_3 = dN1_3*iJ1_1 + dN2_3*iJ1_2;
D1_4 = dN1_4*iJ1_1 + dN2_4*iJ1_2;
D1_5 = dN1_5*iJ1_1 + dN2_5*iJ1_2;
D1_6 = dN1_6*iJ1_1 + dN2_6*iJ1_2;
D1_7 = dN1_7*iJ1_1 + dN2_7*iJ1_2;
D1_8 = dN1_8*iJ1_1 + dN2_8*iJ1_2;
D2_1 = dN1_1*iJ2_1 + dN2_1*iJ2_2;
D2_2 = dN1_2*iJ2_1 + dN2_2*iJ2_2;
D2_3 = dN1_3*iJ2_1 + dN2_3*iJ2_2;
D2_4 = dN1_4*iJ2_1 + dN2_4*iJ2_2;
D2_5 = dN1_5*iJ2_1 + dN2_5*iJ2_2;
D2_6 = dN1_6*iJ2_1 + dN2_6*iJ2_2;
D2_7 = dN1_7*iJ2_1 + dN2_7*iJ2_2;
D2_8 = dN1_8*iJ2_1 + dN2_8*iJ2_2;
end
function [b] = strainDisplacementTensor(D)
b = zeros(3,16);
b(1,1:2:16) = D(1,:);
b(2,2:2:16) = D(2,:);
b(3,2:2:16) = D(1,:);
b(3,1:2:16) = D(2,:);
end
function [N] = shapeTensor4(xi,yi)
xim = 1 - xi;
xip = 1 + xi;
yim = 1 - yi;
yip = 1 + yi;
N1 = 0.25*xim*yim;
N2 = 0.25*xim*yip;
N3 = 0.25*xip*yip;
N4 = 0.25*xip*yim;
N = [N1 N2 N3 N4];
end
function [dN] = localDerivative4(xi,yi)
xim = 1 - xi;
xip = 1 + xi;
yim = 1 - yi;
yip = 1 + yi;
dN1_1 = -0.25*yim;
dN1_2 = -0.25*yip;
dN1_3 = 0.25*yip;
dN1_4 = 0.25*yim;
dN2_1 = -0.25*xim;
dN2_2 = 0.25*xim;
dN2_3 = 0.25*xip;
dN2_4 = -0.25*xip;
dN = [dN1_1 dN1_2 dN1_3 dN1_4;dN2_1 dN2_2 dN2_3 dN2_4;];
end
function [N] = shapeTensor8(xi,yi)
xim = 1 - xi;
xip = 1 + xi;
yim = 1 - yi;
yip = 1 + yi;
N1 = 0.25*xim*yim*(-xi - yi - 1);
N2 = 0.25*xim*yip*(-xi + yi - 1);
N3 = 0.25*xip*yip*(xi + yi - 1);
N4 = 0.25*xip*yim*(xi - yi - 1);
N5 = 0.5*xim*yim*yip;
N6 = 0.5*xim*xip*yip;
N7 = 0.5*xip*yim*yip;
N8 = 0.5*xim*xip*yim;
N = [N1 N2 N3 N4 N5 N6 N7 N8];
end
function [dN,dN1_1,dN1_2,dN1_3,dN1_4,dN1_5,dN1_6,dN1_7,dN1_8,dN2_1,dN2_2,dN2_3,dN2_4,dN2_5,dN2_6,dN2_7,dN2_8] = localDerivative8(xi,yi)
xim = 1 - xi;
xip = 1 + xi;
yim = 1 - yi;
yip = 1 + yi;
dN1_1 = 0.25*yim*(2*xi + yi);
dN1_2 = 0.25*yip*(2*xi - yi);
dN1_3 = 0.25*yip*(2*xi + yi);
dN1_4 = 0.25*yim*(2*xi - yi);
dN1_5 = -0.5*yim*yip;
dN1_6 = -xi*yip;
dN1_7 = 0.5*yim*yip;
dN1_8 = -xi*yim;
dN2_1 = 0.25*xim*(2*yi + xi);
dN2_2 = 0.25*xim*(2*yi - xi);
dN2_3 = 0.25*xip*(xi + 2*yi);
dN2_4 = 0.25*xip*(2*yi - xi);
dN2_5 = -xim*yi;
dN2_6 = 0.5*xim*xip;
dN2_7 = -xip*yi;
dN2_8 = -0.5*xim*xip;
dN = [dN1_1 dN1_2 dN1_3 dN1_4 dN1_5 dN1_6 dN1_7 dN1_8;dN2_1 dN2_2 dN2_3 dN2_4 dN2_5 dN2_6 dN2_7 dN2_8;];
end
function [equivalentPileThickness] = equivalentThicknessRectangle(I)
equivalentPileThickness = (12*I)^(1/3);
end
function [Km_d] = displacementStiffness(no,nc,el,ec,df,ad,E,nu,anchorId,anchorSpringX,anchorSpringY)
[Lambda,Mu] = elasticityParameters(E,nu);
[xi,yi] = integrationPoint();
a1 = zeros(16*16*ec,1);
a2 = a1;
a3 = a1;
for i = 1:ec
    co = no(el(i,:),:);
    km = zeros(16,16);
    for j = 1:4
        id = (i - 1)*4 + j;
        [dN] = localDerivative8(xi(j),yi(j));        
        [D,dA] = globalDerivative(dN,co);
        [b] = strainDisplacementTensor(D);
        Ce = [Lambda(id) + 2*Mu(id) Lambda(id) 0.0;Lambda(id) Lambda(id) + 2*Mu(id) 0.0;0.0 0.0 Mu(id);];
        km = km + b'*Ce*b*dA;
    end
    a1((i - 1)*256 + 1:i*256,1) = repelem(df(i,:),16);
    a2((i - 1)*256 + 1:i*256,1) = repmat(df(i,:),[1 16]);
    a3((i - 1)*256 + 1:i*256,1) = reshape(km',[],1);
end
a1 = cat(1,a1,anchorId*2 - 1);
a1 = cat(1,a1,anchorId*2);
a2 = cat(1,a2,anchorId*2 - 1);
a2 = cat(1,a2,anchorId*2);
a3 = cat(1,a3,anchorSpringX);
a3 = cat(1,a3,anchorSpringY);
Km = sparse(a1,a2,a3,2*nc,2*nc);
Km_d = decomposition(Km(ad,ad));
end
function [Kc] = conductivityStiffness(no,el,ec,ks)
[xi,yi] = integrationPoint();
a1 = zeros(4*4*ec,1);
a2 = a1;
a3 = a1;
for i = 1:ec
    co = no(el(i,1:4),:);
    kc = zeros(4,4);
    for j = 1:4
        id = (i - 1)*4 + j;
        [dN] = localDerivative4(xi(j),yi(j));        
        [D,dA] = globalDerivative(dN,co);
        kc = kc + D'*[ks(id) 0.0;0.0 ks(id);]*D*dA;
    end
    a1((i - 1)*16 + 1:i*16,1) = repelem(el(i,1:4),4);
    a2((i - 1)*16 + 1:i*16,1) = repmat(el(i,1:4),[1 4]);
    a3((i - 1)*16 + 1:i*16,1) = reshape(kc',[],1);
end
Kc = sparse(a1,a2,a3,max(max(el(:,1:4))),max(max(el(:,1:4))));
end
function [nodalSelfWeight] = bodyWeight(no,nc,el,ec,gamma)
[xi,yi] = integrationPoint();
nodalSelfWeight = zeros(2*nc,1);
for i = 1:ec
    co = no(el(i,:),:);
    for j = 1:4
        id = (i - 1)*4 + j;
        [N] = shapeTensor8(xi(j),yi(j));
        [dN] = localDerivative8(xi(j),yi(j));
        [~,dA] = globalDerivative(dN,co);
        nodalSelfWeight(el(i,:)*2) = nodalSelfWeight(el(i,:)*2) - N'*gamma(id)*dA;
    end
end
end
function [df] = degreeFreedom(el,ec)
df = zeros(ec,16);
for i = 1:8
    df(:,i*2 - 1) = el(:,i)*2 - 1;
    df(:,i*2) = el(:,i)*2;
end
end
function [ad,fh,fs] = boundaryDirichlet(wl,X,no,nc,el,freeSurfaceNode)
x0 = find(no(:,1) == min(X));
x1 = find(no(:,1) == max(X));
y0 = find(no(:,2) == 0.0);
rd = unique([x0*2 - 1;x1*2 - 1;y0*2 - 1;y0*2;]);
nd = (1:2*nc)';
ad = nd(ismember(nd,rd) == 0);
freeSurfacePressureNode = freeSurfaceNode(1:2:end);
fsh = interp1(wl(:,1),wl(:,2),no(freeSurfacePressureNode,1));
fh1 = freeSurfacePressureNode(fsh >= no(freeSurfacePressureNode,2));
fh2 = find(no(1:max(max(el(:,1:4))),1) == min(X) & no(1:max(max(el(:,1:4))),2) <= fsh(1));
fh3 = find(no(1:max(max(el(:,1:4))),1) == max(X) & no(1:max(max(el(:,1:4))),2) <= fsh(end));
fs1 = freeSurfacePressureNode(fsh < no(freeSurfacePressureNode,2));
fs2 = find(no(1:max(max(el(:,1:4))),1) == min(X) & no(1:max(max(el(:,1:4))),2) > fsh(1));
fs3 = find(no(1:max(max(el(:,1:4))),1) == max(X) & no(1:max(max(el(:,1:4))),2) > fsh(end));
fh = [fh1;fh2;fh3;];
fs = [fs1;fs2;fs3;];
end
function [nodalWaterPressure] = boundaryNeumann(wl,no,nc,freeSurfaceNode)
freeSurfaceWaterLevel = interp1(wl(:,1),wl(:,2),no(freeSurfaceNode,1));
waterDepth = freeSurfaceWaterLevel - no(freeSurfaceNode,2);
meanWaterDepth = (waterDepth(1:2:end - 2) + waterDepth(3:2:end))/2;
[dL,nxx,nyy] = normalBoundary(no(freeSurfaceNode(1:2:end - 2),:),no(freeSurfaceNode(3:2:end),:));
id1 = 1:2:numel(waterDepth) - 2;
id2 = 2:2:numel(waterDepth) - 1;
id3 = 3:2:numel(waterDepth);
waterPressure = 10.0*meanWaterDepth;
waterPressure(waterPressure <= 0.0) = 0.0;
nodalWaterPressure = zeros(2*nc,1);
nodalWaterPressure(freeSurfaceNode(id1)*2 - 1) = nodalWaterPressure(freeSurfaceNode(id1)*2 - 1) - waterPressure/6.*dL.*nxx;
nodalWaterPressure(freeSurfaceNode(id1)*2) = nodalWaterPressure(freeSurfaceNode(id1)*2) - waterPressure/6.*dL.*nyy;
nodalWaterPressure(freeSurfaceNode(id2)*2 - 1) = nodalWaterPressure(freeSurfaceNode(id2)*2 - 1) - 2*waterPressure/3.*dL.*nxx;
nodalWaterPressure(freeSurfaceNode(id2)*2) = nodalWaterPressure(freeSurfaceNode(id2)*2) - 2*waterPressure/3.*dL.*nyy;
nodalWaterPressure(freeSurfaceNode(id3)*2 - 1) = nodalWaterPressure(freeSurfaceNode(id3)*2 - 1) - waterPressure/6.*dL.*nxx;
nodalWaterPressure(freeSurfaceNode(id3)*2) = nodalWaterPressure(freeSurfaceNode(id3)*2) - waterPressure/6.*dL.*nyy;
end
function [dL,nxx,nyy] = normalBoundary(n1,n2)
dx = n2(:,1) - n1(:,1);
dy = n2(:,2) - n1(:,2);
dL = sqrt(dx.^2 + dy.^2);
nxx = -dy./dL;
nyy = dx./dL;
end
function [id] = selectPileIntegrationPoints(free,pile,equivalentPileThickness,lip)
    id = lip(:,1) >= free(pile(1),1) & lip(:,1) <= free(pile(1),1) + equivalentPileThickness & lip(:,2) >= free(pile(1),2) - pile(2);
end
function [id] = selectInterIntegrationPoints(free,inter,lip)
    id = lip(:,1) >= free(inter(1),1) & lip(:,1) <= free(inter(2),1) & lip(:,2) >= (free(inter(1),2) + free(inter(2),2))/2 - inter(3);
end
function [id] = selectAnchorNode(free,pile,anchor,no)
    [~,id] = min(sqrt((no(:,1) - free(pile(1),1)).^2 + (no(:,2) - anchor(1)).^2));
end
function [id] = selectPileElements(free,pile,no,el)
    pn = find(no(:,1) == free(pile(1),1) & no(:,2) >= free(pile(1),2) - pile(2));
    id = find(ismember(el(:,1),pn));
end
function [snph,csph,snps,cohr,dt] = getShearStrengthParameters(FS,phi,psi,coh,E,nu)
phir = atan(tan(phi*pi/180.0)/FS);
psir = atan(tan(psi*pi/180.0)/FS);
snph = sin(phir);
csph = cos(phir);
snps = sin(psir);
cohr = coh/FS;
dt = (4*(1 + nu).*(1 - 2*nu))./(E.*(1 - 2*nu + snph.^2));
end
function [vId] = vectorDegreeFreedom(df)
vId = [df(:,1);df(:,2);df(:,3);df(:,4);df(:,5);df(:,6);df(:,7);df(:,8);df(:,9);df(:,10);df(:,11);df(:,12);df(:,13);df(:,14);df(:,15);df(:,16);];
end
function [c1_1,c1_2,c2_1,c2_2,c3_1,c3_2,c4_1,c4_2,c5_1,c5_2,c6_1,c6_2,c7_1,c7_2,c8_1,c8_2] = coordinateArrayVectorization(no,el)
c1_1 = no(el(:,1),1);
c1_2 = no(el(:,1),2);
c2_1 = no(el(:,2),1);
c2_2 = no(el(:,2),2);
c3_1 = no(el(:,3),1);
c3_2 = no(el(:,3),2);
c4_1 = no(el(:,4),1);
c4_2 = no(el(:,4),2);
c5_1 = no(el(:,5),1);
c5_2 = no(el(:,5),2);
c6_1 = no(el(:,6),1);
c6_2 = no(el(:,6),2);
c7_1 = no(el(:,7),1);
c7_2 = no(el(:,7),2);
c8_1 = no(el(:,8),1);
c8_2 = no(el(:,8),2);
end
function [U,evpt1,evpt2,evpt3,evpt4,sigma,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16] = arrayDeclarationVectorization(nc,ec)
U = zeros(2*nc,1);
evpt1 = zeros(4*ec,1);
evpt2 = evpt1;
evpt3 = evpt1;
evpt4 = evpt1;
sigma = evpt1;
f1 = zeros(ec,1);
f2 = f1;
f3 = f1;
f4 = f1;
f5 = f1;
f6 = f1;
f7 = f1;
f8 = f1;
f9 = f1;
f10 = f1;
f11 = f1;
f12 = f1;
f13 = f1;
f14 = f1;
f15 = f1;
f16 = f1;
end
function [u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,u15,u16] = displacementVectorization(U,df)
u1 = U(df(:,1));
u2 = U(df(:,2));
u3 = U(df(:,3));
u4 = U(df(:,4));
u5 = U(df(:,5));
u6 = U(df(:,6));
u7 = U(df(:,7));
u8 = U(df(:,8));
u9 = U(df(:,9));
u10 = U(df(:,10));
u11 = U(df(:,11));
u12 = U(df(:,12));
u13 = U(df(:,13));
u14 = U(df(:,14));
u15 = U(df(:,15));
u16 = U(df(:,16));
end
function [e1,e2,e3,e4] = strainDisplacementVectorized(u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,u15,u16,evpt1,evpt2,evpt3,evpt4,D1_1,D1_2,D1_3,D1_4,D1_5,D1_6,D1_7,D1_8,D2_1,D2_2,D2_3,D2_4,D2_5,D2_6,D2_7,D2_8)
e1 = D1_1.*u1 + D1_2.*u3 + D1_3.*u5 + D1_4.*u7 + D1_5.*u9 + D1_6.*u11 + D1_7.*u13 + D1_8.*u15;
e2 = D2_1.*u2 + D2_2.*u4 + D2_3.*u6 + D2_4.*u8 + D2_5.*u10 + D2_6.*u12 + D2_7.*u14 + D2_8.*u16;
e3 = D1_1.*u2 + D1_2.*u4 + D1_3.*u6 + D1_4.*u8 + D2_1.*u1 + D1_5.*u10 + D2_2.*u3 + D1_6.*u12 + D2_3.*u5 + D1_7.*u14 + D2_4.*u7 + D1_8.*u16 + D2_5.*u9 + D2_6.*u11 + D2_7.*u13 + D2_8.*u15;
e1 = e1 - evpt1;
e2 = e2 - evpt2;
e3 = e3 - evpt3;
e4 = -evpt4;
end
function [s1,s2,s3,s4] = stressStrainVectorized(Lambda,Mu,e1,e2,e3,e4)
s1 = (Lambda + 2*Mu).*e1 + Lambda.*(e2 + e4);
s2 = (Lambda + 2*Mu).*e2 + Lambda.*(e1 + e4);
s3 = Mu.*e3;
s4 = (Lambda + 2*Mu).*e4 + Lambda.*(e1 + e2);
end
function [anchorId,anchorSpringX,anchorSpringY] = anchorSpring(free,pile,anchor,E_anchor,d_anchor,HoH_anchor,no,activeAnchor)
if activeAnchor == true
    [anchorId] = selectAnchorNode(free,pile,anchor,no);
    A_anchor = pi*(d_anchor/2)^2;
    dx = no(anchorId,1) - anchor(2);
    dy = no(anchorId,2) - anchor(3);
    L_anchor = sqrt(dx^2 + dy^2);
    angle_anchor = asind(dy/L_anchor);
    equivalentAnchorSpring = E_anchor*A_anchor/(L_anchor*HoH_anchor);
    anchorSpringX = equivalentAnchorSpring*abs(cosd(angle_anchor));
    anchorSpringY = equivalentAnchorSpring*abs(sind(angle_anchor));
else
    anchorId = [];
    anchorSpringX = [];
    anchorSpringY = [];    
end
end
function [dQ1,dQ2,dQ3,dQ4] = dervativeMohrCoulomb(s1,s2,s3,s4,snps)
s1_c = s1 + 1i*eps;
s2_c = s2 + 1i*eps;
s3_c = s3 + 1i*eps;
s4_c = s4 + 1i*eps;
[p,q,t] = stressInvariants(s1_c,s2,s3,s4);
Q1 = p.*snps + q.*(cos(t)/sqrt(3.0) - snps.*sin(t)/3.0);
[p,q,t] = stressInvariants(s1,s2_c,s3,s4);
Q2 = p.*snps + q.*(cos(t)/sqrt(3.0) - snps.*sin(t)/3.0);
[p,q,t] = stressInvariants(s1,s2,s3_c,s4);
Q3 = p.*snps + q.*(cos(t)/sqrt(3.0) - snps.*sin(t)/3.0);
[p,q,t] = stressInvariants(s1,s2,s3,s4_c);
Q4 = p.*snps + q.*(cos(t)/sqrt(3.0) - snps.*sin(t)/3.0);
dQ1 = imag(Q1)/eps;
dQ2 = imag(Q2)/eps;
dQ3 = imag(Q3)/eps;
dQ4 = imag(Q4)/eps;
end
function [F,sigma] = yieldMohrCoulomb(s1,s2,s3,s4,snph,csph,coh)
[p,q,t] = stressInvariants(s1,s2,s3,s4);
a1 = p.*snph;
a2 = coh.*csph;
a3 = cos(t)/sqrt(3.0) - snph.*sin(t)/3.0;
F = a1 + q.*a3 - a2;
sigma = (a2 - a1)./a3;
end
function [CONV,U,evpt1,evpt2,evpt3,evpt4,sigma] = displacementSolver(SF,no,nc,el,ec,df,phi,psi,coh,E,nu,Fext,Km_d,ad,lim,tol)
    CONV = true;
    [Lambda,Mu] = elasticityParameters(E,nu);
    [xi,yi] = integrationPoint();
    [c1_1,c1_2,c2_1,c2_2,c3_1,c3_2,c4_1,c4_2,c5_1,c5_2,c6_1,c6_2,c7_1,c7_2,c8_1,c8_2] = coordinateArrayVectorization(no,el);
    [vId] = vectorDegreeFreedom(df);
    [snph,csph,snps,cohr,dt] = getShearStrengthParameters(SF,phi,psi,coh,E,nu);
    [U,evpt1,evpt2,evpt3,evpt4,sigma,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16] = arrayDeclarationVectorization(nc,ec);
    U(ad) = Km_d\Fext(ad);
    for k = 1:lim
        [u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,u15,u16] = displacementVectorization(U,df);
        for j = 1:4
            id = j:4:4*ec;
            [~,dN1_1,dN1_2,dN1_3,dN1_4,dN1_5,dN1_6,dN1_7,dN1_8,dN2_1,dN2_2,dN2_3,dN2_4,dN2_5,dN2_6,dN2_7,dN2_8] = localDerivative8(xi(j),yi(j));
            [D1_1,D1_2,D1_3,D1_4,D1_5,D1_6,D1_7,D1_8,D2_1,D2_2,D2_3,D2_4,D2_5,D2_6,D2_7,D2_8,dA] = globalDerivativeVectorized(dN1_1,dN1_2,dN1_3,dN1_4,dN1_5,dN1_6,dN1_7,dN1_8,dN2_1,dN2_2,dN2_3,dN2_4,dN2_5,dN2_6,dN2_7,dN2_8,c1_1,c1_2,c2_1,c2_2,c3_1,c3_2,c4_1,c4_2,c5_1,c5_2,c6_1,c6_2,c7_1,c7_2,c8_1,c8_2);
            [e1,e2,e3,e4] = strainDisplacementVectorized(u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,u15,u16,evpt1(id),evpt2(id),evpt3(id),evpt4(id),D1_1,D1_2,D1_3,D1_4,D1_5,D1_6,D1_7,D1_8,D2_1,D2_2,D2_3,D2_4,D2_5,D2_6,D2_7,D2_8);
            [s1,s2,s3,s4] = stressStrainVectorized(Lambda(id),Mu(id),e1,e2,e3,e4);
            [F,sigma(id)] = yieldMohrCoulomb(s1,s2,s3,s4,snph(id),csph(id),cohr(id));
            F(F < 0.0) = 0.0;
            [dQ1,dQ2,dQ3,dQ4] = dervativeMohrCoulomb(s1,s2,s3,s4,snps(id));
            evp1 = dt(id).*F.*dQ1;
            evp2 = dt(id).*F.*dQ2;
            evp3 = dt(id).*F.*dQ3;
            evp4 = dt(id).*F.*dQ4;
            evpt1(id) = evpt1(id) + evp1;
            evpt2(id) = evpt2(id) + evp2;
            evpt3(id) = evpt3(id) + evp3;
            evpt4(id) = evpt4(id) + evp4;
            [sp1,sp2,sp3] = stressStrainVectorized(Lambda(id),Mu(id),evp1,evp2,evp3,evp4);
            f1 = f1 + dA.*(D1_1.*sp1 + D2_1.*sp3);
            f2 = f2 + dA.*(D1_1.*sp3 + D2_1.*sp2);
            f3 = f3 + dA.*(D1_2.*sp1 + D2_2.*sp3);
            f4 = f4 + dA.*(D1_2.*sp3 + D2_2.*sp2);
            f5 = f5 + dA.*(D1_3.*sp1 + D2_3.*sp3);
            f6 = f6 + dA.*(D1_3.*sp3 + D2_3.*sp2);
            f7 = f7 + dA.*(D1_4.*sp1 + D2_4.*sp3);
            f8 = f8 + dA.*(D1_4.*sp3 + D2_4.*sp2);
            f9 = f9 + dA.*(D1_5.*sp1 + D2_5.*sp3);
            f10 = f10 + dA.*(D1_5.*sp3 + D2_5.*sp2);
            f11 = f11 + dA.*(D1_6.*sp1 + D2_6.*sp3);
            f12 = f12 + dA.*(D1_6.*sp3 + D2_6.*sp2);
            f13 = f13 + dA.*(D1_7.*sp1 + D2_7.*sp3);
            f14 = f14 + dA.*(D1_7.*sp3 + D2_7.*sp2);
            f15 = f15 + dA.*(D1_8.*sp1 + D2_8.*sp3);
            f16 = f16 + dA.*(D1_8.*sp3 + D2_8.*sp2);            
        end
        Fint = accumarray(vId,[f1;f2;f3;f4;f5;f6;f7;f8;f9;f10;f11;f12;f13;f14;f15;f16;],[2*nc 1]);
        loads = Fext + Fint;
        U0 = U;
        U(ad) = Km_d\loads(ad);
        er = max(abs(U - U0))/max(abs(U0));
        if er < tol
            break
        end
    end
    if k == lim
        CONV = false;
    end
end
function [H] = hydraulicSolver(H,Kc,no,el,fh,fs,lim,tol)
nh = (1:max(max(el(:,1:4))))';
afs = [];
ah = nh(ismember(nh,[fh;afs;]) == 0);
for k = 1:lim
    Qint = Kc*H;
    H0 = H;
    H(ah) = H(ah) - Kc(ah,ah)\Qint(ah);
    er = max(abs(H - H0))/max(abs(H0));
    if er < tol
        break
    elseif k == lim
        error('no conv (hydraulic phase)')
    end
    [Hs,fsId] = max(H(fs) - no(fs,2));
    if Hs > 0.0
        afs = cat(1,afs,fs(fsId));
    end
    H(afs) = no(afs,2);
    ah = nh(ismember(nh,[fh;afs;]) == 0);
end
end
function [intersectionPoint] = crossCheck(p1,p2,p3,p4)
v1 = p2 - p1;
v2 = p4 - p3;
denom = crossProduct2D(v1, v2);
if denom == 0
    intersectionPoint = [];
    return;
end
v3 = p3 - p1;
t1 = crossProduct2D(v3, v2) / denom;
t2 = crossProduct2D(v3, v1) / denom;
if t1 > 0 && t1 < 1 && t2 > 0 && t2 < 1
    intersectionPoint = p1 + t1 * v1;
else
    intersectionPoint = [];
end
end
function cp = crossProduct2D(v1, v2)
    cp = v1(1) * v2(2) - v1(2) * v2(1);
end
function [intersectionPointList] = findIntersectionPoints(free,wl)
intersectionPointList = [];
for i = 1:size(free,1) - 1
    for j = 1:size(wl,1) - 1
        [intersectionPoint] = crossCheck(free(i,:),free(i + 1,:),wl(j,:),wl(j + 1,:));
        intersectionPointList = cat(1,intersectionPointList,intersectionPoint);
    end
end
end
function [U,evpt1,evpt2,evpt3,evpt4,sigma,FoS] = findSafetyFactor(no,nc,el,ec,df,phi,psi,coh,E,nu,Fext,Km_d,ad,lim,tol)
SF = 0.1:0.1:5.0;
for i = 1:numel(SF)
    [CONV,U,evpt1,evpt2,evpt3,evpt4,sigma] = displacementSolver(SF(i),no,nc,el,ec,df,phi,psi,coh,E,nu,Fext,Km_d,ad,lim,tol);
    if CONV == false
        break
    end
end
FoS = SF(i);
end
function [H] = initialHead(wl,X,no,el,freeSurfaceNode)
fsh = interp1(wl(:,1),wl(:,2),no(freeSurfaceNode(1:2:end),1));
nx = numel(X);
H = zeros(max(max(el(:,1:4))),1);
for i = 1:nx
    H(i:nx:max(max(el(:,1:4))),1) = fsh(i);
end
end
function [nodalPorePressure] = hydroMechanicalCoupling(H,no,nc,el,ec,df)
[xi,yi] = integrationPoint();
P = (H - no(1:max(max(el(:,1:4))),2))*10.0;
P(P < 0.0) = 0.0;
nodalPorePressure = zeros(2*nc,1);
for i = 1:ec
    co = no(el(i,:),:);
    for j = 1:4
        [N] = shapeTensor4(xi(j),yi(j));
        [dN] = localDerivative8(xi(j),yi(j));        
        [D,dA] = globalDerivative(dN,co);
        [b] = strainDisplacementTensor(D);
        nodalPorePressure(df(i,:)) = nodalPorePressure(df(i,:)) + b'*[1;1;0;]*(N*P(el(i,1:4)))*dA;
    end
end
end
function [xi,yi] = integrationPoint()
xi = 1/sqrt(3)*[-1 1 -1 1];
yi = 1/sqrt(3)*[1 1 -1 -1];
end
function [vee,normalFlowRate] = flowDarcy(H,no,el,ec,ks)
[xi,yi] = integrationPoint();
vee = zeros(2,4*ec);
flowRate = zeros(4*ec,1);
for i = 1:ec
    co = no(el(i,1:4),:);
    for j = 1:4
        id = (i - 1)*4 + j;
        [N] = shapeTensor4(xi(j),yi(j));
        [dN] = localDerivative4(xi(j),yi(j));        
        [D,dA] = globalDerivative(dN,co);
        pw = N*(H(el(i,1:4)) - co(:,2))*10.0;
        if pw > 0.0
            vee(:,id) = -[ks(id) 0.0;0.0 ks(id);]*D*H(el(i,1:4));
        end
        flowRate(el(i,1:4)) = flowRate(el(i,1:4)) - D'*vee(:,id)*dA;
    end
end
normalFlowRate = norm(flowRate);
end
function [L] = randomField(x,y,clx,cly)
cfun = @(x,y) exp(-sqrt((x - x').^2/clx^2 + (y - y').^2/cly^2));
L = chol(cfun(x,y),'lower');
end
function [varField] = logNormDistribution(var,z)
log_norm_sig = sqrt(log(1 + var(2)^2/var(1)^2));
log_norm_mu = log(var(1)) - 1/2*log_norm_sig^2;
varField = exp(log_norm_mu + log_norm_sig*z);
end
function [ks,phi,coh] = randomFieldInstance(ks_soil,phi_soil,coh_soil,ec,L)
z = L*randn(4*ec,1);
[ks] = logNormDistribution(ks_soil,z);
phi = phi_soil(1) + (phi_soil(2) - phi_soil(1))*((1 + tanh(z))/2);
[coh] = logNormDistribution(coh_soil,z);
end
function displayResults(free,wl,pile,def_scaling,equivalentPileThickness,anchorId,no,nc,el,ec,nodalWaterPressure,lip,normalFlowRate,FoS,ks,H,U,evpt1,evpt2,evpt3,evpt4,sigma,activePile)
P = (H - no(1:max(max(el(:,1:4))),2))*10.0;
P(P < 0.0) = 0.0;
if activePile == true
    id = selectPileIntegrationPoints(free,pile,equivalentPileThickness,lip);
    sigma(id) = NaN;
    H(id) = NaN;
    P(id) = NaN;
    [pileId] = selectPileElements(free,pile,no,el);
else
    pileId = [];
end
H_mean = 0.25*(H(el(:,1)) + H(el(:,2)) + H(el(:,3)) + H(el(:,4)));
P_mean = 0.25*(P(el(:,1)) + P(el(:,2)) + P(el(:,3)) + P(el(:,4)));
[sigma_mean] = meanElement(sigma);
[evpt1_mean] = meanElement(evpt1);
[evpt2_mean] = meanElement(evpt2);
[evpt3_mean] = meanElement(evpt3);
[evpt4_mean] = meanElement(evpt4);
[ks_mean] = meanElement(ks);
evpt_mean = sqrt(evpt1_mean.^2 + evpt2_mean.^2 + evpt3_mean.^2 + evpt4_mean.^2);
scaling_factor = def_scaling*max(max(no))/max(abs(U));
U_scaling = scaling_factor*U;
no_d = zeros(size(no));
no_d(:,1) = no(:,1) + U_scaling(1:2:2*nc);
no_d(:,2) = no(:,2) + U_scaling(2:2:2*nc);
[vee] = flowDarcy(H,no,el,ec,ks);
%1.0 hydraulic conductivity
figure('Name',['Model Conditions, |Flow Rate| ',num2str(normalFlowRate(end)),' [m3/s]'])
for i = 1:ec
    fill(no(el(i,[1,5,2,6,3,7,4,8]),1),no(el(i,[1,5,2,6,3,7,4,8]),2),ks_mean(i),'EdgeColor','none');
    hold on
end
hold on
for i = 1:size(wl,1) - 1
    plot([wl(i,1) wl(i + 1,1)],[wl(i,2) wl(i + 1,2)],'-red','LineWidth',2)
end
hold on
quiver(no(:,1),no(:,2),-nodalWaterPressure(1:2:end),-nodalWaterPressure(2:2:end))
hold on
scatter(no(anchorId,1),no(anchorId,2),40,'green','s')
colorbar
axis equal tight
xlabel('Width [m]','FontSize',14)
ylabel('Height [m]','FontSize',14)
title('Hydraulic Conductivity [m/s]','FontSize',14)
%2.0 total head & Darcy flow
figure('Name',['Model Conditions, |Flow Rate| ',num2str(normalFlowRate(end)),' [m3/d]'])
for i = 1:ec
    fill(no(el(i,[1,5,2,6,3,7,4,8]),1),no(el(i,[1,5,2,6,3,7,4,8]),2),H_mean(i),'EdgeColor','none');
    hold on
end
quiver(lip(:,1),lip(:,2),vee(1,:)',vee(2,:)')
colorbar
axis equal tight
xlabel('Width [m]','FontSize',14)
ylabel('Height [m]','FontSize',14)
title('Total Head [m] & Darcy flow [m/s]','FontSize',14)
%3.0 porepressure & Darcy flow
figure('Name',['Model Conditions, |Flow Rate| ',num2str(normalFlowRate(end)),' [m3/s]'])
for i = 1:ec
    fill(no(el(i,[1,5,2,6,3,7,4,8]),1),no(el(i,[1,5,2,6,3,7,4,8]),2),P_mean(i),'EdgeColor','none');
    hold on
end
quiver(lip(:,1),lip(:,2),vee(1,:)',vee(2,:)')
colorbar
axis equal tight
xlabel('Width [m]','FontSize',14)
ylabel('Height [m]','FontSize',14)
title('Porepressure [m] & Darcy flow [m/s]','FontSize',14)
%4.0 shear strength
figure('Name',['Model Conditions, Factor of Safety: ',num2str(FoS(end)),' [-]'])
for i = 1:ec
    fill(no(el(i,[1,5,2,6,3,7,4,8]),1),no(el(i,[1,5,2,6,3,7,4,8]),2),sigma_mean(i),'EdgeColor','none');
    hold on
end
hold on
for i = 1:size(wl,1) - 1
    plot([wl(i,1) wl(i + 1,1)],[wl(i,2) wl(i + 1,2)],'-red','LineWidth',2)
end
hold on
quiver(no(:,1),no(:,2),-nodalWaterPressure(1:2:end),-nodalWaterPressure(2:2:end))
hold on
scatter(no(anchorId,1),no(anchorId,2),40,'green','s')
colorbar
axis equal tight
xlabel('Width [m]','FontSize',14)
ylabel('Height [m]','FontSize',14)
title('Shear Strength [kPa]','FontSize',14)
%5.0 viscoplastic strain
figure('Name',['Viscoplastic Strain, Factor of Safety: ',num2str(FoS(end)),' [-]'])
for i = 1:ec
    fill(no(el(i,[1,5,2,6,3,7,4,8]),1),no(el(i,[1,5,2,6,3,7,4,8]),2),evpt_mean(i),'EdgeColor','none');
    hold on
end
colorbar
axis equal tight
xlabel('Width [m]','FontSize',14)
ylabel('Height [m]','FontSize',14)
title('Viscoplastic Strain [-]','FontSize',14)
%6.0 deformed mesh
figure('Name',['Viscoplastic Strain, Factor of Safety: ',num2str(FoS(end)),' [-]'])
for i = 1:ec
    fill(no_d(el(i,[1,5,2,6,3,7,4,8]),1),no_d(el(i,[1,5,2,6,3,7,4,8]),2),'white')
    hold on
end
for i = 1:numel(pileId)
    fill(no_d(el(pileId(i),[1,5,2,6,3,7,4,8]),1),no_d(el(pileId(i),[1,5,2,6,3,7,4,8]),2),'red')
    hold on
end
scatter(no_d(anchorId,1),no_d(anchorId,2),40,'green','s')
axis equal tight
xlabel('Width [m]','FontSize',14)
ylabel('Height [m]','FontSize',14)
title(['Deformed Mesh Scaling: ',num2str(scaling_factor),' [-]'],'FontSize',14)
%7.0 histogram
if numel(normalFlowRate) > 1
    figure('name','Monte Carlo results')
    subplot(1,2,1)
    histogram(normalFlowRate)
    xlabel('Bin','FontSize',14)
    ylabel('Instances','FontSize',14)
    title('|Flow Rate| [m3/s]','FontSize',14)
    subplot(1,2,2)
    histogram(FoS)
    xlabel('Bin','FontSize',14)
    ylabel('Instances','FontSize',14)
    title('Factor of Safety [-]','FontSize',14)
end
end
