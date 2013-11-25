close all; clear all; clc;

unit=1e-3;
load data_oxi_sfaira;

pt={'f_central in GHz:','dipole axis 1 for x 2 for y 3 for z:','human head axis 1 for x 2 for y 3 for z:','distance from dipole antenna (cm)'};
ans=inputdlg(pt,'Input data 1:',1,{'1','3','2','2'});

f.f_kentrikh=str2num(ans{1})*1e9;
aksonas=str2num(ans{2});
aksonas_kef=str2num(ans{3});
r=str2num(ans{4});

if aksonas_kef==1
    prosanatolismos_kef=[1 0 0];
elseif aksonas_kef==2
    prosanatolismos_kef=[0 1 0];
elseif aksonas_kef==3
    prosanatolismos_kef=[0 0 1];
else
    disp('wrong axis');
end;

if aksonas==1
    prosanatolismos=[1 0 0];
elseif aksonas==2
    prosanatolismos=[0 1 0];
elseif aksonas==3
    prosanatolismos=[0 0 1];
else
    disp('wrong axis');
end;

d_apoakrosfairas_se_cm=r;

Sim_Path='save';
Sim_Path=[Sim_Path '_display_version'];
Sim_CSX='save.xml';

f.euros=f.f_kentrikh/2;
f.f_arxikh=f.f_kentrikh-f.euros;
f.f_telikh=f.f_kentrikh+f.euros;

l=3*1e8/f.f_kentrikh;
lmin=3*1e8/f.f_telikh;

dip_mhkos=0.475*l/unit;

space_step=lmin/20/unit;

enhlikas=9;
paidi=3.5*enhlikas/5;

pt={'adult or child?'};
ans=inputdlg(pt,'Input data 2:',1,{'adult'});

if strcmp(ans,'adult')==1    aktina_sfairas=enhlikas;
else    aktina_sfairas=paidi; 
end;

pt={'Give number of layers'};
ans=inputdlg(pt,'Input data 3:',1,{'1'});

layers=str2num(ans{1});

c=1;

while c<layers+1
    
    pt={'layer name from outside to inside','layer permitivity','layer conductivity (S/m)','layer density mass (kg/m^_3)','layer thickness (cm)'};
    s=['layer data ' num2str(c) ':'];
    ans=inputdlg(pt,s,1,{'skin','40.7','0.65','1100','0.1'});
    strwma{c}.onoma=ans{1};
    strwma{c}.epitr=str2num(ans{2});
    strwma{c}.agwg=str2num(ans{3});
    strwma{c}.puknot=str2num(ans{4});
    strwma{c}.paxos=str2num(ans{5});
    
    if c==1
        strwma{c}.aktina_se_cm=aktina_sfairas;
    else
        strwma{c}.aktina_se_cm=strwma{c-1}.aktina_se_cm-strwma{c-1}.paxos;
    end;
    
    strwma{c}.aktina=strwma{c}.aktina_se_cm*10;
    strwma{c}.kentro_se_cm=strwma{1}.aktina_se_cm+d_apoakrosfairas_se_cm;
    strwma{c}.kentro=strwma{c}.kentro_se_cm*10.*prosanatolismos_kef;
    
    c=c+1;
    
end;

feed_imp=70;

CSX=InitCSX();

mesh.x=SmoothMeshLines(2.0*[-l/unit 0 l/unit],space_step);
mesh.y=SmoothMeshLines(2.0*[-l/unit 0 l/unit],space_step);
mesh.z=SmoothMeshLines(2.0*[-1.25*l/unit -dip_mhkos/2-[2*space_step -space_step]/3 ...
    -space_step 0 space_step dip_mhkos/2+[2*space_step -space_step]/3 1.25*l/unit],space_step);

arxh=prosanatolismos.*[-dip_mhkos/2 -dip_mhkos/2 -dip_mhkos/2];
telos=prosanatolismos.*[dip_mhkos/2 dip_mhkos/2 dip_mhkos/2];

CSX=AddMetal(CSX,'Dipole');
CSX=AddBox(CSX,'Dipole',1,arxh,telos);

for n=1:numel(strwma)
  name = ['strwma_' num2str(n)];
  CSX = AddMaterial( CSX, name );
  CSX = SetMaterialProperty( CSX, name, 'Epsilon', strwma{n}.epitr, 'Kappa', strwma{n}.agwg, 'Density', strwma{n}.puknot);
  CSX = AddSphere( CSX, name, 10+n, [0 0 0], 1,'Transform',{'Scale',strwma{n}.aktina, 'Translate', strwma{1}.kentro} ); 
end

epsilon_anom=60;
k_anom=1.5;
puknothta_mazas_anom=1200;
aktina_anom_se_cm=1;
aktina_anom=aktina_anom_se_cm*10;
kentro_anom=prosanatolismos_kef.*(strwma{1}.kentro(find(prosanatolismos_kef==1))-strwma{layers}.aktina_se_cm*10+1.02*aktina_anom);

CSX=AddMaterial(CSX,'anomoiogeneia');
CSX=SetMaterialProperty(CSX,'anomoiogeneia','Epsilon',epsilon_anom,'Kappa',k_anom,'Density',puknothta_mazas_anom);
CSX=AddSphere(CSX,'anomoiogeneia',10+n+1,[0 0 0],1,'Transform',{'Scale',aktina_anom,'Translate',kentro_anom} ); 

arxh=[-space_step/3 -space_step/3 -space_step];
telos=[space_step/3  space_step/3  space_step];
CSX=AddLumpedPort(CSX,10,1,feed_imp,arxh,telos,prosanatolismos,'excite');

start=1.1*[strwma{1}.kentro(1)-strwma{1}.aktina strwma{1}.kentro(2)-strwma{1}.aktina strwma{1}.kentro(3)-strwma{1}.aktina];
stop=1.1*[strwma{1}.kentro(1)+strwma{1}.aktina strwma{1}.kentro(2)+strwma{1}.aktina strwma{1}.kentro(3)+strwma{1}.aktina];
CSX = AddDump( CSX, 'SAR', 'DumpType', 20, 'Frequency', f.f_kentrikh,'FileType',1,'DumpMode',2);
CSX = AddBox( CSX, 'SAR', 0, start, stop);

arxh=[mesh.x(1)   mesh.y(1)   mesh.z(1) ];
telos=[mesh.x(end) mesh.y(end) mesh.z(end) ];
[CSX nf2ff]=CreateNF2FFBox(CSX,'nf2ff',arxh,telos);

mesh = AddPML( mesh, [1 1 1 1 1 1]*8 );
BC = {'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8'};

CSX = DefineRectGrid( CSX, unit, mesh );

FDTD = InitFDTD();  
FDTD = SetGaussExcite( FDTD, f.f_kentrikh, f.euros );  
FDTD = SetBoundaryCond( FDTD, BC );

[~,~,~] = rmdir(Sim_Path,'s');
[~,~,~] = mkdir(Sim_Path);
WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);
CSXGeomPlot([Sim_Path '/' Sim_CSX]);
RunOpenEMS( Sim_Path, Sim_CSX,'');

freq=linspace(f.f_arxikh,f.f_telikh,1001);
U=ReadUI({'port_ut1','et'},Sim_Path,freq);
I=ReadUI('port_it1',Sim_Path,freq);

Zin=U.FD{1}.val./I.FD{1}.val;

uf_inc = 0.5*(U.FD{1}.val + I.FD{1}.val * feed_imp);
if_inc = 0.5*(I.FD{1}.val + U.FD{1}.val / feed_imp);
uf_ref = U.FD{1}.val - uf_inc;
if_ref = if_inc - I.FD{1}.val;
s11 = uf_ref ./ uf_inc;

P_in = 0.5*U.FD{1}.val .* conj( I.FD{1}.val );

Pin_fk=interp1(freq, P_in, f.f_kentrikh);
s11_fk=interp1(freq, s11, f.f_kentrikh);

P_accepted=Pin_fk*(1-abs(s11_fk)^2);

analush=3;
thetaRange=0:analush:359;
phiRange=0:analush:359;
 
r1=1;

nf2ff = CalcNF2FF(nf2ff, Sim_Path, f.f_kentrikh, thetaRange*pi/180, [0 90]*pi/180, 'Outfile', 'theta_Sweep.h5');

h=figure();
polar( thetaRange/180*pi, pedio_E_oxisfaira_phi_90 );
hold on
polar( thetaRange/180*pi, nf2ff.E_norm{1}(:,2)' ,'r');
ylabel( 'theta / deg' );
title( ['e farfield (V/m) r=' num2str(r1) ' m  phi=90 deg'] );
legend( 'e, no sphere', 'e, with sphere',2, 'Location', 'BestOutside' );
grid;
hold off;

nf2ff = CalcNF2FF(nf2ff, Sim_Path, f.f_kentrikh, 90*pi/180, phiRange*pi/180, 'Outfile1', 'phi_Sweep.h51');

SAR_field = ReadHDF5Dump([Sim_Path '/SAR.h5']);

SAR = SAR_field.FD.values{1};

phiRange = 0:15:360;
thetaRange = 0:10:180;
nf2ff = CalcNF2FF(nf2ff, Sim_Path, f.f_kentrikh, thetaRange*pi/180, phiRange*pi/180, 'Outfile', '3D_Sweep.h5');
E_far_normalized = nf2ff.E_norm{1} / max(nf2ff.E_norm{1}(:));
[theta,phi] = ndgrid(thetaRange/180*pi,phiRange/180*pi);
x = E_far_normalized .* sin(theta) .* cos(phi);
y = E_far_normalized .* sin(theta) .* sin(phi);
z = E_far_normalized .* cos(theta);

tomh=0;
[SAR_field SAR_mesh]=ReadHDF5Dump([Sim_Path '/SAR.h5'],'Range',{[],[],tomh});
figure
[X Y]=ndgrid(SAR_mesh.lines{1},SAR_mesh.lines{2});
h=pcolor(X,Y,SAR_field.FD.values{1}/abs(P_accepted)),colorbar;
set(h,'EdgeColor','none');
title(sprintf('SAR slice z=%s',num2str(tomh)));


figure
tomh=[-0.085:0.005:0.08];

for i=1:length(tomh)
hold on;axis on;  
[SAR_field SAR_mesh]=ReadHDF5Dump([Sim_Path '/SAR.h5'],'Range',{[],[],tomh(i)});
[X Y]=ndgrid(SAR_mesh.lines{1},SAR_mesh.lines{2});
bv=SAR_field.FD.values{1}/abs(P_accepted);
h=pcolor(X,Y,log10(squeeze(bv))),title(sprintf('SAR slice z=%s',num2str(tomh(i))));
pause(0.9);

F(i)=getframe;

clf 
  
end
