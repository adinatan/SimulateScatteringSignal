function GS_FormFactor = SimulateScatteringSignal(GS_StructureFile)

XrayEnergy_in_keV=9.5;

IntensityTable=kron([XrayEnergy_in_keV,1],[1;1]);
OutputPath='./out';

AngleList=[0.01:0.01:90];


%GS_StructureFile=GS_StructureFile
[Name, X, Y, Z]=textread(GS_StructureFile, '%s%f%f%f', 'headerlines', 0);

%%%%%%%%%%%%%%%%%%%%%%%
for i=1:numel(X)
    GSatomList(i).name=Name{i};
    GSatomList(i).x=X(i);
    GSatomList(i).y=Y(i);
    GSatomList(i).z=Z(i);
    GSatomList(i).number=i;
    [GSatomList(i).number, GSatomList(i).CMa, GSatomList(i).CMb, GSatomList(i).CMc, GSatomList(i).Z]=CromerMann(GSatomList(i).name);
end

GS_FormFactor=MolscatPolychromCorrect3(IntensityTable,AngleList, GSatomList);
Qlist=GS_FormFactor(:,2);

GS_FormFactor=[Qlist GS_FormFactor(:,3)];

end


function [PolyForm]=MolscatPolychromCorrect3(IntensityTable,AngleRange, atom)

for E_value=1:length(IntensityTable)
    Qvector=(4*pi / (12.398/IntensityTable(E_value,1))) .* sind(AngleRange./2);
    for n=1:numel(atom)
        %test=atom(n);
        atom(n).f0=atomic_f0(atom(n),Qvector); %calculate atomic form factor for atom n
    end
    Favg2(:,E_value)=F_orient_av(Qvector,atom);
    Qvec(:,E_value)=Qvector;
end


for E_value=1:length(IntensityTable)
    %intensity=IntensityTable(E_value,2);
    Favg2Weighted(:,E_value)=Favg2(:,E_value).*IntensityTable(E_value,2);
    QvecWeighted(:,E_value)=Qvec(:,E_value).*IntensityTable(E_value,2);
end

Favg2Corrected=sum(Favg2Weighted,2)/sum(IntensityTable(:,2));
QvecCorrected=sum(QvecWeighted,2)/sum(IntensityTable(:,2));
WeightedMeanEnergy=sum(IntensityTable(:,1).*IntensityTable(:,2))/sum(IntensityTable(:,2));

Angles=AngleRange;
MeanQs=rot90(QvecCorrected);
CorrectedFavg2=rot90(Favg2Corrected);
PolyForm=[Angles;MeanQs;CorrectedFavg2]';

fclose('all');
end

function [Favg2]=F_orient_av(Qvector,AtomList)
% Martin Meedom Nielsen, March 2005
% Modified, K. Haldrup, 2008
NumberOfAtoms=numel(AtomList); %number of atoms in molecule
r=[[AtomList.x]' [AtomList.y]' [AtomList.z]'];
f0=reshape([AtomList.f0],[],NumberOfAtoms);
Favg2=0;%f0(:,NumberOfAtoms).*conj(f0(:,NumberOfAtoms)); % self term
Qvector=Qvector';
AtomNumber=0;
for n=1:NumberOfAtoms
    Favg2=Favg2+f0(:,n).*conj(f0(:,n));
    for m=(n+1):NumberOfAtoms
        Qrnm=Qvector*norm(r(n,:)-r(m,:));
        val=2*f0(:,n).*f0(:,m).*(sin(Qrnm) ./Qrnm); %"both ways"
        Favg2=Favg2+val;%2*f0(:,n).*f0(:,m).*(sin(Qrnm)./Qrnm);
    end
end
end

function [atno,a,b,c,Z]=CromerMann(name)

fid=fopen('f0_CromerMann.txt');
if fid==-1
    error('Cannot find the f0_CromerMann.txt file')
end
CurrentLine=fgetl(fid);
found=0;
while ~strncmpi(CurrentLine,'#S',2) %While the first two characters of CurrentLine are NOT "#S", continue reading next line
    if CurrentLine==-1
        error('Element entry not found1');
    end
    CurrentLine=fgetl(fid);
end

% Now we are in the data-section of CromerMann_f0.txt
name=name;
dums=CurrentLine;
Z=0;
while ~found
    if dums==-1
        error('Element entry not found2');
    end
    ss=sscanf(dums,'%2c %i %c');    %scans for 2 characters, an integer and a number of characters          %% ss(1:2) ='#S*, ss(3)=atomic number, ss(4:end)= name
    bb=ss(4:end)'; %Last three characters
    element=char(bb(bb~=32)); %removes whitespaces
    if strcmp(element,deblank(name)) %if this element is the one we are looking for, exit while-loop next time
        found=1;
    end
    dumn   = fgetl(fid);
    duml   = fgetl(fid);
    dumacb = fgetl(fid);
    dums   = fgetl(fid);
end

if found
    atno=ss(3);
    abc=sscanf(dumacb,'%f %f %f %f %f %f %f %f %f');
    a=abc(1:4);
    c=abc(5);
    b=abc(6:9);
else
    error('Element entry not found');
end
fclose('all');
end

function [f0]=atomic_f0(atom,Q)
f0=atom.CMc;
q2=(Q/(4*pi)).^2;
a=atom.CMa;
b=-atom.CMb;
for n=1:4
    f0=f0+a(n)*exp(b(n)*q2);
end
end