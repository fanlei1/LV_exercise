function [long_dist] = VesselParams1(network_matrix,n)
% -------------------------------------------------------------------------
%
% VesselParams1
%
% Subroutine to retrieve the order specific vessel parameters.
%
% Written by: Ravi Namani
% Date: July 1, 2013.

Unit_Conversion;

% Initialize vessel parameters
long_dist.ap = zeros(n,1);
long_dist.bp=zeros(n,1);
long_dist.php = zeros(n,1);
long_dist.cp=zeros(n,1);
long_dist.rhoa=zeros(n,1);
long_dist.pha=zeros(n,1);
long_dist.ca=zeros(n,1);
long_dist.ma=zeros(n,1);
long_dist.fmax=zeros(n,1);
long_dist.ktau=zeros(n,1);
long_dist.ka=zeros(n,1);
long_dist.a=zeros(n,1);

% Passive vessel property fit data:
%Column data
% R @ 80mmHg, Ap, Bp, Phip, Cp
% data_p = [3.1	3.20	2.43	0	19.28
% 4.7	4.51	3.59	0	17.24
data_p = [3.1	3.25	1.0	0	19.28
4.7	4.55	1.5	0	17.24
33.12	35.02	10.9	0.612	20.11
49.37	51.38	16.41	1.88	14.23
81.04	85.53	32.17	1.64	21.24
125.97	133.52	51.6+5	0.96	23.54
730.5	829.6	401.52	30.39	3.35
];
% data_p = [3.1	3.1	3.0	0	19.28
% 4.7	4.7	4.5	0	17.24
% 33.12	33.12	32.12	0.0	20.11
% 49.37	49.37	48.37	0.0	14.23
% 81.04	81.04	80.04	0.0	21.24
% 125.97	125.97	125.97	0.0	23.54
% 730.5	730.5	729.4	0.0	3.35
% ];

rp = data_p(:,1);
ap = data_p(:,2);
bp = data_p(:,3);
php = data_p(:,4);
cp = data_p(:,5);

% Active vessel property fit data:
% See excel workbook long_distrib_active_params.xls: 
% R @80mmHg, Rhoa, Phia, ca, fmax, ktau, ka, a
data = [3.1     0           0          10         0        200    1.0e-10     0.3333
            4.7     2.75        20        10.0     0.28    150     1.98e-09    0.3333    
            33.1	26.47	69.39  57.21-10	0.43	156     1.98E-09	0.3333
            49.4	50.96	125.4	87.33-10	0.83	199.5	4.11E-08	0.5714
            81      77.86	148.87	135.4   	1          67.5    3.55E-07    0.5714
            126     74.14	77.26	51.31      	0.62     117     2.15E-07	0.5714
            730.5	0           0       20      0.0         200     1.00E-20	0.01
];

ra = data(:,1);
rhoa = data(:,2);
pha = data(:,3);
ca = data(:,4);
fmax = data(:,5);
ktau = data(:,6);
ka = data(:,7);
a = data(:,8);

% Interpolate passive and active data with a smooth spline
Sap = spline(rp,ap);
Sbp = spline(rp,bp);
Sphp = spline(rp,php);
Scp = spline(rp,cp);
Srhoa = pchip(ra,rhoa);
Spha = spline(ra,pha);
Sca = spline(ra,ca);
Sfmax = spline(ra,fmax);
Sktau = spline(ra,ktau);
Ska = spline(ra,ka);
Sa = spline(ra,a);

% r1_4 = [4.7 6.6 9.5 17.2]';
% ppval(Prhom,r1_4)

% -------------------------------------------------------------------------
% Initialize vessel constants
% -------------------------------------------------------------------------
% Initialize parameter matrix for each order (row).  Note that order 0 is
% row 1:

% Radius of each vessel
r = network_matrix(:,3)/2;
%index1 = find(network_matrix(:,4)==1);
%r(index1) = 4.7;
long_dist.ap = ppval(Sap,r).*um2mm;
long_dist.bp = ppval(Sbp,r).*um2mm;
long_dist.php = ppval(Sphp,r).*mmHg2MPa;
long_dist.cp = ppval(Scp,r).*mmHg2MPa;
long_dist.rhoa = fnval(Srhoa,r).*um2mm;
long_dist.pha = ppval(Spha,r).*mmHg2MPa;
long_dist.ca = ppval(Sca,r).*mmHg2MPa;
long_dist.ma = ones(n,1).*2;
long_dist.fmax = ppval(Sfmax,r);
long_dist.ktau = ppval(Sktau,r).*dynpcm22MPa;
long_dist.ka = ppval(Ska,r);
long_dist.a = ppval(Sa,r);

for i = 1:n
    if (network_matrix(i,4) < 1 )
        long_dist.rhoa(i) = 0;
        long_dist.pha(i) = 0;
        long_dist.ca(i) = 20*mmHg2MPa;
        long_dist.fmax(i) = 0;
        long_dist.ktau(i) = 1000*dynpcm22MPa;
        long_dist.ka(i) = 1e-18;
        long_dist.a(i) = 1e-10;
    end
    
   % Precautionary measure to bound the parameters by making them positive
    if(long_dist.bp(i)<0)
        long_dist.bp(i)=0;
    end
    if(long_dist.cp(i)<0)
        long_dist.cp(i)=1e-4;
    end
    
    if(long_dist.php(i)<0)
        long_dist.php(i)=0;
    end
    
    if(long_dist.rhoa(i)<0)
        long_dist.rhoa(i)=1e-4;
    end
    
    if(long_dist.ca(i)<0)
        long_dist.ca(i)=1e-4;
    end
    
    if(long_dist.fmax(i)<0)
        long_dist.fmax(i)=0;
    end
    
    if(long_dist.ka(i)<0)
        long_dist.ka(i)=1e-10;
    end
    
    if(long_dist.a(i)<=0)
        long_dist.a(i)=0.01;
    end
end
% Update properties to MPa, mm, N, s:
%                   1,          2,          3,                  4,              5,              6,              7,               8,    9,       10,                 11, 12      
%                   ap,         bp,         php,          cp,               rhoa,       pha,            ca,               ma, Fmax, Ktau,               Ka, a,       
% dimscale = [um2mm, um2mm, mmHg2MPa, mmHg2MPa, um2mm, mmHg2MPa, mmHg2MPa, 1,   1,       dynpcm22MPa, 1,  1];
% vparams = vparams.*repmat(dimscale, n, 1);