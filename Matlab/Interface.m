%% Clear everything and load test data

% 2) Try economy version
% 3) Integrate to pinv in C

% Clear workspace and figures
clear all
close all
clc

% Load DLL
dllname    = 'libCEA.dll';
headername = 'mydll.h';
funcname   = 'mypinv';
funcsql    = 'sqltest';
funcobs    = 'matobserver';
if ~libisloaded( 'MYDLL' ) 
    loadlibrary( dllname, headername, 'alias', 'MYDLL' );      
end

%% Test DLL
% C function:
% mypinv(integer *m, integer *n, doublereal *a, 
%        doublereal *mytol, doublereal *ainv);

A = [1    2
     3    4
     5    6
     7    8];
A = A';
 
tsmat = tic;
Ai_M = pinv(A);
tmat = toc(tsmat);

m     = size(A,1);
n     = size(A,2);
mytol = 0.0;
Ai    = zeros(n,m);

p_A       = libpointer('doublePtr',A);
p_Ai      = libpointer('doublePtr',Ai);
p_mytol   = libpointer('doublePtr',mytol);
p_m       = libpointer('int64Ptr',m);
p_n       = libpointer('int64Ptr',n);

tslap = tic;
calllib( 'MYDLL', funcname, p_m, p_n, p_A, ...
    p_mytol, p_Ai);
tlap = toc(tslap);
Ai_L   = get(p_Ai,'Value');

clc
disp(['A: ',mat2str(A)]);
disp('MATLAB Result: ');
disp(['A^-1: ',mat2str(Ai_M)]);
disp(['MATLAB Time:',num2str(tmat),' s']);
disp('LAPACK Result: ');
disp(['A^-1: ',mat2str(Ai_L)]);
disp(['LAPACK Time :',num2str(tlap),' s']);

Ai_M
Ai_L

%% Test DLL
% C function:
% sqltest(integer *m, integer *n, doublereal *mymat);

% % Aj = zeros(15,15);
% % m     = size(Aj,1);
% % n     = size(Aj,2);
% % p_Aj      = libpointer('doublePtr',Aj);
% % p_m       = libpointer('int64Ptr',m);
% % p_n       = libpointer('int64Ptr',n);
% % calllib( 'MYDLL', funcsql, p_m, p_n, p_Aj );
% % Aj_L   = get(p_Aj,'Value');
% % 
load('QU_Controller_Parameters.mat');
% % Aj
% % Aj_L
% % difference = max(svd(Aj-Aj_L))

%% Test DLL
% C function:
% matobserver(doublereal *u, doublereal *ym, doublereal *y);

u  = zeros(size(Ap));
ym = zeros(size(Aj));
y  = zeros(size(rollparams));
p_u       = libpointer('doublePtr',u);
p_ym      = libpointer('doublePtr',ym);
p_y       = libpointer('doublePtr',y);
calllib( 'MYDLL', funcobs, p_u, p_ym, p_y );
Ap_L           = get(p_u,'Value');
Aj_L           = get(p_ym,'Value');
rollparams_L   = get(p_y,'Value');

diff1 = max(svd(Ap-Ap_L))
diff2 = max(svd(Aj-Aj_L))
diff3 = max(svd(rollparams-rollparams_L))

%% Unload DLL

unloadlibrary('MYDLL');