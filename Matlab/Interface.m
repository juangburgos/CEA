%% Clear everything and load test data

% 0) Check the contest interface if it uses float or double
%    (change doublePtr to singlePtr accordingly and update DLL)
% 1) Try eliminating lwork
% 2) Try economy version
% 3) Integrate to pinv in C

% Clear workspace and figures
clear all
close all
clc

% Load DLL
dllname    = 'libCEA.dll';
headername = 'mydll.h';
funcname   = 'dgesvd';
if ~libisloaded( 'MYDLL' ) 
    loadlibrary( dllname, headername, 'alias', 'MYDLL' );      
end

%% Test DLL
% C function:
% dgesvd(integer *m, integer *n, doublereal *a, 
% doublereal *s, doublereal *u, doublereal *vt, 
% integer *info);

A = [1    2    4
     3    4    7
     5    6    1
     7    8    6];

tsmat = tic;
[U_M,S_M,V_M] = svd(A);
tmat = toc(tsmat);

tslap = tic;
m = size(A,1);
n = size(A,2);
S = zeros(m,n);
U  = zeros(m,m);
Vt = zeros(n,n);

p_A   = libpointer('doublePtr',A);
p_S   = libpointer('doublePtr',S);
p_U   = libpointer('doublePtr',U);
p_Vt   = libpointer('doublePtr',Vt);
p_m   = libpointer('int64Ptr',m);
p_n   = libpointer('int64Ptr',n);
tslap = tic;
calllib( 'MYDLL', funcname, p_m, p_n, p_A, ...
    p_S, p_U, p_Vt);
tlap = toc(tslap);
U_L   = get(p_U,'Value');
S_L   = get(p_S,'Value');
Vt_L  = get(p_Vt,'Value');
V_L   = (Vt_L)';

clc
disp(['A: ',mat2str(A)]);
disp('MATLAB Result: ');
disp(['U: ',mat2str(U_M)]);
disp(['S: ',mat2str(S_M)]);
disp(['V: ',mat2str(V_M)]);
disp(['MATLAB Time:',num2str(tmat),' s']);
disp('LAPACK Result: ');
disp(['U: ',mat2str(U_L)]);
disp(['S: ',mat2str(S_L)]);
disp(['V: ',mat2str(V_L)]);
disp(['LAPACK Time :',num2str(tlap),' s']);

A
U_M*S_M*V_M'
U_L*S_L*Vt_L
%% Unload DLL

unloadlibrary('MYDLL');