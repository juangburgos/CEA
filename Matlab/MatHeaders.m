clear all
close all
clc

load('QU_Controller_Parameters.mat');

% change names from matlab to c version

sysA = Ap;           
sysB = Bp;           
sysC = Cp;           
disA = Aj;           
disB = Bj;           
disC = Cj;           
disD = Dj;           
sRoll  = rollparams;  
sPitch = pitchparams;
                       
sysCpos  = Cpos;      
pPhi     = phi;          
pJota    = jota;        
pTheta   = theta;      
pGamma   = gama;       
pOmega   = omega;      
pPsi     = pssi ;

save('QU_Header_Matrices.mat', ...
     'sysA', ...
     'sysB', ...
     'sysC', ...
     'disA', ...
     'disB', ...
     'disC', ...
     'disD', ...
     'sRoll', ...
     'sPitch', ...
     'sysCpos', ...
     'pPhi', ...
     'pJota', ...
     'pTheta', ...
     'pGamma', ...
     'pOmega', ...
     'pPsi');
 
 mat2head('QU_Header_Matrices.mat','matrices.h');