param Nnodos; #Nodos
param Nlines; #Lineas
param Npor;

set I := {1..Nnodos}; #Conjunto de nodos (numeracion)
set L := {1..Nlines}; #Conjunto de lineas (numeracion)
set LINES within (I cross I);
param Reacl {(i,j) in LINES}; #Transmission line reactance
param B {(i,j) in {I,I}}; #:=
#						   if (i=j) then (sum{(i,k) in LINES}(1/Reacl[i,k])+sum{(k,i) in LINES}(1/Reacl[k,i])) 
#						   else if(i<>j) then sum{(m,l) in LINES:m=i&&l=j||l=i&&m=j}(-1/Reacl[m,l]);
param CO{i in I, j in 1..3}; #Costos por nodo
param Pg{i in I, j in 1..2}; #Potencia (Min/Max) por nodo
param Pd{i in I}; #Demanda por nodo
param FL{i in L,j in 1..2}; #Flujos de linea (Min/Max)
param KL{i in I,j in 1..2}; #Limites de perturbacion (Min/Max)
param FD{i in I, j in I}; #Numeración de linea para flujo
param D{i in L, j in 1..3}; #Direccion de flujos
param Bp{(i,j) in {I,I}}; #Ybus modificada 
param M{i in I, j in L};
var teta{i in I}; #Angulos de teta
var P{i in I}>= 0; #Potencia por nodo
var F{i in L}; #Flujo de cada linea
var K{i in I}>= 0; #Perturbación por nodo

var K2{i in I}>= 0; #Perturbación por nodo
var muflmin{i in L}>= 0; #multiplicador de lagrange inferior de flujo de linea
var muflmax{i in L}>= 0; #multiplicador de lagrange superior de flujo de linea
var mupmin{i in I}>= 0; #multiplicador de lagrange inferior de potencia
var mupmax{i in I}>= 0; #multiplicador de lagrange superior de potencia
var mukmin{i in I}>= 0; #multiplicador de lagrange inferior de k
var mukmax{i in I}>= 0; #multiplicador de lagrange superior de k
var lambda{i in I}; #multiplicador de lagrange de balance de potencia


param PML{i in I};
param Ri{i in I};
param BGENCO{i in 1..3};

param KO{i in 1..221,j in I}; #Tabla que contiene el porcentaje a evaluar
param PMLT{i in 1..221, j in I}; #Tabla que almacena el PML
param Pot{i in 1..221, j in I}; #Tabla que almacena la potencia
param Flujo{i in 1..221, j in L}; #Tabla que almacena los flujos de lineas
param BEN{i in 1..221, j in I}; #Tabla de beneficio por nodo

param K3{i in I}; #Demanda por nodo

#param K{i in I};

#########################################
#FUNCION OBJETIVO
#minimize f: sum{i in I} K[i]*(2*C[i,1]*P[i]^2+C[i,2]*P[i]);
maximize f: sum{i in I} ((lambda[i]*P[i])-(CO[i,1]*P[i]^2+CO[i,2]*P[i]+CO[i,3]));

s.t. restriccionK{i in I}: KL[i,1]<=K2[i]<=KL[i,2];   ###

#########################################
#RESTRICCIONES DEL PROBLEMA DE MAXIMIZACIÓN


#########################################
#CONDICIONES ESTACIONARIAS
#LAGRANGIANO CON RESPECTO A K
s.t. REK{i in I}: (2*CO[i,1]*P[i]^2+CO[i,2]*P[i])-(mukmin[i])+(mukmax[i])=0;          ###
#LAGRANGIANO CON RESPECTO A LA POTENCIA
#s.t. REP4: K2[1]*(4*CO[1,1]*P[1]+CO[1,2])+(lambda[1])-(mupmin[1])+(mupmax[1])=0;    ###
s.t. REP{i in I}: K2[i]*(4*CO[i,1]*P[i]+CO[i,2])+(lambda[i])-(mupmin[i])+(mupmax[i])=0;  ###  #{i in I: i <> 4}
#LAGRANGIANO CON RESPECTO A TETHA
s.t. RET{i in I}: sum{j in I}(lambda[i]*B[i,j])
	   			-sum{j in I}(lambda[j]*B[i,j])
	   			+sum{j in L}(M[i,j]*(muflmax[j]-muflmin[j]))=0;
	   			
	   			#(lambda[i])*(sum{j in I: j<>i}B[i,j])+(-B[D[i,1],D[i,2]])*(muflmax[i]-muflmin[i])=0;

#########################################
#CONDICIONES DE HOLGURA
s.t. RHFMIN{k in L}: muflmin[k]*(FL[k,1]-(F[D[k,3]]))=0; #FLUJOMIN
s.t. RHFMAX{k in L}: muflmax[k]*((F[D[k,3]]-FL[k,2]))=0;

#s.t. RHFMIN{k in L}: muflmin[k]*(F[D[k,3]]+FL[k,2])=0;
#s.t. RHFMAX{k in L}: muflmax[k]*(F[D[k,3]]+FL[k,1])=0; #FLUJOMAX

s.t. RHPMIN{i in I}: mupmin[i]*(Pg[i,1]-P[i])=0; #POTMMIN
s.t. RHPMAX{i in I}: mupmax[i]*(P[i]-Pg[i,2])=0; #POTMAX

#s.t. RHPMIN{i in I}: mupmin[i]*(Pg[i,1]-P[i])=0; #POTMMIN
#s.t. RHPMAX{i in I}: mupmax[i]*(P[i]-Pg[i,2])=0; #POTMAX

s.t. RHKMIN{i in I}: mukmin[i]*(KL[i,1]-K[i])=0; #KLMMIN  ###
s.t. RHKMAX{i in I}: mukmax[i]*(K[i]-K2[i])=0; #KLMAX     ### 	#########################################################

#########################################
#CONDICIONES DE FACTIBILIDAD
#BALANCE DE POTENCIA
s.t. R1{i in I}: sum{j in I: i<>j} B[i,j]*(teta[i]-teta[j])+P[i]=Pd[i]; 


#FLUJO MIN/MAX DE LÍNEA
#s.t. R2{k in L}: FL[k,1]<=B[D[k,1],D[k,2]]*(teta[D[k,1]]-teta[D[k,2]])<= FL[k,2];
s.t. R21{k in L}: FL[k,1]<=B[D[k,1],D[k,2]]*(teta[D[k,1]]-teta[D[k,2]]);
s.t. R22{k in L}: B[D[k,1],D[k,2]]*(teta[D[k,1]]-teta[D[k,2]])<= FL[k,2];

#POTENCIA MAXIMA DE GENERACION
#s.t. R3{i in I}: Pg[i,1]<=P[i]<=Pg[i,2];
s.t. R31{i in I}: Pg[i,1]<=P[i];
s.t. R32{i in I}: P[i]<=Pg[i,2];

#FLUJOS DE LINEA
s.t. R4{i in L}: F[D[i,3]]=B[D[i,1],D[i,2]]*(teta[D[i,1]]-teta[D[i,2]]);

#ANGULO DE REFERENCIA
s.t. R5: teta[1]=0;

#K FIJAS
#s.t. RK1: K[1]=1.0;   #G1
#s.t. RK2: K[2]=1.0;
#s.t. RK3: K[3]=1.0;   #G2
#s.t. RK4: K[4]=1.0;	  #G3
#s.t. RK5: K[5]=1.0;   #G4
#s.t. RK6: K[6]=1.0;   #G5
#s.t. RK7: K[7]=1.0;   #G6
#s.t. RK8: K[8]=1.0;   


s.t. RKT{i in I}: K2[i]=K3[i];

#s.t. R2K1: K2[1]=0.80;   #G1
#s.t. R2K2: K2[2]=1.0;
#s.t. R2K3: K2[3]=1.0;   #G2
#s.t. R2K4: K2[4]=1.0;	#G3
#s.t. R2K5: K2[5]=1.0;   #G4
#s.t. R2K6: K2[6]=1.0;   #G5
#s.t. R2K7: K2[7]=1.0;   #G6
#s.t. R2K8: K2[8]=1.0;  

#LIMITES DE PERTURBACIÓN
s.t. RKI{i in I}: KL[i,1]<=K[i];    ###
s.t. RKS{i in I}: K[i]<=K2[i];      ###