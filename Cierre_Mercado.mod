option solve cplex;

#DECLARACION DE DATOS CONOCIDOS
param Nnodos; #Numero de nodos
param Nlines; #Numero de lineas
param Npor; #Numero de porcnetaje de k

set I := {1..Nnodos};
set L := {1..Nlines};

param CO{i in I, j in 1..3}; #Costos por nodo
param B{i in I, j in I}; #yBUS GENERADA POR OTRO SOFTWARE
param Pg{i in I, j in 1..2}; #Potencia (Min/Max) por nodo
param Pd{i in I}; #Demanda por nodo
param FL{i in L,j in 1..2}; #Flujos de linea (Min/Max)
param D{i in L, j in 1..3}; #Direccion de flujos 
param K{i in I}; #Valor de k obtenido de K1
param KO{i in 1..21,j in I}; #Tabla que contiene el porcentaje a evaluar
param PML{i in 1..21, j in I}; #Tabla que almacena el PML
param Flujo{i in 1..21, j in L}; #Tabla que almacena los flujos de lineas
param Pot{i in 1..21, j in I}; #Tabla que almacena la potencia

#DECLARACION DE VARIABLES
var teta{i in I}; #Angulos de teta
var P{i in I}; #Potencia por nodo
var F{i in L}; #Flujo de cada linea

#INICIA ALGORITMO

#FUNCION OBJETIVO
minimize f: sum{i in I} (K[i]*(2*CO[i,1]*P[i]^2+CO[i,2]*P[i]));

#RESTRICCIONES

#BALANCE DE POTENCIA
s.t. R1{i in I}: sum{j in I: i<>j} B[i,j]*(teta[i]-teta[j])+P[i]=Pd[i];

#FLUJO MIN/MAX DE LÍNEA
s.t. R2{k in L}: FL[k,1]<=B[D[k,1],D[k,2]]*(teta[D[k,1]]-teta[D[k,2]])<= FL[k,2];

#POTENCIA MAXIMA DE GENERACION
s.t. R3{i in I}: Pg[i,1]<=P[i]<=Pg[i,2];

#FLUJOS DE LINEA
s.t. R4{i in L}: F[D[i,3]]=B[D[i,1],D[i,2]]*(-1)*(teta[D[i,1]]-teta[D[i,2]]);

#ANGULO DE REFERENCIA
s.t. R5: teta[1]=0;
