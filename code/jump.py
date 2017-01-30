#write the jump condition

import numpy as np

fw=open('jump.txt','w')
n_zone=3; #number of nodes in the network
S=np.zeros((n_zone,n_zone),dtype=int);
E=np.zeros((n_zone,n_zone),dtype=int);
I=np.zeros((n_zone,n_zone),dtype=int);
R=np.zeros((n_zone,n_zone),dtype=int); #S,E,I,R initialization
index=0;
V=[S,E,I,R]
for i in range(len(V)):
	for m in range(n_zone):
		for n in range(n_zone):
			V[i][m,n]=index
			index+=1
var=[]
for k in range(len(V)):
	var+=['x'+str(int(V[k][i,j])+1) for i in range(n_zone) for j in range(n_zone)];

exp=" "
for i in range(n_zone):
	for j in range(n_zone):
		exp+=var[I[i][j]]+"+"
exp+="0<=0.6\n" 
fw.writelines(exp)

exp=" "
for i in range(len(var)):
	exp+=var[i]+"':="+var[i]+" "
exp+="\n"
fw.writelines(exp)
fw.close()