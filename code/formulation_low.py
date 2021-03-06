##author Xinwu Qian
##Generate the input for flow* automatically
import numpy as np

fw=open('setting.txt','w')
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


fw.writelines('----------------Variables------------------'+'\n')
#State var x0,x1,x2,.....
fw.writelines('state var '+",".join(var)+',t'+'\n')

#list of variables
sigma=0.1 #latent time 
gamma=0.1 #recover time
betainn=0.1 #contact rate inn zone
n_mode=3; #number of travel modes
betatra=[0.15,0.1,0.001] # within traffic contact rate
g=[0.3,0.3,0.3] #departure rate
Nr=[1,1,1];
Np=[1,1,1]; #residents per zone
N=[[0.5,0.3,0.2],[0.3,0.4,0.3],[0.2,0.3,0.5]]
d=[1,1,1] #control level for mode i 
t=[0.2,0.3,0.3]
c=[0.4,0.3,0.3] #ratio of choosing mode 1, 2,3. mode 1: large capacity, mode 2: medium capactiy, mode 3: small capacity
#calculate the travel ratio
for i in range(n_mode):
	t[i]=t[i]*1.05**(30*(1-d[i]))

for i in range(n_mode):
	c[i]=np.exp(-t[i])/sum([np.exp(-t[j]) for j in range(n_mode)])


m=[[0,0.6,0.4],[0.5,0,0.5],[0.4,0.6,0]] #movement ratio from patch i to patch j
r=[[0,0.6,0.3],[0.6,0,0.3],[0.3,0.6,0]] #return ratio from patch i to patch j

for i in range(n_zone):
	for j in range(n_zone):
		if i==j:
			N[i][j]=1/(1+g[i]*sum([m[i][k]/r[i][k] for k in range(n_zone) if k!=i]))*Nr[i]

for i in range(n_zone):
	for j in range(n_zone):
		if i!=j:
			N[i][j]=g[i]*m[i][j]/r[i][j]*N[i][i]
print N

for i in range(n_zone):
	Np[i]=sum([N[j][i] for j in range(n_zone)])



#load the incident matrix
f_i=open('1.csv')
path_matrix=np.zeros((n_zone**2,n_zone**2))
row=0;
for lines in f_i:
	line=lines.strip().split(',')
	line=np.array([int(i) for i in line]);
	path_matrix[row,:]=line
	row+=1


fw.writelines('----------------ODE formulations------------------'+'\n')

#start writing equations for S
for i in range(n_zone):
	for j in range(n_zone):
		exp=var[S[i][j]]+"'=";
		#Same zone
		if i==j:
			exp+=var[S[i][j]]+'*'+str(-g[i])+'+'+var[S[i][j]]+'*'+var[I[i][j]]+'*'+str(-betainn/Np[i])
			for k in range(n_zone):
				if(k!=j):
					exp+='+'+var[S[i][k]]+'*'+str(r[i][k])
					coef=sum([betainn/Np[i]*c[l]*d[l] for l in range(n_mode)])
					exp+='-'+var[S[i][j]]+'*'+var[I[k][i]]+'*'+str(coef)
		else:
			exp+=var[S[i][i]]+'*'+str(g[i]*m[i][j])+'-'+var[S[i][j]]+'*'+str(r[i][j])+'-'+var[S[i][j]]+'*'+var[I[j][j]]+'*'+str(betainn/Np[j])
			#destination contagion
			for k in range(n_zone):
				if(k!=j):
					coef=sum([betainn/Np[j]*c[l]*d[l] for l in range(n_mode)])
					exp+='-'+var[S[i][j]]+'*'+var[I[k][j]]+'*'+str(coef)
			#calculate within mode contagion
			for l in range(n_mode):
				denom=sum([c[l]*path_matrix[S[i][j],S[k][u]]*N[k][u] for k in range(n_zone) for u in range(n_zone)])
				for k in range(n_zone):
					for u in range(n_zone):
						if path_matrix[S[i][j],S[k][u]]>0:
							exp+='-'+var[S[i][j]]+'*'+var[I[k][u]]+'*'+str(betatra[l]*c[l]*c[l]*d[l]/denom)
		#write the expression to the fire
		fw.writelines(exp+'\n')

#start writing equations for E
for i in range(n_zone):
	for j in range(n_zone):
		exp=var[E[i][j]]+"'="
		#same zone
		if(i==j):
			exp+=var[E[i][j]]+'*'+str(-g[i])+'+'+var[S[i][j]]+'*'+var[I[i][j]]+'*'+str(betainn/Np[i])+'+'+var[E[i][j]]+'*'+str(-sigma)
			for k in range(n_zone):
				if(k!=j):
					exp+='+'+var[E[i][k]]+'*'+str(r[i][k])	
					coef=sum([betainn/Np[i]*c[l]*d[l] for l in range(n_mode)])
					exp+='+'+var[S[i][j]]+'*'+var[I[k][i]]+'*'+str(coef)
		else:
			#departure-return
			exp+=var[E[i][i]]+'*'+str(g[i]*m[i][j])+'-'+var[E[i][j]]+'*'+str(r[i][j])+'+'+var[E[i][j]]+'*'+str(-sigma)+'+'+var[S[i][j]]+'*'+var[I[j][j]]+'*'+str(betainn/Np[j])
			for k in range(n_zone):
				if(k!=j):
					coef=sum([betainn/Np[j]*c[l]*d[l] for l in range(n_mode)])
					exp+='+'+var[S[i][j]]+'*'+var[I[k][j]]+'*'+str(coef)
			#calculate within mode contagion
			for l in range(n_mode):
				denom=sum([c[l]*path_matrix[S[i][j],S[k][u]]*N[k][u] for k in range(n_zone) for u in range(n_zone)])
				for k in range(n_zone):
					for u in range(n_zone):
						if path_matrix[S[i][j],S[k][u]]>0:
							exp+='+'+var[S[i][j]]+'*'+var[I[k][u]]+'*'+str(betatra[l]*c[l]*c[l]*d[l]/denom)
		#write the expression to the fire
		fw.writelines(exp+'\n')

#start writing equations for I
for i in range(n_zone):
	for j in range(n_zone):
		exp=var[I[i][j]]+"'="
		#same zone
		if i==j:
			# - depart + latent - recover
			exp+=var[I[i][j]]+'*'+str(-g[i])+'+'+var[E[i][j]]+'*'+str(sigma)+'+'+var[I[i][j]]+'*'+str(-gamma)
			for k in range(n_zone):
				if(k!=i):
					coef=sum([r[i][k]*d[l]*c[l] for l in range(n_mode)])
					exp+='+'+var[I[i][k]]+'*'+str(coef)
		else:
			# - return + latent
			exp+=var[I[i][j]]+'*'+str(-r[i][j])+'+'+var[E[i][j]]+'*'+str(sigma)
			coef1=sum([g[i]*m[i][j]*c[l]*d[l] for l in range(n_mode)])
			coef2=gamma*sum([c[l]*d[l] for l in range(n_mode)])
			coef3=sum([c[l]*(1-d[l]) for l in range(n_mode)])
			exp+='+'+var[I[i][i]]+'*'+str(coef1)+'+'+var[I[i][j]]+'*'+str(-coef2)+'+'+var[I[i][j]]+'*'+str(-coef3)
		fw.writelines(exp+'\n')

#start writing equations for R
for i in range(n_zone):
	for j in range(n_zone):
		exp=var[R[i][j]]+"'="
		#same zone
		if i==j:
			# -departure + recover
			exp+=var[R[i][i]]+'*'+str(-g[i])+'+'+var[I[i][i]]+'*'+str(gamma)
			for k in range(n_zone):
				if(k!=i):
					exp+='+'+var[R[i][k]]+'*'+str(r[i][k])
					coef=sum([r[i][k]*(1-d[l])*c[l] for l in range(n_mode)])
					exp+='+'+var[I[i][k]]+'*'+str(coef)
		else:
			# - return + departure
			exp+=var[R[i][j]]+'*'+str(-r[i][j])+'+'+var[R[i][i]]+'*'+str(g[i]*m[i][j])
			coef1=gamma*sum([c[l]*d[l] for l in range(n_mode)])
			coef2=sum([c[l]*(1-d[l]) for l in range(n_mode)])
			exp+='+'+var[I[i][j]]+'*'+str(coef1)+'+'+var[I[i][j]]+'*'+str(coef2)
		fw.writelines(exp+'\n')

fw.writelines('----------------Invariant intervals------------------'+'\n')
for item in var:
	fw.writelines(item+' in [0,1]'+'\n')

fw.writelines('----------------Initial conditions------------------'+'\n')

#S init
for i in range(n_zone):
	for j in range(n_zone):
		fw.writelines(var[S[i][j]]+' in [0.8,0.8]'+'\n')
#E init
for i in range(n_zone):
	for j in range(n_zone):
		fw.writelines(var[E[i][j]]+' in [0,0]'+'\n')
#I init
for i in range(n_zone):
	for j in range(n_zone):
		fw.writelines(var[I[i][j]]+' in [0.2,0.2]'+'\n')
#R init
for i in range(n_zone):
	for j in range(n_zone):
		fw.writelines(var[R[i][j]]+' in [0,0]'+'\n')

fw.close()
