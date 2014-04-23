# coding: UTF-8
from pylab import *

from numpy import *

import scipy as sci
import scipy.linalg
import scipy.sparse as sp

import scipy.sparse.linalg as spla
from scipy.sparse import lil_matrix,hstack,vstack
from numpy.random import rand

from numpy import genfromtxt


#データの読み込み
datanamebox=[["A_NBI11-10153A.csv","A_NBI11-10200A.csv","A_NBI11-10200B.csv","A_NBI11-10230.csv"],
			["B_NBI11-10153B.csv","B_NBI11-10155.csv","B_NBI11-10161.csv","B_NBI11-10179.csv","B_NBI11-10201B.csv"],
			["C_NBI11_lla+llc_1.csv","C_NBI11_lla+llc_2.csv","C_NBI11-10244.csv"]]
dataname=datanamebox[0][0]
print dataname
data=genfromtxt("data/"+dataname, delimiter=",")

data = data[1:,:]

# data = log(data)
Adata=data[:,0]
Bdata=data[:,1]
C3data=data[:,2]



data = log(data)


p=0.8
L=3
K=3
I=200
J=1


f = ones([I-1,L,K])*log((1-p)/2)
for i in range(0, I-1):
		for l in range(0,L):
			for k in range(0,K):
				if l==k:
					f[i][l][k]=log(p)


dil=data.flatten()

fijlm= f.flatten()


c = -hstack((dil,fijlm))


c=lil_matrix(c).T



m=I+(I-1)*L+(I-1)*(L-1)
n=I*L+(I-1)*L*K
print m
print n


A = lil_matrix((m,n))

for i in range(0, I):
		for l in range(0,L):
			A[i,i*L+l]=1.0
			if(i<I-1):
				A[I+i*L+l,i*L+l]=-1.0
				if(l<L-1):
					A[I+(I-1)*L+i*(L-1)+l,i*L+l+L]=-1.0
				for k in range(0,K):
					A[I+i*L+l,I*L+i*L*K+l*K+k]=1.0
					if(k<K-1):
						A[I+(I-1)*L+i*(L-1)+k,I*L+i*L*K+l*K+k]=1.0


savetxt("a.csv", A.todense(), fmt="%.10f",delimiter=",")

b = hstack((ones(I),zeros((I-1)*L+(I-1)*(L-1)) ))
b = lil_matrix(b).T


# xx=zeros((I-1)*L*K)
# for i in range(0, I-1):
# 	xx[L*K*i]=1
# xxx=zeros(I*L)
# for i in range(0,I):
# 	xxx[L*i]=1
# x = hstack((xxx,xx))

# x=lil_matrix(x)
# x = lil_matrix(ones(n))
# y = lil_matrix(ones(m))
# z = lil_matrix(ones(n))

x = lil_matrix(rand(n))
y = lil_matrix(rand(m))
z = lil_matrix(rand(n))

chart=x.todense()

chart=chart.reshape((I+(I-1)*K,L)).T

chartA = lil_matrix(chart[0,:I])
chartB = lil_matrix(chart[1,:I])
chartC3 = lil_matrix(chart[2,:I])




X = lil_matrix((n,n))
Z = lil_matrix((n,n))

deltaX1 = lil_matrix((n,n))
deltaZ1 = lil_matrix((n,n))

#初期値やパラメータ，定数行列の設定

E = 0.00000001
Ep = 0.000000001
Ed = 0.000000001


x=x.T
y=y.tocsr().T
z=z.T

print x.T*z

i=0
AA=lil_matrix((m,n+m))
A00 = hstack([A,AA])


ATI = hstack([lil_matrix((n,n)),A.T,lil_matrix(eye(n))])

print "start"

A=A.tocsr()
b=b.tocsr()
c=c.tocsr()



#内点法@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
while x.T*z>E or (A*x-b).multiply(A*x-b).sum()>Ep or (A.T*y+z-c).multiply(A.T*y+z-c).sum()>Ed :
	X.setdiag(x.todense())
	Z.setdiag(z.todense())


	Z0X = hstack([Z,lil_matrix(zeros((n,m))),X])

	A_X = vstack([A00,ATI,Z0X])


	A_Xz = vstack((A*x-b,A.T*y+z-c,X*z))

	i+=1
	print i

	A_X = A_X.tocsr()


	delta1 = -spla.spsolve(A_X,A_Xz)
	

	deltax1=delta1[:n]
	deltay1=delta1[n:n+m]
	deltaz1=delta1[n+m:]



	deltax1 = lil_matrix(deltax1).T
	deltay1 = lil_matrix(deltay1).T
	deltaz1 = lil_matrix(deltaz1).T

	ap1=float("inf")

	ap1x=-(x/deltax1).todense()
	for j in range(0, n):
	  	if 0<ap1x[j][0]<ap1 :
			ap1=ap1x[j][0]


	ad1=float("inf")
	ad1z=-(z/deltaz1).todense()
	for j in range(0, n):
	  	if 0<ad1z[j][0]<ad1 :
			ad1=ad1z[j][0]


	gamma=((x+deltax1*ap1).T*(z+deltaz1*ad1)/(x.T*z).todense())**3



	myu = (x.T*z).todense()/n*gamma



	deltaX1.setdiag(deltax1.todense())


	DeltaXz = lil_matrix(deltaX1*deltaz1-lil_matrix(ones(n)).T*myu)
	


	zero_X = vstack([lil_matrix((n+m,1)),DeltaXz])


	delta2 = -spla.spsolve(A_X,zero_X)

	delta = delta1+delta2



	deltax = delta[:n]
	deltay = delta[n:n+m]
	deltaz = delta[n+m:]

	deltax = lil_matrix(deltax).T
	deltay = lil_matrix(deltay).T
	deltaz = lil_matrix(deltaz).T

	ap2=float("inf")
	ap2x=-(x/deltax).todense()
	for j in range(0, n):
	  	if 0<ap2x[j][0]<ap2:
			ap2=ap2x[j][0]

	if ap2==inf:
		ap2=lil_matrix((1,1))

	ad2=float("inf")
	ad2z=-(z/deltaz).todense()
	for j in range(0, n):
	  	if 0<ad2z[j][0]<ad2 :
			ad2=ad2z[j][0]

	if ad2==inf:
		ad2=lil_matrix((1,1))


	ap = ap2[0,0]*0.99 if ap2[0,0]*0.99<1  else 1.0

	ad = ad2[0,0]*0.99 if ad2[0,0]*0.99<1  else 1.0

	x = x+deltax*ap
	y = y+deltay*ad
	z = z+deltaz*ad


	print x.T*z,

	print (A*x-b).multiply(A*x-b).sum() ,

	print (A.T*y+z-c).multiply(A.T*y+z-c).sum()

	chart=x.todense()

	chart=chart.reshape((I+(I-1)*K,L)).T

	chartA=vstack([chartA,chart[0,:I]])
	chartB=vstack([chartB,chart[1,:I]])
	chartC3=vstack([chartC3,chart[2,:I]])
print "end"


chartA=asarray(chartA.todense())
chartB=asarray(chartB.todense())
chartC3=asarray(chartC3.todense())

x=x.todense()

x=x.reshape((I+(I-1)*K,L))


savetxt("x.csv", x, fmt="%.30f",delimiter=",")

#プロット@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
listx=range(0,I)


# 棒グラフのy軸のラベルを非表示
for j in range(0,i):

	figure(num=None, figsize=(12, 4), dpi=40, facecolor='w', edgecolor='k')

	axes((0.05, 0.06, 0.9, 0.9))

	plot(listx, chartA[j],'b',lw=2)
	plot(listx, chartB[j],'r',lw=2)
	plot(listx, chartC3[j],'g',lw=2)

	axis([0, I, -0.1, 2])

	savefig("result/"+dataname+"_p="+str(p)+"_"+str(j)+".png")

	

close('all')

#プロット@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

x=asarray(x)

listy1=x[:I,0].T
listy2=x[:I,1].T
listy3=x[:I,2].T


#データの読み込み
data2=genfromtxt("data/estimated_label/result_"+str(p)+"_"+dataname, delimiter=",")


Adata2=zeros(I)

Bdata2=zeros(I)

C3data2=zeros(I)

for i in range(0, I):
	if data2[i]==1:
		Adata2[i]=1
	elif data2[i]==2:
		Bdata2[i]=1
	elif data2[i]==3:
		C3data2[i]=1
	else:
		print "data error"

width = 1

x = np.random.randn(I)
y = np.sin(x) + np.random.randn(I) * 0.3

fig = plt.figure(num=None, figsize=(12, 8), dpi=40, facecolor='w', edgecolor='k')

# サブプロットで分割
spase=0.05
height=(1-spase*4)/3
ax1 = fig.add_axes((0.05, height*2+spase*3, 0.93, height))
ax2 = fig.add_axes((0.05,  height+spase*2, 0.93, height), sharex=ax1)
ax3 = fig.add_axes((0.05, spase, 0.93, height), sharex=ax1)

# 棒グラフのy軸のラベルを非表示
ax2.tick_params(labelleft="off")

width=1
ax1.plot(listx, Adata,'b',lw=2)
ax1.plot(listx, Bdata,'r',lw=2)
ax1.plot(listx, C3data,'g',lw=2)

ax2.bar(listx, Adata2, width, color='b', edgecolor='b')
ax2.bar(listx, Bdata2, width, color='r', edgecolor='r')
ax2.bar(listx, C3data2, width, color='g', edgecolor='g')

ax3.axis([0, I, -0.1, 1.1])

ax3.plot(listx,listy1,'b',lw=3)
ax3.plot(listx,listy2,'r',lw=3)
ax3.plot(listx,listy3,'g',lw=3)

legend(('A', 'B','C3'), 'lower right')

fig.savefig("result/"+dataname+"_p="+str(p)+'.png')

show()