# coding: UTF-8
from pylab import *

from numpy import *

import scipy.ndimage as sn

import scipy as sci
import scipy.linalg
import scipy.sparse as sp

import scipy.sparse.linalg as spla
from numpy.random import rand

from numpy import genfromtxt
import time
import Image


imgL=Image.open("data/tsukuba-imL.png")
imgR=Image.open("data/tsukuba-imR.png")

pixL=imgL.load()
pixR=imgR.load()

nColors=len(imgL.mode)

nD=16	#ディスパリティの数＝ラベル数

T,S= imgL.size

depth = Image.new('L', (T, S))

depthpix = depth.load()

dsi=ones([S,T,nD])*inf

for s in range(0,S):
	for t in range(0,T):
		for d in range(0,nD):
			t2=t-d
			if t2>=0:
				sumdiff=0
				for b in range(0,nColors):
					ILm = pixL[t,s][b] if t == 0 else (pixL[t,s][b]+pixL[t-1,s][b])/2
					ILp = pixL[t,s][b] if t == T-1 else (pixL[t,s][b]+pixL[t+1,s][b])/2
					IRm = pixR[t2,s][b] if t2 == 0 else (pixR[t2,s][b]+pixR[t2-1,s][b])/2
					IRp = pixR[t2,s][b] if t2 == T-1 else (pixR[t2,s][b]+pixR[t2+1,s][b])/2
					ILmin=min(ILm,ILp,pixL[t,s][b])
					ILmax=max(ILm,ILp,pixL[t,s][b])
					IRmin=min(IRm,IRp,pixR[t2,s][b])
					IRmax=max(IRm,IRp,pixR[t2,s][b])
					dL=max(0,pixL[t,s][b]-IRmax,IRmin-pixL[t,s][b])
					dR=max(0,pixR[t2,s][b]-ILmax,ILmin-pixR[t2,s][b])
					diff = min(dL,dR)
					sumdiff+=diff*diff
				dsi[s,t,d]=sumdiff
print "fin dsi"

dsi=dsi.flatten()

gradThresh=24
gradPenalty=2.0

Cst=ones([S*(T-1)*nD*nD+(S-1)*T*nD*nD])*2
for s in range(0,S):
	for t in range(0,T):
		sx=sy=0
		for b in range(0,nColors):
			dx = pixL[t,s][b] - pixL[t+1,s][b] if t<T-1 else 1.0
			dy = pixL[t,s][b] - pixL[t,s+1][b] if s<S-1 else 1.0 
			sx += dx*dx
			sy += dy*dy
		wx = gradPenalty if sx < gradThresh else 1.0
		wy = gradPenalty if sy < gradThresh else 1.0
		T_width=t*nD*nD+s*(T-1)*nD*nD
		Twidth=t*nD*nD+s*T*nD*nD
		S_width=S*(T-1)*nD*nD
		for l in range(0,nD):
			for m in range(0,nD):
				if fabs(l-m)<2:
					if t<T-1:
						Cst[m+l*nD+T_width] *= fabs(l-m)*wx
					if s<S-1:
						Cst[m+l*nD+Twidth+S_width] *= fabs(l-m)*wy
				else:
					if t<T-1:
						Cst[m+l*nD+T_width] *= wx
					if s<S-1:
						Cst[m+l*nD+Twidth+S_width] *= wy
print "fin 平滑化項"

c = hstack((dsi,Cst))

print "fin c hstack"

c=sp.lil_matrix(c).T


L=nD

h1=S*T
h2=h1+S*(T-1)*L
h3=h2+(S-1)*T*L
h4=h3+S*(T-1)*(L-1)
h5=h4+(S-1)*T*(L-1)

w1=S*T*L
w2=w1+S*(T-1)*L*L
w3=w2+(S-1)*T*L*L

num1=S*T*L
num2=num1+S*(T-1)*L
num_2=num2+S*(T-1)*L*L
num3=num_2+(S-1)*T*L
num_3=num3+(S-1)*T*L*L
num4=num_3+S*(T-1)*(L-1)
num_4=num4+S*(T-1)*L*(L-1)
num5=num_4+(S-1)*T*(L-1)
num_5=num5+(S-1)*T*L*(L-1)


data=ones(num_5)
data[:num1]=-data[:num1]
data[num1:num_5:L+1]=-data[num1:num_5:L+1]
indices=zeros(num_5,int16)
indptr= zeros(h5+1,int16)
indptr[:h1]=arange(0,num1,L)
indptr[h1:h5+1]=arange(num1,num_5+L+1,L+1)

ia=ib=ic=idd=ie=0

for s in range(0,S):
	print s
	for t in range(0,T):
		for l in range(0,L):
			indices[ia]=ia
			ia+=1
			ltLsTL=l+t*L+s*T*L
			if t<T-1:
				indices[ib+num1]=ltLsTL
				ib+=1
				xb=l*L+t*L*L+s*(T-1)*L*L+w1
				for m in range(0,L):
					indices[ib+num1]=m+xb
					ib+=1
			if s<S-1:
				indices[ic+num_2]=ltLsTL
				ic+=1
				xc=ltLsTL*L+w2
				for m in range(0,L):
					indices[ic+num_2]=m+xc
					ic+=1
			if l<L-1:
				if t<T-1:
					indices[idd+num_3]=ltLsTL+L
					idd+=1
					xd=l+t*L*L+s*(T-1)*L*L+w1
					for m in range(0,L):
						indices[idd+num_3]=m*L+xd
						idd+=1
				if s<S-1:
					indices[ie+num_4]=ltLsTL+T*L 
					ie+=1
					xe=l+t*L*L+s*T*L*L+w2
					for m in range(0,L):
						indices[ie+num_4]=m*L+xe
						ie+=1

print "fin A"


m=h1
n=w3


b = hstack((ones(h1),zeros(h5-h1)))
b = sp.lil_matrix(b).T





x = sp.lil_matrix(ones(n))
y = sp.lil_matrix(ones(m))
z = sp.lil_matrix(ones(n))




X = sp.lil_matrix((n,n))
Z = sp.lil_matrix((n,n))

deltaX1 = sp.lil_matrix((n,n))
deltaZ1 = sp.lil_matrix((n,n))

#初期値やパラメータ，定数行列の設定

E = 0.00000001
Ep = 0.000000001
Ed = 0.000000001


x=x.T
y=y.tocsr().T
z=z.T

print x.T*z

i=0
AA=sp.lil_matrix((m,n+m))
A00 = hstack([A,AA])


ATI = hstack([sp.lil_matrix((n,n)),A.T,sp.lil_matrix(eye(n))])

print "start"

A=A.tocsr()
b=b.tocsr()
c=c.tocsr()



#内点法@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
while x.T*z>E or (A*x-b).multiply(A*x-b).sum()>Ep or (A.T*y+z-c).multiply(A.T*y+z-c).sum()>Ed :
	X.setdiag(x.todense())
	Z.setdiag(z.todense())


	Z0X = hstack([Z,sp.lil_matrix(zeros((n,m))),X])

	A_X = vstack([A00,ATI,Z0X])


	A_Xz = vstack((A*x-b,A.T*y+z-c,X*z))

	i+=1
	print i

	A_X = A_X.tocsr()


	delta1 = -spla.spsolve(A_X,A_Xz)
	

	deltax1=delta1[:n]
	deltay1=delta1[n:n+m]
	deltaz1=delta1[n+m:]



	deltax1 = sp.lil_matrix(deltax1).T
	deltay1 = sp.lil_matrix(deltay1).T
	deltaz1 = sp.lil_matrix(deltaz1).T

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


	DeltaXz = sp.lil_matrix(deltaX1*deltaz1-sp.lil_matrix(ones(n)).T*myu)
	


	zero_X = vstack([sp.lil_matrix((n+m,1)),DeltaXz])


	delta2 = -spla.spsolve(A_X,zero_X)

	delta = delta1+delta2



	deltax = delta[:n]
	deltay = delta[n:n+m]
	deltaz = delta[n+m:]

	deltax = sp.lil_matrix(deltax).T
	deltay = sp.lil_matrix(deltay).T
	deltaz = sp.lil_matrix(deltaz).T

	ap2=float("inf")
	ap2x=-(x/deltax).todense()
	for j in range(0, n):
	  	if 0<ap2x[j][0]<ap2:
			ap2=ap2x[j][0]

	if ap2==inf:
		ap2=sp.lil_matrix((1,1))

	ad2=float("inf")
	ad2z=-(z/deltaz).todense()
	for j in range(0, n):
	  	if 0<ad2z[j][0]<ad2 :
			ad2=ad2z[j][0]

	if ad2==inf:
		ad2=sp.lil_matrix((1,1))


	ap = ap2[0,0]*0.99 if ap2[0,0]*0.99<1  else 1.0

	ad = ad2[0,0]*0.99 if ad2[0,0]*0.99<1  else 1.0

	x = x+deltax*ap
	y = y+deltay*ad
	z = z+deltaz*ad


	print x.T*z,

	print (A*x-b).multiply(A*x-b).sum() ,

	print (A.T*y+z-c).multiply(A.T*y+z-c).sum()


print "end"
