# coding: UTF-8

from pylab import *

from numpy import *

import scipy as sci
import scipy.linalg
import scipy.sparse as sp

import scipy.sparse.linalg as spla
from scipy.sparse import lil_matrix,hstack,vstack
from numpy.random import rand



#データの読み込み
data = loadtxt("data.csv", delimiter=',')


data = data[:,1:]

#data = data[:4]

data = log(data)


p=0.99
L=3
K=3
I=201
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
x = lil_matrix(ones(n))
y = lil_matrix(ones(m))
z = lil_matrix(ones(n))

# x = lil_matrix(rand(n))
# y = lil_matrix(rand(m))
# z = lil_matrix(rand(n))

X = lil_matrix((n,n))
Z = lil_matrix((n,n))

#初期値やパラメータ，定数行列の設定
ganma = 0.95
E = 0.001
Ep = 0.001
Ed = 0.001




# listx = [x[0,0]]
# listy = [x[1,0]]


#配列から行列へ
#A = mat(A)
# b = mat(b).T
# c = mat(c).T

# x = mat(x).T
# y = mat(y).T
# z = mat(z).T

x=x.tocsr().T
y=y.tocsr().T
z=z.tocsr().T

print x.T*z

i=0
AA=lil_matrix((m,n+m))
A00 = hstack([A,AA])


ATI = hstack([lil_matrix((n,n)),A.T,lil_matrix(eye(n))])

print "start"

A=A.tocsr()
b=b.tocsr()
c=c.tocsr()

#内点法
while x.T*z>E or (A*x-b).multiply(A*x-b).sum()>Ep or (A.T*y+z-c).multiply(A.T*y+z-c).sum()>Ed :
	X.setdiag(x.todense())
	Z.setdiag(z.todense())

	myu = ganma * (x.T*z) / n

	Z0X = hstack([Z,lil_matrix(zeros((n,m))),X])

	A_X = vstack([A00,ATI,Z0X])


	A_myue = vstack((A*x-b,A.T*y+z-c,X*z-ones((n,1))*myu))

	i+=1
	print i

	A_X = A_X.tocsr()

	A_myue =A_myue.tocsr()

	delta = - spla.spsolve(A_X,A_myue)*0.4

	# print delta

	deltax=delta[:n]
	deltay=delta[n:n+m]
	deltaz=delta[n+m:]

	deltax = lil_matrix(deltax).tocsc().T
	deltay = lil_matrix(deltay).tocsc().T
	deltaz = lil_matrix(deltaz).tocsc().T

	x = x+deltax
	y = y+deltay
	z = z+deltaz

	print x.T*z

	print (A*x-b).multiply(A*x-b).sum()

	print (A.T*y+z-c).multiply(A.T*y+z-c).sum()
	# njnj=x.todense()[:30]

	# print njnj.reshape((10,3))

print "end"





x=x.todense()

x=x.reshape((I+(I-1)*K,L))



l
savetxt("x.csv", x, fmt="%.10f",delimiter=",")

