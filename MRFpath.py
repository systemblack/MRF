# coding: UTF-8

# from pylab import *
import matplotlib.pylab as pylab

from numpy import *

import scipy as sci
import scipy.linalg
import scipy.sparse as sp

import scipy.sparse.linalg as spla

#solve関数では普通に解けない
#from numpy.linalg import solve


#データの読み込み
data = loadtxt("data.csv", delimiter=',')


data = data[:,1:]



data = log(data)


p=0.9
L=3
M=3
I=201
J=1


f = ones([I-1,L,M])*log((1-p)/2)
for i in range(0, I-1):
		for l in range(0,L):
			for m in range(0,M):
				if l==m:
					f[i][l][m]=log(p)


dil=data.flatten()

fijlm= f.flatten()



c = hstack((dil,fijlm))


A = zeros((I+(I-1)*L+(I-1)*(L-1),I*L+(I-1)*L*M))
for i in range(0, I):
		for l in range(0,L):
			A[i,i*L+l]=1
			if(i<I-1):
				A[I+i*L+l,i*L+l]=-1
				if(l<L-1):
					A[I+(I-1)*L+i*(L-1)+l,i*L+l+L]=-1
				for m in range(0,M):
					A[I+i*L+l,I*L+i*L*M+l*M+m]=1
					if(m<M-1):
						A[I+(I-1)*L+i*(L-1)+m,I*L+i*L*M+l+m*L]=1



savetxt("a.csv", A, fmt="%.0f",delimiter=",")


uarr, spec, vharr = linalg.svd(A)
print spec


m=I+(I-1)*L+(I-1)*(L-1)
n=I*L+(I-1)*L*M
print n
print m


b = hstack((ones(I),zeros((I-1)*L+(I-1)*(L-1))))


# xx=zeros((I-1)*L*M)
# for i in range(0, I-1):
# 	xx[L*M*i]=1
# xxx=zeros(I*L)
# for i in range(0,I):
# 	xxx[L*i]=1
# x = hstack((xxx,xx)) 

# print x

x = ones(n)
y = ones(m)
z = ones(n)


#初期値やパラメータ，定数行列の設定
ganma = 0.9
E = 0.00001
Ep = 0.00001
Ed = 0.00001


# A = matrix([[2.0,1.0,1.0,0.0],[1.0,3.0,0.0,1.0]])
# b = matrix([[4.0],[5.0]])
# c = matrix([[-1.0],[-1.0],[0.0],[0.0]])


# x = matrix([[1.0],[1.0],[1.0],[1.0]])
# y = matrix([[-1.0],[-1.0]])
# z = matrix([[2.0],[3.0],[1.0],[1.0]])



# listx = [x[0,0]]
# listy = [x[1,0]]


#配列から行列へ
A = mat(A)
b = mat(b).T
c = mat(c).T

x = mat(x).T
y = mat(y).T
z = mat(z).T

print "shape"
print A.shape
print b.shape
print c.shape
print x.shape
print y.shape
print z.shape


print x.T*z

i=0


A00 = hstack((A,zeros((m,m)),zeros((m,n))))

ATI = hstack((zeros((n,n)),A.T,eye(n)))


#内点法
while x.T*z>E or linalg.norm(A*x-b)>Ep or linalg.norm(A.T*y+z-c)>Ed :
	X = diag(asarray(x.T)[0])
	Z = diag(asarray(z.T)[0])

	myu = ganma * (x.T*z) / n

	Z0X = hstack((Z,zeros((n,m)),X))

	A_X = vstack((A00,ATI,Z0X))

	savetxt("ax.csv", A_X, fmt="%.0f",delimiter=",")

	A_myue = vstack((A*x-b,A.T*y+z-c,X*z-ones((n,1))*myu))


	#savetxt("a.csv", A, fmt="%.0f",delimiter=",")

	i+=1
	print i


	delta = - linalg.inv(A_X)*A_myue*0.6


	x += delta[:n]
	y += delta[n:n+m]
	z += delta[n+m:]

	print x.T*z

print x
print "end"
