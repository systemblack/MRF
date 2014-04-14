# coding: UTF-8

from pylab import *

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
L=2
M=2
I=2
J=2


f = ones([I,J,L,M])*log((1-p)/2)
for i in range(0, I):
	for j in range(0,J):
		for l in range(0,L):
			for m in range(0,M):
				if l==m:
					f[i][j][l][m]=log(p)


dil=data.flatten()
fijlm= f.flatten()



c = hstack((dil,fijlm))


A = zeros((I+I*J*L*2-J*L,I*L+I*J*L*M))
for i in range(0, I):
	for j in range(0,J):
		for l in range(0,L):
			if j==0:
				A[i,l+L*i]=1.0
			for m in range(0,M):
				A[l+j*L+i*L*J+I,l+L*i] = -1.0

				A[I+l+L*j+L*J*i,I*L+m+M*l+M*L*j+M*L*J*i]=1.0

				if i!=0 and j==0 :
					A[I+I*J*L+l+j*L+i*L*J-L,l+L*i-M] = -1.0
					A[I+I*J*L+m+L*j+L*J*i-L,I*L+m+M*l+M*L*j+M*L*J*i]=1.0
				if i != I-1 and j==1:
					A[I+I*J*L+l+j*L+i*L*J-L,l+L*i+M] = -1.0
					A[I+I*J*L+m+M*j+M*J*i-L,I*L+m+M*l+M*L*j+M*L*J*i]=1.0


savetxt("a.csv", A, fmt="%.0f",delimiter=",")



# print "Ashape"
# print A.shape

b = hstack((ones(I),zeros(I*J*L*2-6)))

x = ones(I*L+I*J*L*M)
y = ones(I+I*J*L*2-6)
z = ones(I*L+I*J*L*M)


#初期値やパラメータ，定数行列の設定
ganma = 0.9
E = 0.00001
Ep = 0.00001
Ed = 0.00001


# A = matrix([[2.0,1.0,1.0,0.0],[1.0,3.0,0.0,1.0]])
# b = matrix([[4.0],[5.0]])
# c = matrix([[-1.0],[-1.0],[0.0],[0.0]])


m=I+I*J*L*2-6
n=I*L+I*J*L*M
print n
print m
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

print "テスト"
X = diag(asarray(x.T)[0])
Z = diag(asarray(z.T)[0])
A00 = hstack((A,zeros((m,m)),zeros((m,n))))


ATI = hstack((zeros((n,n)),A.T,eye(n)))
Z0X = hstack((Z,zeros((n,m)),X))

A_X = vstack((A00,ATI,Z0X))

print linalg.det(A_X)




i=0
#内点法
# while x.T*z>E or linalg.norm(A*x-b)>Ep or linalg.norm(A.T*y+z-c)>Ed:

# 	X = diag(asarray(x.T)[0])
# 	Z = diag(asarray(z.T)[0])


# 	myu = ganma * (x.T*z) / n



# 	A00 = hstack((A,zeros((m,m)),zeros((m,n))))


# 	ATI = hstack((zeros((n,n)),A.T,eye(n)))
# 	Z0X = hstack((Z,zeros((n,m)),X))



# 	A_X = vstack((A00,ATI,Z0X))


# 	A_myue = vstack((A*x-b,A.T*y+z-c,X*z-ones((n,1))*myu))


# 	i+=1
# 	print i


# 	print "success"

# 	delta = - linalg.pinv(A_X)*A_myue


# 	print delta.shape

# 	print "success2"
# 	x += delta[:n]
# 	y += delta[n:n+m]
# 	z += delta[n+m:]



# print "end"
	#print x.T

	# listx.append(x[0,0])
	# listy.append(x[1,0])

#プロット

# def f(x,y):
#     return -x-y

# x = arange(0.1, 2.1, 0.1)
# y = arange(0.1, 2.1, 0.1)
# X,Y = meshgrid(x, y)


# figure(num=None, figsize=(8, 8), dpi=80, facecolor='w', edgecolor='k')

# axes([0.05, 0.05, 0.9, 0.9])

# contourf(X, Y, f(X, Y), 8, alpha=.65, cmap=cm.gray)
# C = contour(X, Y, f(X, Y), 8, colors='black', linewidth=.5)
# clabel(C, inline=1, fontsize=16)


# plot(listx,listy,'ro')
# plot(listx[:],listy[:])

# show()



# f = ones([100,2,3,3])

# print f.shape


# ii=0
# for i in range(0, 10):
# 	for j in range(0,2):
# 		for l in range(0,3):
# 			for m in range(0,3):
# 				f[i][j][l][m]=j
# 				ii+=1

# #print f.flatten()
# print (f.flatten()).shape
