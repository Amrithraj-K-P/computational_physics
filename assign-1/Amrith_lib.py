import numpy as np
class mycomplex():
    def __init__(self,real,imag):
        self.r=real
        self.i=imag

    def display(self):
        print(f"{self.r} + {self.i}i")

    def add(self,c1,c2):
        self.r=c1.r + c2.r
        self.i=c1.i+ c2.i
        return
    def sub(self,c1,c2):
        self.r=c1.r - c2.r
        self.i=c1.i - c2.i
        return
    def mult(self,c1,c2):
        self.r=c1.r*c2.r - c1.i*c2.i
        self.i=c1.i*c2.r + c1.r*c2.i
        return
    def modulus(self):
        return np.sqrt(self.r**2 + self.i**2)
    
def read_matrix(filename):
    with open(filename, 'r') as f:
        matrix = []
        for line in f:
            row = [float(num) for num in line.strip().split()]
            matrix.append(row)
    return matrix

def matrix_mult(m1,m2):
    n_rows_1=len(m1)
    n_cols_1=len(m1[0])
    n_rows_2=len(m2)
    n_cols_2=len(m2[0])
    if n_cols_1 != n_rows_2:
        print("error in input")
        return
    out=[]
    i=0
    while i< n_rows_1:
        j=0
        row=[]
        while j< n_cols_2:
            a_ij=0
            k=0
            while k<n_rows_2:
                a_ij += m1[i][k]*m2[k][j]
                k+=1
            row.append(a_ij)
            j+=1
        out.append(row)
        i+=1
    return(out)

def vec_prod(v1,v2):
    temp=0
    i=0
    while i<len(v1):
        temp+= v1[i][0]*v2[i][0]
        i+=1
    return temp

def list_sum(l):
    out=0
    for i in range(len(l)):
        out+=l[i]
    return out

def lcg(x,a,m,c):
    return (a*x + c)%m 

def lcg_rng(n,x0,a = 1103515245,c = 12345,m = 32768 ,range=None):
    out=[]
    i=0
    tempx=x0
    while i<=n:
        temp=lcg(tempx,a=a,m=m,c=c)
        tempx=temp
        out.append(tempx)
        i+=1
    if range==None:
        return out
    else:
        d=range[1]-range[0]
        templ=[i/m for i in out]
        out2=[range[0] + d*x for x in templ]
        return out2
        

def map(x,c):
    return c*x*(1-x)

def rng(c,length,x0):
    out=[]
    i=0
    tempx=x0
    while i<=length:
        temp=map(tempx,c)
        tempx=temp
        out.append(tempx)
        i+=1
    return out
