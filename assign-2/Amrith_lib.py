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


def augment(A,B):
    out=[]
    temp=[]
    for i in range(len(A)):
        temp=A[i].copy()
        for j in range(len(B[i])):
            temp.append(B[i][j]) 
        out.append(temp)
    return out

# print(augment(A,B))

def pivot(l,n):
    return abs(l[n])

def row_exchange(m,r1,r2): #r1<-->r2
    r1_temp=m[r1]
    r2_temp=m[r2]
    m_temp=m.copy()
    m_temp[r1]=r2_temp
    m_temp[r2]=r1_temp
    return m_temp

def row_mult(m,r,c): # r-->c*r
    m_temp=m.copy()
    r_temp=[c*j for j in m[r]]
    m_temp[r] = r_temp
    return m_temp

def row_mult_and_add(m,r1,r2,c): # r1-->r1 + c.r2
    m_temp = m.copy()
    r_temp = []
    for i in range(len(m[r1])):
        r_temp.append(m[r1][i] + c*m[r2][i])  
    m_temp[r1] = r_temp
    return m_temp

def gauss_jordan(m1,m2):
    aug = augment(m1,m2)
    out=aug
    # print(out)
    for i in range(len(aug)):
        piv_list=[]
        for j in range(len(aug)):
            piv_list.append(pivot(out[j],i))
        piv_copy=piv_list.copy()
        max_piv=max(piv_list)
        max_index=piv_list.index(max_piv)
        if max_index<i:
            temp_max=max_piv
            k=0
            while k<=len(piv_list):
                piv_copy.remove(temp_max)
                temp_max=max(piv_copy)
                if piv_list.index(temp_max)>=i:
                    max_piv=temp_max
                    max_index=piv_list.index(temp_max)
                    break
                k+=1

        out=row_exchange(out,i,max_index)
        # print(out)
        out=row_mult(out,i,1/out[i][i])
        # print(out)
        for j in range(len(aug)):
            if (j!=i and out[j][i]!=0):
                out=row_mult_and_add(out,j,i,-out[j][i])
                # print(out)
    out2=[[out[l][-1]] for l in range(len(out))]

    return out2