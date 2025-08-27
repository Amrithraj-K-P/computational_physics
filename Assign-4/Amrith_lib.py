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
def L_U_decomposition(matrix): #Doolittle LU factorisation
    m=[row.copy() for row in matrix]

    for j in range(len(m)):
        for i in range(1,j+1):
            temp=0
            for k in range(i):
                temp += m[i][k]*m[k][j]
            # print(temp)
            
            m[i][j] = m[i][j] - temp
            # print(U)

        for i in range(j+1,len(m)):
            temp=0
            for k in range(j):
                temp += m[i][k]*m[k][j]
            
            m[i][j] = ( m[i][j] - temp )/m[j][j]
            # print(L)
    return m
def L_U_Solve(A,B):
    m = L_U_decomposition(A)
    n=len(m)
    U=make_zeros(n)
    L=make_zeros(n)

    for i in range(n):
        for j in range(n):
            if i > j:
                L[i][j] = m[i][j]
                U[i][j] = 0
            elif i == j:
                L[i][j] = 1
                U[i][j] = m[i][j]
            else:  # i < j
                L[i][j] = 0
                U[i][j] = m[i][j]


    y=[]
    x=[0]*len(B)
    y.append(B[0][0])
    for i in range(1,len(B)):
        temp=0
        for j in range(i):
            temp+=L[i][j]*y[j]
        y.append(B[i][0]-temp)
    x[-1]=y[-1]/U[-1][-1]
    for i in range(len(B)-1 , -1,-1):
        temp=0
        for j in range(i+1,len(B)):
            temp+=U[i][j]*x[j]
        temp2=y[i]-temp
        x[i] = (temp2/U[i][i])
    return x
def make_zeros(n):
    out=[]
    for i in range(n):
        row=[]
        for j in range(n):
            row.append(0)
        out.append(row)
    return out

def transpose(m):
    n=len(m)
    out=make_zeros(n)
    for i in range(n):
        for j in range(n):
            out[i][j] = m[j][i]
    return out

def cholesky_decomp(m):#only for Symmetric, positive definite
    n=len(m)
    L=make_zeros(n)
    for i in range(n):
        temp=0
        for j in range(i):
            temp+= (L[j][i])**2
        L[i][i] = np.sqrt(m[i][i] - temp)

        for j in range(i+1,n):
            temp=0
            for k in range(i):
                temp += L[k][i]*L[k][j] 
            L[i][j] = (1/L[i][i])*(m[i][j] - temp)

    return L

def cholesky_solve(A,B):
    U=cholesky_decomp(A)
    y=[]
    x=[0]*len(B)
    y.append(B[0][0]/U[0][0])
    for i in range(1,len(B)):
        temp=0
        for j in range(i):
            temp+=U[j][i]*y[j]
        y.append((B[i][0]-temp)/U[i][i])
    x[-1]=y[-1]/U[-1][-1]
    for i in range(len(B)-1 , -1,-1):
        temp=0
        for j in range(i+1,len(B)):
            temp+=U[i][j]*x[j]
        temp2=y[i]-temp
        x[i] = (temp2/U[i][i])
    return x

def D_L_U_sep(m):
    n=len(m)
    D=make_zeros(n)
    L=make_zeros(n)
    U=make_zeros(n)
    for i in range(n):
        D[i][i]=m[i][i]
        for j in range(i+1,n):
            L[j][i]=m[j][i]
            U[i][j]=m[i][j]
    return D,L,U

def dist(v1,v2):
    out=[]
    for i in range(len(v1)):
        out.append(abs(v1[i]-v2[i]))
    return max(out)

def jacobi_iter(m,b,guess,epsilon):
    n=len(m)
    D,L,U= D_L_U_sep(m)
    count=0
    x_k=guess
    x_k2=[0]*n
    while count<=1000:
        count+=1
        for i in range(n):
            temp=0
            for j in range(n):
                if j!=i:
                    temp+=m[i][j]*x_k[j]
            x_k2[i] = (b[i][0] - temp)/m[i][i]

        d=dist(x_k,x_k2)
        if d<=epsilon:
            return x_k2,count
        else:
            x_k = x_k2.copy()

    print('does not converge')
    return
