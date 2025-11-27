import numpy as np
from numpy.polynomial import legendre
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

def matrix_vec_mult(m1,v):
    n_rows_1=len(m1)
    n_cols_1=len(m1[0])
    n_v=len(v)
    if n_cols_1 != n_v:
        print("error in input")
        return

    j=0
    out=[]
    while j< n_cols_1:
        v_j=0
        k=0
        while k<n_v:
            v_j += m1[j][k]*v[k]
            k+=1
        out.append(v_j)
        j+=1

    return(out)

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
    r1_temp=[j for j in m[r1]]
    r2_temp=[j for j in m[r2]]
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

def matrix_inverse(m1):
    m2=make_zeros(len(m1))
    for i in range(len(m1)):
        m2[i][i] = 1
    
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
    out2=[out[l][len(aug)::] for l in range(len(out))]
    return out2

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
    # out=[]
    # for i in range(len(v1)):
    #     out.append(abs(v1[i]-v2[i]))
    # return max(out)
    out=0
    for i in range(len(v1)):
        out += (v1[i]-v2[i])**2
    return np.sqrt(out)

def jacobi_sol(m,b,guess,epsilon):
    n=len(m)
    D,L,U= D_L_U_sep(m)
    count=0
    x_k=guess
    x_k2=[0]*n
    while count<=100000:
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
    return False, count

def check_symmetry(m):
    n=len(m)
    for i in range(n):
        for j in range(n):
            if m[i][j]!=m[j][i]:
                return False
    return True

def L_star_U_sep(m):
    n=len(m)
    L=make_zeros(n)
    U=make_zeros(n)
    for i in range(n):
        L[i][i]=m[i][i]
        for j in range(i+1,n):
            L[j][i]=m[j][i]
            U[i][j]=m[i][j]
    return L,U

def gauss_siedel_solve(a,b,guess,eps):
    n=len(b)
    x=guess
    count=0
    while count<=10000:
        count+=1
        dist=0
        for i in range(n):
            tempx=x[i]
            temp=0
            temp2=0
            for j in range(i):
                temp+=a[i][j]*x[j]
            for j in range(i+1,n):
                temp2+=a[i][j]*x[j]
            
            x[i] = (b[i][0] - temp - temp2)/a[i][i]
            dist += (tempx-x[i])**2
            # if dist<=np.abs(tempx-x[i]):
            #     dist=np.abs(tempx-x[i])
        if np.sqrt(dist)<=eps:
            return x,count
    print('does not converge')
    return False,count


bs_count=0

def bisection_root(f,bracket,accur):
    global bs_count
    bs_count +=1
    a=bracket[0]
    b=bracket[1]
    if f(a)*f(b)>=0:
        print('Wrong bracket input')
        return
    if abs(b-a)<=accur:
        return (a+b)/2 , bs_count
    else:
        c=(a+b)/2
        if f(a)*f(c)<=0:
            return bisection_root(f,[a,c],accur)
        else:
            return bisection_root(f,[c,b],accur)


rf_count=0
def regula_falsi_root(f,bracket,accur):
    global rf_count
    rf_count+=1
    a=bracket[0]
    b=bracket[1]
    if f(a)*f(b)>=0:
        print(a,b)
        print('Wrong bracket input')
        return
    c=b - (b-a)*f(b)/(f(b)-f(a))
    if abs(b-a)<=accur:
        return (a+b)/2 , rf_count
    if abs(f(c))<=accur:
        return c, rf_count  
    else:
        # c1=a-(b-a)*f(a)/(f(b)-f(a))
        
        # if a<=c1<=b:
        #     c=c1
        # else:
        #     c=c2
        
        if f(a)*f(c)<=0:
            return regula_falsi_root(f,[a,c],accur)
        else:
            return regula_falsi_root(f,[c,b],accur)

def find_bracket(f,bracket,beta):
    a=bracket[0]
    b=bracket[1]
    if f(a)*f(b)<=0:
        return [a,b]
    else:
        if abs(f(a))<=abs(f(b)):
            return find_bracket(f,[a-beta*(b-a),b],beta)
        else:
            return find_bracket(f,[a,b+beta*(b-a)],beta)

def newton_raphson_root(f,df,accur,x0_guess):
    count=0
    x0=x0_guess
    x1 = x0 - f(x0)/df(x0)
    while abs(x1-x0)>=accur:
        count+=1
        x0=x1
        x1 = x1 - f(x1)/df(x1)
    return float(x1), count

def fixed_pont_root(g,guess,accur):
    count=0
    x0=guess
    x1=g(x0)
    while abs(x1-x0)>accur:
        count+=1
        x0=x1
        x1=g(x1)
    return x1,count

def multivar_fixed_point_root(g,accur,guess):
    count=0
    x0=guess
    x1=g(x0)
    while dist(x0,x1)/dist(x1,[0]*len(x1))>accur:
        count+=1
        x0=x1
        x1=g(x1)
    return x1,count

def multivar_newton_raphson_root(f,J,accur,x0_guess):
    count=0
    x0=x0_guess
    x1 = [x0[i] - matrix_vec_mult(matrix_inverse(J(x0)),f(x0))[i]  for i in range(len(x0))]
    while dist(x0,x1)/dist(x1,[0]*len(x1))>accur:
        count+=1
        x0=x1
        x1 = [x1[i] - matrix_vec_mult(matrix_inverse(J(x1)),f(x1))[i] for i in range(len(x0))]
    return x1, count

def poly_diff(p):
    n=len(p)
    out=[]
    for i in range(n):
        out.append(p[i]*(n-i-1))

    return out[:-1]

def poly_val(p,x):
    n=len(p)
    out=0
    for i in range(n):
        out+= p[i]*(x**(n-i))
    return out

def poly_divide(p,d):
    out=[p[0]]
    for i in range(len(p)-1):
        out.append(p[i+1]+out[-1]*d)

    return out[:-1], out[-1]
    

def laguerre_step(p,guess,accur):
    count=0
    n=len(p)-1
    x0=guess
    
    if abs(poly_val(p,x0))<=accur:
        return x0,count

    G=poly_val(poly_diff(p),x0)/poly_val(p,x0)
    H=G**2 - poly_val(poly_diff(poly_diff(p)),x0)/poly_val(p,x0)
    if G>=0:
        a = n/(G + np.sqrt(abs((n-1)*(n*H - G**2))))
    else:
        a = n/(G - np.sqrt(abs((n-1)*(n*H - G**2))))
    x1=x0-a

    if abs(poly_val(p,x1))<=accur:
        return x1,count


    while abs(x1-x0)> accur or abs(poly_val(p,x1))>accur:
        if abs(poly_val(p,x0))<=accur:
            return x0,count
        x0=x1
        G=poly_val(poly_diff(p),x0)/poly_val(p,x0)
        H=G**2 - poly_val(poly_diff(poly_diff(p)),x0)/poly_val(p,x0)
        
        if G>=0:
            a = n/(G + np.sqrt(abs((n-1)*(n*H - G**2))))
        else:
            a = n/(G - np.sqrt(abs((n-1)*(n*H - G**2))))
        x1=x0-a 
        count+=1
    return x1,count
    

def laguerre_root(p,guess,accur):
    n=len(p)-1
    p_temp=p.copy()
    root_out=[]
    rem_out=[]
    count=0
    for i in range(n):
        x0=guess[i]
        alpha_i,c = laguerre_step(p_temp,x0,accur)
        root_out.append(float(alpha_i))
        q,r= poly_divide(p_temp,alpha_i)
        rem_out.append(float(r))
        count+=c
        p_temp = q

    return root_out,rem_out



def midpoint_integ(f,lim,n): #lim : [a,b]
    a,b=lim
    h=abs((b-a)/n)
    x_i=a #dummy variable
    out=0
    while x_i<b:
        out += f((x_i + x_i + h)/2)
        x_i += h
    return out*h


def trapezoidal_integ(f,lim,n): #lim : [a,b]
    a,b=lim
    h=abs((b-a)/n)
    out = f(a) + f(b)
    x_i = a+h
    while x_i<b: #summation is done over terms with weight = 2
        out += 2*f(x_i)
        x_i+=h
    return out*h/2

def midpoint_integ_accur(f,lim,err,d2f_max): #lim : [a,b]
    a,b=lim
    n=np.ceil(np.sqrt(abs(d2f_max) * (b-a)**3 / (24 * err))) #finding value of n
    h=abs((b-a)/n)
    x_i=a #dummy variable
    out=0
    while x_i<=b:
        out += f((x_i + x_i + h)/2)
        x_i += h
    return float(out*h), float(n)


def trapezoidal_integ_accur(f,lim,err,d2f_max): #lim : [a,b]
    a,b=lim
    n=np.ceil(np.sqrt(abs(d2f_max) * (b-a)**3 / (12 * err))) #finding value of n
    h=abs((b-a)/n)
    out = f(a) + f(b)
    x_i = a+h
    while x_i<b: #summation is done over terms with weight = 2
        out += 2*f(x_i)
        x_i+=h
    return float(out*h/2), float(n)

def simpson_integ_accur(f,lim,err,d4f_max):#lim : [a,b]
    a,b=lim
    n=np.ceil((abs(d4f_max) * (b-a)**5 / (180 * err))**0.25) #finding value of n
    print(n)
    if n%2!=0:
        n+=1
    h=abs((b-a)/n)
    out = f(a) + f(b)
    i=1
    x_i = a+h
    while x_i<b: #weight is 1 for a and b, 2 for odd and 4 for even
        if i%2!=0:
            out += 4*f(x_i)
        else:
            out += 2*f(x_i)
        i+=1
        x_i+=h
    return float(out*h/3), float(n)

def simpson_integ(f,lim,err,n):#lim : [a,b]
    a,b=lim
    if n%2!=0:
        n+=1
    h=abs((b-a)/n)
    out = f(a) + f(b)
    i=1
    x_i = a+h
    while x_i<b: #weight is 1 for a and b, 2 for odd and 4 for even
        if i%2!=0:
            out += 4*f(x_i)
        else:
            out += 2*f(x_i)
        i+=1
        x_i+=h
    return float(out*h/3)

def monte_carlo_integ(f,lim,n):
    a,b=lim
    ran=lcg_rng(n,x0=0.1,range=lim) #rng with seed=0.1 and range =[a,b]
    i=0
    sum_f=0 #sum f(xi)
    sum_f2=0 #sum f(xi)^2
    while i<n:
        sum_f += f(ran[i])
        sum_f2 += f(ran[i])**2
        i+=1
    integral=sum_f*(b-a)/n
    sigma2 = sum_f2/n - (sum_f/n)**2
    sigma=np.sqrt(sigma2)
    return float(integral), float(sigma)

def gaussian_quadrature_integ(f,limits,n):
    a,b=limits[0],limits[1]
    points, weights = legendre.leggauss(n) # points, weights = legendre.leggauss(degree)
    out=0
    for i in range(n):
        out+=weights[i]*f(((b-a)/2)*points[i] + (b+a)/2) #for when the interval is not [-1,1]
    return out*(b-a)/2

def forward_euler_solve(f,init,t_range,step_size=0.1):
    #  here f is taken as a vector, i.e, f=[f1,f2,...,fn]
    #  initial condition is also taken as a vector, i.e, init=[y1_0,y2_0,...,yn_0]
    # t_range = [t0,t_fin]
    temp_sol=init
    t_init,t_fin=t_range
    temp_t=t_init
    sol_out=[init]
    t_out=[t_init]
    while temp_t<=t_fin:
        k_1 = vect_scalar_mult(step_size,f(temp_sol,temp_t))#finding k1
        temp_sol=[temp_sol[i]+k_1[i] for i in range(len(temp_sol))] #updating solution for each time step
        temp_t+=step_size
        sol_out.append(temp_sol)
        t_out.append(temp_t)
    #sol_out=[[y1_0,y2_0,...,yn_0],...,[y1_tfin,y2_tfin,...,yn_tfin]]
    return t_out, sol_out 

def predictor_corrector_solve(f,init,t_range,step_size=0.1):
    #  here f is taken as a vector, i.e, f=[f1,f2,...,fn]
    #  initial condition is also taken as a vector, i.e, init=[y1_0,y2_0,...,yn_0]
    # t_range = [t0,t_fin]
    temp_sol=init
    t_init,t_fin=t_range
    temp_t=t_init
    sol_out=[init]
    t_out=[t_init]
    while temp_t<=t_fin:
        k_1 = vect_scalar_mult(step_size,f(temp_sol,temp_t))#finding k1
        k_2 = vect_scalar_mult(step_size,f(vect_add(temp_sol,k_1), temp_t + step_size))
        temp_sol = [temp_sol[i] + 0.5* (k_1[i] + k_2[i]) for i in range(len(k_1))]#updating solution for each time step
        temp_t+=step_size
        sol_out.append(temp_sol)
        t_out.append(temp_t)
    #sol_out=[[y1_0,y2_0,...,yn_0],...,[y1_tfin,y2_tfin,...,yn_tfin]]
    return t_out, sol_out

def vect_add(v1,v2): #addition of two vectors
    return [v1[i] + v2[i] for i in range(len(v1))]
def vect_scalar_mult(s,v): #scalar multiplication with a vector
    return [s*i for i in v]
def RK4_solve(f,init,t_range,step_size=0.1):
    #  here f is taken as a vector, i.e, f=[f1,f2,...,fn]
    #  initial condition is also taken as a vector, i.e, init=[y1_0,y2_0,...,yn_0]
    # t_range = [t0,t_fin]
    temp_sol=init
    t_init,t_fin=t_range
    temp_t=t_init
    sol_out=[init]
    t_out=[t_init]
    while temp_t<=t_fin:
        k_1 = vect_scalar_mult(step_size,f(temp_sol,temp_t))#finding k1
        k_2 = vect_scalar_mult(step_size,f(vect_add(temp_sol,vect_scalar_mult(0.5,k_1)), temp_t + step_size*0.5))#finding k2
        k_3 = vect_scalar_mult(step_size, f(vect_add(temp_sol,vect_scalar_mult(0.5,k_2)), temp_t + step_size*0.5))#finding k3
        k_4 = vect_scalar_mult(step_size,f(vect_add(temp_sol,k_3), temp_t + step_size)) #finding k4
        temp_sol = [temp_sol[i] + (k_1[i] + 2*k_2[i] + 2* k_3[i] + k_4[i])/6 for i in range(len(k_1))]#updating solution for each time step
        temp_t+=step_size
        sol_out.append(temp_sol)
        t_out.append(temp_t)
    #sol_out=[[y1_0,y2_0,...,yn_0],...,[y1_tfin,y2_tfin,...,yn_tfin]]
    return t_out, sol_out

def lagrange_interpolate(x,x_vals,y_vals):
    n=len(x_vals)
    out = 0
    i=0
    while i<n:
        k=0
        out2=1
        while k<n:
            if k==i:
                k+=1
            else:
                out2 = out2 * (x - x_vals[k])/(x_vals[i] - x_vals[k])
                k+=1
        out += out2 * y_vals[i]
        i+=1
    return out

def lsqfit(x_data,y_data,sigma=None):
    n = len(x_data)
    if sigma==None:
        sigmas = [1]*n #choosing sigmas to be 1 by default
    else:
        sigmas = sigma

    Sy = sum(y_data[i]/sigmas[i]**2 for i in range(n))
    Sx = sum(x_data[i]/sigmas[i]**2 for i in range(n))
    Sxy = sum(y_data[i]*x_data[i]/sigmas[i]**2 for i in range(n))
    Sxx = sum(x_data[i]**2/sigmas[i]**2 for i in range(n))
    Syy =  sum(y_data[i]**2/sigmas[i]**2 for i in range(n))
    S = sum(1/sigmas[i]**2 for i in range(n))
    Delta = S*Sxx - Sx**2

    #slope and intercept
    intercept = (Sxx * Sy - Sx * Sxy)/Delta
    slope = (Sxy - Sx*Sy)/Delta

    ## Error in slope and intercept
    # sigmasqr_intercept = Sxx/Delta 
    # sigmasqr_slope = S/Delta

    #Pearson’s correlation coefficient
    r_sqr = Sxy**2 / (Sxx * Syy)

    return intercept,slope,r_sqr
