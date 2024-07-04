import sys

ARG = sys.argv[1]



if ARG == 'test':
    M_VALUES = [1,2,4,16]
    K_MAX = 8
    N_MIN = 1
    N_MAX = 180
elif ARG == 'search_small_N_k':
    M_VALUES = [1,2,4,16]
    K_MAX = 200
    N_MIN = 1
    N_MAX = 2000
elif ARG == 'rangeA':
    M_VALUES = [1,2,4]
    N_MIN = int(sys.argv[2])
    N_MAX = int(sys.argv[3])
elif ARG[0] == 'A':
    M_VALUES = [1,2,4]
    idx = int(ARG[1:])
    N_MIN = max(10^8 * idx, 1)
    N_MAX = 10^8 * (idx+1)
elif ARG == 'B61':
    M_VALUES = [1,2,4,16]
    N_MIN = 61
    N_MAX = 10300000
elif ARG[0] == 'B':
    idx = int(ARG[1:])
    assert 1 <= idx <= 60 
    M_VALUES = [1,2,4,16]
    N_MIN = idx
    N_MAX = idx
else:
    assert False, 'invalid argument: ' + ARG






## A1 ####################################################################



def psi_new(N_fact):
    ret = 1
    for p,r in N_fact:
        if r == 1:
            ret *= p-1
        elif r==2:
            ret *= p^2-p-1
        else:
            ret *= p^r - p^(r-1) - p^(r-2) + p^(r-3)
    return ret


def A1_new(m,N_fact, k):
    if round(sqrt(m))^2 == m:
        return (k-1)/12 * psi_new(N_fact) * m^(k//2-1)
    else:
        return 0


## A2 #########################################


def vp(n, p):
    cnt = 0
    while n % p == 0:
        cnt += 1
        n //= p
    return cnt





def get_t_wt(m):
    ret = []
    for t in range(4*m+1):
        if t^2 >= 4*m: 
            break
        wt = 1/2 if (t == 0) else 1
        ret.append((t, wt))
    return ret



class U_seq:
    def __init__(self, t, m):
        self.t = t
        self.m = m
        self.d = [0,1]
        self.cur = 1

    def get(self,k):
        if k != self.cur: 
            if k < self.cur:
                self.d = [0,1]
                self.cur = 1
            # build up to k
            for i in range(self.cur+1, k+1):
                self.d[i%2] = self.t * self.d[(i-1)%2]  -  self.m * self.d[(i-2)%2]
                self.cur += 1
        return self.d[self.cur%2]


U = {}
for m in M_VALUES:
    for t,_ in get_t_wt(m):
        U[(t,m)] = U_seq(t,m)





def get_n(m,t):
    ret = []
    for n in range(1, 4*m - t^2 + 1):
        if (t^2 - 4*m) % (n^2) == 0 and ((t^2 - 4*m)//(n^2))%4 in [0,1]:
            ret.append(n)
    return ret



# iterate over all t,n,m, p<70
MU_SUM_LOCAL = {}
for m in M_VALUES:
    for t,_ in get_t_wt(m):
        for p in prime_range(70):
            # find all the solutions mod p^r
            SOLS = [(0,)]
            for r in range(1, 50):
                if p^r > 10^12:
                    break
                # go from sols mod p^(r-1) to sols mod p^r
                sln_lifts = []
                for sln in SOLS[r-1]:
                    for i in range(p):
                        sln_lft = sln + i*p^(r-1)
                        if (sln_lft^2 - t*sln_lft + m) % p^r == 0:
                            sln_lifts.append(sln_lft)
                SOLS.append(tuple(sln_lifts))
            for n in get_n(m,t):
                vpn = vp(n, p)
                mu_sum_lst = [1]
                for r in range(1,50):
                    if p^r > 12*10^9:
                        break
                    num_sls = len({sln % p^r for sln in SOLS[r+vpn]})
                    mu_sum_lst.append(num_sls)
                MU_SUM_LOCAL[(t,n,m,p)] = mu_sum_lst[:]
            
                




def mu_sum_local(t,n,m,p,r):
    D = t^2 - 4*m
    # if p > 70:
    if D % p != 0 and m % p != 0 and p != 2:
        # from Hensel's Lemma ; just consider it mod p
        return kronecker(D,p) + 1
    else:
        return MU_SUM_LOCAL[(t,n,m,p)][r]


def mu_local(t,n,m,p,r):
    vpn = vp(n,p)
    if r <= vpn:
        ret = (p+1)*p^(r-1)
    else:
        ret = p^vpn
    return ret * mu_sum_local(t,n,m,p,r)



def mu_new(t,n,m, N_fact):
    ret = 1
    for p,r in N_fact:
        if r == 1:
            ret *= mu_local(t,n,m,p,1) - 2
        elif r == 2:
            ret *= mu_local(t,n,m,p,2) - 2*mu_local(t,n,m,p,1) + 1
        else:
            ret *= mu_local(t,n,m,p,r) - 2*mu_local(t,n,m,p,r-1) + mu_local(t,n,m,p,r-2)
    return ret




_6_h_w = {-3: 2, -4: 3, -7: 6, -8: 6, -11: 6, -12: 6, -15: 12, -16: 6, -19: 6, -20: 12, 
        -23: 18, -24: 12, -27: 6, -28: 6, -31: 18, -32: 12, -35: 12, -36: 12, -39: 24,
        -40: 12, -43: 6, -44: 18, -47: 30, -48: 12, -51: 12, -52: 12, -55: 24, -56: 24,
        -59: 18, -60: 12, -63: 24, -64: 12, -67: 6}



def A2_new(m,N_fact,k):
    ret = 0 
    for t, wt in get_t_wt(m):
        for n in get_n(m,t):
            ret += wt * U[(t,m)].get(k-1) * _6_h_w[(t^2 - 4*m)//(n^2)]/6 * mu_new(t,n,m,N_fact)
    return ret





## A3 ######################################

def sigma_local(m,d, p, r):
    h = abs(d-m//d)
    vph = 100 if h == 0 else vp(h,p)
    if h == 0 or r <= 2*vph:
        if r % 2 == 1:
            return 2 * p^((r-1)//2)
        else:
            return (p+1) * p^((r-2)//2)
    else:
        return 2 * p^vph





def sigma_new(m,d, N_fact):
    ret = 1
    for p,r in N_fact:
        if r == 1:
            ret *= sigma_local(m,d,p,1) - 2
        elif r == 2:
            ret *= sigma_local(m,d,p,2) - 2*sigma_local(m,d,p,1) + 1
        else:
            ret *= sigma_local(m,d,p,r) - 2*sigma_local(m,d,p,r-1) + sigma_local(m,d,p,r-2)
    return ret



# d <= sqrt(m) with weight 1/2 if d == sqrt(m)
def get_d_wt(m):
    ret = []
    for d in divisors(m):
        if d^2 < m:
            ret.append((d,1))
        elif d^2 == m:
            ret.append((d,1/2))
    return ret



def A3_new(m,N_fact,k):
    ret = 0
    for d,wt in get_d_wt(m):
        ret += wt * d^(k-1) * sigma_new(m,d,N_fact)
    return ret
    
 
## A4 #################################################

def A4_new(m,N,k):
    if k > 2:
        return 0
    else:
        return sum(divisors(m)) * moebius(N)


## trace and a2 ###########################################


def get_trace(m,N,k):
    assert gcd(m,N)==1
    N_fact = factor(N)
    A1 = A1_new(m,N_fact,k) 
    A2 = A2_new(m,N_fact,k) 
    A3 = A3_new(m,N_fact,k) 
    A4 = A4_new(m,N,k)
    ret = A1 - A2 - A3 + A4
    assert int(ret) == ret
    return int(ret) 


def get_a2_coeff(m, N, k):
    assert gcd(m,N)==1
    ret = get_trace(m, N, k)^2
    for d in divisors(m):
        ret -= d^(k-1) * get_trace((m//d)^2, N, k)
    assert ret % 2 == 0
    return ret // 2





def test_trace_a2(largest_m):
    for N in range(N_MIN, N_MAX+1):
        for m in divisors(largest_m):
            if gcd(m,N) != 1:
                continue
            for k in range(2, K_MAX+1, 2):
                tr = get_trace(m,N,k)
                a2_coeff = get_a2_coeff(m,N,k)
                print(f'{m} ({N},{k}) : \t\t Tr^new = {tr} \t\t a2^new = {a2_coeff}')

                MSnew = ModularSymbols(Gamma0(N),k,sign=1).cuspidal_submodule().new_submodule()
                hecke_op = MSnew.hecke_operator(m)
                assert tr == hecke_op.trace()
                cply = [0,0,0] + hecke_op.charpoly().list()
                assert a2_coeff == cply[-3]





## E_ub_2 ################################################

def pi1(N_fact):
    ret = 1
    for (p,r) in N_fact:
        ret *= 1 + (p+1)/(p^p-p-1)
    return ret

def pi2(N_fact):
    ret = 1
    for (p,r) in N_fact:
        ret *= 1 + 1/(p-1)
    return ret

def omega(N_fact):
    return len(N_fact)


def theta1(N, N_fact):
    return sqrt(N) / (psi_new(N_fact) * pi2(N_fact)^2)
def theta2(N, N_fact):
    return 4^omega(N_fact) / psi_new(N_fact)
def theta3(N, N_fact):
    return 2^omega(N_fact) / psi_new(N_fact)
def theta4(N, N_fact):
    return 1 / psi_new(N_fact)




# we want E_ub_2(N) < (k-1)/16
def E_ub_2(N, N_fact):
    ret = 0.0
    ret += (1/2) * theta1(N, N_fact) 
    ret += 32 * theta2(N, N_fact) 
    ret += (16*sqrt(2) + 65/3) * theta3(N, N_fact)
    ret += 37/4 * theta4(N, N_fact)
    return ret.n()



# return kUB where    E_ub_2(N) < (k-1)/16  for all k > kUB
# i.e. we only need to check k <= kUB
def get_kUB_2(N):
    E_ub_2_N = E_ub_2(N, factor(N))
    for k in range(2,10000,2):
        if E_ub_2_N < (k-1)/16:
            return k-2 
    assert False



## E_ub_4 ####################################################

# we want E_ub_4(N) < (k-1)/192
def E_ub_4(N, N_fact):
    th1 = theta1(N, N_fact) 
    th3 = theta3(N, N_fact) 
    th4 = theta4(N, N_fact)
    ret = 0.0
    ret += 1/8 * th1 + 41/4 * th3 + 37/8 * th4 
    ret += 12*( 1/4 * th1 + 41/2 * th3 + 19/4 * th4 )^2
    ret += th4 * ( 9/2 * th1 + 1258 * th3 + 1515/4 * th4 )
    return ret.n()


# return kUB where    E_ub_4(N) < (k-1)/192  for all k > kUB
# i.e. we only need to check k <= kUB
def get_kUB_4(N):
    E_ub_4_N = E_ub_4(N, factor(N))
    for k in range(2,6000000,2):
        if E_ub_4_N < (k-1)/192:
            return k-2 
    assert False




## find all vanishing a2 ########################################




def find_all_nonneg_a2():
    m = 2
    ret = []
    for N in range(N_MIN, N_MAX+1):
        if gcd(N,m) != 1:
            continue
        kUB = get_kUB_2(N)
        for k in range(2, kUB+1, 2):
            a2_coeff = get_a2_coeff(m,N,k)
            if a2_coeff >= 0:
                print('non-neg:', (m,N,k), a2_coeff)


def find_all_nonpos_a2():
    m = 4
    ret = []
    for N in range(N_MIN, N_MAX+1):
        if gcd(N,m) != 1:
            continue
        kUB = get_kUB_4(N)
        for k in range(2, kUB+1, 2):
            a2_coeff = get_a2_coeff(m,N,k)
            if a2_coeff <= 0:
                print('non-pos:', (m,N,k), a2_coeff)

def search_nonneg_small_N_k(m):
    ret = []
    for N in range(N_MIN, N_MAX+1):
        if gcd(N,m) != 1:
            continue
        for k in range(2, K_MAX+1, 2):
            a2_coeff = get_a2_coeff(m,N,k)
            if a2_coeff >= 0:
                print('non-neg:', (m,N,k), a2_coeff)


def search_nonpos_small_N_k(m):
    ret = []
    for N in range(N_MIN, N_MAX+1):
        if gcd(N,m) != 1:
            continue
        for k in range(2, K_MAX+1, 2):
            a2_coeff = get_a2_coeff(m,N,k)
            if a2_coeff <= 0:
                print('non-pos:', (m,N,k), a2_coeff)









if ARG == 'test':
    test_trace_a2(4)
elif ARG == 'search_small_N_k':
    search_nonpos_small_N_k(4)
    search_nonneg_small_N_k(2)
elif ARG == 'rangeA':
    find_all_nonneg_a2()
elif ARG[0] == 'A':
    find_all_nonneg_a2()
elif ARG[0] == 'B':
    find_all_nonpos_a2()
else:
    assert False


