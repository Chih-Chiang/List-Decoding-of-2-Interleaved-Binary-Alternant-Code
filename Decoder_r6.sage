# cd D:/Course/master\ thesis/Sage/ 

reset()
import random
import time
import sys

# testing parameters
debug = False
times = 500

@cached_function
def mybinomial(a,b):
     # this is just to speed up binomial coefficient computation 
     return binomial(a,b)

@cached_function
def weighted_degree(leading_monomial):
    lm_weighted_degree = 1*leading_monomial.degree(x) + (k_GRS-1)*leading_monomial.degree(y) + (k_GRS-1)*leading_monomial.degree(z)
    return lm_weighted_degree

# testing counters
count_success = 0
count_total = 0
count_size1 = 0

# initialization
random.seed(time.time())

# input some argument
print('input t or input 0 to use default setting')
input_t = int(input('t : '))
if input_t > 0:
    input_m = int(input('m : '))
    input_m_bar = int(input('m_bar : '))
    print("We choose t =", input_t, ", m =", input_m, ", m_bar =", input_m_bar)

# start testing
before = time.time()
for testing_times in range(times):
    M = 2
    
    q_GRS = 2^5
    q_Goppa = 2
    
    n = 32
    
    F_GRS = GF(q_GRS,'a')
    F_Goppa = GF(q_Goppa, 'a')
    Rx.<x> = F_GRS[]
    
    g_x = (x^6 + x^1 + 1)^2
    r = g_x.degree()  # 0 < r < n
    L = [a for a in F_GRS.list() if g_x(a) != 0]
    C_Goppa = codes.GoppaCode(g_x, L[0:n])
    column_multiplier_list = []
    for ii in range(n):
        multiplier = g_x(L[ii])
        for jj in range(n):
            if ii != jj:
                multiplier = multiplier / (L[ii] - L[jj])
        column_multiplier_list.append(multiplier)
    C_GRS = codes.GeneralizedReedSolomonCode(F_GRS.list()[0:n], n-r, column_multiplier_list)
    
    print(C_Goppa)
    print(C_GRS)
    assert C_Goppa.dual_code().generator_matrix().echelon_form() == codes.SubfieldSubcode(C_GRS, F_Goppa).parity_check_matrix().echelon_form(), "not alternant codes"

    k_GRS = C_GRS.dimension()
    k_Goppa = C_Goppa.dimension()
    
    d_GRS = C_GRS.minimum_distance()
    d_Goppa = r+1 #C_Goppa.minimum_distance()
    delta_GRS = d_GRS/n
    delta_Goppa = d_Goppa/n
    
    R_GRS = (k_GRS-1)/n
    K_q = RR(1 - q_Goppa/(q_Goppa-1)*(1-R_GRS))
 
    tau = RR((q_Goppa-1)/q_Goppa*(1-sqrt(1-q_Goppa/(q_Goppa-1)*delta_Goppa)))
    epsilon = 0.2
    gamma = RR((1-epsilon) * tau)
    # t = floor(gamma * n)
    t = 11 # sync. error
    # m_total = RR(1/(2*epsilon*tau*sqrt(K_q)))
    m_total = 6
    # m = ceil(m_total*(1-gamma))
    # m_bar = ceil(m_total*gamma/(q_Goppa^2-1))
    m = 5
    m_bar = 1
    
    if input_t > 0:
        t = input_t
        m = input_m
        m_bar = input_m_bar
        m_total = m+m_bar
    
    print("\n")
    print("We choose t =", t, ", m =", m, ", m_bar =", m_bar)
    
    #DELTA = m_total*n*((1-gamma)^2+gamma^2/(q_Goppa^2-1))
    DELTA = RR(pow(n*(k_GRS-1)^2*m*(m+1)*(m+2) + n*(q_Goppa^2-1)*(k_GRS-1)^2*m_bar*(m_bar+1)*(m_bar+2), 1/3))

    E_Goppa = codes.encoders.GoppaCodeEncoder(C_Goppa)
    print("\n")
    print(E_Goppa)
    print(E_Goppa.message_space())

    VS = F_Goppa^k_Goppa
    info1 = VS.random_element()
    info2 = VS.random_element()
    c1 = E_Goppa.encode(info1)
    c2 = E_Goppa.encode(info2)

    print("\n")
    print("c1:")
    print(c1)
    print("c2:")
    print(c2)
    
    E_GRS = C_GRS.encoder("EvaluationPolynomial")
    p1 = E_GRS.unencode(c1)
    p2 = E_GRS.unencode(c2)
    
    # Error
    error_count = t      #number of synchronized error
    error_position = sample(range(n), error_count)

    error1 = vector(F_Goppa, n)
    error2 = vector(F_Goppa, n)
    for ii in range(error_count):
        error1[error_position[ii]] = F_Goppa.list()[randint(0, q_Goppa-1)]
        error2[error_position[ii]] = F_Goppa.list()[randint(0, q_Goppa-1)]
        while error1[error_position[ii]] == 0 and error2[error_position[ii]] == 0:
            error1[error_position[ii]] = F_Goppa.list()[randint(0, q_Goppa-1)]
            error2[error_position[ii]] = F_Goppa.list()[randint(0, q_Goppa-1)]
    error_position.sort()
    print("\n")
    print (error_count, "error @", error_position)
    error_matrix = MatrixSpace(F_Goppa, 2, n).matrix([*error1, *error2])
    print(error_matrix)
    
    r1 = c1 + error1
    r2 = c2 + error2
    
    print("\n")
    print("r1:")
    print(r1)
    print("r2:")
    print(r2)
    
    #Groebner basis
    print("\n")
    print("Groebner basis")
    mu = 0
    while True:
        if (k_GRS-1)*(mybinomial(mu+2,3)) > n*(mybinomial(m+2,3)) + n*(q_Goppa^2-1)*(mybinomial(m_bar+2,3)):
            break
        mu = mu + 1
    while (k_GRS-1)*mu < DELTA: # when m and m_bar are small. (k_GRS-1)*mu is about but less than DELTA
        mu = mu + 1
    l = mybinomial(mu+1,2)

    Rxyz.<x,y,z> = PolynomialRing(F_GRS, 3, order = TermOrder('wdeglex',(1,k_GRS-1,k_GRS-1)))
    G_set = []
    for ii in range(mu):
        for jj in range(mu-ii):
            G_set.append(y^ii*z^jj)

    delta = [0] * l
    iteration_count = RR(1)
    total_iteration = RR(n*len(F_Goppa)^M)
    print ("iterations:", total_iteration)
    for ii in range(n):
        for elements_in_based_field in [(element1, element2) for element1 in F_Goppa for element2 in F_Goppa]:
            # print progress here
            print("progress: {:.2%}".format(iteration_count/total_iteration), end = '\r')
            iteration_count = iteration_count + 1
            if elements_in_based_field[0] == r1[ii] and elements_in_based_field[1] == r2[ii]:
                for jj in range(m):
                    for kk in range(m-jj):
                        for ll in range(m-jj-kk):
                            min_weighted_degree = sys.maxsize
                            mm_min_weighted_degree = l 
                            for mm in range(l):
                                monomials = G_set[mm].monomials()
                                delta[mm] = 0
                                for monomial in monomials:
                                    degx = monomial.degree(x)
                                    degy = monomial.degree(y)
                                    degz = monomial.degree(z)
                                    if degx<jj or degy<kk or degz<ll:
                                        continue
                                    g_mm = G_set[mm].monomial_coefficient(monomial)
                                    delta[mm] = delta[mm] + mybinomial(degx,jj) * \
                                                            mybinomial(degy,kk) * \
                                                            mybinomial(degz,ll) * \
                                                            g_mm * \
                                                            L[ii]                                                        ^ (degx-jj) * \
                                                            (elements_in_based_field[0]*column_multiplier_list[ii]^(-1)) ^ (degy-kk) * \
                                                            (elements_in_based_field[1]*column_multiplier_list[ii]^(-1)) ^ (degz-ll)
                                if delta[mm] != 0 and weighted_degree(G_set[mm].lm()) < min_weighted_degree:
                                    mm_min_weighted_degree = mm
                                    min_weighted_degree = weighted_degree(G_set[mm].lm())
                            if mm_min_weighted_degree == l:
                                continue
                            for mm in range(l):
                                if mm != mm_min_weighted_degree:
                                    G_set[mm] = G_set[mm] - delta[mm]/delta[mm_min_weighted_degree]*G_set[mm_min_weighted_degree]
                            G_set[mm_min_weighted_degree] = (x-L[ii]) * G_set[mm_min_weighted_degree]
            else: #elements_in_based_field[0] != r1[ii] or elements_in_based_field[1] != r2[ii]
                for jj in range(m_bar):
                    for kk in range(m_bar-jj):
                        for ll in range(m_bar-jj-kk):
                            min_weighted_degree = sys.maxsize
                            mm_min_weighted_degree = l
                            # parallelize here
                            for mm in range(l):
                                monomials = G_set[mm].monomials()
                                delta[mm] = 0
                                for monomial in monomials:
                                    degx = monomial.degree(x)
                                    degy = monomial.degree(y)
                                    degz = monomial.degree(z)
                                    if degx<jj or degy<kk or degz<ll:
                                        continue
                                    g_mm = G_set[mm].monomial_coefficient(monomial)
                                    delta[mm] = delta[mm] + mybinomial(degx,jj) * \
                                                            mybinomial(degy,kk) * \
                                                            mybinomial(degz,ll) * \
                                                            g_mm * \
                                                            L[ii]                                                        ^ (degx-jj) * \
                                                            (elements_in_based_field[0]*column_multiplier_list[ii]^(-1)) ^ (degy-kk) * \
                                                            (elements_in_based_field[1]*column_multiplier_list[ii]^(-1)) ^ (degz-ll)
                                if delta[mm] != 0 and weighted_degree(G_set[mm].lm()) < min_weighted_degree:
                                    mm_min_weighted_degree = mm
                                    min_weighted_degree = weighted_degree(G_set[mm].lm())
                            if mm_min_weighted_degree == l:
                                continue
                            for mm in range(l):
                                if mm != mm_min_weighted_degree:
                                    G_set[mm] = G_set[mm] - delta[mm]/delta[mm_min_weighted_degree]*G_set[mm_min_weighted_degree]
                            G_set[mm_min_weighted_degree] = (x-L[ii]) * G_set[mm_min_weighted_degree]
    G_set.sort(key = lambda G: weighted_degree(G.lm()))
    
    # assert RR(m*(k-1)) < weighted_degree(G_set[0].lm()) and weighted_degree(G_set[0].lm()) < RR(pow(n^3*R^2*m*(m+1)*(m+2), 1/3)), "DELTA1 degree constraint"
    # assert weighted_degree(G_set[1].lm()) <= RR(weighted_degree(G_set[0].lm())/2*(1+sqrt(-1/3+4/3*R^2*n^3*m*(m+1)*(m+2)/weighted_degree(G_set[0].lm())^3))), "DELTA2 degree constraint 1"
    # assert weighted_degree(G_set[1].lm()) < RR(m*n*sqrt(R)/2*(1+sqrt(-1/3+4/3/R*m*(m+1)*(m+2)/m^3))), "DELTA2 degree constraint 2"
    # assert RR(weighted_degree(G_set[0].lm())/2*(1+sqrt(4/3*R^2*(n*m/weighted_degree(G_set[0].lm()))^3-1/3))) <= RR(n*m*R/2*(1+sqrt(4/3/R-1/3))), "Theorem"
    # if debug:
        # xx = var("xx")
        # fig = plot(xx/2*(1+sqrt(4/3*R^2*(n*m/xx)^3-1/3)), (xx, n*m*R-2, weighted_degree(G_set[0].lm())+2))
        # show(fig)

    #check
    if debug:
        print("\n")
        print("checking")
        iteration_count = RR(1)
        total_iteration = RR(l*n)
        for Gxyz in G_set:
            print(Gxyz)
            print("wdeg:", weighted_degree(Gxyz.lm()))
            if weighted_degree(Gxyz.lm()) < DELTA:
                #assert Gxyz(x,p1,p2) == 0, "Gxyz(x,p1,p2) != 0"
                if Gxyz(x,p1,p2) != 0:
                    print("alert!!! Gxyz(x,p1,p2) != 0")
            else:
                print("wdeg > polynomial weighted degree constraint")
            for ii in range(n):
                print("progress: {:.2%}".format(iteration_count/total_iteration), end = '\r')
                iteration_count = iteration_count + 1
                #print(Gxyz(L[ii],r1[ii],r2[ii]))
                for elements_in_based_field in [(element1, element2) for element1 in F_Goppa for element2 in F_Goppa]:
                    if elements_in_based_field[0] == r1[ii] and elements_in_based_field[1] == r2[ii]:
                        for jj in range(m):
                            for kk in range(m-jj):
                                for ll in range(m-jj-kk):
                                    hasse_derivative = 0
                                    for monomial in Gxyz.monomials():
                                        degx = monomial.degree(x)
                                        degy = monomial.degree(y)
                                        degz = monomial.degree(z)
                                        if degx<jj or degy<kk or degz<ll:
                                            continue
                                        hasse_derivative = hasse_derivative + mybinomial(degx,jj) * \
                                                                              mybinomial(degy,kk) * \
                                                                              mybinomial(degz,ll) * \
                                                                              Gxyz.monomial_coefficient(monomial) * \
                                                                              L[ii]                                                        ^ (degx-jj) * \
                                                                              (elements_in_based_field[0]*column_multiplier_list[ii]^(-1)) ^ (degy-kk) * \
                                                                              (elements_in_based_field[1]*column_multiplier_list[ii]^(-1)) ^ (degz-ll)
                                    if hasse_derivative != 0:
                                        print(elements_in_based_field, m, ii, jj, kk, ll)
                                        print(hasse_derivative)
                                    assert hasse_derivative == 0, "not fulfill interpolation constraints" 
                    else:
                        for jj in range(m_bar):
                            for kk in range(m_bar-jj):
                                for ll in range(m_bar-jj-kk):
                                    hasse_derivative = 0
                                    for monomial in Gxyz.monomials():
                                        degx = monomial.degree(x)
                                        degy = monomial.degree(y)
                                        degz = monomial.degree(z)
                                        if degx<jj or degy<kk or degz<ll:
                                            continue
                                        hasse_derivative = hasse_derivative + mybinomial(degx,jj) * \
                                                                              mybinomial(degy,kk) * \
                                                                              mybinomial(degz,ll) * \
                                                                              Gxyz.monomial_coefficient(monomial) * \
                                                                              L[ii]                                                        ^ (degx-jj) * \
                                                                              (elements_in_based_field[0]*column_multiplier_list[ii]^(-1)) ^ (degy-kk) * \
                                                                              (elements_in_based_field[1]*column_multiplier_list[ii]^(-1)) ^ (degz-ll)
                                    if hasse_derivative != 0:
                                        print(elements_in_based_field, m_bar, ii, jj, kk, ll)
                                        print(hasse_derivative)
                                    assert hasse_derivative == 0, "not fulfill interpolation constraints" 
            print("\n")
      
    #Resultant
    output = []
    Qxyz = G_set[0]
    Pxyz = G_set[1]
    PHIxyz = Rxyz(0)
    U = Rxyz(0)
    V = Rxyz(0)
    Hxy = Rxyz(0)
    Hxz = Rxyz(0)
    Hxy_factors = []
    Hxz_factors = []
    Hxy_factor = Rxyz(0)
    Hxz_factor = Rxyz(0)
    fx = Rxyz(0)
    gx = Rxyz(0)

    for ii in range(2,l):
        PHIxyz = Qxyz.gcd(Pxyz)
        if PHIxyz.degree(y) == 0 and PHIxyz.degree(z) == 0:
            Hxy = Qxyz.resultant(Pxyz,z)
            Hxy_factors = list(Hxy.factor())
            for Hxy_factor in Hxy_factors:
                if Hxy_factor[0].degree(y)==1:
                    fx = Hxy_factor[0] - Rxyz(y)
                    if fx.degree(x) < k_GRS:
                        Hxz = Qxyz.resultant(Pxyz,y)
                        Hxz_factors = list(Hxz.factor())
                        for Hxz_factor in Hxz_factors:
                            if Hxz_factor[0].degree(z)==1:
                                gx = Hxz_factor[0] - Rxyz(z)
                                if gx.degree(x) < k_GRS:
                                    if (fx, gx) not in output:
                                        output.append((fx, gx))
            break
            # else:
                # Pxyz = G_set[ii]
                # continue
        else:
            U = Rxyz(Qxyz/PHIxyz)
            V = Rxyz(Pxyz/PHIxyz)
            Hxy = U.resultant(V,z)
            Hxy_factors = list(Hxy.factor())
            for Hxy_factor in Hxy_factors:
                if Hxy_factor[0].degree(y)==1:
                    fx = Hxy_factor[0] - Rxyz(y)
                    if fx.degree(x) < k_GRS:
                        Hxz = U.resultant(V,y)
                        Hxz_factors = list(Hxz.factor())
                        for Hxz_factor in Hxz_factors:
                            if Hxz_factor[0].degree(z)==1:
                                gx = Hxz_factor[0] - Rxyz(z)
                                if gx.degree(x) < k_GRS:
                                    if (fx, gx) not in output:
                                        output.append((fx, gx))
        Qxyz = PHIxyz
        Pxyz = G_set[ii]

    print("\n")
    print("Returned list:")
    print(output)
    print("\n")
    print("(p1,p2):")
    print([(Rxyz(p1), Rxyz(p2))])
    print("(f(X), g(x)) in output:", (Rxyz(p1), Rxyz(p2)) in output)
    if (Rxyz(p1), Rxyz(p2)) in output:
        count_success = count_success + 1
        if len(output) == 1:
            count_size1 = count_size1 + 1
    # else:
        # break
    count_total = count_total + 1
    print("successful decoding rate:", count_success, '/', count_total)
    print("size1 rate:", count_size1, '/', count_total)
    print("error_count:", error_count, "<=", RR((n-n*R_GRS/2*(1+sqrt(4/3/R_GRS+4/3/R_GRS*(q_Goppa^2-1)*(m_bar/m)^3-1/3)))/(1-m_bar/m)))

    deltaset = []
    for Gxyz in G_set[0:2]:
        terms = list(Gxyz)
        for term in terms:
            if term[1] not in deltaset:
                deltaset.append(term[1])
    print("deltaset size:", len(deltaset))
    print("size of ideal:", n*m*(m+1)*(m+2)/6 + n*(q_Goppa^2-1)*m_bar*(m_bar+1)*(m_bar+2)/6)
    print("wdeg G_set[0](X,Y,Z):", weighted_degree(G_set[0].lm()))
    print("G_set[0](X,p1,p2) =", G_set[0](x,p1,p2))
    print("wdeg G_set[1](X,Y,Z):", weighted_degree(G_set[1].lm()))
    print("G_set[1](X,p1,p2) =", G_set[1](x,p1,p2))
    for ii in range(l):
        if G_set[ii+1](x,p1,p2) != 0:
            break
    print("ii =", ii)
    print("wdeg G_set[ii](X,Y,Z):", weighted_degree(G_set[ii].lm()))
    print("G_set[ii](X,p1,p2) =", G_set[ii](x,p1,p2))
    print("wdeg G_set[l-1](X,Y,Z):", weighted_degree(G_set[l-1].lm()))
    print("G_set[l-1](X,p1,p2) =", G_set[l-1](x,p1,p2))
    print("DELTA =", DELTA)
    print("m*(n-t) + m_bar*t =", m*(n-t)+m_bar*t)
    print("m*n*R =", m*n*R_GRS , "==", m*(k_GRS-1))
    print("m*n*R^(2/3)", RR(m*n*R_GRS^(2/3)))
    print("m*n*R^(1/2)", RR(m*n*R_GRS^(1/2)))
    print("t<=(n-(n*(k_GRS-1)^2+n*(q_Goppa^2-1)*(k_GRS-1)^2*(m_bar/m)^3)^(1/3))/(1-(m_bar/m)) =", (n-(n*(k_GRS-1)^2+n*(q_Goppa^2-1)*(k_GRS-1)^2*(m_bar/m)^3)^(1/3))/(1-(m_bar/m)))
    
    print("G_1(X,Y,Z) =",G_set[0])
    print("G_2(X,Y,Z) =",G_set[1])

after = time.time()

if input_t>0:
    f = open('output_ABC_'+str(t)+'_'+str(input_m)+'_'+str(input_m_bar)+'.txt', 'a')
    f.write("We choose n = 32, r = 6, t = " + str(t) + ", m = " +  str(m) + ", m_bar = " + str(m_bar) + "   " + str(count_success) + "/" + str(count_total) + ", time = " + str(int(after-before)) + "\n")
    f.close()
exit