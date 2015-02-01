#jaharvey@gmail.com

## Given the q-expansion of p(tau), returns the q-expansion of p(n tau)
def shift_powers(p,n):	
	precision = p.prec()
	coef = p.list()
	exps = p.exponents()
	if exps[0] > 0:
		exps = range(0,exps[-1] + 1)
	else:
		exps = range(exps[0],exps[-1] + 1)
	pnew = 0
	for i in range(0,len(exps)):
		pnew = pnew + coef[i]*q^(n*exps[i])
	return pnew.add_bigoh(precision*n)


## Produce Eisenstein series with specified weight and precision of q-expansion
def eis(weight,precision):
	t = EisensteinForms(1,weight)
	t.set_precision(precision)
	t = t.eisenstein_series()
	e = t[0].q_expansion()
	return e*(1/e.list()[0])

## Take the 'q-derivative' of q-expansion, in other words take q d/dq p(q)
def q_deriv_(p):
	coef = p.coefficients()
	exps = p.exponents()
	for i in range(0,len(coef)):
		coef[i] = exps[i]*coef[i]
		pnew = coef[i]*q^(exps[i])
	return pnew.add_bigoh(p.prec())

def q_deriv(p,r):
	for i in range(r):
		p = q_deriv_(p)
	return p

## Takes the order n RC bracket of functions f and g of weight k and l respectively
def RC_bracket(f,g,k,l,n):
	summing_list = [(r,s) for r in range(0,n+1) for s in range(0,n+1) if r+s==n]
	rc = 0
	for i in range(0,len(summing_list)):
		r = summing_list[i][0]
		s = summing_list[i][1]
		rc = rc + (-1)^r*binomial(n+k-1,s)*binomial(n+l-1,r)*q_deriv(f,r)*q_deriv(g,s)
	return rc.add_bigoh(f.common_prec(g))

## Given some form p of weight k, computes RC bracket with the theta function with the 
## appropriate order so that the resulting form after diving by cusp form c of weight m is d, and then 
## eliminates the constant by subtracting off the appropriate multiple of theta	
def gen_moonshine_please(p,c,k,m,d):
	n = (d-k-1/2+m)/2
	u = RC_bracket(theta,p,1/2,k,n)
	u = u/c
	if len(u.exponents()) > 0:
		exp = u.exponents()[0]
		u = u - u.list()[-exp]*theta
		if u.exponents()[0] < 0:
			return u/(u.list()[0])
		else:
			return u
	else:
		return u

## Retrieves the coefficients in the q-expansion of p between n and m inclusive (must satisfy n <= m)
def get_coeffs(p,n,m):
	assert(n<=m)
	if len(p.exponents()) > 0:
		expi = p.exponents()[0]
		expf = p.exponents()[-1]
		if n < expi:
			if m < expi:
				coeffs = [0 for i in range(0,m-n+1)]
				return coeffs
			else:
				coeffs = [0 for i in range(0,expi-n)]
				if m <= expf:
					coeffs.extend(p.list()[0:(m-expi+1)])
					return coeffs
				else:
					coeffs.extend(p.list())
					coeffs.extend([0 for i in range(0,m-expf)])
					return coeffs
		else:
			if m <= expf:
				coeffs = p.list()[n-expi:m-expi+1]
				return coeffs
			else:
				coeffs = p.list()[n-expi:]
				coeffs.extend([0 for i in range(0, m-expf)])
				return coeffs
	else:
		return [0 for i in range(0, m-n+1)]

	
## Creates m grid elements of weight 1/2 weakly holomorphic functions on Gamma_0(4N)
## This function only intended to work for N=1,2,3,5,6
## Still haven't performed construction for N=6 and N=5 has a bug
def construct_grid(m,N,p):
	assert(N in [1,2,3,5,6])
	leading_coeffs = [[3],[4,7],[3,8,11],[4,11,15,16,19],[8,12,15,20,23]]
	mapping = {1:0, 2:1, 3:2, 5:3, 6:4}
	seeds = leading_coeffs[mapping[N]]
	
	# Constructs theta function, f_0
	t = theta_fun(p)
	
	# Sets f_0 = theta
	f = [t]
	
	# Constructs delta function
	d = CuspForms(1,12).0
	d = d.q_expansion(p)
	ds = shift_powers(d, 4*N)
	
	# Constructs shifted Eisenstein series of various weights
	e = []
	w = [4,6,8,10,12]
	for n in [1,2,3,4,5]:
		t = eis(2+2*n, p)
		e.append(shift_powers(t,4*N))
	ed = {4:e[0], 6:e[1], 8:e[2], 10:e[3], 12:e[4]}
	rcb = []
	for n in [0,1,2,3,4]:
		t = RC_bracket(theta,ed[12-2*n],1/2,12-2*n,n)/ds
		exp = t.exponents()[0]
		rcb.append(t - t.list()[-exp]*f[0])
	
	# Computes j function
	g2 = 60*eis(4,p)/240
	j = 1728*g2^3/(27*delta) 
	
	## Creates seeds manually for different N
	if N == 1:
		t = (-1/20)*(RC_bracket(f[0],ed[10],1/2,10,1)/ds + 608*f[0])
		f.append(t)
	
	if N == 2:
		f.append((1/17280)*(72*rcb[0]+20*rcb[1]))
		f.append((1/-17280)*(1152*rcb[0]+80*rcb[1]))
	
	if N == 3:
		f.append((-1/1800)*rcb[0]-(1/5184)*rcb[1] - (1/40320)*rcb[2])
		f.append((3/400)*rcb[0] + (1/432)*rcb[1] + (1/6720)*rcb[2])
		f.append((-3/40)*rcb[0] - (13/1728)*rcb[1] - (1/2688)*rcb[2])
		
	if N == 5:
		f.append((1/11200)*rcb[0] + (7/207360)*rcb[1] + (1/161280)*rcb[2] + (1/1411200)*rcb[3])
		f.append((-2/1575)*rcb[0] - (1/2160)*rcb[1] - (1/13440)*rcb[2] - (1/176400)*rcb[3])
		for n in [1,2,3,4]:
			t = RC_bracket(f[1],ed[12-2*n],1/2,12-2*n,n)/ds
			exp = t.exponents()[0]
			rcb.append(t-t.list()[-exp]*f[0])
		f.append((39/1600)*rcb[4] + (293/138240)*rcb[5] + (1/7168)*rcb[6] + (1/134400)*rcb[7])
		f.append((1/100)*rcb[0] + (169/51840)*rcb[1] + (13/40320)*rcb[2] + (1/50400)*rcb[3])
		f.append((-2/25)*rcb[0] - (61/6480)*rcb[1] - (29/40320)*rcb[2] - (1/25200)*rcb[3])
	
	if N == 6:
		f.append((1/11200)*rcb[0] + (7/207360)*rcb[1] + (1/161280)*rcb[2] + (1/1411200)*rcb[3])
		#f.append((-17/960)*rcb[4] + (-673/622080)*rcb[5] + (-23/483840)*rcb[6] + (-1/604800)*rcb[7])
		f.append((-2/1575)*rcb[0] + (-1/2160)*rcb[1] + (-1/13440)*rcb[2] + (-1/176400)*rcb[3])
		f.append((1/100)*rcb[0] + (169/51840)*rcb[1] + (13/40320)*rcb[2] + (1/50400)*rcb[3])
		f.append((-2/25)*rcb[0] + (-61/6480)*rcb[1] + (-29/40320)*rcb[2] + (-1/25200)*rcb[3])
		
		
		
	
	## Creates remaining grid elements from seeds
	if len(f) > m:
		return f[0:m]
	else:
		j = shift_powers(j,4*N)
		l = len(f)
		leads = []
		for i in range(0,len(f)):
			leads.append(f[i].exponents()[0])
		while l < m:
			t = f[l-N-1]*j
			for i in range(0,len(leads)):
				t = t - f[i]*get_coeffs(t,leads[i],leads[i])[0]
			f.append(t)
			leads.append(f[l].exponents()[0])
			l = l + 1
		return f	
		
##########################################################################################		
## Starts Laurent series in q
R.<q> = LaurentSeriesRing(QQ)

precision = 300

## Computes basis for modular forms of weight 10 on Gamma0(12)
m = ModularForms(Gamma0(12),10,prec = precision)
p = m.basis()

## Basis for cusp forms on Gamma0(12)
l = CuspForms(Gamma0(12),12,prec = precision)
c = l.basis()

i = 1
squares_list =[1]
while i^2 < precision:
	i = i+1
	squares_list.append(i^2)

i = 1
alt_squares_list = [-1]
while i^2 < precision:
	i = i+1
	alt_squares_list.append((-1)^i*i^2)
	
def gen_theta(i,squares_list):
	if i in squares_list:
		return 2
	else:
		return 0
		
def gen_alt_theta(i,alt_squares_list):
	if i in map(abs,alt_squares_list):
		j = map(abs,alt_squares_list).index(i)
		return 2*sign(alt_squares_list[j])
	else:
		return 0
		
	

## Gives you back q_expansion of theta function to precision p
def theta_fun(p):
	thetalist = [1]
	thetalist.extend([gen_theta(i,squares_list) for i in range(1,p)])
	theta = R(thetalist)
	theta = theta.add_bigoh(p)
	return theta
	
def alt_theta_fun(p):
	thetalist = [1]
	thetalist.extend([gen_alt_theta(i,alt_squares_list) for i in range(1,p)])
	theta = R(thetalist)
	theta = theta.add_bigoh(p)
	return theta

## Defining theta function
theta = theta_fun(precision)

## Defining some Eisenstein series
e10 = eis(10,precision)
e4 = eis(4,precision)
e6 = eis(6,precision)

## Defining delta cusp form
delta = CuspForms(1,12).0
delta = delta.q_expansion(precision)

## Verifying construction of f_3,3
u1 = RC_bracket(theta,shift_powers(e10,12),1/2,10,1)
u2 = u1/(shift_powers(delta,12))
u3 = u2 + 1584*theta
u = (-1/20)*u3

v1 = RC_bracket(theta,shift_powers(eis(8,precision),12),1/2,8,2)
v2 = v1/(shift_powers(delta,12))
v3 = v2 - 25920*theta
v = (1/72)*v3

w1 = RC_bracket(theta,shift_powers(e6,12),1/2,6,3)
w2 = w1/(shift_powers(delta,12))
w3 = w2 + 272160*theta
w = (-1/112)*w3

f33 = (4*u - 5*v + w)/360


## gen_moonshine_please dividing by delta(12 tau) as cusp form
#g = []
#for i in range(0,len(p)):
#	g.append(gen_moonshine_please(p[i].q_expansion(),shift_powers(theta,12), 10, 12, 1/2))

## gen_moonshine_please taking all combinations of basis elements of cusp forms on Gamma0(12) of weight 12 
## and basis elements of weight 10 forms on Gamma0(12)
#h = []
#for i in range(0,len(p)):
#	for j in range(0,len(c)):
#		h.append(gen_moonshine_please(p[i].q_expansion(),c[j].q_expansion(),10,12,1/2))
		
## Initializes variables for use in system of equation solving
#var_string = 'a0'
#num_forms = len(h)
#for i in range(1,num_forms):	
#	var_string = var_string + ', a' + str(i)
#a = var(var_string)

#coeffs = []	
#polys = range(0,num_forms)
#for i in polys:
#	coeffs.append(get_coeffs(h[i],-11,len(h)))

#num_constraints = 398
#eq = []
#f33_coeffs = get_coeffs(f33,-11,len(h))
#for i in range(0, num_constraints):
#	eq.append(0)
#	for j in range(0,len(polys)):
#		eq[i] = eq[i] + coeffs[j][i]*a[j]
#	eq[i] = eq[i] == f33_coeffs[i]
		
			

#answer = solve(eq, a)
#result = [0 for i in range(0,len(polys))]
##modular_form = 0
##for i in range(0,len(polys)):
##	modular_form = modular_form + result[i]*s[polys[i]]