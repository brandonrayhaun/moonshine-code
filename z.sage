## implementing corollary 1.3
precision = 100
I = CC.0 # imaginary unit

def pair(z):
	w = floor((sqrt(8*z+1) - 1)/2)
	t = (w^2 + w)/2
	y = z-t
	x = w-y
	return (x,y)
	
def eval_form(q_,p):
	return q_[0]*p[0]^2 + q_[1]*p[0]*p[1] + q_[2]*p[1]^2
	
## Starts Laurent series in q
R.<q> = LaurentSeriesRing(QQ)

## Take the 'q-derivative' of q-expansion, in other words take q d/dq p(q)
def q_deriv(p):
	coef = p.coefficients()
	exps = p.exponents()
	pnew = 0
	for i in range(len(coef)):
		coef[i] = exps[i]*coef[i]
		pnew = pnew + coef[i]*q^(exps[i])
	return pnew.add_bigoh(p.prec())
	
def in_quads(q_,qs):
	if len(qs)==0:
		return False
	for k in qs:
		if k.is_equivalent(q_):
			return True
	return False
	
def uniquify(qs):
	new_qs = []
	for q_ in qs:
		if not in_quads(q_,new_qs):
			new_qs.append(q_)
	return new_qs

def tpose(g):
	l=[[0,0],[0,0]]
	l[0][0] = g[0][0]
	l[1][1] = g[1][1]
	l[0][1] = g[0][1]
	l[1][0] = g[1][0]
	return l

def ac(gamma, tau):
	return ((gamma[0,0]*tau + gamma[0,1])/(gamma[1,0]*tau + gamma[1,1]))
	
# generates representatives of Q_D/Gamma_0(N)
# potential issue-- gen_quad_reps(-3,3) has two *equal* representatives, fixed by taking unique values
def gen_quad_reps(D,N):
	reps = uniquify(BinaryQF_reduced_representatives(-abs(D)))
	if N==1:
		return reps
	else:
		cos_reps = list(Gamma0(N).coset_reps())
		cos_reps = map(lambda g:~g,cos_reps)
		#cos_reps = map(lambda g:tpose(g), cos_reps)
		cos_reps = map(lambda x: matrix(ZZ,[[x[0][0],x[0][1]],[x[1][0],x[1][1]]]), cos_reps)
		new_reps = []
		for q_ in reps:
			for g in cos_reps:
				new_reps.append(q_.matrix_action_right(g))
		return new_reps

def omega(q_):
	if q_.is_equivalent( BinaryQF([q_[0],0,q_[0]]) ):
		return 2
	if q_.is_equivalent( BinaryQF([q_[0],q_[0],q_[0]]) ):
		return 3
	else:
		return 1
		
def alpha(q_):
	return q_.complex_point()
	
def chi(q_,D1):
	D = q_.discriminant()
	if not floor(D/D1)==D/D1:
		raise NameError('not well defined!')
	if gcd(D1,gcd(q_[0],gcd(q_[1],q_[2]))) != 1:
		return 0
	else:
		z = 2
		while True:
			x = eval_form(q_, pair(z))
			if gcd(x,D1) == 1:
				return kronecker(D1,x)
			z = z + 1

def qt(tau):
	return exp(2*pi*I*tau)

def evaluate(f,z):
	result = 0
	coeffs = f.coefficients()
	exps   = f.exponents()
	for i in range(len(coeffs)):
		result = result + coeffs[i]*z^(exps[i])
	return result

def eta(tau,p):
	result = 1
	q_ = exp(2*pi*I*tau).n()
	for n in range(1,p+1):
		result = result*(1-q_^n)
	return ((exp((pi*I*tau)/12)).n())*result
	
def eta_2(tau,p):
	result = 0
	q_ = exp(2*pi*I*tau).n()
	for n in range(-p,p+1):
		result = result + ((-1)^n)*q_^((n*(3*n-1))/2)
	return ((q_^(1/24)).n())*result.n()
	
def example_2(p):
	q_inf = 1
	q_inf_25 = 1
	for i in range(1,p+1):
		q_inf = q_inf*(1-q^i)
		q_inf_25 = q_inf_25*(1-q^(25*i))
	return q^(-1)*(q_inf/q_inf_25) + 1
	
def tr_eta(N,D1,D2):
	qs = gen_quad_reps(D1*D2, N)
	f = lambda tau: (eta(tau,1000)/eta(25*tau,1000) + 1).n()
	trace = 0
	for q_ in qs:
		if mod(q_[0],N) == 0:
			print(q_)
			print(chi(q_,D1))
			print(alpha(q_))
			print(f(alpha(q_)))
			print(omega(q_))
			trace = trace + (chi(q_,D1)*f(alpha(q_)))/(omega(q_))
	return trace
	
def ep(d):
	if mod(d,4) == 1:
		return 1
	if mod(d,4) == 3:
		return I
	else:
		raise NameError('delta must be 1 or 3 mod 4')
	
def klooster(lamb, m, n, c):
	result = 0
	for v in range(1,c): #0 and c don't appear in this sum?
		if gcd(v,c) == 1:
			vbar = inverse_mod(v,c)
			x = ((m*vbar+n*v)/c).n()
			result = result + ((kronecker(c,v)*(ep(v)^(2*lamb + 1)))*exp(2*pi*I*x)).n()
	return result.n()
	
def kr_delta(n,m):
	if n == m:
		return 1
	else:
		return 0
	
def tr_j(d):
	qs = gen_quad_reps(d, 1)
	trace = 0
	for q_ in qs:
		trace = trace + (chi(q_,1))*((evaluate(j,exp(2*pi*I*alpha(q_))))/(omega(q_)))
	return trace.n()

def dodd(c):
	if mod(c,2) == 0:
		return 0
	else:
		return 1
	
def b0(N,m,n):
	if n == 0:
		factor = ((pi^(3/2))*2*(m^(1/2))*(1-I)*(1/((1/2)*sqrt(pi)))).n()
		result = 0
		for c in map(lambda g: 4*N*g, range(1,150)):
			result = result + (((1+ dodd(c/4))*klooster(0,-m, 0, c))/(c^(3/2))).n()
		return (result*factor).n()
	if n > 0:
		factor = pi*sqrt(2)*(n/m)^(-(1/4))*(1-I).n()
		result = 0
		for c in map(lambda g:4*N*g, range(1,150)):
			result = result + ((1 + dodd(c/4))*(klooster(0,-m,n,c)/c)*bessel_I(1/2,4*pi*sqrt(m*n)/c)).n()
		return (result*factor).n()
			

def b_minus(lamb,N,m,n):
	m = abs(m)
	summation = 0
	for l in range(1,1000):
		c = 4*N*l
		summation = summation + (1+dodd(c/4))*(klooster(lamb,-m,-n,c)/c)*bessel_J(1/2 - lamb, (4*pi*sqrt(m*n))/c)
	result = ((-1)^(floor((lamb+1)/2)))*pi*sqrt(2)*((n/m)^((lamb/2)-(1/4)))*(1-(-1)^(lamb)*I)
	result = result*summation
	result = result - kr_delta(n,m)
	return result

def tr(N,D1,f,D2):
	qs = gen_quad_reps(D1*D2, N)
	trace = 0
	for q_ in qs:
		if mod(q_[0],N) == 0:
			trace = trace + (chi(q_,D1)*evaluate(f,exp(2*pi*I*alpha(q_))))/(omega(q_.reduced_form()))
	return trace
	
def tr2(N,D1,f,D2):
	qs = gen_quad_reps(D1*D2, N)
	trace = 0
	for q_ in qs:
		if mod(q_[0],N) == 0:
			trace = trace + (chi(q_,D1)*f(alpha(q_)))/(omega(q_.reduced_form()))
	return trace
		
	
t_2b = lambda tau: (eta(tau,10000)^24)/(eta(2*tau,10000)^24) + 24
t_3b = lambda tau: (eta(tau,10000)^(12))/(eta(3*tau,10000)^12) + 12
t_4c = lambda tau: 1/((eta(4*tau,10000)^8)/(eta(tau,10000)^8))
t_6e = lambda tau: (((eta(tau,10000)^5)*eta(3*tau,10000)))/(((eta(2*tau,10000)*eta(6*tau,10000)^5)))
t_8e = lambda tau: ((eta(tau,10000)^4)*eta(4*tau,10000)^2)/((eta(2*tau,10000)^2)*(eta(8*tau,10000)^4))
t_9b = lambda tau: ((eta(tau,10000)^3)/(eta(9*tau,10000)^3))
t_25z  = lambda tau: (eta(tau,10000))/(eta(25*tau,10000)) + 1
t_25z2 = lambda tau: (eta_2(tau,1000))/(eta_2(25*tau,1000)) + 1
t_5b = lambda tau: (eta(tau,10000)^6)/(eta(5*tau,10000)^6) + 6
t_7b = lambda tau: (eta(tau,10000)^4)/(eta(7*tau,10000)^4) + 4
t_13b = lambda tau: (eta(tau,10000)^2)/(eta(13*tau,10000)^2) + 2
t_10e = lambda tau: 1/((eta(2*tau,10000)*(eta(10*tau,10000)^3))/((eta(tau,10000)^3)*(eta(5*tau,10000))))
t_12i = lambda tau: ( (eta(tau,10000)^3)*(eta(4*tau,10000))*(eta(6*tau,10000)^2) ) / (((eta(2*tau,10000)^2)*(eta(3*tau,10000))*(eta(12*tau,10000)^3)) )
t_18d = lambda tau: ( (eta(tau,10000)^2)*(eta(6*tau,10000))*(eta(9*tau,10000)) ) / ( (eta(2*tau,10000))*(eta(3*tau,10000))*(eta(18*tau,10000)^2) )
t_39 = lambda tau: ( (eta(3*tau,10000)*eta(13*tau,10000))/(eta(tau,10000)*eta(39*tau,10000)) )

#f2a = q^-3 - 4 + 120*q^4 + 21*q^5 - 768*q^8 + 3584*q^12 - 42*q^13 - 13320*q^16 + 43008*q^20 + 155*q^21 - 125440*q^24 + 337920*q^28 - 405*q^29 - 854016*q^32 + 2048120*q^36 + 933*q^37 - 4698624*q^40 + 10375680*q^44 - 1925*q^45 - 22163456*q^48 + 45975552*q^52 + 3849*q^53 - 92912640*q^56 + 183416832*q^60 - 7455*q^61 - 354476040*q^64 + 671956992*q^68 + 13825*q^69 - 1251443968*q^72 + 2293025280*q^76 - 24606*q^77 - 4138745856*q^80 + 7366512640*q^84 + 42813*q^85 - 12942094848*q^88 + 22462973952*q^92 - 73116*q^93 - 38546206720*q^96 + 65440383096*q^100

def evaluate(f,z):
	result = 0
	coeffs = f.coefficients()
	exps   = f.exponents()
	for i in range(0,len(coeffs)):
		result = result + coeffs[i]*z^(exps[i])
	return result
	
## Defining delta cusp form
delta = CuspForms(1,12).0
delta = delta.q_expansion(precision)

## Produce Eisenstein series with specified weight and precision of q-expansion
def eis(weight,precision):
	t = EisensteinForms(1,weight)
	t.set_precision(precision)
	t = t.eisenstein_series()
	e = t[0].q_expansion()
	return e*(1/e.list()[0])
	
# Computes j function
g2 = 60*eis(4,precision)/240
j = 1728*g2^3/(27*delta) 
j = j - 744

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

def coeff(m,N,n,f):
	if m >= 1:
		result = tr(N,-n,f,m)
		result = result/(sqrt(m))
		return result.n()
	else:
		raise NameError('havent implemented singular coefficients')
	
def coeff2(m,N,n,f):
	if m>= 1:
		result = tr2(N,-n,f,m)
		result = result/(sqrt(m))
		return result
	else:
		raise NameError('havent implemented singular coefficients')

#precision=200
#q_inf = 1
#for i in range(1,precision+1):
#	q_inf = q_inf + ((-1)^i)*(q^((i*(3*i+1))/2)) + ((-1)^(-i))*(q^(((-i)*(3*(-i)+1))/2))

#q_inf_25 = shift_powers(q_inf,25)

#eta_f = q^(-1)*(q_inf/q_inf_25) + 1


# things to look into: why do different prescriptions for eta differ?
# right vs left action?		

#prints out latex formatting for a modular form
def print_to_latex(N,n,f,f_name, level,order):
	string = f_name + ' = q^{-3} + c_{' + str(level) + '} '
	for i in range(1,order+1):
		if mod(i,4) == 0 or mod(i,4) == 1:
			co = real_part(coeff2(i,N,n,f).n())
			if abs(imag_part(co)) > 10^(-6):
				raise NameError('careful! nonvanishing complex part on coefficient of q^' + str(i))
			if abs(round(co)-co) < 10^(-8):
				string = string + ' + ' + str(round(co)) + 'q^{' + str(i) + '} '
			else:
				string = string + ' + ' + ("%.4f" % co) + 'q^{' + str(i) + '} '
	print string
	
#prints out latex formatting for a modular form
def print_to_latex_(N,n,f,f_name, level,order):
	string = f_name + ' = q^{-3} + c_{' + str(level) + '} '
	for i in range(1,order+1):
		if mod(i,4) == 0 or mod(i,4) == 1:
			co = real_part(coeff(i,N,n,f).n())
			if abs(imag_part(co)) > 10^(-6):
				raise NameError('careful! nonvanishing complex part on coefficient of q^' + str(i))
			if abs(round(co)-co) < 10^(-8):
				string = string + ' + ' + str(round(co)) + 'q^{' + str(i) + '} '
			else:
				string = string + ' + ' + ("%.4f" % co) + 'q^{' + str(i) + '} '
	print string
	
def print_to_latex2(N,n,f_name, level,order):
	string = f_name + ' = q^{-3}'
	for i in range(0,order+1):
		if mod(i,4) == 0 or mod(i,4) == 1:
			co = b0(N,n,i).n()
			if abs(imag_part(co)) > 10^(-6):
				raise NameError('careful! nonvanishing complex part on coefficient of q^' + str(i))
			co = real_part(co)
			if abs(round(co)-co) < (.15):
				string = string + ' + ' + str(round(co)) + 'q^{' + str(i) + '} '
			else:
				string = string + ' + ' + ("%.4f" % co) + 'q^{' + str(i) + '} '
	print string
	

def p_series(N,n,f,order):
	result = q^(-3) + round(real_part(b0(N,n)))
	for i in range(1,order+1):
		if mod(i,4) == 0 or mod(i,4) == 1:
			co = round(real_part(coeff2(i,N,n,f).n()))
			result = result + co*q^(i)
	return result
	
def qt(tau):
	return exp(2*pi*I*tau)	

def action(gamma,tau):
	return ((gamma[0]*tau + gamma[1])/(gamma[2]*tau + gamma[3]))
	
def moonshine_3(func,gamma,weight):
	num_terms = len(func.exponents())
	tau = -(gamma[3]/gamma[2]) + I*(1/gamma[2]).n()
	tau_p = action(gamma,tau).n()
	q_tau = qt(tau).n()
	q_tau_p = qt(tau_p).n()
	print(abs(q_tau).n())
	print(abs(q_tau_p).n())
	w = []
	for i in range(1,num_terms):
		g = func.add_bigoh(func.exponents()[i])
		print(func.exponents()[i-1])
		print(evaluate(g,q_tau_p).n(),((gamma[2]*tau+gamma[3])^(weight))*evaluate(g,q_tau).n())
		temp = abs( evaluate(g,q_tau_p) ).n() - abs(((gamma[2]*tau+gamma[3])^(weight))*evaluate(g,q_tau) ).n()
		w.append(temp.n())
	g = func.add_bigoh(func.exponents()[num_terms-1]+1)
	print(func.exponents()[num_terms-1])
	print(evaluate(g,q_tau_p).n(),((gamma[2]*tau+gamma[3])^(weight))*evaluate(g,q_tau).n() )
	temp = abs( evaluate(g,q_tau_p) ).n() - abs(((gamma[2]*tau+gamma[3])^(weight))*evaluate(g,q_tau) ).n()
	w.append(temp.n())
	return w

i = 1
result = 0 
po = BinaryQF_reduced_representatives(-i*3)
for p in po:
	result = result + chi(p,-3)*evaluate(j,qt(p.complex_point()))/omega(p)
print(round(real_part(result.n())))