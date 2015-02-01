I = CC.0 # imaginary unit

def over_prim_residues(N,li,g_l,g_u,lamb,m,n):
	if g_l != 0:
		di = g_u/g_l
	c = 4*N*g_u
	result = 0
	newli = []
	for v in li:
		if gcd(di,v) == 1:
			newli.append(v)
			vbar = inverse_mod(v,c)
			x = ((m*vbar+n*v)/c).n()
			result = result + ((kronecker(c,v)*(ep(v)^(2*lamb+1)))*exp(2*pi*I*x)).n()
	for v in range(4*N*g_l+1,c):
		if gcd(v,c) == 1:
			newli.append(v)
			vbar = inverse_mod(v,c)
			x = ((m*vbar+n*v)/c).n()
			result = result + ((kronecker(c,v)*(ep(v)^(2*lamb+1)))*exp(2*pi*I*x)).n()
	return (result,newli)
			
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
	
def klooster2(lamb, m, n, c):
	result = 0
	for v in range(1,c): #0 and c don't appear in this sum?
		if is_relatively_prime(v,c):
			vbar = inverse_mod(v,c)
			x = ((m*vbar+n*v)/c).n()
			result = result + ((kronecker(c,v)*(ep(v)^(2*lamb + 1)))*exp(2*pi*I*x)).n()
	return result.n()
	
def kr_delta(n,m):
	if n == m:
		return 1
	else:
		return 0
	
def dodd(c):
	if mod(c,2) == 0:
		return 0
	else:
		return 1
	
def is_relatively_prime(a,b):
	P = Primes()
	p = P.first()
	c = min(a,b)
	while p <= c:
		if p.divides(a) and p.divides(b):
			return False
		p = P.next(p)
	return True

def lowest_nontrivial_prime(N):
	P = Primes()
	p = P.first()
	while True:
		if p.divides(N):
			return p
		p = P.next(p)
	
def b0(N,m,n):
	if n == 0:
		factor = ((pi^(3/2))*2*(m^(1/2))*(1-I)*(1/((1/2)*sqrt(pi)))).n()
		result = 0
		for c in map(lambda g: 4*N*g, range(1,500)):
			result = result + (((1+ dodd(c/4))*klooster(0,-m, 0, c))/(c^(3/2))).n()
		return (result*factor).n()
	if n > 0:
		factor = pi*sqrt(2)*(n/m)^(-(1/4))*(1-I).n()
		result = 0
		for c in map(lambda g:4*N*g, range(1,50)):
			result = result + ((1 + dodd(c/4))*(klooster(0,-m,n,c)/c)*bessel_I(1/2,4*pi*sqrt(m*n)/c)).n()
		return (result*factor).n()
		
def b0_klooster2(N,m,n):
	if n == 0:
		factor = ((pi^(3/2))*2*(m^(1/2))*(1-I)*(1/((1/2)*sqrt(pi)))).n()
		result = 0
		for c in map(lambda g: 4*N*g, range(1,20)):
			result = result + (((1+ dodd(c/4))*klooster2(0,-m, 0, c))/(c^(3/2))).n()
		return (result*factor).n()
	if n > 0:
		factor = pi*sqrt(2)*(n/m)^(-(1/4))*(1-I).n()
		result = 0
		for c in map(lambda g:4*N*g, range(1,20)):
			result = result + ((1 + dodd(c/4))*(klooster2(0,-m,n,c)/c)*bessel_I(1/2,4*pi*sqrt(m*n)/c)).n()
		return (result*factor).n()

def b0_3(N,m,n):
	if n == 0:
		factor = ((pi^(3/2))*2*(m^(1/2))*(1-I)*(1/((1/2)*sqrt(pi)))).n()
		prims = {}
		(kloos, li) = over_prim_residues(N,[],0,1,0,-m,0)
		prims[1] = li
		result = (((1+ dodd(N))*kloos)/((4*N)^(3/2))).n()
		for g in range(2,500):
			c = 4*N*g
			if g not in Primes():	
				di = lowest_nontrivial_prime(g)
				g_l = g/di
				(kloos,li) = over_prim_residues(N,prims[g_l],g_l,g,0,-m,0)
				prims[g] = li
				result = result + (((1+ dodd(c/4))*kloos)/(c^(3/2))).n()
			else:
				g_l = 0
				(kloos,li) = over_prim_residues(N,[],g_l,g,0,-m,0)
				prims[g] = li
				result = result + (((1+dodd(c/4))*kloos)/(c^(3/2))).n()
		return (result*factor).n()
	if n > 0:
		factor = pi*sqrt(2)*(n/m)^(-(1/4))*(1-I).n()
		result = 0
		for c in map(lambda g:4*N*g, range(1,20)):
			result = result + ((1 + dodd(c/4))*(klooster2(0,-m,n,c)/c)*bessel_I(1/2,4*pi*sqrt(m*n)/c)).n()
		return (result*factor).n()
