I = CC.0 # imaginary unit
R.<q> = LaurentSeriesRing(QQ)

def pair(z):
	w = floor((sqrt(8*z+1) - 1)/2)
	t = (w^2 + w)/2
	y = z-t
	x = w-y
	return (x,y)

T = matrix(2,[1,1,0,1])
S = matrix(2,[0,-1,1,0])
id = matrix(2,[1,0,0,1])

def reduce_form(Q,N):
	a = Q[0]
	b = Q[1]
	c = Q[2]
	if (b^2-4*a*c>=0 or a<0):
		raise Exception("Must be positive definite!")
	M = matrix(2,[1,0,0,1])
	while not Q.is_reduced():
		a = Q[0]
		b = Q[1]
		c = Q[2]
		if c<a:
			M = M*S
			Q = Q.matrix_action_right(S)
		elif abs(b)>a or -b == a:
			k = floor((a-b)/(2*a))
			M = M*(T^k)
			Q = Q.matrix_action_right(T^k)
	if N == 1:
		return Q,M
	else:
		reps = map(lambda g:~g, list(Gamma0(N).coset_reps()))
		for r in reps:
			r = matrix(2,[ r[0,0],r[0,1],r[1,0],r[1,1] ])
			M_ = M*r
			if mod(M_[1,0],N) == 0:
				return Q.matrix_action_right(r),M*r
		else:
			print("ERROR!")

def is_equiv(Q1,Q2,N):
	Q1,M1 = reduce_form(Q1,N)
	Q2,M2 = reduce_form(Q2,N)
	return Q2==Q1

def removeDups(duptest,iterable):
	res = []
	for e in iterable:
		if not any(duptest(e,r) for r in res):
			res.append(e)
	return res

def gen_quad_reps(D,N):
	reps = BinaryQF_reduced_representatives(-abs(D))
	if N == 1:
		return reps
	else:
		cos_reps = map(lambda g:~g,list(Gamma0(N).coset_reps()))
		cos_reps = map(lambda x: matrix(2,[[x[0][0],x[0][1]],[x[1][0],x[1][1]]]), cos_reps)
		new_reps = []
		for q_ in reps:
			for g in cos_reps:
				new_reps.append(q_.matrix_action_right(g))
		return removeDups(lambda Q1,Q2: is_equiv(Q1,Q2,N), new_reps)

def eval_form(q_,p):
	return q_[0]*p[0]^2 + q_[1]*p[0]*p[1] + q_[2]*p[1]^2

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

safety_set = [id,T,T^(-1)]

elliptic_points = {}
n3o3 = ((1/2) + (1/6)*sqrt(-3)).n()
n5o2_1 = ((2/5) + (1/5)*sqrt(-1)).n()
n5o2_2 = ((3/5) + (1/5)*sqrt(-1)).n()
n31o3_1 = ((11/62) + (1/62)*sqrt(-3)).n()
n31o3_2 = ((17/62) + (1/186)*sqrt(-3)).n()
elliptic_points[3] = ([],[n3o3,n3o3+1,n3o3-1])
elliptic_points[5] = ([ n5o2_1,n5o2_1+1,n5o2_1-1,n5o2_2,n5o2_2+1,n5o2_2-1],[])
elliptic_points[31] = ([], [n31o3_1,n31o3_1+1,n31o3_1-1,n31o3_2,n31o3_2+1,n31o3_2-1])


def omega(Q,N):
	tau = Q.complex_point()
	if tau in elliptic_points[N][0]:
		return 2
	if tau in elliptic_points[N][1]:
		return 3
	else:
		return 1

## should check these are all right
t_2b = lambda tau: (eta(tau)^24)/(eta(2*tau)^24) + 24
t_3b = lambda tau: (eta(tau)^(12))/(eta(3*tau)^12) + 12
t_4c = lambda tau: 1/((eta(4*tau)^8)/(eta(tau)^8))
t_6e = lambda tau: (((eta(tau )^5)*eta(3*tau)))/(((eta(2*tau )*eta(6*tau )^5)))
t_8e = lambda tau: ((eta(tau )^4)*eta(4*tau )^2)/((eta(2*tau )^2)*(eta(8*tau )^4))
t_9b = lambda tau: ((eta(tau )^3)/(eta(9*tau )^3))
t_25z  = lambda tau: (eta(tau ))/(eta(25*tau )) + 1
t_25z2 = lambda tau: (eta_2(tau))/(eta_2(25*tau)) + 1
t_5b = lambda tau: (eta(tau )^6)/(eta(5*tau )^6) + 6
t_7b = lambda tau: (eta(tau )^4)/(eta(7*tau )^4) + 4
t_13b = lambda tau: (eta(tau )^2)/(eta(13*tau )^2) + 2
t_10e = lambda tau: 1/((eta(2*tau )*(eta(10*tau )^3))/((eta(tau )^3)*(eta(5*tau ))))
t_12i = lambda tau: ( (eta(tau )^3)*(eta(4*tau ))*(eta(6*tau )^2) ) / (((eta(2*tau )^2)*(eta(3*tau ))*(eta(12*tau )^3)) )
t_18d = lambda tau: ( (eta(tau )^2)*(eta(6*tau ))*(eta(9*tau )) ) / ( (eta(2*tau ))*(eta(3*tau ))*(eta(18*tau )^2) )
t_39 = lambda tau: ( (eta(3*tau )*eta(13*tau ))/(eta(tau )*eta(39*tau )) )


def tr(N,D1,f,D2):
	qs = gen_quad_reps(D1*D2, N)
	trace = 0 
	for q_ in qs:
		if mod(q_[0],N) == 0:
			tau = q_.complex_point()
			trace = trace + (chi(q_,D1)*f(tau))/omega(q_,N)
	return trace
