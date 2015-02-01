## Starts Laurent series in q
R.<q> = LaurentSeriesRing(QQ)

I = CC.0 # imaginary unit

precision = 75

def parse_characters():
	fid = open('thompson_characters.txt','r')
	rep1 = []
	rep2 = []
	line = fid.readline()
	while(line != ''):
		xy = line.split(' ')
		rep1.append(int(xy[0]))
		rep2.append(int(xy[1]))
		line = fid.readline()
	return (rep1,rep2)
# 	
# def split_nonzeros(l):
# 	r = [[]]
# 	for i in l:
# 		if i > 0:
# 			r[-1].append(i)
# 			r.append([])
# 		else:
# 			r[-1].append(i)
# 	if not r[-1]:
# 		return r[:-1]
# 	else:
# 		return r
# 
# def generate_cross_field(n):
# 	if n == 1:
# 		result = [[1],[-1]]
# 	if n > 1:
# 		temp = generate_cross_field(n-1)
# 		result = []
# 		for l in temp:
# 			result.append(l + [1])
# 			result.append(l +[-1])
# 	return result
# 	 	
# def action_cross_field(cf, split_v):
# 	result = []
# 	for i in range(len(cf)):
# 		result.extend([x*cf[i] for x in split_v[i]])
# 	return result
# 	
# def action_cross_field2(cf, v):
# 	le = len(v)
# 	i = 0
# 	j = 0
# 	result = v
# 	while i < len(v):
# 		if v[i] != 0:
# 			result[i] = result[i]*cf[j]
# 			j = j + 1
# 		i = i + 1
# 	return result
# 			
# 		
# def return_integer_vectors(summing, length):
# 	vecs = IntegerVectors(summing, max_length=length).list()
# 	cross_fields = []
# 	result = []
# 	for i in range(1,summing+1):
# 		cross_fields.append(generate_cross_field(i))
# 	for v in vecs:
# 		split_v = split_nonzeros(v)
# 		for i in cross_fields[len(split_v)-1]:
# 			result.append(action_cross_field(i, split_v))
# 	return result
# 
# def return_integer_vectors3(summing, length):
# 	if length == 0:
# 		if summing == 0: return [[]]
# 		else: return []
# 	else:
# 		r = []
# 		for x in range(-summing,summing+1):
# 			for vec in return_integer_vectors3(summing-abs(x), length-1):
# 				vec.append(x)
# 				r.append(vec)
# 		return r
# 
# def integer_vectors_gen(summing, length):
# 	if length == 0:
# 		if summing == 0:
# 			yield []
# 	else:
# 		for x in range(-summing,summing+1):
# 			for vec in return_integer_vectors3(summing-abs(x), length-1):
# 				vec.append(x)
# 				yield vec
# 
# def action_vector_characters(int_vector,char_list):
# 	return sum([int_vector[i]*char_list[i] for i in range(len(int_vector))])
# 	
# def search_decomposition(int_vectors,char_list,sum_to):
# 	result = []
# 	for v in int_vectors:
# 		if action_vector_characters(v,char_list) == sum_to:
# 			result.append(v)
# 	return result
# 
# def search_decomposition(int_vectors,char_list,sum_to):
# 	for v in int_vectors:
# 		if action_vector_characters(v,char_list) == sum_to:
# 			yield v
	
def evaluate(f,z):
	result = 0
	coeffs = f.coefficients()
	exps   = f.exponents()
	for i in range(0,len(coeffs)):
		result = result + coeffs[i]*z^(exps[i])
	return result
	
def action(gamma,tau):
	return ((gamma[0]*tau + gamma[1])/(gamma[2]*tau + gamma[3]))
	
def moonshine_1(f,gamma):
	I = CC.0
	epsilon = .0001
	num_terms = len(f.exponents())
	tau = -gamma[3]/gamma[2] + I/gamma[2]
	tau_p = action(gamma,tau)
	q_tau = exp(I*2*pi*tau)
	q_tau_p = exp(I*2*pi*tau_p)
	theta = tau + epsilon + epsilon*I
	theta_p = action(gamma,theta)
	q_theta = exp(I*2*pi*theta)
	q_theta_p = exp(I*2*pi*theta_p)
	print(abs(q_theta).n())
	print(abs(q_theta_p).n())
	print(abs(q_tau).n())
	print(abs(q_tau_p).n())
	w = []
	for i in range(1,num_terms):
		g = f.add_bigoh(f.exponents()[i])
		numerator = log(abs( (evaluate(g,q_tau_p)*evaluate(g,q_theta))/(evaluate(g,q_tau)*evaluate(g,q_theta_p)) ))
		denominator = log(abs( (gamma[2]*tau+gamma[3])/(gamma[2]*theta+gamma[3]) ))
		w = (numerator/denominator).n()
		print(w)
	return w
	
def moonshine_2(f,gamma):
	I = CC.0
	epsilon = .001
	num_terms = len(f.exponents())
	tau = -gamma[3]/gamma[2] + I/gamma[2] + I*epsilon
	tau_p = action(gamma,tau)
	q_tau = exp(I*2*pi*tau)
	q_tau_p = exp(I*2*pi*tau_p)
	print(abs(q_tau).n())
	print(abs(q_tau_p).n())
	w = []
	for i in range(1,num_terms):
		g = f.add_bigoh(f.exponents()[i])
		numerator = log(abs( (evaluate(g,q_tau_p))/(evaluate(g,q_tau)) ))
		denominator = log(abs( (gamma[2]*tau+gamma[3]) ))
		print(numerator.n()/denominator.n())
	return w
	
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
	print(evaluate(g,q_tau_p).n())
	print( ((gamma[2]*tau+gamma[3])^(weight))*evaluate(g,q_tau).n() )
	temp = abs( evaluate(g,q_tau_p) ).n() - abs(((gamma[2]*tau+gamma[3])^(weight))*evaluate(g,q_tau) ).n()
	w.append(temp.n())
	return w

def pp_c(c):
	string = "$ "
	if real_part(c) != 0:
		if abs(real_part(c)) > .1:
			string = string + ("%.3f" % real_part(c))
		else:
			flo = "%.3e" % real_part(c)
			flo = flo.split('e')
			string = string + "\\left(" + flo[0] + "\\times 10^{" + flo[1] + "}\\right)"
	if imag_part(c) != 0:
		if (real_part(c) != 0) and imag_part(c) > 0:
			string = string + " + "
		if abs(imag_part(c)) > .1:
			string = string + ("%.3f" % imag_part(c)) + "i"
		else:
			flo = "%.3e" % imag_part(c)
			flo = flo.split('e')
			string = string + "\\left(" + flo[0] + "\\times 10^{" + flo[1] + "}\\right)i"
	return string + " $"
	
	
def moonshine_3_latex(func,gamma,weight,latex_header,dorder):
	if dorder == 0:
		num_terms = len(func.exponents())
	else:
		num_terms = min(dorder+1, len(func.exponents()))
	tau = -(gamma[3]/gamma[2]) + I*(1/gamma[2]).n()
	tau_p = action(gamma,tau).n()
	q_tau = qt(tau).n()
	q_tau_p = qt(tau_p).n()
	string = "\\begin{tabular}{ | c | c | c | } \\hline \\multicolumn{3}{ |c| }{ $|q(\\tau)| \\approx $ " + pp_c(abs(q_tau).n()) + " $, |q(\\gamma_{" + str(gamma[2]) + "}\\tau)| \\approx $ " + pp_c(abs(q_tau_p).n()) + "} \\\\ \\hline\\hline Order & " + latex_header + "\\\\ \\hline\\hline "
	for i in range(1,num_terms):
		g = func.add_bigoh(func.exponents()[i])
		string = string + str(func.exponents()[i-1]) + " & " + pp_c(evaluate(g,q_tau_p).n()) + " & " + pp_c((gamma[2]*tau+gamma[3])^(weight)*evaluate(g,q_tau).n()) + "\\\\ "
	g = func.add_bigoh(func.exponents()[num_terms-1] + 1)
	string = string + str(func.exponents()[num_terms-1]) + " & " + pp_c(evaluate(g,q_tau_p).n()) + " & " + pp_c((gamma[2]*tau+gamma[3])^(weight)*evaluate(g,q_tau).n()) + "\\\\ "	
	string = string + "\\hline \\end{tabular}"
	print(string)	
	return(string)

## gives you q which corresponds to tau	
def qt(tau):
	return exp(2*pi*I*tau)
	
	
def moonshine_4(func,gamma,weight):
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
		print(evaluate(g,q_tau_p).n())
		print( ((gamma[2]*tau+gamma[3])^(weight))*evaluate(g,q_tau).n() )
		temp = (( evaluate(g,q_tau_p) ).n())/((((gamma[2]*tau+gamma[3])^(weight))*evaluate(g,q_tau) ).n())
		w.append(temp.n())
	g = func.add_bigoh(func.exponents()[num_terms-1]+1)
	print(func.exponents()[num_terms-1])
	print(evaluate(g,q_tau_p).n())
	print( ((gamma[2]*tau+gamma[3])^(weight))*evaluate(g,q_tau).n() )
	temp = (( evaluate(g,q_tau_p) ).n())/((((gamma[2]*tau+gamma[3])^(weight))*evaluate(g,q_tau) ).n())
	w.append(temp.n())
	return w
	
	
def temp_func(f,gamma,tau):
	tau_p = action(gamma,tau)
	q_tau = exp(I*2*pi*tau)
	q_tau_p = exp(I*2*pi*tau_p)
	return abs(evaluate(f,q_tau_p)*((gamma[2]*tau+gamma[3])^(-1/2)) - evaluate(f,q_tau)).n()
	
def construct_grid(m,N,p):
	assert(N in [1,2,3,5,6])
	leading_coeffs = [[3],[4,7],[3,8,11],[4,11,15,16,19],[8,12,15,20,23]]
	mapping = {1:0, 2:1, 3:2, 5:3, 6:4}
	seeds = leading_coeffs[mapping[N]]
	
	# Constructs theta func, f_0
	t = theta_fun(p)
	
	# Sets f_0 = theta
	f = [t]
	
	# Constructs delta function
	d = CuspForms(1,12).0
	d = d.q_expansion(p)
	ds = shift_powers(d, 4*N)
	
	# Constructs shifted Eisenstein series of various weights
	e = []
	w = [4,6,8,10]
	for n in [1,2,3,4]:
		t = eis(2+2*n, p)
		e.append(shift_powers(t,4*N))
	ed = {4:e[0], 6:e[1], 8:e[2], 10:e[3]}
	rcb = []
	for n in [1,2,3,4]:
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

p = construct_grid(3,1,200)
f3 = p[1]
p = construct_grid(3,3,500)
f33 = p[1]		


# f2A = q^(-3) + 8*q + 128*q^4 +21*q^5 - 768*q^8 + 8*q^9 + 3584*q^12 -13312*q^16
# f2A = f2A.add_bigoh(17)
# 
# f3B = q^(-3) - 5*q + 22*q^4 + 27*q^5 - 54*q^8 + 3*q^9 + 6*q^12 +346*q^16
# f3B = f3B.add_bigoh(17)
# 
# f3C = q^(-3) - 14*q + 4*q^4 - 27*q^5 - 54*q^8 
# f3C = f3C.add_bigoh(9)
# 
# f5A = q^(-3) + 2*q + 2*q^4 + 5*q^5 + 14*q^8 + 2*q^9 -4*q^12
# 
# f7A = q^(-3) - 3*q - 2*q^4 + 0*q^5 + 6*q^8 + 5*q^9 + 0*q^12




#int_vectors = []
#for i in range(1,5):
#	int_vectors.extend(return_integer_vectors(i,48))
	
	
char_lists = parse_characters()
char_list = char_lists[0]

#potentials = search_decomposition(int_vectors,char_list,708938752)



# IntegerPartitions[91951146, {8}, 
#  List[1, 248, 4123, 27000, 30628, 30875, 61256, 85995, 147250, 767637,
#    767637, 779247, 957125, 1707264, 2450240, 2572752, 3376737, 
#   4096000, 4123000, 4881384, 4936750, 6669000, 10822875, 11577384, 
#   16539120, 18154500, 21326760, 28861000, 30507008, 40199250, 
#   44330496, 51684750, 72925515, 76271625, 77376000, 81153009, 
#   91171899, 111321000, 
#   190373976, -1, -248, -4123, -27000, -30628, -30875, -61256, -85995, \
# -147250, -767637, -779247, -957125, -1707264, -2450240, -2572752, \
# -3376737, -4096000, -4123000, -4881384, -4936750, -6669000, \
# -10822875, -11577384, -16539120, -18154500, -21326760, -28861000, \
# -30507008, -40199250, -44330496, -51684750, -72925515, -76271625, \
# -77376000, -81153009, -91171899, -111321000, -190373976]]



