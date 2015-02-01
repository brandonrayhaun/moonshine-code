## Starts Laurent series in q
R.<q> = LaurentSeriesRing(QQ)

I = CC.0 # imaginary unit

precision = 1000

##evaluates a function using its q-expansion
def evaluate(f,z):
	result = 0
	coeffs = f.coefficients()
	exps   = f.exponents()
	for i in range(0,len(coeffs)):
		result = result + coeffs[i]*z^(exps[i])
	return result
	
## computes the action of a member of the modular group on tau in the upper half plane
def action(gamma,tau):
	return ((gamma[0]*tau + gamma[1])/(gamma[2]*tau + gamma[3]))
	
## Produce Eisenstein series with specified weight and precision of q-expansion
def eis(weight,precision):
	t = EisensteinForms(1,weight)
	t.set_precision(precision)
	t = t.eisenstein_series()
	e = t[0].q_expansion()
	return e*(1/e.list()[0])

## gives you q which corresponds to tau	
def qt(tau):
	return exp(2*pi*I*tau)
	
## Defining delta cusp form
delta = CuspForms(1,12).0
delta = delta.q_expansion(precision)

# Computes j function
g2 = 60*eis(4,precision)/240
j = 1728*g2^3/(27*delta) 

tau = (I+1)/11
gamma = [0,-1,1,0]
tau_p = action(gamma,tau)
f_tau = evaluate(j,qt(tau)).n()
f_tau_p = evaluate(j,qt(tau_p)).n()