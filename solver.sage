def get_coeff(f,exp):
	if exp in f.exponents():
		idx = f.exponents().index(exp)
		return f.coefficients()[idx]
	else:
		return 0
		
def dott(vec1,vec2):
	result = 0
	for i in range(len(vec1)):
		result = result + vec1[i]*vec2[i]
	return result
		
def get_equation(objective,fss,variables,exp):
	coeffs = []
	for i in range(len(fss)):
		coeffs.append(get_coeff(fss[i],exp))
	return get_coeff(objective,exp) == dott(coeffs,variables)
		

def moonshinesolver(objective, fss):
	variables = []
	test_exps = [-3,0,1,4,5,8,9,12,13,16,17]
	for i in range(len(fss)):
		variables.append(var('x' + str(i)))
	eqns = []
	for i in range(len(fss)):
		eqns.append(get_equation(objective,fss,variables,test_exps[i]))
	return (eqns,variables)
	
theta = 1 + 2*q + 2*q^4 + 2*q^9 + 2*q^(16)
theta4 = 1 + 2*q^4 + 2*q^(16)
theta9 = 1 + 2*q^9 + 2*q^(36)
f3a = (1/2)*(get_twist(f_3decomps,f_3orders,'3a') - 14*theta)
f3b = (1/2)*(get_twist(f_3decomps,f_3orders,'3b') - 5*theta)
f3c = (1/2)*(get_twist(f_3decomps,f_3orders,'3c') + 4*theta)
f3 = (1/2)*(get_twist(f_3decomps,f_3orders,'1a') - 248*theta)
f6a = 2*q^(-3) + 4 - 6*q^5  - 12*q^8 -8*q^12
f6b = 2*q^(-3) - 2 - 6*q^4 + 16*q^(12)-6*q^(13)+18*q^(16)
f6c = 2*q^(-3) + 1 + 6*q^4 + 6*q^5 + 12*q^8+4*q^(12)+6*q^(13) - 18*q^(16)
f9a = 2*q^(-3) + 5 + 16*q^9 + 12*q^(12)
f9b = 2*q^(-3) -4 -2*q^9 + 12*q^(12)
f9c = 2*q^(-3) + 2 -8*q^9-6*q^(12)
f33 = q^(-3) - 14*q + 40*q^4 - 78*q^9 + 168*q^(12) - 378*q^(13) + 688*q^(16)-897*q^(21)
F6 = q^(-3) + (1/2) + 3*q^4 + 3*q^5 + 6*q^8 + 2*q^(12) + 3*q^(13)-9*q^(16)
F9 = q^(-3) - (1/2) + 2*q^9 + 6*q^(12) - 6*q^(21) -2*q^(24)
F12 = q^(-3) + 3*q^5 + 3*q^(13) +2*q^(21)
F18 = q^(-3) + (1/2) + 2*q^(12)
F21 = q^(-3)
F27 = q^(-3)

f4a = 2*q^(-3) + 8 + 16*q^4 + 42*q^5-84*q^(13)+ 16*q^(16)
f4b = 2*q^(-3) - 22*q^5+108*q^(13)

F4 = q^(-3) + 21*q^5 - 42*q^(13)+155*q^(21)-405*q^(29)
F8 = q^(-3) + 5*q^5 + 6*q^(13) - 5*q^(21) - 5*q^(29)





		