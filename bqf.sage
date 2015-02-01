##
I = CC.0 # imaginary unit
R.<q> = LaurentSeriesRing(QQ)
precision = 200

def q_inf(c,prec):
	q_ = 1
	for i in range(1,prec):
		q_ = (1-q^(c*i))*q_
	return q_
	
	
num = q^(-1)*(q_inf(1,precision)^12)
denom = q_inf(3,precision)^12

# q_inf_ = 1
# for i in range(1,precision):
# 	q_inf_ = q_inf_ + ((-1)^i)*q^((3*i^2 - i)/2)
# 	i = -i
# 	q_inf_ = q_inf_ + ((-1)^i)*q^((3*i^2 - i)/2)

#t_3 = q^(-1)*((q_inf(1,precision)^12)/(q_inf(3,precision)^12)) + 8
#t_6 = q^(-1)*( (q_inf(1,precision)^5)*(q_inf(3,precision)) )/( (q_inf(2,precision))*(q_inf(6,precision)^5) ) + 5

#f3b = q^(-3) -5*q + 22*q^4 + 27*q^5 - 54*q^8 + 3*q^9 + 6*q^12 - 189*q^(13)+346*q^(16)
#f3  = q^(-3) - 248*q +  26752*q^4 - 85995*q^5 + 1707264*q^8 - 4096248*q^9 + 44330496*q^(12) - 91951146*q^(13) + 708938752*q^(16)
#f33 = q^(-3) - 14*q + 40*q^4 - 78*q^9 + 168*q^(12) - 378*q^(13) + 688*q^(16)
#f6c = q^(-3) -q + 2*q^4 + 3*q^5 + 6*q^8-q^9+ 2*q^(12) + 3*q^(13)
#a,b,c,d = var('a b c d')

#eqn = []
#eqn.append(a+b+c == 1)
#eqn.append(27*a-85995*b == -27)
#eqn.append(-54*a+1707264*b == 54 )