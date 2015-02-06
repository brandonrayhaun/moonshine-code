## starting fresh with the trace formulation

S = matrix(2,[0,-1,1,0])
T = matrix(2,[1, 1,0,1])

def action_matrix(M):
	alpha = M[0,0]   
	beta  = M[0,1]
	gamma = M[1,0]
	delta = M[1,1]
	return (matrix(3,[alpha^2, alpha*gamma, gamma^2, 2*alpha*beta,  
		alpha*delta + beta*gamma, 2*gamma*delta, beta^2, beta*delta, delta^2]))

def b_action(gamma,Q):
	M = action_matrix(gamma)
	Q_ = (M^(-1))*matrix(3,1,[Q[0],Q[1],Q[2]])
	return BinaryQF([Q_[0,0],Q_[1,0],Q_[2,0]])

def trans(Q,t):
	return BinaryQF([Q[0],-2*t*Q[0]+Q[1], t^2*Q[0]-t*Q[1] + Q[2]])

def normal_form(Q):
	t = -floor((Q[0]-Q[1])/(2*Q[0]))
	M = T^t
	return b_action(M,Q), M

def is_reduced(Q):
	return Q[0]<Q[2] or (Q[0]==Q[2] and -Q[1] <= 0)

def reduce_form(Q,N):
	(Q, M) = normal_form(Q)
	while not is_reduced(Q):
		(Q,A) = normal_form(b_action(S,Q))
		M = A*S*M
	if N == 1:
		return Q,M
	else:
		RN = range(1, (N-1)/2 + 1)
		RN = RN + map(lambda g:-g, RN)
		reps = [S] + map(lambda g: T^g, RN)
		for r in reps:
			M_ = r*M
			if mod(M_[0,1] ,N) == 0:
				return b_action(M,Q), M
		error("doesn't give congruence subgroup!")



