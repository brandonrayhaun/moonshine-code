import ast
## v1

def parse_number_list(string):
	return map(int,string.split(','))
	
def parse_intofint_list(string):
	string = ("".join(string.split()))[1:-1] + ','
	string = string.split('],')[:-1]
	return map(lambda x: parse_number_list(x[1:]), string)

def parse_decompositions():
	fid = open('potential_decomps.txt')
	l = fid.read()
	return parse_intofint_list(l)
	
def parse_characters():
	fid = open('thompson_characters.txt','r')
	rep1 = []
	rep2 = []
	rep3 = []
	rep5 = []
	rep7 = []
	rep13 = []
	line = fid.readline()
	line_num = 1
	while(line != ''):
		xy = filter(lambda s: s.strip() != '',line.split(' '))
		if not len(xy) == 6:
			print line_num
			print len(xy)
			print xy
			raise NameError('character table incorrect!')
		rep1.append(int(xy[0]))
		rep3.append(int(xy[1]))
		rep2.append(int(xy[2]))
		rep5.append(int(xy[3]))
		rep7.append(int(xy[4]))
		rep13.append(int(xy[5]))
		line = fid.readline()
		line_num = line_num+1
	return (rep1,rep2,rep3,rep5,rep7,rep13)
	
def sum_decomp(decomp, n):
	result = 0
	if n == 2:
		for de in decomp:
			result = result + repDict2[de]
	if n == 3:
		for de in decomp:
			result = result + repDict3[de]
	if n == 5:
		for de in decomp:
			result = result + repDict5[de]
	if n == 7:
		for de in decomp:
			result = result + repDict7[de]
	if n == 13:
		for de in decomp:
			result = result + repDict13[de]
	return result
# 	
# cl1 = parse_characters()[0]
# cl2 = parse_characters()[1]
# cl3 = parse_characters()[2]
# cl5 = parse_characters()[3]
# cl7 = parse_characters()[4]
# cl13 = parse_characters()[5]
# cln1 = [cl1[0]]
# cln2 = [cl2[0]]
# cln3 = [cl3[0]]
# cln5 = [cl5[0]]
# cln7 = [cl7[0]]
# cln13 = [cl13[0]]
# 
# for i in range(1,len(cl1)):
# 	if (cl1[i] != cl1[i-1]):
# 		cln1.append(cl1[i])
# 		cln2.append(cl2[i])
# 		cln3.append(cl3[i])
# 		cln5.append(cl5[i])
# 		cln7.append(cl7[i])
# 		cln13.append(cl13[i])
# 		
# repDict2 = {}
# repDict3 = {}
# repDict5 = {}
# repDict7 = {}
# repDict13 = {}
# for i in range(len(cln1)):
# 	repDict2[cln1[i]] = cln2[i]
# 	repDict2[-cln1[i]] = -cln2[i]
# 	repDict3[cln1[i]] = cln3[i]
# 	repDict3[-cln1[i]] = -cln3[i]
# 	repDict5[cln1[i]] = cln5[i]
# 	repDict5[-cln1[i]] = -cln5[i]
# 	repDict7[cln1[i]] = cln7[i]
# 	repDict7[-cln1[i]] = -cln7[i]
# 	repDict13[cln1[i]] = cln13[i]
# 	repDict13[-cln1[i]] = -cln13[i]
# 	
# test_coefficient3 = 378
# test_coefficient2 = 42
# test_coefficient5 = 21
# test_coefficient7 = 0
# test_coefficient13 = 1
# 
# decomps = parse_decompositions()
# 
# t_decomps = []
# for d in decomps:
# 	cond2 = sum_decomp(d,2) == test_coefficient2
# 	cond3 = sum_decomp(d,3) == test_coefficient3 
# 	cond5 = sum_decomp(d,5) == test_coefficient5
# 	cond7 = sum_decomp(d,7) == test_coefficient7
# 	cond13 =sum_decomp(d,13)== test_coefficient13
#  	if cond2 and cond3 and cond5 and cond7 and cond13:
#  		t_decomps.append(d)
#  		
# t_decomps = filter(lambda x: (779247 in x) or (-779247 in x), t_decomps)
# t_decomps = filter(lambda x: not x.count(779247) == x.count(-779247),t_decomps)

##v2

test_coefficients = {}
test_coefficients['3b'] = -756
# test_coefficients['2a'] = 43008
# test_coefficients['5a'] = -30
# test_coefficients['7a'] = 0
# test_coefficients['10a'] = -2
test_coefficients['13a'] = 0

def E(N):
	return exp(2*pi*I/N)
	
def parse_characters2():
	a = 3*E(3)-E(3)^2
	b = 6*E(3) - 2*E(3)^2
	c = -E(15)^7-E(15)^11 - E(15)^13 - E(15)^14
	d = -E(3)+E(3)^2
	e = E(24)+E(24)^11 - E(24)^17 - E(24)^19
	f = 2*E(3) - E(3)^2
	g = E(31)+E(31)^2+E(31)^4+E(31)^5+E(31)^7+E(31)^8+E(31)^9+E(31)^10+E(31)^14+E(31)^16+E(31)^18+E(31)^19+E(31)^20+E(31)^25+E(31)^28
	h = 2*E(3)
	i = -E(39)^7-E(39)^14-E(39)^17-E(39)^19-E(39)^23-E(39)^28-E(39)^29-E(39)^31-E(39)^34-E(39)^35-E(39)^37-E(39)^38
	notation = {'a':a, 'b':b, 'c':c, 'd':d, 'e':e, 'f':f, 'g':g, 'h':h, 'i':i, '/a':conjugate(a), '/b':conjugate(b), '/c':conjugate(c), '/d':conjugate(d), '/e':conjugate(e), '/f':conjugate(f), '/g':conjugate(g), '/h':conjugate(h), '/i':conjugate(i) }
	fid = open('thompson_characters2.txt','r')
	rep_table =[]
	for i in range(48):
		rep_table.append([])
	classes = filter(lambda s: s.strip() != '',fid.readline().split(' '))
	line = fid.readline()
	while (line != ''):
		xy = filter(lambda s: s.strip() != '',map(lambda s: s.strip(),line.split(' ')))
		row = map(lambda x: parse_entries(x,notation),xy[1:])
		for i in range(len(row)):
			rep_table[i].append(row[i])
		line = fid.readline()
	classes[-1] = '39a'
	return (classes[1:],rep_table[:])

def parse_entries(s,notation):
	if '-' in s:
		s_ = s[1:]
		if s_.isdigit():
			return -int(s_)
		else:
			return (-notation[s_.lower()]).n()
	if s.isdigit():
		return int(s)
	else:
		return (notation[s.lower()]).n()

	

(classes,table)  = parse_characters2()

def conj_dict(classi):
	return classes.index(classi)


k = 1
while k < len(table[0]):
	if table[0][k] == table[0][k-1]:
		for j in range(len(table)):
			del table[j][k]
	else:
		k=k+1
		
repDict = []
for j in range(len(table)):
	repDict.append({})

# def negate_entry(e):
# 	if isinstance(e, (int, long, float, complex)):
# 		return -e
# 	if type(e) is str:
# 		return '-' + e
			
for i in range(len(table)):
	for j in range(len(table[i])):
		repDict[i][table[0][j]] = table[i][j]
		repDict[i][-table[0][j]] = -(table[i][j])
	
#decomps = parse_decompositions()

def sum_decomp2(decomp, classi):
	result = 0
	for de in decomp:
		result = result + real_part(repDict[conj_dict(classi)][de])
	return result
	
#t_decomps2 = []
#cond2 = []
#for key in test_coefficients:
#	cond2.append(0)
#for d in decomps:
#	k = 0
#	for key in test_coefficients:
#		cond2[k] = (sum_decomp2(d,key) == test_coefficients[key])
#		k = k + 1
# 	if all(cond2):
# 		t_decomps2.append(d)
 		
# t_decomps2 = filter(lambda x: (779247 in x) or (-779247 in x), t_decomps)
# t_decomps2 = filter(lambda x: not x.count(779247) == x.count(-779247),t_decomps)

# def compute_twist(classi):
# 	curr_decomps = ([1,4,5,8,9,12],[[-248],[27000,-248],[-85995],[1707264],[-4096000,-248],[44330496]])
# 	decomps = curr_decomps[1]
# 	li = map(lambda d: sum_decomp2(d,classi), decomps)
# 	for i in range(len(curr_decomps[1])):
# 		print(curr_decomps[0][i],li[i])

o16 = [190373976,190373976,91171899,81153009,76271625,72925515,6669000]
o17 = [-190373976,-190373976,-190373976,-111321000,-111321000,-91171899,-81153009,-77376000,-76271625,-72925515,-51684750,-40199250,-30507008,-28861000]
o16.extend(o16)
o17.extend(o17)
f_3orders = [-3,0,4,5,8,9,12,13,16,17]
f_3decomps = [[1,1],[248],[27000,27000],[-85995,-85995],[1707264,1707264],[-4096000,-4096000],[44330496,44330496],[-91171899,-91171899,-779247,-779247],o16,o17]
theta = 1 + 2*q + 2*q^4 + 2*q^9 + 2*q^16  
## Starts Laurent series in q
R.<q> = LaurentSeriesRing(QQ)
def get_twist(f_decomps,orders,classi):
	result = 0
	for i in range(len(orders)):
		co = sum_decomp2(f_decomps[i],classi)
		result = result + co*q^(orders[i])
	return result
	
  
def print_to_latex_twist(f_decomps,orders,classi):
	string = "$\mathcal{F}_{" + classi.upper() + "}(\\tau)=$ &$2q^{-3}"
	for i in range(len(orders)):
		co = sum_decomp2(f_decomps[i],classi)
		if co != 0:
			if orders[i] != 0:
				if co >=0:
					string = string + " + " + str(co) + "q^{" + str(orders[i]) + "} " 
				else:
					string = string + " " + str(co) + "q^{" + str(orders[i]) + "} "
			else:
				if co >= 0:
					string = string + " + " + str(co) 
				else:
					string = string + " " + str(co)
	print(string + "+ \mathcal{O}(q^{20})$ & \\\\" )		