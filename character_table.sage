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

def E(N):
	return exp(2*pi*I/N)

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

def print_table(classes, rep_table, classrange, characterrange):
	strng = """\\begin{sidewaystable}
\\begin{center}
\caption{Character table of the Thompson group, part }\label{thchartab}
\smallskip
\\begin{small}
\\begin{tabular}{c@{ }r@{ }r@{ }r@{ }r@{ }r@{ }r@{ }r@{ }r@{ }r@{ }r@{ }r@{ }r@{ }r@{ }r@{ }r@{ }r@{ }r@{ }r@{ }r@{ }r@{ }r@{ }r@{ }r@{ }r} \\toprule
$[g]$"""
	for i in classrange:
		strng = strng + " & " + classes[i].upper() + " "
	strng = strng + """ \\\\ 
	 \\midrule 
	 """
	for i in characterrange:
		strng = strng + "$\\chi_{" + str(i+1) + "}$ "
		for j in classrange:
			strng = strng + " & $" + str(rep_table[j][i]) + "$ "
		strng = strng + """\\\\
		 """
	strng = strng + """\\bottomrule
	\end{tabular}
	\end{small}
	\end{center}
	\end{sidewaystable}
	"""
	return strng

