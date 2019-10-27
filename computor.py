import argparse
import re
import collections
import cmath
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt

def parsing():
	parser = argparse.ArgumentParser()
	parser.add_argument('equation', help='equation to be solved')
	parser.add_argument('-p', '--plot', help='plot for polynmial with degree 1 or 2', default=False, action='store_true')
	parser.add_argument('-v', '--verbose', help='show intermediate steps', default=False, action='store_true')
	args = parser.parse_args()

	eq = args.equation
	plot_arg = args.plot
	verbose = args.verbose

	if (len(eq) == 0):
		print ('No equation.')
		exit(0)

	pattern = re.compile("(^\\s*|\\+\\s*|-\\s*|=\\s*)((-?\\d+)(\\.\\d+)?)?\\s*(\\*?\\s*X\\s*(\\^\\s*(-?\\d+)(\\.\\d+)?)?)?(\\s*|$)")
	result = pattern.findall(eq)

	check_eq = re.sub(pattern, '', eq)
	check_eq = re.sub('\\s+', '', check_eq)
	computor_symbols = '+-*=.^X0123456789'
	if (check_eq):
		for ch in check_eq:
			if ch not in computor_symbols:
				print ('Lexical error: unrecognizable symbol ' + ch + '. Available symbols: [' + computor_symbols + ']')
				exit(0)
		print ('Lexical error: incorrect input.')
		exit(0)

	return result, plot_arg, verbose

'''
Terms:
0 - sign of number (+, -, = or space before)
1 - number, if empty, number equals to 1
2 - whole part of number
3 - fractional part of number
4 - if empty, that means that pow of term is equal 0
5 - if empty, that means that pow of term is equal 1
6 - whole part of degree
7 - fractional part of degree
8 - space or end of line, check ending
'''

def get_degrees(result):
	degrees = {0:0, 1:0, 2:0}

	counter = 0
	right_side_flag = 1
	error = ''
	for term in result:

		if (len(term[1]) == 0 and len(term[4]) == 0):
			print ('Syntax error: missing number.') #write number
			exit(0)

		#get degree
		if (len(term[4]) == 0):
			degree = 0
		elif (len(term[5]) == 0):
			degree = 1
		else:
			if (len(term[7])):
				degree_float_part = float(term[6])
			else:
				degree_float_part = 0

			degree = int(term[6]) + degree_float_part
			if (degree > 2 or degree < 0 or degree_float_part != 0):
				error = 'Syntax error: wrong degree ' + str(degree) + '. The degree should be integer from 0 to 2, I can\'t solve.'

		#get sign of number
		if (len(term[0]) == 0):
			if (counter != 0):
				print ('Syntax error: missing sign before number.') #write number
				exit(0)
			else:
				sign = 1
		else:
			sign_term = term[0][0]
			if (sign_term == '+'):
				sign = 1
			elif (sign_term == '-'):
				sign = -1
			elif (sign_term == '='):
				if (right_side_flag == -1):
					print ('Syntax error: more than two equal signs')
					exit(0)
				if (counter == 0):
					print ('Syntax error: equal sign in begin of equation')
					exit(0)

				right_side_flag = -1
				sign = 1

		#get number
		if (len(term[1]) == 0):
			number = 1
		else:
			number = int(term[2])
			if (len(term[3])):
				number += float(term[3])

		if degree in degrees:
			degrees[degree] += number * sign * right_side_flag
		else:
			degrees[degree] = number * sign * right_side_flag

		counter += 1

	degrees = collections.OrderedDict(sorted(degrees.items()))
	return degrees, error

def get_max_degree(degrees):
	max_degree = 0
	for degree in reversed(list(degrees.keys())):
		if (degrees[degree] != 0):
			max_degree = degree
			break
	return max_degree

def print_reduced_form(degrees):
	max_degree = get_max_degree(degrees)
	counter = 0
	print ('Reduced form:', end='')
	for p in sorted(degrees):
		if (degrees[p] == 0 and (p != 0 or max_degree != 0)):
			continue
		if (counter == 0 and degrees[p] >= 0):
			sign = ''
		elif (counter == 0 and degrees[p] < 0):
			sign = '-'
		else:
			sign = '+ ' if degrees[p] > 0 else '- '
		print (' %s' % sign, end='')
		if (abs(degrees[p]) != 1 or p == 0):
			print ('%g' %abs(degrees[p]), end='')
		if (abs(degrees[p]) != 1 and p != 0):
			print (' * ', end='')
		if (p == 1):
			print ('X', end='')
		if (p != 1 and p != 0):
			print ('X^%g' % (p), end='')
		counter += 1
	print (' = 0')
	print ('Polynomial degree: %g' % max_degree)

def results_output(degrees, verbose):
	if (degrees[2] != 0):
		D = degrees[1] * degrees[1] - 4 * degrees[0] * degrees[2]
		if (verbose):
			print ('D = ' + str(degrees[1]) + ' * ' + str(degrees[1]) \
				+ ' - 4 * ' + str(degrees[0]) + ' * ' + str(degrees[2]) + ' = ' \
				+ str(degrees[1] * degrees[1]) + ' - 4 * ' + str(degrees[0] * degrees[2]) + ' = ' + str(D))
		if (D > 0):
			x1 = (-degrees[1] - sqrt(D)) / (2 * degrees[2])
			x2 = (-degrees[1] + sqrt(D)) / (2 * degrees[2])
			print ('Discriminant is strictly positive, the two solutions are:')
			if (verbose):
				print ('x1 = ' + '(-1 * ' + str(degrees[1]) + ' - ' + str(sqrt(D)) + ') / (2 * ' + str(degrees[2]) + ') = '\
						+ str(degrees[1] - sqrt(D)) + ' / ' + str(2 * degrees[2]) + ' = ' + str(x1))
				print ('x2 = ' + '(-1 * ' + str(degrees[1]) + ' + ' + str(sqrt(D)) + ') / (2 * ' + str(degrees[2]) + ') = '\
						+ str(degrees[1] + sqrt(D)) + ' / ' + str(2 * degrees[2]) + ' = ' + str(x2))
			else:
				print (x1)
				print (x2)
		if (D == 0):
			x = (-degrees[1]) / (2 * degrees[2])
			print ('Discriminant is equal 0, the solution is:')
			if (verbose):
				print ('x = -1 * ' + str(degrees[1]) + ' / 2 * ' + str(degrees[2]) + ' = ' + str(-degrees[1]) \
				 + ' / ' + str(2 * degrees[2]) + ' = ' + str(x))
			else:
				print (x)
		if (D < 0):
			x1 = (-degrees[1] - cmath.sqrt(D)) / (2 * degrees[2])
			x2 = (-degrees[1] + cmath.sqrt(D)) / (2 * degrees[2])
			print ('Discriminant is strictly negative, the two solutions are:')
			if (verbose):
				print ('x1 = ' + '(-1 * ' + str(degrees[1]) + ' - ' + str(cmath.sqrt(D)) + ') / (2 * ' + str(degrees[2]) + ') = '\
						+ str(-degrees[1] - cmath.sqrt(D)) + ' / ' + str(2 * degrees[2]) + ' = ' + str(x1))
				print ('x2 = ' + '(-1 * ' + str(degrees[1]) + ' + ' + str(cmath.sqrt(D)) + ') / (2 * ' + str(degrees[2]) + ') = '\
						+ str(-degrees[1] + cmath.sqrt(D)) + ' / ' + str(2 * degrees[2]) + ' = ' + str(x2))
			else:
				print (x1)
				print (x2)
	elif (degrees[1] != 0):
		x = -degrees[0] / degrees[1]
		print ('The solution is:')
		if (verbose):
			print ('x = -1 * ' + str(degrees[0]) + ' / ' + str(degrees[1]) + ' = ' + str(x))
		else:
			print (x)
	elif (degrees[0] != 0):
		print ('No solution.')
	elif (degrees[0] == 0):
		print ('The solution are all the real numbers')

def get_plot(degrees):
	if (degrees[2] != 0):
		top = (-degrees[1]) / (2 * degrees[2])
		x = np.linspace(top - 10, top + 10, 200)
		y = degrees[0] + degrees[1] * x + degrees[2] * x * x
	elif (degrees[1] != 0):
		top = degrees[0] / degrees[1]
		x = np.linspace(top - 10, top + 10, 200)
		y = degrees[0] + degrees[1] * x

	plt.plot(x, y)
	plt.title('Polynomial\'s plot')
	plt.xlabel('x')
	plt.ylabel('y')
	plt.grid()
	plt.show()


result, plot_arg, verbose = parsing()
degrees, error = get_degrees(result)
print_reduced_form(degrees)
if (error):
	print (error)
	exit(0)
results_output(degrees, verbose)
if (plot_arg == True and (degrees[2] != 0 or degrees[1] != 0)):
	get_plot(degrees)
