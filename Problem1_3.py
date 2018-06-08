#solution to problem 3, which as I understand also works for problem 1 for the deterministic case
#assumes the automaton is deterministic and complete (transitions specified for all letters)

from sympy import *
import collections

#reads the description of an automaton / markov chain into a 2D dictionary - trans[start_vertex][letter] = (end_vertex, weight / probability)
def read_aut():
    n, m = map(int, raw_input().split(" "))
    trans = collections.defaultdict(dict)
    for i in range(m):
        ln = raw_input().split(" ")
        trans[int(ln[0])][ln[1]] = int(ln[3]), Rational(ln[2])
        #print(int(ln[0]), ln[1], trans[int(ln[0])][ln[1]])
    return n, m, trans

aut_n, aut_m, aut_trans = read_aut()
mar_n, mar_m, mar_trans = read_aut()


#finds a solution by considering a product of the automaton and the markov chain where the edges are labeled by both weights and probabilities (but no letters)
#the expected value of the original automaton under the distribution is the expected value obtained in this product from the starting state (corresponding to the pair (starting state of automaton, starting state of markov chain)
#the expected value is found using a system of linear equations of the form: expected_value_x = sum_over_y (probability_of_transition_to_y * (expected_value_y + weight_of_transition_to_y))
nn = aut_n * mar_n
eq_mat = eye(nn)
rhs = zeros(nn, 1)

for frm_aut in aut_trans:
    for l in aut_trans[frm_aut]:
        for frm_mar in mar_trans:
            frm = frm_aut * mar_n + frm_mar
            to = aut_trans[frm_aut][l][0] * mar_n + mar_trans[frm_mar][l][0]
            w = aut_trans[frm_aut][l][1]
            p = mar_trans[frm_mar][l][1]
            #print(frm, l, to, w, p)
            eq_mat[frm,to] -= p
            rhs[frm] += p * w

#print(eq_mat)
#print(rhs)
x = symbols('x')
print(next(iter(linsolve((eq_mat, rhs), x)))[0])
