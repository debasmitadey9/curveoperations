#E1:= y^2 = x^3 + ax + b mod p
#The P-224 weierstrass curve is provided in E1 form, we need to convert it into E2
#E2 = y^2 + a1xy + a3y = x^3

import math
import sys

#~~~ Changing E1 to E2 ~~~
#parameters of E1
#b=18958286285566608000408668544493926415504680968679321075787234672564
#p = 2**224 - 2**96 + 1
#a=-3

#x = 19277929113566293071110308034699488026831934219452440156649784352033
#y = 19926808758034470970197974370888749184205991990603949537637343198772
def modInverse(a, m):
     
    for x in range(1, m):
        if (((a%m) * (x%m)) % m == 1):
            return x
    return -1
    
def cuberoot(a, p):
    if p == 2:
        return a
    if p == 3:
        return a
    if (p%3) == 2:
        return pow(a,(2*p - 1)//3, p)
    if (p%9) == 4:
        root = pow(a,(2*p + 1)//9, p)
        if pow(root,3,p) == a%p:
            return root
        else:
            return None
    if (p%9) == 7:
        root = pow(a,(p + 2)//9, p)
        if pow(root,3,p) == a%p:
            return root
        else:
            return None
    else:
        print("Not implemented yet. See the second paper")
 
 
#legendre's symbol

def isPrime(a):
    return all(a % i for i in range(2, a))
    
def factorize(n):
    factors = []

    p = 2
    while True:
        while n % p == 0 and n > 0:  # while we can divide by smaller number, do so
            factors.append(p)
            n = n / p
        p += 1  # p is not necessary prime, but n%p == 0 only for prime numbers
        if p > n / p:
            break
    if n > 1:
        factors.append(n)
    return factors

def calculateLegendre(a, p):
    if a >= p or a < 0:
        return calculateLegendre(a % p, p)
    elif a == 0 or a == 1:
        return a
    elif a == 2:
        if p % 8 == 1 or p % 8 == 7:
            return 1
        else:
            return -1
    elif a == p - 1:
        if p % 4 == 1:
            return 1
        else:
            return -1
    elif not isPrime(a):
        factors = factorize(a)
        product = 1
        for pi in factors:
            product *= calculateLegendre(pi, p)
            return product
    else:
        if ((p - 1) / 2) % 2 == 0 or ((a - 1) / 2) % 2 == 0:
            return calculateLegendre(p, a)
        else:
            return (-1) * calculateLegendre(p, a)



def modular_sqrt(a, p):

    def legendre_symbol(a, p):
        """ Compute the Legendre symbol a|p using
            Euler's criterion. p is a prime, a is
            relatively prime to p (if p divides
            a, then a|p = 0)
            Returns 1 if a has a square root modulo
            p, -1 otherwise.
        """
        ls = pow(a, (p - 1) // 2, p)
        return -1 if ls == p - 1 else ls

    """ Find a quadratic residue (mod p) of 'a'. p
        must be an odd prime.
        Solve the congruence of the form:
            x^2 = a (mod p)
        And returns x. Note that p - x is also a root.
        0 is returned is no square root exists for
        these a and p.
        The Tonelli-Shanks algorithm is used (except
        for some simple cases in which the solution
        is known from an identity). This algorithm
        runs in polynomial time (unless the
        generalized Riemann hypothesis is false).
    """
    # Simple cases
    #
    if legendre_symbol(a, p) != 1:
        return 0
    elif a == 0:
        return 0
    elif p == 2:
        return p
    elif p % 4 == 3:
        return pow(a, (p + 1) // 4, p)

    # Partition p-1 to s * 2^e for an odd s (i.e.
    # reduce all the powers of 2 from p-1)
    #
    s = p - 1
    e = 0
    while s % 2 == 0:
        s //= 2
        e += 1

    # Find some 'n' with a legendre symbol n|p = -1.
    # Shouldn't take long.
    #
    n = 2
    while legendre_symbol(n, p) != -1:
        n += 1

    # Here be dragons!
    # Read the paper "Square roots from 1; 24, 51,
    # 10 to Dan Shanks" by Ezra Brown for more
    # information
    #

    # x is a guess of the square root that gets better
    # with each iteration.
    # b is the "fudge factor" - by how much we're off
    # with the guess. The invariant x^2 = ab (mod p)
    # is maintained throughout the loop.
    # g is used for successive powers of n to update
    # both a and b
    # r is the exponent - decreases with each update
    #
    x = pow(a, (s + 1) // 2, p)
    b = pow(a, s, p)
    g = pow(n, s, p)
    r = e

    while True:
        t = b
        m = 0
        for m in range(r):
            if t == 1:
                break
            t = pow(t, 2, p)

        if m == 0:
            return x

        gs = pow(g, 2 ** (r - m - 1), p)
        g = (gs * gs) % p
        x = (x * gs) % p
        b = (b * g) % p
        r = m
	   
p = 2003
x3 = 522
y3 = 1914
a =1132
b = 278

#Going over the points starting from x=1
#calculating x^3 + ax + b

for i in range(2,2003):
        send = (i**3)%p
        send1 = (a*i)%p
        send2 = (send + send1 + b)%p
        sq1 = modular_sqrt(send2,p)
        sq2 = (-sq1)%p
        #check if point is order 3
        if i == 598:
            print("Value of y for x=598 is=",sq1,", or ",sq2)
        temp1 = (3*(i**4))%p
        temp2 = (6*a*(i**2))%p
        temp3 = (12*b*i)%p
        temp4 = (a**2)%p

        temp1 = (temp1 + temp2 + temp3)%p
        temp2 = temp1 - temp4
        if temp2 == 0:
            print("result of order:",temp2)
            print("Order 3 point is:(",i,sq1,")")
            print("Order 3 point is:(",i,sq2,")")


temp1 = ((3*((x3)**2)) + a)%p
temp2 = (2*1914)%p
temp2 = modInverse(temp2,p)
slope = (temp1*temp2)%p
#choosing a point of order 3 = (522,1914)
#moving the point to the origin to obtain the equation
# of form y^2 + a1xy + a3y = x^3


#a1xy

a1 = (slope*2)%p
print("a1:",a1)
#a3y

a3 = (2*y3)%p
print("a3:",a3)

 # derived from equation 2
delta = ((a1**3) - (27*a3))%p
print("delta",delta)
temp1 = ((-1)*(27*a3*(delta**2)))%p
temp2 = ((delta**3))%p
temp1 = ((temp1 - temp2))%p
temp1 = cuberoot(temp1,p) #cube root
temp1 = (temp1 + delta)%p
miu = (temp1 * modInverse(3,p))%p
print("miu",miu)

#Equation of Hessian curve: X^3 + Y^3 + Z^3 = cXYZ
#calculating points on the curve

#(x1,y1) is a point on E2

x1 = 598
y1 = 85

#calculating other points of E2

for i in range(2,p):
    for j in range(2,p):
        temp1 = (i**3)%p
        temp2 = (j**2)%p
        temp2 = (temp2 + (a1 * i * j)%p)%p
        temp2 = (temp2 + (a3 * j)%p)%p
        if temp1 == temp2:
            x1 = i
            y1 = j

#X1
            temp1 = (a1*((2*miu)-delta)*x1)%p #(x1,y1) is a point on E2
            temp2 = ((3*miu)-delta)%p
            temp2 = modInverse(temp2,p)
            temp1 = (temp1*temp2)%p
            X1 = (temp1 + y1 + a3)%p

#Y1
            temp1 = ((-1)*a1*miu*x1)%p
            temp2 = ((3*miu)-delta)%p
            temp2 = modInverse(temp2,p)
            temp1 = (temp1*temp2)%p
            Y1 = (temp1 - y1)%p

#Z1
            temp1= ((-1)*a1*miu*x1)%p
            temp2 = ((3*miu)-delta)%p
            temp2 = modInverse(temp2,p)
            temp1 = (temp1*temp2)%p
            Z1 = (temp1 - a3)%p

#c
            miuinv = modInverse(miu,p)
            c = ((3*(miu-delta))*miuinv)%p

            print("X1:",X1,"Y1:",Y1,"Z1:",Z1)
