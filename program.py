
lam=3
mu=-4
u0=3
u1=1

import math as maths

def factorial(n):
    
    if n < 0:
        return maths.inf
    
    if n==0:
        return 1
    
    total = 1
    
    for i in range(1, n+1):
        total *= i
        
    return total

def F(n,k):
    return  ( factorial(k) * ( lam ** (2*k-n) ) * ( mu ** (n-k) ) ) / ( factorial(2*k-n) * factorial(n-k) ) 

def recurrence(n):
    total =0
    
    for k in range(1,n-1 +1):
        total += ( F(n-2,k-1) * mu * u0 + F(n-1,k) * u1)
        
    return total


for i in range(0,25):
    print(i , recurrence(i))