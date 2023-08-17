
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

def binomial(n,k):
    return factorial(n) / ( factorial(k) * factorial(n-k) ) 


def gregory_newton(n, numbers, kthdifference):
    
    if kthdifference > len(numbers):
        return None
    
    difference_lists = [numbers[0:kthdifference+1] ]
    
    
    for ithdifference in range(kthdifference):
        
        next_list = []
        last_list = difference_lists[-1]
        
        
        for i in range(0, len(last_list) - 2 + 1):
            
            this_diff = last_list[i+1] - last_list[i]
            
            next_list.append(this_diff)
            
        difference_lists.append(next_list)    
            
    total = 0            

    for i in range(0,len(difference_lists)):
        
        total += difference_lists[i][0] * binomial(n, i)
        
    return total
            
        
thenumbers = [x**2 - 3*x + 1 for x in range(0,5)]

print(gregory_newton(35, thenumbers, 3))
            
            
            
    