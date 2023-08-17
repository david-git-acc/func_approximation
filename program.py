
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
            
        


def first_n_numbers_of_recurrence(n):
    
    if n == 0:
        return [u0]
    
    if n== 1:
        return [u0,u1]
    
        
    nums = [u0,u1]
    
    for _ in range(n-2):
        
        nextnum = lam * nums[-1] + mu * nums[-2]
        
        nums.append(nextnum)
        
    return nums


from matplotlib import pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

extent = 25
sign = "-" if mu < 0 else "+"

fig, (ax1,ax2,ax3) = plt.subplots(ncols=3,nrows=1)

recurrence_from_0_to_extent = first_n_numbers_of_recurrence(extent) 
my_recurrence_version = [recurrence(i) for i in range(0, extent)]
greg_numbers = [gregory_newton(n, recurrence_from_0_to_extent, extent-1) for n in range(0,24+1)]


ax1.grid(True)
ax2.grid(True)
ax3.grid(True)




x = np.arange(0,extent)

growth_rate_limiter = lambda x : x ** (1/3)

ax1plot_original = ax1.plot(x, [ growth_rate_limiter(x) for x in recurrence_from_0_to_extent] , color="green")
my_rec_values = ax2.plot(x, [growth_rate_limiter(x) for x in my_recurrence_version], color="red")
greg_formula = ax3.plot(x, [ growth_rate_limiter(x) for x in greg_numbers], color = "blue")

ax1.set_ylabel("u(n)", fontsize=15)
ax1.set_xlabel("n (original recurrence)" , fontsize=15)

ax2.set_xlabel("Finite series approximation")
ax3.set_xlabel("Gregory-Newton approximation")

ax3.yaxis.set_label_position("right")
ax3.set_ylabel(f"(data coded: f(x) = cbrt(x))")

plt.suptitle(f"Approximating the recurrence: u(n) = {lam}u(n-1) {sign} {abs(mu)}u(n-2)")

plt.show()
            
            
            
    