
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
    
    bottomleft = 2*k-n
    bottomright = n-k
    
    if bottomleft < 0 or bottomright < 0:
        return 0
    
    minimum = min(bottomleft,bottomright)
    maximum = max(bottomleft,bottomright)
    
    total = 1
    
    divide = []
    mult = []
    
    for i in range(1, minimum + 1):
        
        divide.append( i )
        

          
    for i in range(maximum+1, k+1):
        mult.append(i)
        
      
    if len(divide) > len(mult):
        
        for i in range(len(mult)):
            
            total *= (mult[i] / divide[i])
            
        for j in range(len(mult), len(divide) ):
            total /= divide[j]
            
    else:
        
        for i in range(len(divide)):
            
            total *= (mult[i] / divide[i])
            
        for j in range(len(divide), len(mult) ):
            total *= mult[j]
        
    
    return total * (lam ** bottomleft) * (mu ** bottomright)
    

    
    return  ( factorial(k) * ( lam ** (2*k-n) ) * ( mu ** (n-k) ) ) / ( factorial(2*k-n) * factorial(n-k) ) 

def recurrence(n,ktimes):
    total =0
    
    for k in range(1,ktimes +1):
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

px =1/96
extent = 30
sign = "-" if mu < 0 else "+"

fig, (ax1,ax2,ax3) = plt.subplots(ncols=3,nrows=1, figsize=(1920*px, 1080*px))

ax1.grid(True)

x = np.arange(0,extent)

growth_rate_limiter = lambda x : x ** (1/27)



title = plt.suptitle(f"Approximating the recurrence: u(n) = {lam}u(n-1) {sign} {abs(mu)}u(n-2)")
    
recurrence_from_0_to_extent = first_n_numbers_of_recurrence(extent) 
ax1plot_original = ax1.plot(x, [ growth_rate_limiter(x) for x in recurrence_from_0_to_extent] , color="green")    

original_recurrence_frames = []
my_recurrence_frames = []
gregs_frames = []

for t in range(extent):
    
    original_rec_frames = first_n_numbers_of_recurrence(t+4)
    my_recurrence_version = [recurrence(n, t) for n in range(0, extent)]
    greg_numbers = [gregory_newton(n, recurrence_from_0_to_extent, t) for n in range(0,extent)]
    
    original_recurrence_frames.append(original_rec_frames)
    my_recurrence_frames.append(my_recurrence_version)
    gregs_frames.append(greg_numbers)
    
    
interpolation_extent = 100

def interpolate(my_recurrence_frames):

    my_interpolated_recurrence_frames = []

    for i in range(len(my_recurrence_frames) - 1):
        
        prev_framelist = my_recurrence_frames[i]
        next_framelist = my_recurrence_frames[i+1]
        
        interpolation_lists = [prev_framelist]
        
        for j in range(1, interpolation_extent-1 + 1):
            
            the_next_new_interpolation_list = []
            
            for k, value in enumerate(prev_framelist):
                
                current_frame_from_prevframelist = prev_framelist[k]
                current_frame_from_nextframelist = next_framelist[k]
                diff = current_frame_from_nextframelist - current_frame_from_prevframelist
                
                the_next_new_frame = current_frame_from_prevframelist + (j/interpolation_extent) * diff
                
                the_next_new_interpolation_list.append(the_next_new_frame)
                
            interpolation_lists.append(the_next_new_interpolation_list)
            
        my_interpolated_recurrence_frames += interpolation_lists
        
    my_interpolated_recurrence_frames.append( my_recurrence_frames[-1] )
    
    return my_interpolated_recurrence_frames

original_recurrence_frames = interpolate(original_recurrence_frames)
my_recurrence_frames = interpolate(my_recurrence_frames)
gregs_frames = interpolate(gregs_frames)
        
top = max( [growth_rate_limiter(x).real for x in original_recurrence_frames[-1]]) 
bottom = min( [growth_rate_limiter(x).real for x in original_recurrence_frames[-1]])

text = plt.text(0.08, -0.2, f"u(0) = {u0}, u(1) = {u1}",fontsize=12    , ha="left", va = "bottom", transform=ax1.transAxes)


def animate(t):

    global my_reccurrence_frames
    global gregs_frames

    ax1.cla()
    ax2.cla()
    ax3.cla()
    
    ax1.grid(True)
    ax2.grid(True)
    ax3.grid(True)
    
    text = plt.text(0.1, -0.1, f"u(0) = {u0}, u(1) = {u1}",fontsize=12    , ha="left", va = "bottom", transform=ax1.transAxes)
    

    
    the_true_t = t // interpolation_extent
    title = plt.suptitle(f"Approximating the recurrence: u(n) = {lam}u(n-1) {sign} {abs(mu)}u(n-2)", fontsize=20)
    
    ax1.set_ylabel("u(n)", fontsize=15)
    ax1.set_xlabel("n (original recurrence)" , fontsize=15)

    ax2.set_xlabel("Finite series approximation", fontsize= 15)
    ax3.set_xlabel("Gregory-Newton approximation",fontsize = 15)

    ax3.yaxis.set_label_position("right")
    ax3.set_ylabel(f"(data coded: f(x) = x^(1/27)")


    
    next_orig = original_recurrence_frames[t]


    original_rec_values = ax1.plot(np.arange(0,len(next_orig)) , [growth_rate_limiter(x) for x in next_orig], color="green")
    my_rec_values = ax2.plot(x, [growth_rate_limiter(x) for x in my_recurrence_frames[t] ], color="red")
    greg_formula = ax3.plot(x, [ growth_rate_limiter(x) for x in gregs_frames[t]], color = "blue")
    
    
    ax1.set_ylim(bottom,top)
    ax2.set_ylim(bottom, top )
    ax3.set_ylim(bottom, top)
    ax1.set_xlim(0,30)
    ax2.set_xlim(0,30)
    ax3.set_xlim(0,30)



anim = FuncAnimation(fig, animate, frames = extent*(interpolation_extent ) - interpolation_extent, interval=1, repeat_delay = 2500)


anim.save("first.gif", fps=1000 )



            
            
            
    