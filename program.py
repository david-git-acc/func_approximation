
''' 

The goal of this program is to compare the approximation made by the finite series (see the LinkedIn post or GitHub) of some arbitrary recurrence with 
that of the Gregory Newton interpolation formula. We'll also show the original recurrence to show how "good" each recurrence is doing.

This will create an animation, where each frame of the animation progressively has 1 more "iteration" of the finite series / GN formula applied to it.

The recurrence is of the form u(n) = (lam) * u(n-1) + (mu) * u(n-2), where lam and mu are constants. u0 and u1 represent the base cases of the recurrence.

'''

# Defining our constants - this was one of the more interesting recurrences
lam=3
mu=-4
u0=3
u1=1

import math as maths

# Factorial is used in binomial coefficients, and I was too lazy to google the maths factorial
def factorial(n):
    
    # Inputs will only be integers, and any the factorial of any negative integer = infinity
    if n < 0:
        return maths.inf
    
    # You can't really capture this with the product from 1 to n
    if n==0:
        return 1
    
    # Calculating the factorial
    total = 1
    for i in range(1, n+1):
        total *= i
        
    return total

# This is a helper function used in the finite series calculation.
# Instead of using the standard factorial calculation which was too slow and encountered overflows for n > 84, I used a different way
def F(n,k):
    
    # These are the bottom left and bottom right parts of the formula as shown by the picture
    bottomleft = 2*k-n
    bottomright = n-k
    
    # If either < 0, then as (-x)! = infinity, 1/(x!) = 0
    if bottomleft < 0 or bottomright < 0:
        return 0
    
    # minimum and maximum refer to the bounds. Because this function is the ratio of a factorial to 2 factorials in the denominator, there will be
    # a point where you have (1x2x3x...) / ((1x2x3x...) * (1x2x3x...)) = 1 / (1x2x3...) which is equivalent to dividing by these numbers up to the minimum
    minimum = min(bottomleft,bottomright)
    
    # Between the minimum and maximum points, you will have the ratio of 2 factorials, which cancel out and equal 1. The growth of the function only begins from 
    # (maximum) to k
    maximum = max(bottomleft,bottomright)
    
    # Our product accumulator - always start with the identity element
    total = 1
    
    # List of numbers to divide and multiply by
    divide = []
    mult = []
    
    # As stated above, from 1 to minimum we need to divide by these numbers
    for i in range(1, minimum + 1):
        
        divide.append( i )
        

    # Again as stated above, during this range we multiply  
    for i in range(maximum+1, k+1):
        mult.append(i)
        
    # This IF statement is just to avoid "index out of range" errors, the idea is that we multiply the total by the RATIO between the next number to multiply by
    # and the next number to divide by, that way we won't get overflows and less likely to encounter floating point errors on small values
    if len(divide) > len(mult):
        
        # See above
        for i in range(len(mult)):
            
            total *= (mult[i] / divide[i])
            
        # As there are more numbers to divide by than to multiply by, for the difference in the number of multipliers/divisors we only multiply
        for j in range(len(mult), len(divide) ):
            total /= divide[j]
            
    # Symmetrically inverse case to above - if we have more multipliers than divisors, do the exact opposite
    else:
        
        # See above
        for i in range(len(divide)):
            
            # Ditto
            total *= (mult[i] / divide[i])
            
        # As more multipliers than divisors, after we run out of divisors we only multiply
        for j in range(len(divide), len(mult) ):
            total *= mult[j]
        
    
    # Don't forget to multiply by the exponents as stated in the formula for F(n,k)
    return total * (lam ** bottomleft) * (mu ** bottomright)
    
# This function performs the finite series recurrence
# First argument is the number, second is how many terms we sum up (minimum:0, maximum: n-1)
def recurrence(n,ktimes):
    total =0
    
    # Performing the formula shown in the recurrence
    for k in range(1,ktimes +1):
        total += ( F(n-2,k-1) * mu * u0 + F(n-1,k) * u1)
        
    return total

# Binomial formula, using in Gregory Newton formula
def binomial(n,k):
    return factorial(n) / ( factorial(k) * factorial(n-k) ) 

# Unfortunately because I didn't have the time to learn how to do partial function application in Python, I have to do this every single time that
# I need to find n, if I had more time then I would've realised I simply could've made it have an array of numbers instead of a single number to avoid complexity
# O(N^2), but since the sample size (30) is low then it's worth the time saved
# kthdifference = how many numbers we want to fit, so if we wanted to interpolate between k consecutive numbers, we'd need to know from 0 up to the (k-1)th difference
def gregory_newton(n, numbers, kthdifference):
    
    # You can't use this formula without at least as many numbers as the kth difference
    if kthdifference > len(numbers):
        return None
    
    # Start with the numbers themselves, which are the 0th difference
    # Limit the number array to only as many numbers as we need to find the kth difference (hence, k numbers)
    difference_lists = [numbers[0:kthdifference+1] ]
    
    
    for ithdifference in range(kthdifference):
        
        # We need to find the 1st, 2nd, ... , kth difference arrays
        # Hence we need to get the difference list we just computed to help us compute the next difference list
        next_list = []
        last_list = difference_lists[-1]
        
        # Calculate each individual difference
        for i in range(0, len(last_list) - 2 + 1):
            
            # The difference itself
            this_diff = last_list[i+1] - last_list[i]
            next_list.append(this_diff)
            
        # Add the now completed ith difference to the list of differences
        difference_lists.append(next_list)    
    
    # accumulator, now we sum up the differences multiplied by the binomial coefficients of n
    total = 0            

    # See above
    for i in range(0,len(difference_lists)):
        
        total += difference_lists[i][0] * binomial(n, i)
        
    return total
            
        

# Get the first n numbers of the recurrence (starting from u0)
# O(n)
def first_n_numbers_of_recurrence(n):
    
    if n == 0:
        return [u0]
    
    if n== 1:
        return [u0,u1]
    
        
    nums = [u0,u1]
    
    # Performing the recurrence O(n) times to get the n numbers
    for _ in range(n-2):
        
        nextnum = lam * nums[-1] + mu * nums[-2]
        
        nums.append(nextnum)
        
    return nums

# It's plotting time

from matplotlib import pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

# Standard data
# px conversion ratio because I don't know what inches are
# extent gives us the x-axis range, from 0 to extent, in this case from 0 to 30
px =1/96
extent = 30

# Seeing +-(number) instead of just -(number) in the title triggers me so I made this to tell the difference
sign = "-" if mu < 0 else "+"

fig, (ax1,ax2,ax3) = plt.subplots(ncols=3,nrows=1, figsize=(1920*px, 1080*px))

# Generating the x-axis for our plot
x = np.arange(0,extent)

# Because the recurrence grows so fast, in order to see a reasonable view of the plot we need to apply some limit to the growth rate
growth_rate_limiter = lambda x : x ** (1/27)

title = plt.suptitle(f"Approximating the recurrence: u(n) = {lam}u(n-1) {sign} {abs(mu)}u(n-2)")

# Because the animation is cleared and reset every frame, we need to initialise them BEFORE we begin the animation to avoid an undefined error (since we'd be
# removing the plots before they were initialised)
recurrence_from_0_to_extent = first_n_numbers_of_recurrence(extent) 
ax1plot_original = ax1.plot(x, [ growth_rate_limiter(x) for x in recurrence_from_0_to_extent] , color="green")    

# These will be lists of lists, with each list representing all of the y-values for a given frame
original_recurrence_frames = []

# Finite series frames
my_recurrence_frames = []

# Greg-Newton formula frames
gregs_frames = []

# Generating the lists of frames for each of the 3 functions
for t in range(extent):
    
    original_rec_frames = first_n_numbers_of_recurrence(t+4)
    my_recurrence_version = [recurrence(n, t) for n in range(0, extent)]
    greg_numbers = [gregory_newton(n, recurrence_from_0_to_extent, t) for n in range(0,extent)]
    
    original_recurrence_frames.append(original_rec_frames)
    my_recurrence_frames.append(my_recurrence_version)
    gregs_frames.append(greg_numbers)
    
    
# Because these sequences are discrete, trying to apply half-time intervals is much more difficult since the finite series only applies to integer n
# Which means if we just use the frame lists themselves
# Therefore we need to interpolate their points between the frames, and since matplotlib plots interpolate linearly by default I thought it a good idea
# The interpolation extent defines the ratio of frames after interpolation to frames before, so e.g an interpolation extent = 10 means we have 10x more frames
# after interpolation than before
interpolation_extent = 10

# This function takes our lists of frames and creates frame lists inbetween consecutive frames using linear interpolation to estimate the y-values
def interpolate(my_recurrence_frames):

    # Stores all frames created, both the original ones and the newly created ones
    my_interpolated_recurrence_frames = []

    # For each interval between 2 existing framelists, we need to create E-1 (where E is the interpolation extent) total framelists using linear interpolation
    for i in range(len(my_recurrence_frames) - 1):
    
        prev_framelist = my_recurrence_frames[i]
        next_framelist = my_recurrence_frames[i+1]
        
        # We will store a total of E framelists - the first is the one already made, the next (E-1) will be generated     
        interpolation_lists = [prev_framelist]
        
        # This loop is for creating the new framelists, with j corresponding to each one
        for j in range(1, interpolation_extent-1 + 1):
            
            # The next interpolation list to create
            the_next_new_interpolation_list = []
            
            # This loop covers every frame in each individual framelist, this is the most atomic level of interpoltion
            for k, value in enumerate(prev_framelist):
                
                current_frame_from_prevframelist = prev_framelist[k]
                current_frame_from_nextframelist = next_framelist[k]
                # Get the difference between the previous frame at this value and the next frame at this value so we can interpolate
                diff = current_frame_from_nextframelist - current_frame_from_prevframelist
                
                # Performs the interpolation
                the_next_new_frame = current_frame_from_prevframelist + (j/interpolation_extent) * diff
                
                the_next_new_interpolation_list.append(the_next_new_frame)
                
            interpolation_lists.append(the_next_new_interpolation_list)#
            
        # Add all of the new lists to the list of interpolation lists    
        my_interpolated_recurrence_frames += interpolation_lists
        
    # Once finished, add the one at the end to complete it
    my_interpolated_recurrence_frames.append( my_recurrence_frames[-1] )
    
    return my_interpolated_recurrence_frames

# Creating new framelists that are interpolated so we can have higher FPS
original_recurrence_frames = interpolate(original_recurrence_frames)
my_recurrence_frames = interpolate(my_recurrence_frames)
gregs_frames = interpolate(gregs_frames)
        
# Y-axis vertical limits, somehow x^(1/27) creates complex values, even though it shouldn't since it's just cbrt 3 times? (hence the .real)
top = max( [growth_rate_limiter(x).real for x in original_recurrence_frames[-1]]) 
bottom = min( [growth_rate_limiter(x).real for x in original_recurrence_frames[-1]])

# Need to specify what the base cases are, I didn't want to put them in the title since it was too crowded
text = plt.text(0.08, -0.2, f"u(0) = {u0}, u(1) = {u1}",fontsize=12    , ha="left", va = "bottom", transform=ax1.transAxes)

def animate(t):

    global my_reccurrence_frames
    global gregs_frames

    # Clear all the plots from the previous frame
    ax1.cla()
    ax2.cla()
    ax3.cla()
    
    # Therefore we need to reinitialise everything, including the grid
    ax1.grid(True)
    ax2.grid(True)
    ax3.grid(True)
    
    # Need to reinitialise our text as well
    text = plt.text(0.1, -0.1, f"u(0) = {u0}, u(1) = {u1}",fontsize=12    , ha="left", va = "bottom", transform=ax1.transAxes)
    

    title = plt.suptitle(f"Approximating the recurrence: u(n) = {lam}u(n-1) {sign} {abs(mu)}u(n-2), {(t/interpolation_extent):.2f} iterations", fontsize=20)
    
    ax1.set_ylabel("u(n)", fontsize=15)
    ax1.set_xlabel("n (original recurrence)" , fontsize=15)

    ax2.set_xlabel("Finite series approximation", fontsize= 15)
    ax3.set_xlabel("Gregory-Newton approximation",fontsize = 15)

    ax3.yaxis.set_label_position("right")
    ax3.set_ylabel(f"(data coded: f(x) = x^(1/27)")


    # Getting the next framelist for the original recurrence
    next_orig = original_recurrence_frames[t]

    # Plotting the graphs themselves, using the x-axis and the growth limiter
    original_rec_values = ax1.plot(np.arange(0,len(next_orig)) , [growth_rate_limiter(x) for x in next_orig], color="green")
    my_rec_values = ax2.plot(x, [growth_rate_limiter(x) for x in my_recurrence_frames[t] ], color="red")
    greg_formula = ax3.plot(x, [ growth_rate_limiter(x) for x in gregs_frames[t]], color = "blue")
    
    # Without the limits, the graphs look like they are becoming distorted over time
    ax1.set_ylim(bottom,top)
    ax2.set_ylim(bottom, top )
    ax3.set_ylim(bottom, top)
    ax1.set_xlim(0,30)
    ax2.set_xlim(0,30)
    ax3.set_xlim(0,30)



anim = FuncAnimation(fig, animate, frames = extent*(interpolation_extent ) - interpolation_extent, interval=1, repeat_delay = 2500)


anim.save("first.gif", fps=100 )



            
            
            
    