#!/usr/bin/env python
# coding: utf-8

# In[88]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
import pandas as pd
import math


# In[89]:


def features(x,y):
    
    '''
    Extracts the features from a dataset defined by (x,y)
    '''
    
    x_orig = x.copy()                 # The original list of x-values without changes
    
    x = strict_increase(x)            # Makes sure that x is a strict increase (necessary if repeats are present)
    
    m = len(x)                        # Number of data points
    
    '''
    First, a smooth curve is generated
    '''
    
    if m == 2:                        # If there are two points, a linear polynomial is used
        sp = polynomial(x,y,1)
    
    elif m == 3:                      # If there are three points, a quadratic polynomial is used
        sp = polynomial(x,y,2)
    
    elif m == 4:                      # If there are four points, a cubic polynomial is used
        sp = polynomial(x,y,3)
        
    else:                             # If m>= 5, a cubic spline is generated (see the function 'smoothcurve' for more information)
        sp = smoothcurve(x,y)
    
    '''
    The shape of this curve is extracted and based on this, it can be seen whether the curve needs to be simplified or not
    '''
    
    shape1, shape2, extremes = primitives(sp,x)      # Lists containing the information on the shape of the curve (see the function 'primitives' for more information)
    p_string = primitivesstring(shape1,shape2)       # List of strings with letters representing the primitives (see the function 'primitivesstring' for more information)
    too_complex = primcomplex(p_string)              # Based on the list of strings, it decided whether the curve is too complex or not (see the function 'too_complex' for more information)
    
    sp = simplifysmoothcurve(x,y,sp,too_complex)     # The curve can now be simplified if necessary (see the function 'simplifysmoothcurve' for more information)
    
    A, B, C, sign_x = log_fit(x_orig, y)
    
    choose_log = log_simplify(x, y, sp, A, B, C, sign_x)
    
    if choose_log == True:
        
        shape1, shape2, extremes = log_primitives(A, B, C, sign_x, x)      # Shapes of the simplified curve
        p_string = primitivesstring(shape1,shape2)                         # Primitives of the simplified cuve
    
        log_primitivesvisual(x, y, A, B, C, sign_x, extremes, p_string)    # The curve is now visualised (see the function 'primitivesvisual' for more information)
    
        extremes = np.array(extremes)
    
        y_values = A*np.log(sign_x*(extremes - C)) + B
        
    else:
    
        shape1, shape2, extremes = primitives(sp,x)      # Shapes of the simplified curve
        p_string = primitivesstring(shape1,shape2)       # Primitives of the simplified cuve
    
        primitivesvisual(x,y,sp,extremes,p_string)       # The curve is now visualised (see the function 'primitivesvisual' for more information)
    
        extremes = np.array(extremes)
    
        y_values = sp(extremes)                          # y-values at the extremes between the primitives
    
    extrema = extremum(p_string)
    
    '''
    Important information is now printed: the primitives, the x-values at their extremes and the y-values at their extremes
    '''
    
    print("Primitives:", p_string)
    print()
    
    print("Extremes of the primitives:", extremes)
    print()
    
    print("Value at extremes of the primitives:", y_values)
    print()
    
    print("Extrema:", extrema)
    print()
    
    return p_string, extremes, y_values, extrema


# In[90]:


def strict_increase(x):
    
    '''
    Python's UnivariateSpline algorithm only works if x is strictly increasing, so it does not work when repeats are present
    
    To avoid this, a very small value is added to x-values if it is equal to another x-value
    
    For example:
    
    x = [1, 1, 1, 2]
    
    x[1] is first compared to x[0]: they are equal, so x[1] is increased by a small value
    x[2] must then be compared to x[0] as well, as x[1] was changed
    
    This means that x[i] must be compared to x[i-1-nr_repeats] with nr_repeats the number of repeats that were already found
    '''
    
    delta_x = x[-1] - x[0]           # x_max - x_min, the range of x-values
    
    nr_repeats = 0                   # Number of repeats of a certain x-value
    
    for i in range(1,len(x)):
        
        if x[i] == x[i-1-nr_repeats]:                   # x[i] is compared to x[i-1-nr_repeats] 
            
            nr_repeats +=1                              # The number of repeats increases by 1
            
            x[i] += nr_repeats * delta_x * 10**-10      # Delta_x * 10^(-10) is added to the value, multiplied with nr_repeats to make sure that every repeat has a different value
        
        else:
            
            nr_repeats = 0           # nr_repeats is set to 0 once a new x-value is encountered
    
    return x


# In[91]:


def polynomial(x,y,k):
    
    '''
    Draw a curve through a dataset defined by (x,y) using a polynomial
    
    If k = 1: Straight line is fitted through data
    If k = 2: Quadratic polynomial is fitted through data
    If k = 3: Cubic polynomial is fitted through data
    '''
    
    m = len(x)                       # Number of points in dataset
    s = m                            # Default smoothing factor is equal to the number of data points m
    done = False                     # Determines if the algorithm is finished or not
    
    '''
    A first spline is calculated
    If it has two knots, the algorithm is finished
    '''
    
    sp = UnivariateSpline(x,y,None,[None,None],k,s)           # A new spline is calculated using the new residual s
    nr_knots = len(sp.get_knots())                            # Number of knots of the new spline
    
    if nr_knots == 2:
        
        return sp
        
    '''
    If the number of knots is higher than 2, a higher smoothing factor is needed
    The smoothing factor is doubled until the number of knots equals 2
    '''
    
    while not done:
        
        s *= 2                         # Residual is doubled
        
        sp = UnivariateSpline(x,y,None,[None,None],k,s)           # A new spline is calculated using the new residual s
        nr_knots = len(sp.get_knots())                            # Number of knots of the new spline
        
        if nr_knots == 2:                                         # Algorithm is finished if nr_knots equals 2
            done = True                                              
    
    return sp


# <font size="4"><center>Smooth curve algorithm:</center></font>
# 
# <br>
# <br>
# 
# <font size="3"><center>GOAL:</center></font>
# 
# <font size="3"><center>Out of all the splines in the figure, find the realistic spline with the lowest smoothing factor possible</center></font>
# 
# ![Considered%20Splines.PNG](attachment:Considered%20Splines.PNG)
#  
# <br>
# <br>
# <br>
# 
# <font size="3"><center>IMPLEMENTATION:</center></font>
# 
# ![SmoothCurve%20Shapes.PNG](attachment:SmoothCurve%20Shapes.PNG)

# In[92]:


def smoothcurve(x,y):
    
    '''
    Draw a smooth curve through a dataset defined by (x,y) using a cubic spline
    
    This is done using two loops:
    - One which halves the smoothing factor of the spline in each step
    - One which doubles the smoothing factor of the spline in each step
    
    The former loop is controlled by the Boolean 'done1' and the latter loop is controlled by the Boolean 'done2'
    '''
    
    k = 3                            # Degree of cubic spline
    
    m = len(x)                       # Number of points in dataset
    s = m                            # Default smoothing factor is equal to the number of data points m
    
    done1, done2 = False, False      # The loops are not done at the beginning
    
    '''
    First of all, a maximal number of knots must be set for the spline
    
    For a cubic spline, this is a function of the number of data points and is equal to or greater than 2
    
    The first two knots of a spline are always placed at the beginning and at the end of the considered interval
    New knots are placed one by one in the middle of the interval with the highest residual sum of squares
    '''
    
    if m <= 9:           
        
        nr_knots_max = m - 2         # If m <= 9 (low number of data points), the maximal number of knots is equal to m - 2
    
    else:                            # If m > 9, the maximal number of knots is set to 7, as after thorough testing,
                                     # it was verified that this is enough to capture possible shapes expected in catalytic datasets
        nr_knots_max = 7             
                                
    nr_knots_interpol = m - 2        # Number of knots needed for interpolation; it is used to make sure that algorithm stops once this number is reached
    
    '''
    A first spline is calculated using the default smoothing factor equal to the number of data points
    '''
    
    
    sp_orig = UnivariateSpline(x,y,None,[None,None],k,s)   # First spline with smoothing factor s = m
    nr_knots = len(sp_orig.get_knots())                    # Number of knots of the spline
    
    shape1, shape2, extremes = primitives(sp_orig,x)       # Shapes of the spline
    p_string = primitivesstring(shape1,shape2)             # Primitives of the spline
    too_complex_orig = primcomplex(p_string)               # Indicates whether the spline is too complex or not
    
    if nr_knots > nr_knots_max:       # If the number of knots surpasses nr_knots_max, the loop which halves the smoothing factor is not entered
        
        done1 = True
        
    if nr_knots == 2:                 # If the number of knots equals 2, the loop which doubles the smoothing factor is not entered
        
        done2 = True
        
    '''
    The first spline has now been calculated using the default value of the smoothing factor
    A while-loop is now started to change this smoothing factor
    
    In the first loop, the number of knots of the first spline is lower than or equal to nr_knots_max and
    the smoothing factor is halved until nr_knots_max is surpassed
    
    In each step, the spline is stored in the list 'splines'
    '''
    
    if done1 == False:                       # The loop is entered if 'done1' is False
        
        splines = [sp_orig]                  # List of splines
        complexity = [too_complex_orig]      # List of Booleans indiciting whether the spline is too complex or not
        
        while not done1:
        
            s /= 2                           # Residual is halved
        
            sp = UnivariateSpline(x,y,None,[None,None],k,s)           # A new spline is calculated using the new residual s
            nr_knots = len(sp.get_knots())                            # Number of knots of the new spline
    
            shape1, shape2, extremes = primitives(sp,x)    # Shapes of the new spline
            p_string = primitivesstring(shape1,shape2)     # Primitives of the new spline
            too_complex = primcomplex(p_string)            # Complexity of the new spline
        
            if nr_knots > nr_knots_max:               # Algorithm is finished if nr_knots_max is surpassed
            
                done1 = True                          # The last calculated spline has a number of knots which is too high, so it is not added to the list of splines

            elif nr_knots == nr_knots_interpol:       # Algorithm is finished if the number of knots equals the number of knots needed for interpolation
            
                done1 = True                          # The last calculated spline has a number of knots which is not too high, so it is added to the list of splines
                
                splines.append(sp)
                complexity.append(too_complex)
        
            else:                                     # Otherwise, the spline is added to the list of splines
            
                splines.append(sp)
                complexity.append(too_complex)
                
        '''
        The list of splines has been generated, so now the algorithm looks for the index of the last spline which is not considered unrealistic,
        the last spline with too_complex == False
        '''
     
        index_optimal_spline = -1                     # The algorithm starts with an index outside the range of indices
    
        for i in range(0, len(splines)):              # Each spline in the list of splines is considered
        
            if complexity[i] == False:                # If the spline has a realistic shape, the index is updated
            
                index_optimal_spline = i
    
        if index_optimal_spline != -1:                # If an index was found, the corresponding spline is returned
    
            return splines[index_optimal_spline]
    
    '''
    In the second loop, the smoothing factor is halved until it has two knots
    A cubic spline with two knots is a cubic polynomial, so further changes are impossible
    
    The first calculated spline which has a realistic shape is returned, as that is the most accurate spline with a realistic shape
    '''
    
    if done2 == False:                       # The loop is entered if 'done2' is False
    
        while not done2:
        
            s *= 2                           # Residual is halved
        
            sp = UnivariateSpline(x,y,None,[None,None],k,s)            # A new spline is calculated using the new residual s
            nr_knots = len(sp.get_knots())                             # Number of knots of the new spline
        
            shape1, shape2, extremes = primitives(sp,x)    # Shapes of the new spline
            p_string = primitivesstring(shape1,shape2)     # Primitives of the new spline
            too_complex = primcomplex(p_string)            # Complexity of the new spline        
        
            if nr_knots <= nr_knots_max and too_complex == False:    # If a spline is found with a number of knots below nr_knots_max and with a realistic shape, it is returned 
            
                return sp
            
            elif nr_knots == 2:                            # Algorithm is finished if the number of knots equals 2
                
                done2 = True
    
    '''
    If no spline was found with a realistic shape, sp_orig (default smoothing factor) is returned.
    '''
    
    return sp_orig


# <font size="3"><center>Simplify smooth curve algorithm:</center></font>
# 
# ![Simplify.png](attachment:Simplify.png)
# 
# ![Simplify%20Num.png](attachment:Simplify%20Num.png)
# 

# In[93]:


def simplifysmoothcurve(x,y,sp,too_complex):
    
    '''
    Generates cubic spline for a single dataset and simplifies it to a linear, quadratic or cubic polynomial if necessary
    '''
    
    m = len(x)                    # Number of data points
    n = len(sp.get_knots())       # Number of knots
    
    shape1_s, shape2_s, extremes_s = primitives(sp,x)    # Shape of the original curve
    Rs = rvalue(sp,x)                                    # R-value of the original curve
    
    p_string_s = primitivesstring(shape1_s, shape2_s)
    nr_ABCD_s = p_string_s.count("A") + p_string_s.count("B") + p_string_s.count("C") + p_string_s.count("D")
    
    if n == 2 and m >= 4:         # If there are only two knots and m >= 4, the curve is a cubic polynomial
        R3 = Rs                   # The R-value of the cubic polynomial is the same as the one calculated before
    
    elif m == 3:                  # If there are only two knots and m = 3, the curve is a quadratic polynomial
        R2 = Rs                   # The R-value of the quadratic polynomial is the same as the one calculated before
    
    '''
    If the number of knots is higher than two, the original curve is a cubic spline, so in the following section, the
    current curve is compared to a cubic polynomial
    
    If the number of knots is equal to two, this part of the code is skipped as the curve is already a cubic polynomial
    or a polynomial of a lower degree
    '''
    
    if n != 2:
        
        sp3 = polynomial(x,y,3)                                  # Cubic polynomial
        shape1_3, shape2_3, extremes_3 = primitives(sp3,x)       # Shape of the cubic polynomial
        R3 = rvalue(sp3,x)                                       # R-value of the cubic polynomial
        
        p_string_3 = primitivesstring(shape1_3, shape2_3)
        nr_ABCD_3 = p_string_3.count("A") + p_string_3.count("B") + p_string_3.count("C") + p_string_3.count("D")
        
        if too_complex == True:                                # If too_complex is True, the simplification is forced
            
            sp = sp3                                           # sp is updated with 'sp3'
            shape1_s, shape2_s = shape1_3, shape2_3
            Rs = R3
            too_complex = primcomplex(primitivesstring(shape1_s,shape2_s))
        
        elif R3/Rs < 0.998 or nr_ABCD_3 > nr_ABCD_s:           # If the difference in R-value is too big, the curve is not simplified and the algorithm is finished
                                                               # (the value of 0,998 was obtained after testing many examples)
            return sp                                          # The algorithm is also stopped if the cubic polynomial is more complex than sp (if it is represented by more primitives A, B, C and D)
        
        elif shape1_3 != shape1_s or shape2_3 != shape2_s:     # If the shape of the cubic polynomial and sp are different, sp is updated
                                                               # If they have the same shape, sp is not updated to keep the more accurate curve
            sp = sp3
            shape1_s, shape2_s = shape1_3, shape2_3
            Rs = R3
            too_complex = primcomplex(primitivesstring(shape1_s,shape2_s))
    
    '''
    The curve sp is now compared to a quadratic polynomial
    
    This step is skipped if the curve is already quadratic or linear (m < 4)
    '''
    
    if m >= 4:
        
        sp2 = polynomial(x,y,2)                                # Quadratic polynomial
        shape1_2, shape2_2, extremes_2 = primitives(sp2,x)     # Shape of the quadratic polynomial
        R2 = rvalue(sp2,x)                                     # R-value of the quadratic polynomial
        
        p_string_2 = primitivesstring(shape1_2, shape2_2)
        nr_ABCD_2 = p_string_2.count("A") + p_string_2.count("B") + p_string_2.count("C") + p_string_2.count("D")
        
        if too_complex == True:                                # If too_complex is True, the simplification is forced
            
            sp = sp2                                           # sp is updated with 'sp2'
            shape1_s, shape2_s = shape1_2, shape2_2
            Rs = R2
            too_complex = primcomplex(primitivesstring(shape1_s,shape2_s))
            
        elif R2/R3 < 0.998 or nr_ABCD_2 > nr_ABCD_s:           # If the difference in R-value is too big, the curve is not simplified and the algorithm is finished
                                                               # (the value of 0,998 was obtained after testing many examples)
            return sp                                          # The algorithm is also stopped if the quadratic polynomial is more complex than sp (if it is represented by more primitives A, B, C and D)
        
        elif shape1_2 != shape1_s or shape2_2 != shape2_s:     # If the shape of the quadratic polynomial and sp are different, sp is updated
                                                               # If they have the same shape, sp is not updated to keep the more accurate curve
            sp = sp2
            shape1_s, shape2_s = shape1_2, shape2_2
            Rs = R2
            too_complex = primcomplex(primitivesstring(shape1_s,shape2_s))
            
    '''
    The curve sp is now compared to a linear polynomial
    
    This step is skipped if the curve is already  linear (m = 2)
    '''
    
    if m >= 3:
        
        sp1 = polynomial(x,y,1)                               # Linear polynomial
        shape1_1, shape2_1, extremes_1 = primitives(sp1,x)    # Shape of the linear polynomial
        R1 = rvalue(sp1,x)                                    # R-value of the linear polynomial
        
        if len(shape1_s) == 1 and shape1_s[0] == shape1_1[0] and shape2_s[0] == shape2_1[0]:      # If sp is already linear, it is simplified to a linear polynomial for visualisation
            return sp1
        
        if R1/R2 < 0.998:         # The curve is simplified if the difference in R-value is small
            return sp             # Otherwise, sp is kept
        
        else:
            return sp1
    
    return sp                     # Nothing is done when m = 2


# In[94]:


def rvalue(sp, x):
   
    '''
    Calculates multiple correlation coefficient R of a spline fit
    '''
    
    ssy = np.sum(np.square(sp(x)))             # Sum of squares of prediction
    sse = sp.get_residual()                    # Sum of squares of error
    R = np.sqrt(ssy/(sse+ssy))                 # Multiple correlation coefficient
    
    return R


# In[95]:


def primitives(sp, x):
    
    '''
    Returns two lists storing the evolution of the signs of the first and second derivative of a spline, and a list with the extremes of the intervals
    
    For example:
    
    shape1 = [1, -1, 0, 1]
    shape2 = [1, 1, 0, 0]
    extremes = [0, 0.5, 1, 2.5, 4]
    
    Between x-values 0 and 0.5: first derivative is positive (1), second dertivative is positive (1)
    Between x-values 0.5 and 1: first derivative is negative (-1), second derivative is positive (1)
    Between x-values 1 and 2.5: first derivative is close to zero (0), second dertivative is close to zero (0)
    Between x-values 2.5 and 4: first derivative is positive (1), second derivative is close to zero (0)
    '''
    
    m = len(x)                                # Number of data points
    n = len(sp.get_knots())                   # Number of knots of the spline
    
    '''
    The first step is to determine in which points the derivatives will be calculated and to determine which values will be
    used to normalise the derivatives.
    
    The derivatives are calculated in 100 equidistant points in the considered interval (100 is enough to extract the primitives)
    
    The first derivative is normalised by multiplying it with (x_max - x_min)/(y_max - y_min)
    The second derivative is normalised by multiplying it with (x_max - x_min)/(y_max - y_min)
    '''
    
    xs = np.linspace(x[0],x[m-1],100)         # Points where derivatives are calculated
    x_m = x[m-1] - x[0]                       # x_max - x_min
    y_m = max(sp(xs)) - min(sp(xs))           # y_max - y_min
    
    '''
    A function exists to calculate the derivatives of the UnivaririateSpline objects, but the k'th derivative can only be
    calculated if the degree of the spline is at least k, so for a linear spline the second derivative cannot be calculated.
    If this is still done, a ValueError will be returned
    
    The following code is used for splines of a degree of at least 2. For linear splines, an error will be returned and
    the code goes to a later section of the code which is meant to handle exceptions (except ValueError)
    '''
    
    try:
        
        der1_sp, der2_sp = sp.derivative(1), sp.derivative(2)           # First and second derivative of spline
        der1, der2 = der1_sp(xs[0]), der2_sp(xs[0])                     # First and second derivative of spline in first point
        
        shape1, shape2 = [], []                         # Lists containing the signs
        extremes = []                                   # List containing the extremes of the intervals
        
        '''
    The signs of the derivatives are now calculated. If the absolute value of the normalised derivative is below a certain value,
    it is considered equal to 0
    
    If the absolute value of the first derivative is lower than 0.5, it is considered 0
    If the absolute value of the second derivative is lower than 1, it is considered 0
    
    If the first derivative equals 0, both derivatives are set to 0
        '''
        
        if np.abs(der1)*x_m/y_m < 0.5:                # If the normalised first derivative is lower than 0.5, it's considered 0
            sign_der1 = 0
            
        else:
            sign_der1 = int(np.sign(der1))
            
        if np.abs(der2)*np.square(x_m)/y_m < 0.5 or sign_der1 == 0:       # If the normalised derivative is lower than 1, it's considered 0
            sign_der2 = 0                                               # If the first derivative equals 0, both derivatives are set to 0
        
        else:
            sign_der2 = int(np.sign(der2))
        
        '''
    A for-loop is now started to calculate the derivative in all points and based on this, the signs are added to the lists of signs
    
    The signs are added to lists of signs when at least two consecutive points have the same sign
    To do this, the value of 'sign_change' becomes equal to True when the sign changed in the previous point
        '''
        
        sign_change = True                                              # If the sign changes, this becomes True. It is used to make sure that two consecutive point have the same primitive
        
        for i in range(1,100):                                          # If one of the signs changes, the lists are updated with the new signs
            
            der1_new, der2_new = der1_sp(xs[i]), der2_sp(xs[i])         # The derivatives in point i      
            
            if np.abs(der1_new)*x_m/y_m < 0.5:                          # The previous equations for the derivatives are repeated for point i
                sign_der1_new = 0
            
            else:
                sign_der1_new = int(np.sign(der1_new))
            
            if np.abs(der2_new)*np.square(x_m)/y_m < 0.5 or sign_der1_new == 0:
                sign_der2_new = 0
            
            else:
                sign_der2_new = int(np.sign(der2_new))
            
            if sign_der1_new != sign_der1 or sign_der2_new != sign_der2:    # If sign changes, sign_change becomes True, or when i == 1 as this is the first point in the loop
                sign_change = True
            
            elif sign_change == True:                                   # Two consecutive points have the same derivative
                                                                        # (this elif-statement can only be reached if the previous once was False,
                if shape1 == []:                                        # this means that there was no sign change this time)
                    extremes.append(x[0])
                    shape1.append(sign_der1_new)                        # If the shape-list is still empty, x[0] (the first point) is added to the list 'extremes'
                    shape2.append(sign_der2_new)
                
                elif shape1[-1] != sign_der1_new or shape2[-1] != sign_der2_new:
                    extremes.append(xs[i-1])                 
                    shape1.append(sign_der1_new)                        # If the shape-list is not empty, point i-1 is added to the list 'extremes'
                    shape2.append(sign_der2_new)                        # Neither of the if-statement are reached if signs are the same as for the signs in the shape-list
                                                                        # as there is no need to have the same primitive twice consecutively
                
                sign_change = False         # 'sign_change' is set to False until a new sign change occurs
            
            sign_der1, sign_der2 = sign_der1_new, sign_der2_new         # Stored as der1 and der2 for the following iteration
        
        extremes.append(x[m-1])               # Final x-value is also added to the list of extremes
        
        '''
    In some cases, a flat section (derivatives both 0) represents a maximum or a minimum and it would be more intuituive
    to remove it and keep the sections next to it
        '''
        
        shape1_new, shape2_new, extremes_new = shape1.copy(), shape2.copy(), extremes.copy()
        nr_removed = 0
        
        if len(shape1) >= 3:                              # The shape-lists need to contain at least three signs for a maximum or a minimum
            
            for i in range(0,len(shape1)-2):              # signs i are compared to signs i+1 and i+2
                
                if shape1[i] == 1 and shape2[i] == -1:            # Maximum is described by first derivative: [1, 0, -1] and second derivative: [-1, 0, 1]
                    
                    if shape1[i+1] == 0 and shape2[i+1] == 0 and shape1[i+2] == -1 and shape2[i+2] == -1:
                        
                        del shape1_new[i+1-nr_removed]            # The flat sections are deleted
                        del shape2_new[i+1-nr_removed]
                        
                        t_new = (extremes[i+1] + extremes[i+2])/2     # The two extremes of the flat section are replaced by one value in the middle
                        del extremes_new[i+2-nr_removed]
                        extremes_new[i+1-nr_removed] = t_new
                        nr_removed += 1
                
                if shape1[i] == -1 and shape2[i] == 1:            # Minimum is described by first derivative: [-1, 0, 1] and second derivative [1, 0, -1]
                    
                    if shape1[i+1] == 0 and shape2[i+1] == 0 and shape1[i+2] == 1 and shape2[i+2] == 1:
                        
                        del shape1_new[i+1-nr_removed]            # The flat sections are deleted
                        del shape2_new[i+1-nr_removed]
                        
                        t_new = (extremes[i+1] + extremes[i+2])/2     # The two extremes of the flat section are replaced by one value in the middle
                        del extremes_new[i+2-nr_removed]
                        extremes_new[i+1-nr_removed] = t_new
                        nr_removed += 1
                        
        shape1, shape2, extremes = shape1_new, shape2_new, extremes_new
    
        '''
    The following code is used for linear polynomials; in this case, the polynomial is a straight line, which can
    be represented by one shape
    
    The sign of the first derivative is calculated, while the sign of the second derivative is always 0
        '''
    
    except ValueError:                             # ValueError occurs when linear spline is used
        
        der1_sp = sp.derivative(1)
        der1 = der1_sp(xs[0])                      # First derivative of spline in first point
        
        shape1, shape2 = [], []                    # List containing the signs
        extremes = []
        
        if np.abs(der1)*x_m/y_m <= 0.5:            # Only the first derivative is calculated
            sign1 = 0
        
        else:
            sign1 = int(np.sign(der1))
        
        shape1.append(sign1)                       # The signs are added to the shape-lists
        shape2.append(0)
        extremes.append(x[0])                      # The extremes of the considered interval are added to the 'extremes'-list
        extremes.append(x[m-1])
    
    return shape1, shape2, extremes


# <font size="3"><center>Possible primitives:</center></font>
# 
# ![PrimitiveS.PNG](attachment:PrimitiveS.PNG)

# In[96]:


def primitivesstring(shape1, shape2):
    
    '''
    Creates list of strings with letters representing the primitives of the curve
    
    '''
    
    primitives = []         # List storing the strings
    n = len(shape1)         # Number of primitives
    
    for i in range(0,n):
        
        if shape1[i] == 1 and shape2[i] == -1:        # Primitive A is represented by the signs 1 and -1
            primitives.append('A')
        
        elif shape1[i] == -1 and shape2[i] == 1:      # Primitive A is represented by the signs -1 and 1
            primitives.append('B')
        
        elif shape1[i] == 1 and shape2[i] == 1:       # Primitive A is represented by the signs 1 and 1
            primitives.append('C')
        
        elif shape1[i] == -1 and shape2[i] == -1:     # Primitive A is represented by the signs -1 and -1
            primitives.append('D')
        
        elif shape1[i] == 1 and shape2[i] == 0:       # Primitive A is represented by the signs 1 and 0
            primitives.append('E')
        
        elif shape1[i] == -1 and shape2[i] == 0:      # Primitive A is represented by the signs -1 and 0
            primitives.append('F')
        
        else:                                         # Primitive G is represented by the signs 0 and 0
            primitives.append('G')
    
    return primitives


# In[97]:


def primitivesvisual(x,y,sp,extremes,p_string):
    
    '''
    Visualises the features
    '''
    
    m = len(x)           # Number of data points
    n = len(extremes)      # Number of extremes
    
    xs = np.linspace(x[0],x[m-1],1000)       # Points to plot
    
    y_min, y_max = min([min(y),min(sp(xs))]), max([max(y),max(sp(xs))])    # Smallest and highest y-values reached by the spline and the data points 
    delta_y = y_max - y_min                                                # The range of y-values
    
    plt.plot(x, y, 'ro', ms=5)               # Data points are plotted in red
    
    '''
    The data points have been plotted; now, the spline and its features need to plotted
    '''
    
    for i in range(0,n-1):
        
        xs = np.linspace(extremes[i],extremes[i+1],1000)       # Points to plot for one primitive
        
        if p_string[i] == 'A':                             # If the primitive is A, the curve is plotted in green
            plt.plot(xs, sp(xs), 'g', lw=3)
        
        elif p_string[i] == 'B':                           # If the primitive is B, the curve is plotted in yellow
            plt.plot(xs, sp(xs), 'y', lw=3)
        
        elif p_string[i] == 'C':                           # If the primitive is C, the curve is plotted in blue
            plt.plot(xs, sp(xs), 'b', lw=3)
        
        elif p_string[i] == 'D':                           # If the primitive is D, the curve is plotted in red
            plt.plot(xs, sp(xs), 'r', lw=3)
        
        elif p_string[i] == 'E':                           # If the primitive is E, the curve is plotted in purple
            plt.plot(xs, sp(xs), color='#9400d3', lw=3)
        
        elif p_string[i] == 'F':                           # If the primitive is F, the curve is plotted in orange
            plt.plot(xs, sp(xs), color='#ff8c00', lw=3)
        
        else:                                              # If the primitive is G, the curve is plotted in grey
            plt.plot(xs, sp(xs), color='#808080', lw=3)
        
        xt = (extremes[i]+extremes[i+1])/2 - (extremes[n-1]-extremes[0])/40    # x-position where letter will be placed on graph
        
        if delta_y != 0:
            
            yt = y_min - 0.2*delta_y                       # y-position where letter will be placed on graph
 
        elif y_min != 0:

            yt = 0.95*y_min
            
        else:
            
            yt = -5
        
        plt.text(xt,yt,p_string[i],fontsize=15)            # Letter representing primitive is printed on graph
    
    '''
    The black vertical lines are now plotted between the primitives
    '''
    
    if delta_y != 0:
        
        yp = np.linspace(y_min - 0.2*delta_y, y_max + 0.2*delta_y, 1000)         # y points to plot in the vertical line
    
    elif y_min != 0:
        
        yp = np.linspace(0.95*y_min, 1.05*y_max, 1000)
        
    else:
        
        yp = np.linspace(-5, 5, 1000)
    
    for i in range(0,n):
        
        xp = [extremes[i]]*1000                # List with 1000 times the x-position of the line
        
        plt.plot(xp, yp, 'k', lw=2)          # Vertical line is plotted


# <font size="3"><center>Smooth curve is considered too complex if one of the following conditions is satisfied:</center></font>
# 
# ![TooComplex.png](attachment:TooComplex.png)

# In[98]:


def primcomplex(p_string):
    
    '''
    Determines if curve is too complex or not
    '''
    
    nA = p_string.count('A')       # Number of A-primitives
    nB = p_string.count('B')       # Number of B-primitives
    nC = p_string.count('C')       # Number of C-primitives
    nD = p_string.count('D')       # Number of D-primitives
    
    if nA + nB + nC + nD > 2:                         # The curve is considered complex if there are more than two primitives of the types A, B, C or D
        return True
    
    if nA > 1 or nB > 1 or nC > 1 or nD > 1:          # The curve is considered complex if there is more than one primitive of same the type of the types A, B, C or D
        return True
    
    if nA + nB == 2 or nC + nD == 2:                  # The curve is considered complex if there are both A and B primitives, or both C and D primitives
        return True
    
    return False


# In[99]:


def extremum(p_string):
    
    '''
    Determines if there is an extremum
    '''
    
    n = len(p_string)
    
    extrema = [" "]
    
    for i in range(1, n):
        
        if p_string[i-1] == "B" and p_string[i] == "C":
            
            extrema.append("minimum")
            
        elif p_string[i-1] == "A" and p_string[i] == "D":
            
            extrema.append("maximum")
        
        else:
            
            extrema.append(" ")
    
    extrema.append(" ")
    
    return extrema


# In[100]:


def features_parametric(list_x,list_y):
    
    '''
    Extracts features from parametric datasets
    
    list_x is a list containing the x-values of all datasets
    list_y is a list containing the y-values of all datasets
    
    Dataset i has x-values list_x[i] and y-values list_y[i]
    '''
    
    n_list = len(list_x)        # Number of datasets
    list_x_orig = []
    
    list_sp = []                # The splines for each dataset will be stored in this list
    list_shape1 = []            # The list shape1 will be stored in this list for each dataset
    list_shape2 = []            # The list shape2 will be stored in this list for each dataset
    list_extremes = []          # The list extremes will be stored in this list for each dataset
    list_too_complex = []       # The list will indicate which splines are too complex
    
    '''
    To access the objects corresponding to dataset i, the index i must be used in the lists
    
    First, a smooth curve is generated
    '''
    
    for i in range(0,n_list):                                # Smooth curve is calculated for each dataset
        
        x_orig = list_x[i]
        list_x_orig.append(x_orig.copy())
        
        list_x[i] = strict_increase(list_x[i])               # Makes sure that x has a strict increase (necessary if repeats are present)
        
        m = len(list_x[i])                                   # Number of data points in dataset i
        
        if m == 2:
            
            sp = polynomial(list_x[i], list_y[i],1)          # If two points: linear polynomial
            
        elif m == 3:
            
            sp = polynomial(list_x[i], list_y[i],2)          # If three points: quadratic polynomial
            
        elif m == 4:
            
            sp = polynomial(list_x[i], list_y[i], 3)         # If four points: cubic polynomial
            
        else:
            
            sp = smoothcurve(list_x[i],list_y[i])                      # In other cases, a cubic spline is used
            
        list_sp.append(sp)                                             # The spline is added to the list of splines
        
        shape1, shape2, extremes = primitives(sp,list_x[i])            # Shape lists of the splines
        
        list_shape1.append(shape1)                                     # shape1 is added to the list of shape1
        list_shape2.append(shape2)                                     # shape2 is added to the list of shape2
        list_extremes.append(extremes)
        too_complex = primcomplex(primitivesstring(shape1,shape2))     # Calculates if the spline is too complex
        list_too_complex.append(too_complex)                           # Stored in list for values of 'too_complex'
        
    '''
    The smooth curve can be simplified if possible
    '''
    
    list_sp = simplifysmoothcurve_parametric(list_x,list_y,list_sp,list_too_complex,list_shape1,list_shape2,list_extremes)    # Simplification algorithm
    
    list_shape1, list_shape2, list_extremes = [], [], []          # The shape-lists will be recalculated for the new curves
    list_p_string = []                                          # The letters representing the primitives will also be generated for each curve
    list_y_extremes = []
    list_extrema = []
    
    for i in range(0,n_list):
        
        sp = list_sp[i]
        
        shape1, shape2, extremes = primitives(sp, list_x[i])      # Shape lists of the splines
        
        list_shape1.append(shape1)                      # shape1 is added to the list of shape1
        list_shape2.append(shape2)                      # shape2 is added to the list of shape2
        list_extremes.append(extremes)                      # extremes is added to the list of extremes
        
        p_string = primitivesstring(shape1,shape2)      # Lists of letters representing primitives
        list_p_string.append(p_string)                  # Added to the list containg all these lists
        
        y_extremes = sp(extremes)
        list_y_extremes.append(y_extremes)
        
        extrema = extremum(p_string)
        list_extrema.append(extrema)
        
    list_A = []
    list_B = []
    list_C = []
    list_sign_x = []
    list_choose_log = []
        
    for i in range(0,n_list):
        
        A, B, C, sign_x = log_fit(list_x_orig[i], list_y[i])
        
        list_A.append(A)
        list_B.append(B)
        list_C.append(C)
        list_sign_x.append(sign_x)
    
        choose_log = log_simplify(list_x[i], list_y[i], list_sp[i], A, B, C, sign_x)
        
        list_choose_log.append(choose_log)
        
        if choose_log == True:
            
            shape1, shape2, extremes = log_primitives(A, B, C, sign_x, list_x[i])      # Shapes of the simplified curve
            p_string = primitivesstring(shape1, shape2)                                # Primitives of the simplified cuve
    
            extremes = np.array(extremes)
    
            y_extremes = A*np.log(sign_x*(extremes - C)) + B
        
            list_shape1[i] = shape1
            list_shape2[i] = shape2
            list_extremes[i] = extremes
            list_p_string[i] = p_string
            list_y_extremes[i] = y_extremes
    
    for i in range(0,n_list):
        
        nr = i+1
        
        if list_choose_log[i] == True:
            
            log_primitivesvisual_parametric(list_x[i], list_y[i], list_A[i], list_B[i], list_C[i], list_sign_x[i], list_extremes[i], list_p_string[i])
        
        else:
            
            primitivesvisual_parametric(list_x[i],list_y[i],list_sp[i],list_extremes[i],list_p_string[i])

        print("DATASET", nr)
        print("-----------")
        print()
        print("Primitives:", list_p_string[i])
        print()
        print("Extremes of the primitives:", list_extremes[i])
        print()
        print("Value at extremes of the primitives of:", list_y_extremes[i])
        print()
        print()
        
    return list_p_string, list_extremes, list_y_extremes, list_extrema


# In[101]:


def simplifysmoothcurve_parametric(list_x,list_y,list_sp,list_too_complex,list_shape1,list_shape2,list_start):
    
    '''
    Generates cubic spline for a single dataset and simplifies it to a linear, quadratic or cubic polynomial if necessary
    '''
    
    n_list = len(list_x)
    
    list_R = []
    list_nr_ABCD = []
    
    if list_too_complex.count(True) == 0:             # Check if trend is present (only if none of the curves are too complex)
        
        if list_shape1.count(list_shape1[0]) == n_list and list_shape2.count(list_shape2[0]) == n_list and list_shape2[0] != [0]:    # Trend if all curves heve the same primitives
            return list_sp
        
    for i in range(0,n_list):                   # R-value is calculated for each curve
        
        R = rvalue(list_sp[i],list_x[i])
        list_R.append(R)
        
        p_string = primitivesstring(list_shape1[i], list_shape2[i])
        nr_ABCD = p_string.count("A") + p_string.count("B") + p_string.count("C") + p_string.count("D")
        list_nr_ABCD.append(nr_ABCD)
        
    '''
    Step 1: The curve is compared to a cubic polynomial
    '''
        
    for i in range(0,n_list):       # Step 1 of the algorithm
        
        sp = list_sp[i]             # Curve i
        
        n = len(sp.get_knots())     # Number of knots of the curve
        
        if n != 2:                  # Step 1 is only performed if the curve is a cubic spline with more than two knots
            
            sp3 = polynomial(list_x[i],list_y[i],3)                     # Cubic polynomial
            
            shape1_3, shape2_3, start_3 = primitives(sp3,list_x[i])      # Shape of cubic polynomial
            
            R3 = rvalue(sp3,list_x[i])                                   # R-value of the cubic polynomial
            
            p_string_3 = primitivesstring(shape1_3, shape2_3)
            nr_ABCD_3 = p_string_3.count("A") + p_string_3.count("B") + p_string_3.count("C") + p_string_3.count("D")
            
            if list_too_complex[i] == True:                              # Simplification is forced if curve is too complex
                
                list_sp[i] = sp3                                         # All lists are updated with the cubic polynomial
                list_shape1[i], list_shape2[i] = shape1_3, shape2_3
                list_R[i] = R3
                list_too_complex[i] = primcomplex(primitivesstring(list_shape1[i],list_shape2[i]))
                list_nr_ABCD[i] = list_nr_ABCD_3
                
            elif R3/list_R[i] < 0.998 or nr_ABCD_3 > list_nr_ABCD[i]:     # the cubic spline is kept if the difference in R-value is large or the cubic polynomial has more primitives than the cubic spline
                
                continue                                                          # Do nothing                        
                
            elif shape1_3 != list_shape1[i] or shape2_3 != list_shape2[i]:        # The cubic polynomial is used unless the cubic spline and the cubic polynomial have the same shape
                
                list_sp[i] = sp3                                          # All lists are updated with the cubic polynomial
                list_shape1[i], list_shape2[i] = shape1_3, shape2_3
                list_R[i] = R3
                list_too_complex[i] = primcomplex(primitivesstring(list_shape1[i],list_shape2[i]))
                list_nr_ABCD[i] = list_nr_ABCD_3
                
    if list_too_complex.count(True) == 0:              # If a trend is found, the algorithm is stopped
        
        if list_shape1.count(list_shape1[0]) == n_list and list_shape2.count(list_shape2[0]) == n_list and list_shape2[0] != [0]:
            return list_sp
        
    '''
    Step 2: The curve is now compared to a quadratic polynomial
    '''
        
    for i in range(0,n_list):
        
        sp = list_sp[i]         # Curve i
        
        m = len(list_x[i])      # Number of data points of the curve
        
        if m >= 4:              # Step 2 is only performed if the curve has at least 4 points (the curve is a cubic spline or a cubic polynomial)
            
            sp2 = polynomial(list_x[i],list_y[i],2)                     # Quadratic polynomial
            
            shape1_2, shape2_2, start_2 = primitives(sp2,list_x[i])      # Shape of quadratic polynomial
            
            R2 = rvalue(sp2,list_x[i])                                   # R-value of quadratic polynomial
            
            p_string_2 = primitivesstring(shape1_2, shape2_2)
            nr_ABCD_2 = p_string_2.count("A") + p_string_2.count("B") + p_string_2.count("C") + p_string_2.count("D")
            
            if list_too_complex[i] == True:                              # Simplification is forced if curve is too complex
                
                list_sp[i] = sp2                                         # All lists are updated with the quadratic polynomial
                list_shape1[i], list_shape2[i] = shape1_2, shape2_2
                list_R[i] = R2
                list_too_complex[i] = primcomplex(primitivesstring(list_shape1[i],list_shape2[i]))
                
            elif R2/list_R[i] < 0.998 or nr_ABCD_2 > list_nr_ABCD[i]:      # the curve is kept if the difference in R-value is large or the quadratic polynomial has more primitives than the curve
                
                continue                                                           # Do nothing 
                
            elif shape1_2 != list_shape1[i] or shape2_2 != list_shape2[i]:         # The quadratic polynomial is used unless the curve and the quadratic polynomial have the same shape
                
                list_sp[i] = sp2                                                   # All lists are updated with the quadratic polynomial
                list_shape1[i], list_shape2[i] = shape1_2, shape2_2
                list_R[i] = R2
                list_too_complex[i] = primcomplex(primitivesstring(list_shape1[i],list_shape2[i]))
                
    if list_too_complex.count(True) == 0:                # If a trend is found, the algorithm is stopped
        
        if list_shape1.count(list_shape1[0]) == n_list and list_shape2.count(list_shape2[0]) == n_list and list_shape2[0] != [0]:
            return list_sp
        
    '''
    Step 3: The curve is now compared to a linear polynomial
    '''
        
    for i in range(0,n_list):
        
        sp = list_sp[i]       # Curve i
        
        m = len(list_x[i])    # Number of data points of the curve
        
        if m >= 3:            # Step 3 is only performed if the curve has at least 3 points (the curve is not already linear)
            
            sp1 = polynomial(list_x[i],list_y[i],1)                    # Linear polynomial
            
            shape1_1, shape2_1, start_1 = primitives(sp1,list_x[i])     # Shape of the linear polynomial
            
            R1 = rvalue(sp1,list_x[i])                                  # R-value of the linear polynomial
            
            if shape1_1 == list_shape1[i] and shape2_1 == list_shape2[i]:     # If the curve has the same shape is the linear polynomial, it is replaced by the linear polynomial for visual reasons
                
                list_sp[i] = sp1
                
            elif R1/list_R[i] < 0.998:    # If the difference in R-value is large, the original curve is kept
                
                continue                  # Do nothing
                
            else:                         # Otherwise, the curve is simplified to a linear polynomial
                
                list_sp[i] = sp1
                
    return list_sp


# In[102]:


def primitivesvisual_parametric(x,y,sp,extremes,p_string):
    
    '''
    Visualisation of the parametric datasets:
    Same method as for non-parametric datasets, but without vertical black lines and letters
    '''
    
    m = len(x)             # Number of data points
    
    n = len(p_string)      # Number of primitives
    
    xs = np.linspace(x[0],x[m-1],1000)       # Points to plot
    
    plt.plot(x, y, 'ro', ms=5)               # Data points are plotted in red
    
    for i in range(0,n):
        
        xs = np.linspace(extremes[i],extremes[i+1],1000)       # Points to plot for one primitive
        
        if p_string[i] == 'A':                             # If the primitive is A, the curve is plotted in green
            plt.plot(xs, sp(xs), 'g', lw=3)
        
        elif p_string[i] == 'B':                           # If the primitive is B, the curve is plotted in yellow
            plt.plot(xs, sp(xs), 'y', lw=3)
        
        elif p_string[i] == 'C':                           # If the primitive is C, the curve is plotted in blue
            plt.plot(xs, sp(xs), 'b', lw=3)
        
        elif p_string[i] == 'D':                           # If the primitive is D, the curve is plotted in red
            plt.plot(xs, sp(xs), 'r', lw=3)
        
        elif p_string[i] == 'E':                           # If the primitive is E, the curve is plotted in purple
            plt.plot(xs, sp(xs), color='#9400d3', lw=3)
        
        elif p_string[i] == 'F':                           # If the primitive is F, the curve is plotted in orange
            plt.plot(xs, sp(xs), color='#ff8c00', lw=3)
        
        else:                                              # If the primitive is G, the curve is plotted in grey
            plt.plot(xs, sp(xs), color='#808080', lw=3)


# In[103]:


def log_fit(x, y):
    
    '''
    Fits logarithm to a dataset defined by the lists x and y
    
    The fitted function is:
    
    y = A*ln(±(x-C)) + B
    
    The algorithm returns the values for the constants A, B and C, and the sign before (x-C)
    '''
    
    x_orig = x.copy()             # Original unchanged list of x-values
    y_orig = y.copy()             # Original unchanged list of y-values
    
    x = strict_increase(x)        # x now consists of a strict increase
    
    x = np.array(x)               # list x is turned into a Numpy-array
    y = np.array(y)               # list y is turned into a Numpy-array
    
    m = len(x)
    
    '''
    First, the function y = A*ln(x-C) + B is fitted to the data
    
    C is a constant which is slightly lower than x_min, the lowest x-value
    '''
    
    C_plus = estimate_C(x_orig, y_orig, 1)        # A value for the constant C is estimated
    
    ln_x = np.log(x - C_plus)                     # ln(x - C) is calculated
    
    sp = polynomial(ln_x,y,1)                     # y = A*ln(x - C) + B is now fitted to obtain the parameters A and B
    
    B_plus = sp([0])                              # y = A*0 + B equals the parameter B
    
    A_plus = sp([1]) - B_plus                     # A = (y - B)/ln(x - C); the value 1 is arbitrarily chosen for ln(x - C) to calculate this

    res_plus = log_residual(x, y, A_plus, B_plus, C_plus, 1)      # The residual is calculated for this logarithmic fit
    
    '''
    Then, the function y = A*ln(-(x-C)) + B is fitted to the data
    
    C is now a constant which is slightly higher than x_max, the highest x-value
    '''
    
    C_min = estimate_C(x_orig, y_orig, -1)        # A value for the constant C is estimated
    
    x = np.array(list(reversed(x)))               # The lists x and y are reversed, as x will be multiplied by -1 later
    y = np.array(list(reversed(y)))               # and x must increase for the UnivariateSpline algorithm
    
    ln_x = np.log(-(x - C_min))                   # ln(-(x - C)) is calculated
    
    sp = polynomial(ln_x,y,1)                     # y = A*ln(-(x - C)) + B is now fitted to obtain the parameters A and B
    
    B_min = sp([0])                               # y = A*0 + B equals the parameter B
    
    A_min = sp([1]) - B_min                       # A = (y - B)/ln(-(x - C)); the value 1 is arbitrarily chosen for ln(-(x - C)) to calculate this
        
    res_min = log_residual(x, y, A_min, B_min, C_min, -1)         # The residual is calculated for this logarithmic fit
    
    if res_min > res_plus:                        # If the logarithm y = A*ln(x-C) + B has a better fit than A*ln(-(x-C)) + B,
                                                  # The results of the former are returned
        return A_plus, B_plus, C_plus, 1
    
    else:                                         # Otherwise, the results for the latter are returned
        
        return A_min, B_min, C_min, -1


# In[104]:


def remove_repeats_at_extremes(x, y):
    
    '''
    If there are points with the x-values x_min or x_max, they are replaced by one points with the average of
    their corresponding y-values
    
    For example:
    
    x = [1,1,1,2,3,4,4]
    y = [1,2,3,4,4,5,7]
    
    New lists:
    
    x = [1,2,3,4]
    y = [2,4,4,6]
    '''
    
    m = len(x)
    
    '''
    First, all points with the value x_min are checked
    
    x_min corresponds to the first x-value x[0]
    '''
    
    y_tot = y[0]              # The y-values of all points with x-value x_min will be added to this sum
    
    i = 1                     # Counts how many points have the x-value x_min
    
    while x[1] == x[0] or abs((x[1] - x[0])/(x[-1] - x[0])) < 10**-6:       # If x[1] equals x[0], the point with index 1 (second point) must be removed
        
        y_tot += y[1]         # y[1] is added to the sum of the y-values
        
        del x[1]              # The point with index 1 is removed
        del y[1]              # In the next iteration, a new point will have index 1
        
        i += 1                # The counter i is updated
        
    y[0] = y_tot/i            # The average of the y-values of all points with x-value x_min
    
    '''
    Next, all points with the value x_max are checked
    
    x_max corresponds to the last x-value x[-1]
    '''
    
    y_tot = y[-1]             # The y-values of all points with x-value x_max will be added to this sum
    
    i = 1                     # Counts how many points have the x-value x_max
    
    while x[-2] == x[-1] or abs((x[-2] - x[-1])/(x[-1] - x[0])) < 10**-6:     # If x[-2] equals x[-1], the point with index -2 (second to last point) must be removed
        
        y_tot += y[-2]        # y[-2] is added to the sum of the y-values
        
        del x[-2]             # The point with index -2 is removed
        del y[-2]             # In the next iteration, a new point will have index -2
        
        i += 1                # The counter i is updated
        
    y[-1] = y_tot/i           # The average of the y-values of all points with x-value x_max
    
    return x, y


# In[105]:


def estimate_C(x, y, sign_x):
    
    '''
    An equation is used to estimate C in y = A*ln(±(x-C)) + B
    
    -------------------------------------------------------------------------
    
    For y = A*ln(x-C) + B, the equation is:
    
    C = ( x3 - x1 * exp(k) ) / ( 1 - exp(k) )   with   k = (y3 - y1) / (y2 - y3) * ln(α)
    
    In this equation, x1 and x3 represents the first and last x-value in the considered interval respectively.
    
    x2 is a point between these two, which can also be written as x2 = x1 + α*(x3 - x1) = α*x3 + (1-α)*x1
    
    yi represents the y-value corresponding to xi.
    
    This equation estimates the C-value necessary to to fit an exponential between three points and was derived
    under the assumption that x3 - C >> x1 - C, which means that x1 is relatively close to C
    
    -------------------------------------------------------------------------
    
    For y = A*ln(-(x-C)) + B, the equation is:
    
    C = ( x1 - x3 * exp(k) ) / ( 1 - exp(k) )   with   k = (y1 - y3) / (y2 - y1) * ln(1 - α)
    
    This equation was derived under the assumption that C - x1 >> C - x3, which means that x3 is relatively close to C
    
    --------------------------------------------------------------------------
    
    To constant C is estimated for all data points (x2,y2) with an x-value between x1 and x3
    
    The average of the different C-values is then taken as the final estimate.
    '''
    
    x, y = remove_repeats_at_extremes(x, y)  # Repeats are removed at the extremes of the considered intervals, see the function for more information
    
    m = len(x)           # Number of data points
    
    C_est = []           # The estimates for C are stored in this list
    
    for i in range(1, m-1):
        
        alpha = (x[i] - x[0])/(x[-1] - x[0])                    # Calculate α: x2 = x1 + α*(x3 - x1) ⇒ α = (x2 - x1)/(x3 - x1)
        
        if sign_x == 1:                                        
        
            k = np.log(alpha)*(y[-1] - y[0])/(y[i] - y[-1])     # The factor k is calculated as an intermediate calculation
        
            C = (x[-1] - x[0]*np.exp(k))/(1 - np.exp(k))        # C is estimated
            
        else:
            
            k = np.log(1 - alpha)*(y[0] - y[-1])/(y[i] - y[0])  # The factor k is calculated as an intermediate calculation
        
            C = (x[0] - x[-1]*np.exp(k))/(1 - np.exp(k))        # C is estimated
        
        if math.isnan(C) == False and C > x[0] - (x[-1] - x[0]) and C < x[-1] + (x[-1] - x[0]):
            
            C_est.append(C)                    # The estimation is added to C_est
        
    if sign_x == 1 and len(C_est) == 0:      
        
        return x[0] - 10**-8*(x[-1] - x[0])
    
    elif sign_x == -1 and len(C_est) == 0:
    
        return x[-1] + 10**-8*(x[-1] - x[0])
    
    else:
        
        C_est_avg = sum(C_est)/len(C_est)      # The average of the estimates is taken
    
    '''
    For y = A*ln(x-C) + B, C cannot be equal to or larger than the smallest x-value x1,
    as the logarithm cannot be calculated if x - C ≤ 0
    
    For this reason, C is estimated as x1 - 10^(-8)*(x3 - x1) if C ≥ x1
    
    ---------------------------------------------------------
    
    For y = A*ln(-(x-C)) + B, C cannot be equal to or smaller than the largest x-value x3,
    as the logarithm cannot be calculated if x - C ≥ 0
    
    For this reason, C is estimated as x3 + 10^(-8)*(x3 - x1) if C ≤ x3
    '''
    
    if sign_x == 1 and C_est_avg > x[0] - 10**-8*(x[-1] - x[0]):      
        
        return x[0] - 10**-8*(x[-1] - x[0])
    
    elif sign_x == -1 and C_est_avg < x[-1] + 10**-8*(x[-1] - x[0]):
    
        return x[-1] + 10**-8*(x[-1] - x[0])
    
    else:
    
        return C_est_avg


# In[106]:


def log_residual(x, y, A, B, C, sign_x):
    
    '''
    Calculates the residual of a logarithmic fit
    '''
    
    m = len(x)
    
    res_tot = 0
    
    for i in range(0, m):
        
        y_pred = A*np.log(sign_x*(x[i] - C)) + B
        
        res = (y[i] - y_pred)**2
        
        res_tot += res
        
    return res_tot


# In[107]:


def log_primitives(A, B, C, sign_x, x):
    
    '''
    Returns two lists storing the evolution of the signs of the first and second derivative of a logarithmic function, and a list with the extremes of the intervals
    
    For example:
    
    shape1 = [1, -1, 0, 1]
    shape2 = [1, 1, 0, 0]
    extremes = [0, 0.5, 1, 2.5, 4]
    
    Between x-values 0 and 0.5: first derivative is positive (1), second dertivative is positive (1)
    Between x-values 0.5 and 1: first derivative is negative (-1), second derivative is positive (1)
    Between x-values 1 and 2.5: first derivative is close to zero (0), second dertivative is close to zero (0)
    Between x-values 2.5 and 4: first derivative is positive (1), second derivative is close to zero (0)
    '''
    
    xs = np.linspace(x[0],x[-1],100)          # Points where derivatives are calculated
    ys = A*np.log(sign_x*(xs - C)) + B
    
    ys_der_1 = A/(xs - C)
    ys_der_2 = -A/np.square(xs - C)
    
    x_m = x[-1] - x[0]                        # x_max - x_min
    y_m = max(ys) - min(ys)                   # y_max - y_min
    
    shape1, shape2 = [], []                         # Lists containing the signs
    extremes = []                                   # List containing the extremes of the intervals
    
    if np.abs(ys_der_1[0])*x_m/y_m < 0.5:                # If the normalised first derivative is lower than 0.5, it's considered 0
        sign_der1 = 0
            
    else:
        sign_der1 = int(np.sign(ys_der_1[0]))
            
    if np.abs(ys_der_2[0])*np.square(x_m)/y_m < 0.5 or sign_der1 == 0:       # If the normalised derivative is lower than 1, it's considered 0
        sign_der2 = 0                                                   # If the first derivative equals 0, both derivatives are set to 0
        
    else:
        sign_der2 = int(np.sign(ys_der_2[0]))
        
    sign_change = True                                              # If the sign changes, this becomes True. It is used to make sure that two consecutive point have the same primitive
        
    for i in range(1,100):                                          # If one of the signs changes, the lists are updated with the new signs
            
        if np.abs(ys_der_1[i])*x_m/y_m < 0.5:                          # The previous equations for the derivatives are repeated for point i
            sign_der1_new = 0
            
        else:
            sign_der1_new = int(np.sign(ys_der_1[i]))
            
        if np.abs(ys_der_2[i])*np.square(x_m)/y_m < 0.5 or sign_der1_new == 0:
            sign_der2_new = 0
            
        else:
            sign_der2_new = int(np.sign(ys_der_2[i]))
            
        if sign_der1_new != sign_der1 or sign_der2_new != sign_der2:    # If sign changes, sign_change becomes True, or when i == 1 as this is the first point in the loop
            sign_change = True
            
        elif sign_change == True:                                   # Two consecutive points have the same derivative
                                                                        # (this elif-statement can only be reached if the previous once was False,
            if shape1 == []:                                        # this means that there was no sign change this time)
                extremes.append(x[0])
                shape1.append(sign_der1_new)                        # If the shape-list is still empty, x[0] (the first point) is added to the list 'extremes'
                shape2.append(sign_der2_new)
                
            elif shape1[-1] != sign_der1_new or shape2[-1] != sign_der2_new:
                extremes.append(xs[i-1])                 
                shape1.append(sign_der1_new)                        # If the shape-list is not empty, point i-1 is added to the list 'extremes'
                shape2.append(sign_der2_new)                        # Neither of the if-statement are reached if signs are the same as for the signs in the shape-list
                                                                        # as there is no need to have the same primitive twice consecutively
                
            sign_change = False         # 'sign_change' is set to False until a new sign change occurs
            
        sign_der1, sign_der2 = sign_der1_new, sign_der2_new         # Stored as der1 and der2 for the following iteration
        
    extremes.append(x[-1])               # Final x-value is also added to the list of extremes
    
    return shape1, shape2, extremes


# In[108]:


def log_primitivesvisual(x, y, A, B, C, sign_x, extremes, p_string):
    
    '''
    Visualises the features for a logarithmic curve
    '''
    
    x = np.array(x)               # list x is turned into a Numpy-array
    y = np.array(y)               # list y is turned into a Numpy-array
    
    m = len(x)             # Number of data points
    n = len(extremes)      # Number of extremes
    
    xs = np.linspace(x[0],x[-1],1000)       # Points to plot
    ys = A*np.log(sign_x*(xs - C)) + B
    
    y_min, y_max = min([min(y),min(ys)]), max([max(y),max(ys)])          # Smallest and highest y-values reached by the spline and the data points 
    delta_y = y_max - y_min                                              # The range of y-values
    
    plt.plot(x, y, 'ro', ms=5)               # Data points are plotted in red
    
    '''
    The data points have been plotted; now, the spline and its features need to plotted
    '''
    
    for i in range(0,n-1):
        
        xs = np.linspace(extremes[i],extremes[i+1],1000)       # Points to plot for one primitive
        ys = A*np.log(sign_x*(xs - C)) + B
        
        if p_string[i] == 'A':                             # If the primitive is A, the curve is plotted in green
            plt.plot(xs, ys, 'g', lw=3)
        
        elif p_string[i] == 'B':                           # If the primitive is B, the curve is plotted in yellow
            plt.plot(xs, ys, 'y', lw=3)
        
        elif p_string[i] == 'C':                           # If the primitive is C, the curve is plotted in blue
            plt.plot(xs, ys, 'b', lw=3)
        
        elif p_string[i] == 'D':                           # If the primitive is D, the curve is plotted in red
            plt.plot(xs, ys, 'r', lw=3)
        
        elif p_string[i] == 'E':                           # If the primitive is E, the curve is plotted in purple
            plt.plot(xs, ys, color='#9400d3', lw=3)
        
        elif p_string[i] == 'F':                           # If the primitive is F, the curve is plotted in orange
            plt.plot(xs, ys, color='#ff8c00', lw=3)
        
        else:                                              # If the primitive is G, the curve is plotted in grey
            plt.plot(xs, ys, color='#808080', lw=3)
        
        xt = (extremes[i]+extremes[i+1])/2 - (extremes[n-1]-extremes[0])/40    # x-position where letter will be placed on graph
        
        if delta_y != 0:
            
            yt = y_min - 0.2*delta_y                       # y-position where letter will be placed on graph
 
        elif y_min != 0:

            yt = 0.95*y_min
            
        else:
            
            yt = -5
        
        plt.text(xt,yt,p_string[i],fontsize=15)            # Letter representing primitive is printed on graph
    
    '''
    The black vertical lines are now plotted between the primitives
    '''
    
    if delta_y != 0:
        
        yp = np.linspace(y_min - 0.2*delta_y, y_max + 0.2*delta_y, 1000)         # y points to plot in the vertical line
    
    elif y_min != 0:
        
        yp = np.linspace(0.95*y_min, 1.05*y_max, 1000)
        
    else:
        
        yp = np.linspace(-5, 5, 1000)
    
    for i in range(0,n):
        
        xp = [extremes[i]]*1000                # List with 1000 times the x-position of the line
        
        plt.plot(xp, yp, 'k', lw=2)          # Vertical line is plotted


# In[109]:


def log_primitivesvisual_parametric(x, y, A, B, C, sign_x, extremes, p_string):
    
    '''
    Visualises the features for a logarithmic curve without the letters and the vertical black lines
    '''
    
    x = np.array(x)               # list x is turned into a Numpy-array
    y = np.array(y)               # list y is turned into a Numpy-array
    
    m = len(x)             # Number of data points
    n = len(extremes)      # Number of extremes
    
    xs = np.linspace(x[0],x[-1],1000)       # Points to plot
    ys = A*np.log(sign_x*(xs - C)) + B
    
    y_min, y_max = min([min(y),min(ys)]), max([max(y),max(ys)])          # Smallest and highest y-values reached by the spline and the data points 
    delta_y = y_max - y_min                                              # The range of y-values
    
    plt.plot(x, y, 'ro', ms=5)               # Data points are plotted in red
    
    '''
    The data points have been plotted; now, the spline and its features need to plotted
    '''
    
    for i in range(0,n-1):
        
        xs = np.linspace(extremes[i],extremes[i+1],1000)       # Points to plot for one primitive
        ys = A*np.log(sign_x*(xs - C)) + B
        
        if p_string[i] == 'A':                             # If the primitive is A, the curve is plotted in green
            plt.plot(xs, ys, 'g', lw=3)
        
        elif p_string[i] == 'B':                           # If the primitive is B, the curve is plotted in yellow
            plt.plot(xs, ys, 'y', lw=3)
        
        elif p_string[i] == 'C':                           # If the primitive is C, the curve is plotted in blue
            plt.plot(xs, ys, 'b', lw=3)
        
        elif p_string[i] == 'D':                           # If the primitive is D, the curve is plotted in red
            plt.plot(xs, ys, 'r', lw=3)
        
        elif p_string[i] == 'E':                           # If the primitive is E, the curve is plotted in purple
            plt.plot(xs, ys, color='#9400d3', lw=3)
        
        elif p_string[i] == 'F':                           # If the primitive is F, the curve is plotted in orange
            plt.plot(xs, ys, color='#ff8c00', lw=3)
        
        else:                                              # If the primitive is G, the curve is plotted in grey
            plt.plot(xs, ys, color='#808080', lw=3)


# In[110]:


def log_rvalue(x, y, A, B, C, sign_x):
   
    '''
    Calculates multiple correlation coefficient R of a spline fit
    '''
    
    x = np.array(x)               # list x is turned into a Numpy-array
    y = np.array(y)               # list y is turned into a Numpy-array
    
    ssy = np.sum(np.square(A*np.log(sign_x*(x - C)) + B))      # Sum of squares of prediction
    sse = log_residual(x, y, A, B, C, sign_x)                  # Sum of squares of error
    R = np.sqrt(ssy/(sse+ssy))                                 # Multiple correlation coefficient
    
    return R


# In[111]:


def log_simplify(x, y, sp, A, B, C, sign_x):
    
    R_sp = rvalue(sp, x)
    R_log = log_rvalue(x, y, A, B, C, sign_x)
    
    k = len(sp.get_coeffs()) - 1                      # Degree of the spline
    
    if R_log/R_sp < 0.998 or k == 1:
        
        return False
    
    else:
        
        return True