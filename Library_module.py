# -*- coding: utf-8 -*-
"""
Created on Mon May 11 14:15:47 2020

@author: Vasco Saltão
"""

import matplotlib.pyplot as plt
import FeatureExtraction as ft
import math
import numpy as np
import pandas as pd
import decimal

Tags_Lib = {'a': ('All species adsorb molecularly', 'Adsorption of A controlling'),
            'b': ('All species adsorb molecularly', 'Adsorption of B controlling'),
            'c': ('All species adsorb molecularly', 'Desorption of R controlling'),
            'd': ('All species adsorb molecularly', 'Surface reaction controlling'),
            'e': ('A adsorbs dissociatively', 'Adsorption of A controlling'),
            'f': ('A adsorbs dissociatively', 'Adsorption of B controlling'),
            'g': ('A adsorbs dissociatively', 'Desorption of R controlling'),
            'h': ('A adsorbs dissociatively', 'Surface reaction controlling'),
            'i': ('B does not adsorb', 'Adsorption of A controlling'),
            'j': ('B does not adsorb', 'Desorption of R controlling'),
            'k': ('B does not adsorb', 'Surface reaction controlling'),
            'l': ('A adsorbs dissociatively', 'B does not adsorb', 'Adsorption of A controlling'),
            'm': ('A adsorbs dissociatively', 'B does not adsorb', 'Desorption of R controlling'),
            'n': ('A adsorbs dissociatively', 'B does not adsorb', 'Surface reaction controlling'),
            'o': ('A does not adsorb', 'Impact of A controlling'),
            'p': ('A does not adsorb', 'Desorption of R controlling'),
            'q': ('A does not adsorb', 'Adsorption of B controlling'),
            'r': ('Uncatalyzed reaction')}

IDs_full = list(Tags_Lib.keys())

IDs_swap = ['e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q'] #non-symmetrical mechanisms that are run with swapped reactants

"""
Lists of possible values for the equilibrium constants
"""
Keq_adsorbed_full = [0.1, 0.3, 0.5, 0.7, 1, 10, 100]  # Possible values for adsorption equilibrium constants of species that meaningfully adsorb

Keq_adsorbed_redux = [0.1, 0.5, 1, 10, 100] # Redux - Used only for reducing runtime of tool

Keq_all_full = [1e-7, 0.01, 0.1, 0.3, 0.5, 0.7, 1, 10, 100]  # Possible values for equilibrium constants

Keq_all_redux = [1e-7, 0.01, 0.1, 0.5, 1, 10, 100] # Redux - Used only for reducing runtime of tool

Keq_impact_full = [0.1, 0.3, 0.5, 0.7, 1, 10, 100]  # Possible values for impact equilibrium constants

Keq_impact_redux = [0.1, 0.5, 1, 10, 100] # Redux - Used only for reducing runtime of tool


def Theoretical_Features(Stoich, yA, yB, pressure, Bad_IDs, Redux, Swap_run = False):
    """
    Creates an archive with the features of all the generated curves for all 
    tested initial rate equations

    Parameters
    ----------
    Stoich : List
        Stoichiometric coefficients of species A, B, R and S, respectively.
    yA : Float
        Molar fraction of reactant A.
    yB : Float
        Molar fraction of reactant B.
    pressure : List
        Pressure values used to generate the theoretical curves.
    Bad_IDs : List
        Contains the IDs of the mechanisms that the user knows 'a priori' that 
        cannot be a valid model for the experimental data.
    Redux : Boolean
        True if the user chose to use the reduced lists of possible values for 
        the equilibrium constants.
    Swap_run : Boolean, optional
        True if the current run is the second one with the reactants swapped. The default is False.

    Returns
    -------
    all_Theor_Features : Dictionary
        Contains all of the possible features for every tested initial rate 
        equation.

    """
    
    all_Theor_Features = {}
    
    if Swap_run == True:
        
        ID_list = IDs_swap
    
    else:
        
        ID_list = IDs_full
    
    for ID in ID_list:
        
        if ID in Bad_IDs:
            
            continue
        
        Tags = Tags_Lib[ID]
        
        if (Stoich[1] == 0) and (('Adsorption of B controlling' in Tags) or ('B does not adsorb' in Tags) or ('A does not adsorb' in Tags)): # These mechanisms are only valid for more than 1 reactant
            
            continue
        
        Theor_Features = Possible_Features(Tags, Stoich, yA, yB, pressure, Redux)
        
        if Theor_Features != 'Unrealistic':
            
            all_Theor_Features[ID] = Theor_Features
        
    return all_Theor_Features


def Lib_Check(all_Theor_Features, Exp_Features, Tolerance, Exp_rate_norm):
    """
    Searches for initial rate equations that can produce curves with features 
    similar to the ones from the experimental data.

    Parameters
    ----------
    all_Theor_Features : Dictionary
        Contains all of the possible features for every tested initial rate 
        equation.
    Exp_Features : List
        Contains the curve features extracted from the experimental data.
    Tolerance : Float
        Maximum relative error allowed when comparing theoretical and 
        experimental extremes. Chosen by the user.
    Exp_rate_norm : List
        Normalized initial rate values from the experimental data.

    Returns
    -------
    Tables : Dictionary
        Contains for each initial rate equation a table with the parameters 
        and matching criteria of each theoretical curve that was considered as 
        a possibility by the tool.

    """
    Exp_Primitives = tuple(Exp_Features[0])
    
    Exp_Extremes = Exp_Features[1]
    
    Tables={}
    
    for ID in all_Theor_Features:
        
        Theor_Features = all_Theor_Features[ID]
            
        if Exp_Primitives in Theor_Features: # Checking if this rate equation produces any curves with the same primitives
            
            Current_Table = Theor_Features[Exp_Primitives]
            
            New_Table = pd.DataFrame(columns=['KA', 'KB', 'KR', 'KS', 'Kr', 'Extremes', 'y_values', 'MSE', 'SRE'])
            
            i=0
            while i<len(Current_Table):
                
                Theor_Extremes = Current_Table['Extremes'][i]
        
                Extremes_in_range = True    # Start assuming that all extremes are in range
                
                SRE = 0
                
                j=1
                while j<len(Theor_Extremes)-1: # Because the first and last extremes are always the start and end of the dataset
                    
                    if not (Exp_Extremes[j]*(1-Tolerance) <= Theor_Extremes[j] <= Exp_Extremes[j]*(1+Tolerance)):
                        
                        Extremes_in_range = False
                        break
                    
                    else:
                        
                        SRE += abs((Theor_Extremes[j] - Exp_Extremes[j])/Exp_Extremes[j])
                    
                    j+=1
                    
                if Extremes_in_range == True: # Means that all extremes are in range
                    
                    rate_norm = Current_Table['rate_norm'][i]
                    
                    Sum_SE = 0
                    n=0
                    while n<len(rate_norm):
                        
                        Sum_SE += (rate_norm[n] - Exp_rate_norm[n])**2 # Adding the squared errors of all datapoints
                        
                        n+=1
                    
                    MSE = Sum_SE/len(rate_norm) # Calculating the normalized mean squared error (MSE_norm)
                    
                    index = len(New_Table)
                    
                    row = list(Current_Table.iloc[i])
                    
                    new_row = row[0:5] + row[6:]
                    
                    new_row += [MSE, SRE]
                    
                    New_Table.loc[index] = new_row
                
                i+=1
                
            if len(New_Table) > 0: # Means that there is at least 1 theoretical curve for which the exp extremes are in range
                
                Tables[ID] = New_Table.sort_values(by = ['SRE', 'MSE'], ignore_index = True)
    
    return Tables


def Build_Initial_Rate(Stoich, Tags):
    """
    Generates a text version of the initial rate equation

    Parameters
    ----------
    Stoich : List
        Stoichiometric coefficients of species A, B, R and S, respectively.
    Tags : Tuple
        Text strings that describe a mechanism.

    Returns
    -------
    Eq_str : String
        Initial rate equation in printable format.

    """
    a=Stoich[0]
    b=Stoich[1]
    r=Stoich[2]
    s=Stoich[3]
    
    power = 1
        
    if 'A adsorbs dissociatively' in Tags:
            
        A_term_str = ' + (KA*yA)^(1/2)'
        
        P_A_str = '*P^(1/2)'
        
        power_a = 2*a
    
    else:
        
        A_term_str = ' + KA*yA'
        
        P_A_str = '*P'
        
        power_a = a
    
    if b == 0 or 'B does not adsorb' in Tags:
        
        B_term_str = ''
        
        P_B_str = ''
        
        power_b = 0
        
    else:
        
        B_term_str = " + KB*yB"
        
        P_B_str = '*P'
        
        power_b = b
        
    if 'Adsorption of A controlling' in Tags:
            
        numerator_str = "k'*yA*P"
        
        if B_term_str == "":
                    
            denominator_str = ""
            
        else:
                    
            denominator_str = "1" + B_term_str + P_B_str
            
        if 'A adsorbs dissociatively' in Tags:
            
            power = 2
        
    elif 'Adsorption of B controlling' in Tags:
            
        numerator_str = "k'*yB*P"
        
        if 'A does not adsorb' in Tags:
            
            denominator_str = ""
        
        else:
        
            denominator_str = "1" + A_term_str + P_A_str
    
    elif 'Desorption of R controlling' in Tags:
        
        if s == 0:
            
            if b == 0:
                
                if r == 1:
                    
                    numerator_str = "k'*Kg*yA^a*P^a"
                    
                else:
                    
                    numerator_str = "k'*(Kg*yA^a)^(1/r)*P^(a/r)"
                
            else:
                
                if r == 1:
                    
                    numerator_str = "k'*Kg*yA^a*yB^b*P^{}".format(a+b)
                
                else:
                    
                    numerator_str = "k'*(Kg*yA^a*yB^b)^(1/r)*P^({}/r)".format(a+b)
            
            if 'A does not adsorb' in Tags:
                
                A_term_str = " + KA*KB*yB*yA^(a/b)"
                
                P_A_str = "*P^(sum/b)"
               
                if a == b:
                   
                    A_term_str =  A_term_str.replace("^(a/b)", "")
                    
                    P_A_str = P_A_str.replace('^(sum/b)', '^2')
                
                elif b == 1:
                    
                    A_term_str =  A_term_str.replace("^(a/b)", "^a") 
                    
                    P_A_str = P_A_str.replace('^(sum/b)', '^' + str(1+a))
                
                else:
                    
                    P_A_str = P_A_str.replace('sum', str(b+a))
                    
            desorption_term_str = numerator_str.replace("k'", "KR")
            
            if 'All species adsorb molecularly' in Tags:
                
                denominator_str = "1 + " + desorption_term_str + ' + (' + A_term_str.replace(' + ', '') + B_term_str + ')*P'
                
            else:
            
                denominator_str = "1 + " + desorption_term_str + A_term_str + P_A_str + B_term_str + P_B_str
            
        else:
            
            numerator_str = "k'"
            
            denominator_str = "KR"
        
    elif 'Surface reaction controlling' in Tags:
            
        if b == 0:
            
            numerator_str = "k'*(KA*yA)^a*P^a"
        
        else:
            
            if 'B does not adsorb' in Tags:
                
                numerator_str = "k'*(KA*yA)^a*yB^b*P^{}".format(a+b)
                
            else:
            
                numerator_str = "k'*(KA*yA)^a*(KB*yB)^b*P^{}".format(a+b)
            
        if 'All species adsorb molecularly' in Tags:
                
            denominator_str = "1" + ' + (' + A_term_str.replace(' + ', '') + B_term_str + ')*P'
                
        else:
                
            denominator_str = "1" + A_term_str + P_A_str + B_term_str + P_B_str
            
        if power_a + power_b >= r + s:
            
            power = power_a + power_b
            
        else:
            
            power = r + s
    
    elif 'Impact of A controlling' in Tags:
        
        numerator_str = "k'*KB*yB*yA^(a/b)*P^(sum/b)"
        
        if a == b:
                   
            numerator_str = numerator_str.replace('^(sum/b)', '^2')
            
        elif b == 1:
            
            numerator_str = numerator_str.replace('^(a/b)', '^a')
            
            numerator_str = numerator_str.replace('^(sum/b)', '^' + str(1+a))
            
        else:
            
            numerator_str = numerator_str.replace('sum', str(b+a))
                
        denominator_str = "1" + B_term_str + '*P'
        
    elif 'Uncatalyzed reaction' in Tags:
            
        if b == 0:
            
            numerator_str = "k*yA^a*P^a"
            
        else:
            
            numerator_str = "k*yA^a*yB^b*P^{}".format(a+b)
                
        denominator_str = ""
        
    numerator_str = numerator_str.replace('a', str(a))
    numerator_str = numerator_str.replace('b', str(b))
    numerator_str = numerator_str.replace('r', str(r))
    
    numerator_str = numerator_str.replace('.0', "")
    
    numerator_str = numerator_str.replace('^1', "")
    numerator_str = numerator_str.replace('^(1/1)', "")
    numerator_str = numerator_str.replace('^(2/2)', "")
    numerator_str = numerator_str.replace('^(4/2)', "^2")
    
    if denominator_str == "":
            
        Eq_str = 'r = ' + numerator_str
    
    else:
        
        if power != 1:
            
            denominator_str = "(" + denominator_str + ")^power"
        
        denominator_str = denominator_str.replace('a', str(a))
        denominator_str = denominator_str.replace('b', str(b))
        denominator_str = denominator_str.replace('power', str(power))
        denominator_str = denominator_str.replace('r', str(r))
        
        denominator_str = denominator_str.replace('.0', "")
        
        denominator_str = denominator_str.replace('^1', "")
        denominator_str = denominator_str.replace('^(2/2)', "")
        denominator_str = denominator_str.replace('^(4/2)', "^2")
            
        fraction_dash = 'r = ---------------------------------'
        
        Eq_str = "      " + numerator_str + "\n" + fraction_dash + "\n" + "      " + denominator_str
        
    return Eq_str


def Exp_and_Theoretical(Stoich, yA, yB, pressure, Exp_rate, Tables, R_names, P_names, Swap = False):
    """
    Selects the best theoretical curve from each initial rate equation, 
    creates the text strigs to appear in the output and prepares the 
    information to create the graphical output.

    Parameters
    ----------
    Stoich : List
        Stoichiometric coefficients of species A, B, R and S, respectively.
    yA : Float
        Molar fraction of reactant A.
    yB : Float
        Molar fraction of reactant B.
    pressure : List
        Pressure values used to generate the theoretical curves.
    Exp_rate : List
        Initial rate values present in the experimental data.
    Tables : Dictionary
        Contains for each initial rate equation a table with the parameters 
        and matching criteria of each theoretical curve that was considered as 
        a possibility by the tool.
    R_names : List
        Contains the names of the reactants.
    P_names : List
        Contains the names of the products.
    Swap : Boolean, optional
        True if the current run is the second one with the reactants swapped. The default is False.

    Returns
    -------
    Final_table : Dataframe
        Table with the best theoretical curve from each initial rate equation 
        that produced curves that follow the experimantal data trends.
    Text_list : Dictionary
        Contains the text output for each proposed initial rate equation.
    Tags_list : Dictionary
        Contains the altered descriptive tags of every proposed initial rate 
        equation.
    Theor_press : List
        Pressure values used to plot the theoretical curves on the output 
        graphs.

    """
    a = Stoich[0]
    b = Stoich[1]
    r = Stoich[2]
    s = Stoich[3]
    
    Theor_press = np.linspace(0, max(pressure), 100)
    
    Final_table = pd.DataFrame(columns=['ID', 'Theor_rate', 'k', 'KA', 'KB', 'KR', 'KS', 'Kr', 'MSE', 'SRE'])
    
    for ID in Tables:
        
        Tags = Tags_Lib[ID]
        
        [KA, KB, KR, KS, Kr, Extremes, y_values, MSE_norm, SRE] = list(Tables[ID].iloc[0])
        
        Eq_Consts = [KA, KB, KR, KS, Kr]
        
        k, MSE = Get_k(Exp_rate, Tags, Stoich, yA, yB, Eq_Consts, pressure)
        
        Theor_rate = Generate_Curve(Tags, Stoich, yA, yB, k, Eq_Consts, Theor_press)
        
        index = len(Final_table)
        
        Final_table.loc[index] = [ID, Theor_rate, k, KA, KB, KR, KS, Kr, MSE, SRE]
    
    Text_list = {}
    
    Tags_list = {}
    
    New_IDs = []
    
    i=0
    while i<len(Final_table):
        
        [ID, Theor_rate, k, KA, KB, KR, KS, Kr, MSE, SRE] = list(Final_table.iloc[i])
        
        Tags = list(Tags_Lib[ID])
        
        Output_tags = []
        
        for Tag in Tags:
            
            New_Tag = Tag
            
            New_Tag = New_Tag.replace('A ', R_names[0] + ' ')
            
            if b > 0:
                
                New_Tag = New_Tag.replace('B ', R_names[1] + ' ')
            
            if s == 0:
                
                New_Tag = New_Tag.replace('R ', P_names[0] + ' ')
            
            else:
                
                New_Tag = New_Tag.replace('R ', 'a product ')
            
            Output_tags.append(New_Tag)
        
        
        if Swap == True:
            
            ID_new = ID + '_2'
        
        else:
            
            ID_new = ID + '_1'
        
        New_IDs.append(ID_new)
        
        Tags_list[ID_new] = Output_tags
        
        
        Eq_str = Build_Initial_Rate(Stoich, Tags)
        
        text = "\n" + "Rate equation:"
        
        text += "\n\n" + Eq_str
        
        text += "\n\n" + "A = {}".format(R_names[0])
        
        if b > 0:
            
            text += "\n" + "B = {}".format(R_names[1])
        
        text += "\n\n" + "Specific constants of this simulated dataset (Units coincide with pressure and rate units):"
        
        text += "\n\n" + "k' = {:.3e}".format(k)
        
        if "KA" in Eq_str:
            
            text += "\n" + "KA = {:.3e}".format(KA)
        
        if "KB" in Eq_str:
            
            text += "\n" + "KB = {:.3e}".format(KB)
        
        if "KR" in Eq_str:
            
            text += "\n" + "KR = {:.3e}".format(KR)
        
        if "KS" in Eq_str:
            
            text += "\n" + "KS = {:.3e}".format(KS)
        
        if "Kg" in Eq_str:
            
            if 'B does not adsorb' in Tags: # j & m
                
                Kg = Kr*KA**a/(KR**r)
            
            elif 'A does not adsorb' in Tags: # p
                
                Kg = Kr*(KA*KB)**b/(KR**r)
                
            else: # c & g
                
                Kg = Kr*KA**a*KB**b/(KR**r) 
            
            text += "\n" + "Kg = {:.3e}".format(Kg)
        
        text += "\n\n" + "Ranking criteria:"
        
        text += "\n\n" + "SRE = {:.3e}".format(SRE)
        
        text += "\n" + "MSE = {:.3e}".format(MSE) + "\n\n"
        
        
        Text_list[ID_new] = text
        
        i+=1
    
    Final_table['ID'] = New_IDs
        
    return Final_table, Text_list, Tags_list, Theor_press


def Create_output(Final_table, Text_list, Tags_list, Theor_press, pressure, Exp_rate, P_units, r0_units):
    """
    Creates the graphical output and prints the full output in the console

    Parameters
    ----------
    Final_table : Dataframe
        Table with the best theoretical curve from each initial rate equation 
        that produced curves that follow the experimantal data trends.
    Text_list : Dictionary
        Contains the text output for each proposed initial rate equation.
    Tags_list : Dictionary
        Contains the altered descriptive tags of every proposed initial rate 
        equation.
    Theor_press : List
        Pressure values used to plot the theoretical curves on the output 
        graphs.
    pressure : List
        Pressure values used to generate the theoretical curves.
    Exp_rate : List
        Initial rate values present in the experimental data.
    P_units : String
        Pressure units.
    r0_units : String
        Initial rate units.

    Returns
    -------
    None.

    """
    i=0
    while i < len(Final_table):
        
        ID = Final_table['ID'][i]
        
        Tags = Tags_list[ID]
        
        text = Text_list[ID]
        
        Theor_rate = Final_table['Theor_rate'][i]
        
        """
        Creating the plot
        """
        plt.figure(i+1)
        
        Title = 'Mechanism ' + str(ID) + ": " + str(Tags)
        
        plt.title(Title)
        
        plt.xlabel('P (' + P_units + ')')
        plt.ylabel('r0 (' + r0_units + ')')
        
        plt.plot(pressure, Exp_rate, 'ro', label='Experimental') # plotting the experimental datapoints
        
        plt.plot(Theor_press, Theor_rate, color='tab:blue', label='Simulated')
        
        plt.legend()
        
        plt.show(block=False)
        
        print(text)
        
        i+=1
    
    if len(Final_table) > 0:
        
        Output_table = Final_table.drop('Theor_rate', 1)
    
        print(Output_table)
    
    else:
        
        print('\n', 'No mechanisms were a good enough match')


def Get_k(Exp_rate, Tags, Stoich, yA, yB, Eq_Consts, pressure):
    """
    Estimates an initial guess for the parameter k' of a theoretical curve

    Parameters
    ----------
    Exp_rate : List
        Initial rate values present in the experimental data.
    Tags : Tuple
        Text strings that describe a mechanism.
    Stoich : List
        Stoichiometric coefficients of species A, B, R and S, respectively.
    yA : Float
        Molar fraction of reactant A.
    yB : Float
        Molar fraction of reactant B.
    Eq_Consts : List
        Contains the equilibrium constants of the theoretical curve.
    pressure : List
        Pressure values used to generate the theoretical curves.

    Returns
    -------
    best_k : Float
        k' value that has the smallest MSE from all the tested ones.
    best_MSE : Float
        Smallest MSE from all the tested k' values.

    """
    decimal.getcontext().prec = 100
    
    a = decimal.Decimal(Stoich[0])
    b = decimal.Decimal(Stoich[1])
    r = decimal.Decimal(Stoich[2])
    s = decimal.Decimal(Stoich[3])
    
    KA = decimal.Decimal(Eq_Consts[0])
    KB = decimal.Decimal(Eq_Consts[1])
    KR = decimal.Decimal(Eq_Consts[2])
    KS = decimal.Decimal(Eq_Consts[3])
    Kr = decimal.Decimal(Eq_Consts[4])
    
    yA = decimal.Decimal(yA)
    yB = decimal.Decimal(yB)
    
    k_list = []
    
    if 'Adsorption of A controlling' in Tags:    # a, e, i & l
        
        if ('B does not adsorb' in Tags) or (b==0):    # i & l (a & e for b=0)
            
            i=0
            while i<len(Exp_rate):
                
                p = decimal.Decimal(pressure[i])
                r = decimal.Decimal(Exp_rate[i])
                
                k = decimal.Decimal(r/(yA*p))
                
                k_list.append(k)
                
                i+=1
        
        else:    # a & e
            
            if 'A adsorbs dissociatively' in Tags:    # e
            
                i=0
                while i<len(Exp_rate):
                    
                    p = decimal.Decimal(pressure[i])
                    r = decimal.Decimal(Exp_rate[i])
                    
                    k = decimal.Decimal(r*(1 + KB*yB*p)**2/(yA*p))
                    
                    k_list.append(k)
                    
                    i+=1
        
            else:    # a
            
                i=0
                while i<len(Exp_rate):
                    
                    p = decimal.Decimal(pressure[i])
                    r = decimal.Decimal(Exp_rate[i])
                    
                    k = decimal.Decimal(r*(1 + KB*yB*p)/(yA*p))
                    
                    k_list.append(k)
                    
                    i+=1
    
    elif 'Adsorption of B controlling' in Tags:    # b, f & q
        
        if 'A does not adsorb' in Tags:    # q
            
            i=0
            while i<len(Exp_rate):
                
                p = decimal.Decimal(pressure[i])
                r = decimal.Decimal(Exp_rate[i])
                
                k = decimal.Decimal(r/(yB*p))
                
                k_list.append(k)
                
                i+=1
        
        elif 'A adsorbs dissociatively' in Tags:    # f
            
            i=0
            while i<len(Exp_rate):
                
                p = decimal.Decimal(pressure[i])
                r = decimal.Decimal(Exp_rate[i])
                
                k = decimal.Decimal(r*(1 + decimal.Decimal(math.sqrt(KA*yA*p)))/(yB*p))
                
                k_list.append(k)
                
                i+=1
        
        else:    # b
            
            i=0
            while i<len(Exp_rate):
                
                p = decimal.Decimal(pressure[i])
                r = decimal.Decimal(Exp_rate[i])
                
                k = decimal.Decimal(r*(1 + KA*yA*p)/(yB*p))
                
                k_list.append(k)
                
                i+=1
            
    elif 'Desorption of R controlling' in Tags: # c, g, j, m & p
        
        if s!=0: # If there is more than 1 product, the initial rate equation is the same for all mechanisms where desorption of 1 product is the RDS
            
            i=0
            while i<len(Exp_rate):
                
                p = decimal.Decimal(pressure[i])
                r = decimal.Decimal(Exp_rate[i])
                
                k = decimal.Decimal(r*KR)
                
                k_list.append(k)
                
                i+=1
            
        else:
            
            if 'B does not adsorb' in Tags: # j & m
                
                Kg = decimal.Decimal(Kr*KA**a/(KR**r))
                
                num = decimal.Decimal((Kg*yA**a*yB**b)**(1/r))
                
                if 'A adsorbs dissociatively' in Tags: # m
                    
                    i=0
                    while i<len(Exp_rate):
                        
                        p = decimal.Decimal(pressure[i])
                        r = decimal.Decimal(Exp_rate[i])
                        
                        k = decimal.Decimal(r*(1 + decimal.Decimal(math.sqrt(KA*yA*p)) + KR*num*p**((a+b)/r))/(num*p**((a+b)/r)))
                        
                        k_list.append(k)
                        
                        i+=1
                
                else: # j
                    
                    i=0
                    while i<len(Exp_rate):
                        
                        p = decimal.Decimal(pressure[i])
                        r = decimal.Decimal(Exp_rate[i])
                        
                        k = decimal.Decimal(r*(1 + KA*yA*p + KR*num*p**((a+b)/r))/(num*p**((a+b)/r)))
                        
                        k_list.append(k)
                        
                        i+=1
            
            elif 'A does not adsorb' in Tags: # p
                
                Kg = decimal.Decimal(Kr*(KA*KB)**b/(KR**r))
                
                num = decimal.Decimal((Kg*yA**a*yB**b)**(1/r))
                
                i=0
                while i<len(Exp_rate):
                    
                    p = decimal.Decimal(pressure[i])
                    r = decimal.Decimal(Exp_rate[i])
                    
                    k = decimal.Decimal(r*(1 + KA*KB*yA**(a/b)*yB*p**(1+a/b) + KR*num*p**((a+b)/r) + KB*yB*p)/(num*p**((a+b)/r)))
                    
                    k_list.append(k)
                    
                    i+=1
                
            else: # c & g
                
                Kg = decimal.Decimal(Kr*KA**a*KB**b/(KR**r))
                
                if b==0:
                    
                    num = decimal.Decimal((Kg*yA**a)**(1/r))
                    
                else:
                    
                    num = decimal.Decimal((Kg*yA**a*yB**b)**(1/r))
                
                if 'All species adsorb molecularly' in Tags: # c
                    
                    i=0
                    while i<len(Exp_rate):
                        
                        p = decimal.Decimal(pressure[i])
                        r = decimal.Decimal(Exp_rate[i])
                        
                        k = decimal.Decimal(r*(1 + (KA*yA + KB*yB)*p + KR*num*p**((a+b)/r))/(num*p**((a+b)/r)))
                        
                        k_list.append(k)
                        
                        i+=1
                
                elif 'A adsorbs dissociatively' in Tags: # g
                    
                    i=0
                    while i<len(Exp_rate):
                        
                        p = decimal.Decimal(pressure[i])
                        r = decimal.Decimal(Exp_rate[i])
                        
                        k = decimal.Decimal(r*(1 + decimal.Decimal(math.sqrt(KA*yA*p)) + KB*yB*p + KR*num*p**((a+b)/r))/(num*p**((a+b)/r)))
                        
                        k_list.append(k)
                        
                        i+=1
    
    elif 'Surface reaction controlling' in Tags:    # d, h, k & n
        
        if 'B does not adsorb' in Tags:    # k & n
            
            num = decimal.Decimal((KA*yA)**a*yB**b)
            
            if 'A adsorbs dissociatively' in Tags:    # n
            
                n = max(2*a, r+s)
                
                i=0
                while i<len(Exp_rate):
                    
                    p = decimal.Decimal(pressure[i])
                    r = decimal.Decimal(Exp_rate[i])
                    
                    k = decimal.Decimal(r*(1 + decimal.Decimal(math.sqrt(KA*yA*p)))**n/(num*p**(a+b)))
                    
                    k_list.append(k)
                    
                    i+=1
            
            else:    # k
                
                n = max(a, r+s)
                
                i=0
                while i<len(Exp_rate):
                    
                    p = decimal.Decimal(pressure[i])
                    r = decimal.Decimal(Exp_rate[i])
                    
                    k = decimal.Decimal(r*(1 + KA*yA*p)**n/(num*p**(a+b)))
                    
                    k_list.append(k)
                    
                    i+=1
        
        else:    # d & h
            
            if b==0:
                
                num = decimal.Decimal((KA*yA)**a)
                
            else:
                
                num = decimal.Decimal((KA*yA)**a*(KB*yB)**b)
        
            if 'All species adsorb molecularly' in Tags:    # d
                
                n = max(a+b, r+s)
                
                i=0
                while i<len(Exp_rate):
                    
                    p = decimal.Decimal(pressure[i])
                    r = decimal.Decimal(Exp_rate[i])
                    
                    k = decimal.Decimal(r*(1 + (KA*yA + KB*yB)*p)**n/(num*p**(a+b)))
                    
                    k_list.append(k)
                    
                    i+=1
            
            elif 'A adsorbs dissociatively' in Tags:    # h
            
                n = max(2*a+b, r+s)
                
                i=0
                while i<len(Exp_rate):
                    
                    p = decimal.Decimal(pressure[i])
                    r = decimal.Decimal(Exp_rate[i])
                    
                    k = decimal.Decimal(r*(1 + decimal.Decimal(math.sqrt(KA*yA*p)) + KB*yB*p)**n/(num*p**(a+b)))
                    
                    k_list.append(k)
                    
                    i+=1
        
    elif 'Impact of A controlling' in Tags:     # o
        
        i=0
        while i<len(Exp_rate):
            
            p = decimal.Decimal(pressure[i])
            r = decimal.Decimal(Exp_rate[i])
            
            k = decimal.Decimal(r*(1 + KB*yB*p)/(KB*yB*yA**(a/b)*p**(1+a/b)))
            
            k_list.append(k)
            
            i+=1
    
    elif 'Uncatalyzed reaction' in Tags:     # r
        
        i=0
        while i<len(Exp_rate):
            
            p = decimal.Decimal(pressure[i])
            r = decimal.Decimal(Exp_rate[i])
            
            k = decimal.Decimal(r/(yA**a*yB**b*p**(a+b)))
            
            k_list.append(k)
            
            i+=1
    
    k = decimal.Decimal(sum(k_list)/len(k_list))
    
    k_list.append(k)
    
    MSE_list = []
    
    yA = float(yA)
    yB = float(yB)
    
    for k in k_list:
        
        k = float(k)
        
        rate = Generate_Curve(Tags, Stoich, yA, yB, k, Eq_Consts, pressure)
        
        Sum_SE = 0
        n=0
        while n<len(rate):
                          
            Sum_SE += (rate[n] - Exp_rate[n])**2 # Adding the squared errors of all datapoints
            
            n+=1
        
        MSE = Sum_SE/len(rate) # Calculating the mean squared error (MSE)
        
        MSE_list.append(MSE)
    
    best_MSE = min(MSE_list)
    
    min_idx = MSE_list.index(best_MSE)
    
    best_k = float(k_list[min_idx])
        
    return best_k, best_MSE


def Possible_Features(Tags, Stoich, yA, yB, pressure, Redux):
    """
    Extracts all possible curve features for a given initial rate equation

    Parameters
    ----------
    Tags : Tuple
        Text strings that describe a mechanism.
    Stoich : List
        Stoichiometric coefficients of species A, B, R and S, respectively.
    yA : Float
        Molar fraction of reactant A.
    yB : Float
        Molar fraction of reactant B.
    pressure : List
        Pressure values used to generate the theoretical curves.
    Redux : Boolean
        True if the user chose to use the reduced lists of possible values for 
        the equilibrium constants.

    Returns
    -------
    Theor_Features : Dictionary
        Contains all the equilibrium constants and the features of curves 
        generated by the given initial rate equation, discriminated by 
        primitives.
        NOTE: if the initial rate equation is considered unrealistic, then 
        this variable becomes of the type String

    """
    
    [a, b, r, s] = Stoich
    
    all_lines = {}
    
    Theor_Features = {}
    
    Unrealistic_mech = False
    
    if (Redux == True) and ('Desorption of R controlling' in Tags) and (b > 0) and (s == 0):
        
        Keq_adsorbed = Keq_adsorbed_redux
        
        Keq_all = Keq_all_redux
        
        Keq_impact = Keq_impact_redux
    
    else:
        
        Keq_adsorbed = Keq_adsorbed_full
        
        Keq_all = Keq_all_full
        
        Keq_impact = Keq_impact_full
        
    j=0
    while j<len(Keq_adsorbed):
                
        if ('A does not adsorb' and 'Desorption of R controlling' in Tags) and (s==0): # A does not adsorb, it impacts B while it is adsorbed
            
            KA = Keq_impact[j]
        
        elif ('Adsorption of A controlling' in Tags) or (('Desorption of R controlling' in Tags) and (s!=0)) or (('A does not adsorb' in Tags) and ('Desorption of R controlling' not in Tags)) or ('Uncatalyzed reaction' in Tags): # No KA on equation
            
            KA = 0
            j = len(Keq_adsorbed)-1 # No point in cycling through different values for KA
                
        else:
            
            KA = Keq_adsorbed[j]
            
        l=0
        while l<len(Keq_adsorbed):
            
            if (b==0) or ('Adsorption of B controlling' in Tags) or ('B does not adsorb' in Tags) or (('Desorption of R controlling' in Tags) and (s!=0)) or ('Uncatalyzed reaction' in Tags): # No KB on equation
                
                KB = 0
                l = len(Keq_adsorbed)-1 # No point in cycling through different values for KB
                    
            else:
                
                KB = Keq_adsorbed[l]
                
            z=0
            while z<len(Keq_adsorbed):
                
                if 'Desorption of R controlling' in Tags:
                    
                    KR = Keq_adsorbed[z]
                
                else: # No KR on equation
                    
                    KR = 0
                    z = len(Keq_adsorbed)-1 # No point in cycling through different values for KR
                
                q=0
                while q<len(Keq_all):
                    
                    if ('Desorption of R controlling' in Tags) and (s==0):
                        
                        Kr = Keq_all[q]
                    
                    else: # No Kr on equation
                        
                        Kr = 0
                        q = len(Keq_all)-1 # No point in cycling through different values for Kr
                    
                    Eq_Consts = [KA, KB, KR, 0, Kr]
                        
                    rate = Generate_Curve(Tags, Stoich, yA, yB, 1, Eq_Consts, pressure)
                    
                    if rate == 'Unrealistic':
                        
                        Unrealistic_mech = True
                        
                        break
                    
                    rate_norm = [x/max(rate) for x in rate] # NORMALIZATION
                    
                    [Primitives, Extremes, y_values, Extrema]=ft.features(pressure, rate_norm)
                    
                    line = [KA, KB, KR, 0, Kr, rate_norm, Extremes, y_values]
                    
                    if tuple(Primitives) in all_lines:
                        
                        all_lines[tuple(Primitives)].append(line)
                        
                    else:
                            
                        all_lines[tuple(Primitives)] = [line]
                    
                    q+=1
                
                if Unrealistic_mech == True:
                    
                    break
                
                z+=1
            
            if Unrealistic_mech == True:
                    
                break
            
            l+=1
        
        if Unrealistic_mech == True:
                    
            break
        
        j+=1
        
           
    if Unrealistic_mech == True:
        
        Theor_Features = 'Unrealistic'
    
    else:
        
        for Primitives in all_lines:
            
            Consts_Table = pd.DataFrame(columns=['KA', 'KB', 'KR', 'KS', 'Kr', 'rate_norm', 'Extremes', 'y_values'])
            
            n=0
            while n<len(all_lines[Primitives]):
                
                Consts_Table.loc[n] = all_lines[Primitives][n]
                
                n+=1
            
            Theor_Features[Primitives] = Consts_Table
    
    return Theor_Features


def Generate_Curve(Tags, Stoich, yA, yB, k, Eq_Consts, pressure):
    """
    Creates coefficients c1, c2, c3, c4, c_extra, power1, power2, power3 & power_extra that describe
    a generic initial rate equation of the type
    
                                     c1*P^power1
    r = --------------------------------------------------------------------
        (1 + c2*P^(1/2) + c3*P + c4*P^power2 + c_extra*P^power_extra)^power3

    Parameters
    ----------
    Tags : Tuple
        Text strings that describe a mechanism.
    Stoich : List
        Stoichiometric coefficients of species A, B, R and S, respectively.
    yA : Float
        Molar fraction of reactant A.
    yB : Float
        Molar fraction of reactant B.
    k : Float
        Value of k' used to generate the curve.
    Eq_Consts : List
        Contains the equilibrium constants of the theoretical curve.
    pressure : List
        Pressure values used to generate the theoretical curves.

    Returns
    -------
    rate : List
        Initial rate values of the generated curve.
        NOTE: if the initial rate equation used to generate this curve is 
        considered unrealistic, then this variable becomes of the type String

    """
   
    [a, b, r, s] = Stoich
    
    [KA, KB, KR, KS, Kr] = Eq_Consts
    
    c_extra = 0     # These are only needed for mechanism p
    power_extra = 0
    
    if 'Adsorption of A controlling' in Tags:    # a, e, i & l
        
        c1 = k*yA
        power1 = 1
        
        c2 = c4 = power2 = 0
        
        if ('B does not adsorb' in Tags) or (b==0):    # i & l
            
            c3 = 0
        
        else:    # a & e
            
            c3 = KB*yB
        
        if 'A adsorbs dissociatively' in Tags:    # e & l
            
            power3 = 2
        
        else:    # a & i
            
            power3 = 1
    
    elif 'Adsorption of B controlling' in Tags:    # b, f & q
        
        c1 = k*yB
        power1 = 1
        
        c4 = power2 = 0
        power3 = 1
        
        if 'A does not adsorb' in Tags:    # q
            
            c2 = c3 = 0
        
        elif 'A adsorbs dissociatively' in Tags:    # f
            
            c2 = math.sqrt(KA*yA)
            c3 = 0
        
        else:    # b
            
            c2 = 0
            c3 = KA*yA
            
    elif 'Desorption of R controlling' in Tags: # c, g, j, m & p
        
        if s!=0: # If there is more than 1 product, the initial rate equation is the same for all mechanisms where desorption of 1 product is the RDS
            
            c1 = k
            power1 = c2 = c3 = power2 = 0
            
            c4 = KR-1 # Subtract 1 because the generic function has a +1 on the denominator and in this case the equation is k/KR
            power3 = 1
            
        else:
            
            if 'B does not adsorb' in Tags: # j & m
                
                Kg = Kr*KA**a/(KR**r)
                
                if 'A adsorbs dissociatively' in Tags: # m
                    
                    c2 = math.sqrt(KA*yA)
                    c3 = 0
                
                else: # j
                    
                    c2 = 0
                    c3 = KA*yA
            
            elif 'A does not adsorb' in Tags: # p
                
                Kg = Kr*(KA*KB)**b/(KR**r)
                
                c2 = 0
                c3 = KB*yB
                
                c_extra = KA*KB*yA**(a/b)*yB
                power_extra = 1 + (a/b)
                
            else: # c & g
                
                Kg = Kr*KA**a*KB**b/(KR**r)
                
                if 'All species adsorb molecularly' in Tags: # c
                    
                    c2 = 0
                    c3 = KA*yA + KB*yB
                
                elif 'A adsorbs dissociatively' in Tags: # g
                    
                    c2 = math.sqrt(KA*yA)
                    c3 = KB*yB
            
            c1 = k*(Kg*yA**a*yB**b)**(1/r)
            
            c4 = KR*(Kg*yA**a*yB**b)**(1/r)
            
            power1 = power2 = (a+b)/r
            power3 = 1
    
    elif 'Surface reaction controlling' in Tags:    # d, h, k & n
        
        c4 = power2 = 0
        power1 = a+b
        
        if 'B does not adsorb' in Tags:    # k & n
            
            c1 = k*(KA*yA)**a*yB**b
            
            if 'A adsorbs dissociatively' in Tags:    # n
            
                c2 = math.sqrt(KA*yA)
                c3 = 0
            
                if 2*a>r+s:
                    power3 = 2*a
                else:
                    power3 = r+s
            
            else:    # k
                
                c2 = 0
                c3 = KA*yA
                
                if a>r+s:
                    power3 = a
                else:
                    power3 = r+s
        
        else:    # d & h
        
            c1 = k*(KA*yA)**a*(KB*yB)**b
        
            if 'All species adsorb molecularly' in Tags:    # d
                
                c2 = 0
                c3 = KA*yA + KB*yB
                
                if a+b>r+s:
                    power3 = a+b
                else:
                    power3 = r+s
            
            elif 'A adsorbs dissociatively' in Tags:    # h
            
                c2 = math.sqrt(KA*yA)
                c3 = KB*yB
            
                if 2*a+b>r+s:
                    power3 = 2*a+b
                else:
                    power3 = r+s
        
    elif 'Impact of A controlling' in Tags:     # o
        
        c1 = k*KB*yB*yA**(a/b)
        power1 = 1+(a/b)
        
        c3 = KB*yB
        power3 = 1
        
        c2 = c4 = power2 = 0
    
    elif 'Uncatalyzed reaction' in Tags:     # r
        
        c1 = k*yA**a*yB**b
        power1 = a+b
        
        c2 = c3 = c4 = power2 = 0
        power3 = 1
        
    
    if power3 > 3:
        
        rate = 'Unrealistic'
    
    else:
        
        rate = [c1*x**power1/((1 + c2*math.sqrt(x) + c3*x + c4*x**power2 + c_extra*x**power_extra)**power3) for x in pressure]
    
    return rate