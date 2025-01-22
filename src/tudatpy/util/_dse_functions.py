import numpy as np
from itertools import combinations as comb

def get_orthogonal_array(no_of_factors : int, 
                no_of_levels : int) -> np.ndarray:
    """ 
    Create orthogonal arrays from Latin Square in 4 successive steps:
    
    0) Take the column from the smaller array to create 2 new
       columns and 2x new rows,
    1) block 1 (1/2 rows): take old values 2x for new columns,
    2) block 2 (1/2 rows): take old values, use Latin-Square for new
       columns,
    3) column 1: divide experiments into groups of 1,2.

    Parameters
    ----------
    no_of_factors : int
        Integer representing the number of design variables.
    no_of_levels : int
        Integer representing the number of levels that will be explored within each bound.

    Returns
    -------
    Lx : np.ndarray
        Orthogonal array containing the FFD experiments.

    """
    # TODO: merge the no_of_levels = 2 and no_of_levels = 3 parts of the code (instead of two separate if statements
                                                                               # with a lot of overlap

    assert no_of_levels >= 2
    assert no_of_factors >= 1

    row_number = lambda icount, no_of_levels : no_of_levels**(icount+1)
    col_number = lambda row_number : row_number-1

    if no_of_levels == 2:

        if no_of_factors == 1:
            icount = 0
        elif no_of_factors >= 2 and no_of_factors <= 3:
            icount = 1
        elif no_of_factors >= 4 and no_of_factors <= 7:
            icount = 2
        elif no_of_factors >= 8 and no_of_factors <= 15:
            icount = 3
        elif no_of_factors >= 16 and no_of_factors <= 31:
            icount = 4
        elif no_of_factors >= 32 and no_of_factors <= 63:
            icount = 5
        elif no_of_factors >= 64 and no_of_factors <= 127:
            icount = 6
        elif no_of_factors >= 128 and no_of_factors <= 255:
            icount = 7
        else:
            raise RuntimeError("The number of factors provided is not implemented. Provide a value between 1 and 255")

        Lxrow = row_number(icount, no_of_levels)
        Lxcol = col_number(Lxrow)
        Lx = np.zeros((Lxrow,Lxcol), dtype='int')
        iaux = Lx.copy()
        
        # Define the 2-level Latin Square
        index_list = [0, 1]
        two_level = [-1, 1]
        LS = np.zeros((2,2), dtype='int')
        LS[0,0] = -1
        LS[0,1] =  1
        LS[1,0] =  1
        LS[1,1] = -1
        
        # In case of only one factor, copy the first Latin Square and leave the subroutine.
        if icount == 0:
            Lx[0,0] = LS[0,0]
            Lx[1,0] = LS[0,1]
            return Lx
        
        iaux[0,0] = -1
        iaux[1,0] =  1
        irow = 2
        icol = 1

        # Algorithm in Matlab starts from index 1
        Lx = np.hstack((np.zeros((len(Lx), 1)), Lx))
        Lx = np.vstack((np.zeros((1, len(Lx[0,:]))), Lx))
        iaux = np.hstack((np.zeros((len(iaux), 1)), iaux))
        iaux = np.vstack((np.zeros((1, len(iaux[0,:]))), iaux))
        
        # Fill in orthogonal array
        for i1 in range(1, icount + 1):
            for i2 in range(1, irow + 1):
                for i3 in range(1, icol + 1):
                    for p in range(2):
                        for q in range(2):
                            for r in range(2):

                                # Block 1.
                                if iaux[i2,i3] == two_level[q] and p == 0:
                                    Lx[i2,i3*2 + index_list[r]] = two_level[q] 

                                # Block 2
                                if iaux[i2,i3] == two_level[q] and p == 1:
                                    Lx[i2 + irow,i3*2 + index_list[r]] = LS[index_list[q], index_list[r]]

                        Lx[i2 + irow*p,1] = two_level[p]

            if i1 == icount:
                # Deleting extra row from Matlab artifact
                Lx = np.delete(Lx, 0, 0)
                Lx = np.delete(Lx, 0, 1)
                return Lx

            irow = 2*irow
            icol = 2*icol+1
            for i2 in range(1, irow + 1):
                for i3 in range(1, icol + 1):
                    iaux[i2,i3] = Lx[i2,i3]

    elif no_of_levels == 3:

        if no_of_factors == 1:
            icount = 0
        elif no_of_factors >= 2 and no_of_factors <= 4:
            icount = 1
        elif no_of_factors >= 5 and no_of_factors <= 13:
            icount = 2
        elif no_of_factors >= 14 and no_of_factors <= 40:
            icount = 3
        elif no_of_factors >= 41 and no_of_factors <= 121:
            icount = 4
        else:
            raise RuntimeError("The number of factors provided is not implemented. Provide a value between 0 and 255")

        Lxrow = row_number(icount, no_of_levels)
        Lxcol = col_number(Lxrow) // 2
        Lx = np.zeros((Lxrow,Lxcol))
        iaux = Lx.copy()
        
        # Define Latin Square 1
        index_list = [0, 1, 2]
        
        LS1 = np.zeros((3,3))
        three_level = [-1, 0, 1]
        for i in range(3):
            for j in range(3):
                LS1[i,index_list[j]] = three_level[(j+i)%3];

        # Define Latin Square 2
        LS2 = np.zeros((3,3))
        three_level_2 = [-1, 1, 0]
        for i in range(3):
            for j in range(3):
                LS2[i, index_list[j]] = three_level_2[j-i]
     
        # In case of only one factor, copy the first Latin Square and leave the subroutine
        if icount == 0:
           Lx[0,0] = LS1[0,0];
           Lx[1,0] = LS1[0,1];
           Lx[2,0] = LS1[0,2];
           return Lx

        # Define iaux for loops
        iaux[0,0] = -1
        iaux[1,0] = 0
        iaux[2,0] =  1
        irow = 3
        icol = 1

        # Algorithm in Matlab starts from index 1
        Lx = np.hstack((np.zeros((len(Lx), 1)), Lx))
        Lx = np.vstack((np.zeros((1, len(Lx[0,:]))), Lx))
        iaux = np.hstack((np.zeros((len(iaux), 1)), iaux))
        iaux = np.vstack((np.zeros((1, len(iaux[0,:]))), iaux))
        
        # Filling in Lx
        for i1 in range(1, icount + 1):
            for i2 in range(1, irow + 1):
                for i3 in range(1, icol + 1):
                    for p in range(3):
                        for q in range(3):
                            for r in range(3):

                                #Block 1.
                                if iaux[i2,i3] == three_level[q] and p == 0:
                                    Lx[i2 + irow*p,i3*3 + three_level[r]] = three_level[q] 

                                #Block 2.
                                if iaux[i2,i3] == three_level[q] and p == 1:
                                    Lx[i2 + irow*p,i3*3 + three_level[r]] = LS1[index_list[q], index_list[r]]

                                #Block 3.
                                if iaux[i2,i3] == three_level[q] and p == 2:
                                    Lx[i2 + irow*p,i3*3 + three_level[r]] = LS2[index_list[q], index_list[r]]

                        Lx[i2 + irow*p,1] = three_level[p]

            if i1 == icount:
                # Deleting extra row from Matlab artifact
                Lx = np.delete(Lx, 0, 0)
                Lx = np.delete(Lx, 0, 1)
                return Lx

            irow = 3*irow
            icol = 3*icol+1
            for i2 in range(1, irow + 1):
                for i3 in range(1, icol + 1):
                        iaux[i2,i3] = Lx[i2,i3]
    else:
        raise RuntimeError("These levels are not implemented yet. (You may wonder whether you need them)")

def get_yates_array(no_of_factors : int, 
                        no_of_levels : int) -> np.ndarray:
    """
    Function that creates a Yates array according to Yates algorithm

    Parameters
    ----------
    no_of_factors : int
        Integer representing the number of design variables.
    no_of_levels : int
        Integer representing the number of levels that will be explored within each bound.

    Returns
    -------
    yates_array : np.ndarray
        Orthogonal array containing the FFD experiments.
    """

    assert no_of_factors >= 1
    assert no_of_levels >= 2

    if no_of_levels == 2: #for the ANOVA analysis with interactions
        levels = [-1, 1]
    elif no_of_levels == 3:
        levels = [-1, 0, 1]
    elif no_of_levels == 7: #for the response surfaces
        levels = [-3, -2, -1, 0, 1, 2, 3]
    else:
        raise RuntimeError("The number of levels you have requested is not implemented")

    n_rows = no_of_levels**no_of_factors
    n_cols = no_of_factors
    yates_array = np.zeros((n_rows, n_cols), dtype='int')

    row_seg = n_rows
    for col in range(n_cols):
        repetition_amount = no_of_levels**col # Number of times to repeat the row segment to fill the array
        row_seg = row_seg // no_of_levels # Get row segment divided by number of levels

        for j in range(repetition_amount):
            for it, level in enumerate(levels):
                # fill in from position i to position i + no_of_levels
                yates_array[(it*row_seg + j*row_seg*no_of_levels):((it+1)*row_seg +
                j*row_seg*no_of_levels), col] = levels[it] 
    return yates_array 

# ANOVA short explanation
# After some definitions, we iterate through yates array, depending on the value being -1 or 1, we add the occurrence
# to a certain container. Later the amount of items in the container determines the contribution of that specific
# variable to an objective value. This iteration is done for individual, linear, and quadratic effects, thereby
# determing the contribution of all interactions.

def anova_analysis(objective_values, 
                   factorial_design_array : np.ndarray, 
                   no_of_factors : int, 
                   no_of_levels : int,
                   level_of_interactions=2) -> tuple:
    """
    Function that performs an Analysis of Variance (ANOVA) for no_of_levels levels and no_of_factors factors. 

    Parameters
    ----------
    objective_values : list 
        List of objective function values following from the factorial design array.
    factorial_design_array : np.array 
        An array compatible with factorial design, generally a Yates array.
    no_of_factors : int
        Integer representing the number of design variables.
    no_of_levels : int
        Integer representing the number of levels that will be explored within each bound.
    level_of_interactions : int
        Integer either equal to 2 or 3, depending on what interactions you want to
        include

    Returns
    -------
    percentage_contribution_i : np.ndarray
        Array of percentage contributions of individual effects for each factor
    percentage_contribution_ii : np.ndarray
        Array of percentage contributions of linear effects for each factor
    percentage_contribution_iii : np.ndarray
        Array of percentage contributions of quadratic effects for each factor
    percentage_contribution_error : np.ndarray
        Float for percentage contribution of error

    """

    assert len(objective_values) == len(factorial_design_array) # Quick check for validity of data

    number_of_simulations = no_of_levels ** no_of_factors

    #########################
    ### Array definitions ###
    #########################

    # Lambda function to determine number of interaction columns
    interactions = lambda iterable, k : [i for i in comb(range(iterable), k)]
    no_of_interactions_2 = interactions(no_of_factors, 2) # number of 2-level interactions
    no_of_interactions_3 = interactions(no_of_factors, 3) # number of 3-level interactions

    # Arrays for individual, 2, and 3 level interactions - sum of squares
    sum_of_squares_i = np.zeros(no_of_factors)
    sum_of_squares_ij = np.zeros(len(no_of_interactions_2))
    sum_of_squares_ijk = np.zeros(len(no_of_interactions_3))

    #Arrays for individual, 2, and 3 level interactions - percentage contribution
    percentage_contribution_i = np.zeros(no_of_factors)
    percentage_contribution_ij = np.zeros(len(no_of_interactions_2))
    percentage_contribution_ijk = np.zeros(len(no_of_interactions_3))

    # Sum of objective values and mean of parameters
    sum_of_objective = np.sum(objective_values)
    mean_of_param = sum_of_objective/number_of_simulations

    # Variance and standard deviation of data
    variance_per_run = np.zeros(number_of_simulations)
    objective_values_squared = np.zeros(number_of_simulations)
    for i in range(len(factorial_design_array)):
        variance_per_run[i] = (objective_values[i] - mean_of_param)**2 / (number_of_simulations-1)
        objective_values_squared[i] = objective_values[i]**2
    variance = np.sum(variance_per_run)
    standard_deviation = np.sqrt(variance)
    
    # Square of sums
    CF = sum_of_objective * mean_of_param
    # Difference square of sums and sum of squares
    sum_of_deviation = np.sum(objective_values_squared) - CF

    #####################################
    ### Iterations through yates array ###
    #####################################

    ### Linear effects ###
    # Container for appearance of minimum/maximum value
    Saux = np.zeros((no_of_factors, no_of_levels))
    number_of_appearanes = number_of_simulations/no_of_levels
    for j in range(number_of_simulations):
        for i in range(no_of_factors):
            if factorial_design_array[j,i] == -1:
                Saux[i, 0] += objective_values[j]
            else:
                Saux[i, 1] += objective_values[j]


    if sum_of_deviation > 1e-6: # If there is a deviation, then there is a contribution.
        for i in range(no_of_factors):
            sum_of_squares_i[i] = (1/2)*(Saux[i, 1] - Saux[i, 0])**2/number_of_appearanes
            percentage_contribution_i[i] = 100 * sum_of_squares_i[i]/sum_of_deviation

    ### 2-level interaction ###
    # Container for appearance of minimum/maximum value
    Saux = np.zeros((len(no_of_interactions_2), no_of_levels))
    for j in range(number_of_simulations):
        for i in range(len(no_of_interactions_2)):
            # Interaction sequence of all possible 2-level interactions is created by multiplying the
            # two respective elements from yates array
            interaction = factorial_design_array[j, no_of_interactions_2[i][0]] * \
                            factorial_design_array[j, no_of_interactions_2[i][1]]
            if interaction == -1:
                Saux[i, 0] += objective_values[j]
            else:
                Saux[i, 1] += objective_values[j]

    if sum_of_deviation > 1e-6: # If there is a deviation, then there is a contribution.
        for i in range(len(no_of_interactions_2)):
            sum_of_squares_ij[i] = (1/2)*(Saux[i, 1] - Saux[i, 0])**2/number_of_appearanes
            percentage_contribution_ij[i] = 100 * sum_of_squares_ij[i]/sum_of_deviation

    ### 3-level interaction ###
    # Container for appearance of minimum/maximum value
    Saux = np.zeros((len(no_of_interactions_3), no_of_levels))
    for j in range(number_of_simulations):
        for i in range(len(no_of_interactions_3)):
            # Interaction sequence of all possible 3-level interactions is created by multiplying the
            # three respective elements from yates array
            interaction = factorial_design_array[j, no_of_interactions_3[i][0]] * \
                            factorial_design_array[j, no_of_interactions_3[i][1]] * \
                              factorial_design_array[j, no_of_interactions_3[i][2]]
            if interaction == -1:
                Saux[i, 0] += objective_values[j]
            else:
                Saux[i, 1] += objective_values[j]

    if sum_of_deviation > 1e-6: # If there is a deviation, then there is a contribution.
        for i in range(len(no_of_interactions_3)):
            sum_of_squares_ijk[i] = (1/2)*(Saux[i, 1] - Saux[i, 0])**2/number_of_appearanes
            percentage_contribution_ijk[i] = 100 * sum_of_squares_ijk[i]/sum_of_deviation

    ### Error contribution ###
    sum_of_squares_error = sum_of_deviation - np.sum(sum_of_squares_i) - \
    np.sum(sum_of_squares_ij) - np.sum(sum_of_squares_ijk)

    percentage_contribution_error = 0
    if sum_of_deviation > 1e-6: # If there is a deviation, then there is a contribution.
        percentage_contribution_error = 100 * sum_of_squares_error / sum_of_deviation

    return percentage_contribution_i, percentage_contribution_ij, \
    percentage_contribution_ijk, percentage_contribution_error

