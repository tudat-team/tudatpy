from tudatpy import util
import numpy as np

def test_dse_functions():
    # print('Running unit tests on DSE functions')

    ##################################
    ## get_orthogonal_array testing ##
    ##################################

    # boundary cases
    assert util.get_orthogonal_array(1, 2).all() == np.array([[-1], [1]]).all() 
    assert util.get_orthogonal_array(1, 3).all() == np.array([[-1], [0], [1]]).all() 

    # from https://www.york.ac.uk/depts/maths/tables/orthogonal.htm
    assert util.get_orthogonal_array(3, 2).all() == np.array([[-1, -1, -1],
                                                                       [-1,  1,  1],
                                                                       [ 1, -1,  1],
                                                                       [ 1,  1, -1]]).all() 

    # separate not independently verified
    assert util.get_orthogonal_array(4, 3).all() == np.array([[-1, -1, -1, -1],
                                                                       [-1,  0,  0,  0],
                                                                       [-1,  1,  1,  1],
                                                                       [ 0, -1,  0,  1],
                                                                       [ 0,  0,  1, -1],
                                                                       [ 0,  1, -1,  0],
                                                                       [ 1, -1,  1,  0],
                                                                       [ 1,  0, -1,  1],
                                                                       [ 1,  1,  0, -1]]).all()

    # TODO: There should maybe be testing for raising errors
    # try:
    #     util.get_orthogonal_array(0, 0)
    #     assert True == False
    # except RuntimeError as e:
    #     print("We got here")
    #     assert e == "These levels are not implemented yet. (You may wonder whether you need them)"


    #############################
    ## get_yates_array testing ##
    #############################

    #boundary case
    assert util.get_yates_array(1, 2).all()  == np.array([[-1, 1]]).all()

    # from https://www.itl.nist.gov/div898/handbook/eda/section3/eda35i.htm
    assert util.get_yates_array(3, 2).all()  == np.array([[-1, -1, -1],
                                                                   [-1, -1,  1],
                                                                   [-1,  1, -1],
                                                                   [-1,  1,  1],
                                                                   [ 1, -1, -1],
                                                                   [ 1, -1,  1],
                                                                   [ 1,  1, -1],
                                                                   [ 1,  1,  1]]).all() 

    # separate not independently verified
    assert util.get_yates_array(5, 2).all() == np.array([[-1, -1, -1, -1, -1],
                                                           [-1, -1, -1, -1,  1],
                                                           [-1, -1, -1,  1, -1],
                                                           [-1, -1, -1,  1,  1],
                                                           [-1, -1,  1, -1, -1],
                                                           [-1, -1,  1, -1,  1],
                                                           [-1, -1,  1,  1, -1],
                                                           [-1, -1,  1,  1,  1],
                                                           [-1,  1, -1, -1, -1],
                                                           [-1,  1, -1, -1,  1],
                                                           [-1,  1, -1,  1, -1],
                                                           [-1,  1, -1,  1,  1],
                                                           [-1,  1,  1, -1, -1],
                                                           [-1,  1,  1, -1,  1],
                                                           [-1,  1,  1,  1, -1],
                                                           [-1,  1,  1,  1,  1],
                                                           [ 1, -1, -1, -1, -1],
                                                           [ 1, -1, -1, -1,  1],
                                                           [ 1, -1, -1,  1, -1],
                                                           [ 1, -1, -1,  1,  1],
                                                           [ 1, -1,  1, -1, -1],
                                                           [ 1, -1,  1, -1,  1],
                                                           [ 1, -1,  1,  1, -1],
                                                           [ 1, -1,  1,  1,  1],
                                                           [ 1,  1, -1, -1, -1],
                                                           [ 1,  1, -1, -1,  1],
                                                           [ 1,  1, -1,  1, -1],
                                                           [ 1,  1, -1,  1,  1],
                                                           [ 1,  1,  1, -1, -1],
                                                           [ 1,  1,  1, -1,  1],
                                                           [ 1,  1,  1,  1, -1],
                                                           [ 1,  1,  1,  1,  1]]).all()

    # separate not independently verified
    assert util.get_yates_array(2, 7).all() == np.array([[-3, -3],
                                                           [-3, -2],
                                                           [-3, -1],
                                                           [-3,  0],
                                                           [-3,  1],
                                                           [-3,  2],
                                                           [-3,  3],
                                                           [-2, -3],
                                                           [-2, -2],
                                                           [-2, -1],
                                                           [-2,  0],
                                                           [-2,  1],
                                                           [-2,  2],
                                                           [-2,  3],
                                                           [-1, -3],
                                                           [-1, -2],
                                                           [-1, -1],
                                                           [-1,  0],
                                                           [-1,  1],
                                                           [-1,  2],
                                                           [-1,  3],
                                                           [ 0, -3],
                                                           [ 0, -2],
                                                           [ 0, -1],
                                                           [ 0,  0],
                                                           [ 0,  1],
                                                           [ 0,  2],
                                                           [ 0,  3],
                                                           [ 1, -3],
                                                           [ 1, -2],
                                                           [ 1, -1],
                                                           [ 1,  0],
                                                           [ 1,  1],
                                                           [ 1,  2],
                                                           [ 1,  3],
                                                           [ 2, -3],
                                                           [ 2, -2],
                                                           [ 2, -1],
                                                           [ 2,  0],
                                                           [ 2,  1],
                                                           [ 2,  2],
                                                           [ 2,  3],
                                                           [ 3, -3],
                                                           [ 3, -2],
                                                           [ 3, -1],
                                                           [ 3,  0],
                                                           [ 3,  1],
                                                           [ 3,  2],
                                                           [ 3,  3]]).all()

    ############################
    ## anova_analysis testing ##
    ############################

    # test from asteroid orbit optimization with PyGMO example.
    # TODO: This testing requires independent verification. Now it only verifies itself.

    #These results also indirectly verify the yates array, as the results are physically plausible

    no_of_levels = 2
    no_of_factors = 4
    yates_array = util.get_yates_array(no_of_factors, no_of_levels)
    mean_distances = np.array([ 299.99999298,  299.99999298,  299.99999298,  299.99999298,
        354.26698789,  354.26698789,  354.26698789,  354.26698789,
       2028.0757654 , 2028.11203823, 1911.74108971, 1911.73835608,
       2035.58331439, 2036.22834897, 2016.1089744 , 2015.59407609])

    i, ij, ijk, ierror = util.anova_analysis(objective_values=mean_distances,
                                        factorial_design_array=yates_array,
                                        no_of_factors=no_of_factors,
                                        no_of_levels=no_of_levels,
                                        level_of_interactions=2)

    i = np.array([round(p, 12) for p in i])
    ij = np.array([round(p, 12) for p in ij])
    ijk = np.array([round(p, 12) for p in ijk])
    ierror = round(ierror, 12)

    i_verify = np.array([9.97668486e+01, 1.08564179e-01, 4.15641300e-02, 1.49603978e-08])
    i_verify = np.array([round(p, 12) for p in i_verify])
    ij_verify = np.array([2.56666572e-05, 4.15641300e-02, 1.49603978e-08, 2.07151148e-02, 5.21079366e-09, 8.02731604e-07])
    ij_verify = np.array([round(p, 12) for p in ij_verify])
    ijk_verify = np.array([2.07151148e-02, 5.21079366e-09, 8.02731604e-07, 7.01665549e-07])
    ijk_verify = np.array([round(p, 12) for p in ijk_verify])
    ierror_verify = round(7.016655796880146e-07, 12)

    assert  (i.all(), ij.all(), ijk.all(), ierror) == (i_verify.all(), ij_verify.all(), ijk_verify.all(), ierror_verify)

#############################
### testing functionality ###
#############################

# test_dse_functions() #for separate testing (commented when used for tudatpy compilation)
