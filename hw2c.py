from copy import deepcopy
from NumericalMethods import GaussSeidel

def main():
    """
    Using Gauss-Seidel method to solve linear systems
    """
    #first system
    Aaug1 = [
        [3, 1, -1, 2], # [A | b]
        [1, 4, 1, 12],
        [2, 1, 2, 10]
    ]
    x_guess1 = [0, 0, 0] # initial guess

    # second system
    Aaug2 = [
        [1, -10, 2, 4, 2], # [A | b]
        [3, 1, 4, 12, 12],
        [9, 2, 3, 4, 21],
        [-1, 2, 7, 3, 37]
    ]
    x_guess2 = [0, 0, 0, 0] # initial guess

    # solving using Gauss-Seidel
    x_solution1 = GaussSeidel(Aaug1, x_guess1, 15)
    x_solution2 = GaussSeidel(Aaug2, x_guess2, 15)

    # print results
    print("Solution for first system:", ["{:.4f}".format(val) for val in x_solution1])
    print("Solution for second system:", ["{:.4f}".format(val) for val in x_solution2])

    pass

if __name__=="__main__":
    main()