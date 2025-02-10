#region imports
from math import sqrt, pi, exp
from NumericalMethods import GPDF, Simpson, Probability
#endregion

#region function definitions
def main():
    """
    I want to integrate the Gaussian probability density function between
    a left hand limit = (mean - 5*stDev) to a right hand limit = (c).  Here
    is my step-by-step plan:
    1. Decide mean, stDev, and c and if I want P(x>c) or P(x<c).
    2. Define args tuple and c to be passed to Probability
    3. Pass args, and a callback function (GPDF) to Probability
    4. In probability, pass along GPDF to Simpson along with the appropriate args tuple
    5. Return the required probability from Probability and print to screen.
    :return: Nothing to return, just print results to screen.
    """
    """ 
    Computes probabilities using the Probability function in NumericalMethods.py 
    """
    #region testing user input
    # The following code solicites user input through the CLI.
    mean = float(input("Population mean? ").strip()) # convert to float
    stDev = float(input("Standard deviation?").strip()) # convert to float
    c = float(input("c value?").strip()) # convert to float
    GT = input("Probability greater than c? (y/n)").strip().lower() in ["y","yes","true"]

    # compute probability
    p = Probability(GPDF, (mean, stDev), c, GT)

    print(f"P(x{'<' if GT else '<'}{c:.2f}|N({mean},{stDev}))={p:.2f}")
    #endregion
#endregion

if __name__ == "__main__":
    main()