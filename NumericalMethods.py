#region imports
import Gauss_Elim as GE  # this is the module from lecture 2 that has usefule matrix manipulation functions
from math import sqrt, pi, exp, cos
from copy import deepcopy as dc
#endregion

#region function definitions
def Probability(PDF, args, c, GT=True):
    """
    This is the function to calculate the probability that x is >c or <c depending
    on the GT boolean.
    Step 1:  unpack args into mu and stDev
    Step 2:  compute lhl and rhl for Simpson
    Step 3:  package new tuple args1=(mu, stDev, lhl, rhl) to be passed to Simpson
    Step 4:  call Simpson with GNPDF and args1
    Step 5:  return probability
    :param PDF: the probability density function to be integrated
    :param args: a tuple with (mean, standard deviation)
    :param c: value for which we ask the probability question
    :param GT: boolean deciding if we want probability x>c (True) or x<c (False)
    :return: probability value
    """
    # Step 1: Unpack args into mu and stDev
    mu, stDev = args

    # Step 2: Compute left-hand limit (lhl) and right-hand limit (rhl) for Simpson
    if GT:
        # probability of x > c --> integrate from c to (mu + 5sigma)
        lhl, rhl = c, mu + 5 * stDev
    else:
        # probability of x < c --> integrate from (mu - 5sigma) to c
        lhl, rhl = mu - 5 * stDev, c

    # Step 3: Package new tuple args1=(mu, stDev, lhl, rhl)
    args1 = (mu, stDev, lhl, rhl)

    # Step 4: Call Simpson with GPDF and args1
    probability = Simpson(GPDF, args1)

    # Step 5: Return probability
    return probability

def GPDF(args):
    """
    Here is where I will define the Gaussian probability density function.
    This requires knowing the population mean and standard deviation.
    To compute the GPDF at any value of x, I just need to compute as stated
    in the homework assignment.
    Step 1:  unpack the args tuple into variables called: x, mu, stDev
    Step 2:  compute GPDF value at x
    Step 3:  return value
    :param args: (x, mean, standard deviation)  tuple in that order
    :return: value of GPDF at the desired x
    """
    # Step 1: unpack args
    x, mu, sig = args
    # step 2: compute GPDF at x
    fx = (1 / (sig * sqrt(2 * pi))) * exp(-0.5 * ((x - mu) / sig) ** 2)
    # step 3: return value
    return fx

def Simpson(fn, args, N=100):
    """
    This executes the Simpson 1/3 rule for numerical integration (see page 832, Table 19.4).
    As I recall:
    1. divide the range from x=lhl to x=rhl into an even number of parts. Perhaps 20?
    2. compute fx at each x value between lhl and rhl
    3. sum the even and odd values of fx as prescribed
    4. return the area beneath the function fx
    :param fx: some function of x to integrate
    :param args: a tuple containing (mean, stDev, lhl, rhl)
    :return: the area beneath the function between lhl and rhl
    """
    mu, sig, a, b = args
    if N % 2 == 1:
        N += 1 # ensure N is even
    h = (b -a) / N
    s= fn((a, mu, sig)) + fn((b, mu, sig)) # first and last terms

    for i in range(1, N, 2):
        s += 4 * fn((a + i * h, mu, sig)) # odd terms
    for i in range(2, N-1, 2):
        s += 2 * fn((a + i * h, mu, sig)) # even terms
    return (h/3) * s

def Secant(fcn, x0, x1, maxiter=10, xtol=1e-5):
    """
    This funciton implements th Secant method to find the root of an equation.  You should write your equation in a form
    fcn = 0 such that when the correct value of x is selected, the fcn actually equals zero (or very close to it).
    :param fcn: the function for which we want to find the root
    :param x0: x value in neighborhood of root (or guess 1)
    :param x1: another x value in neighborhood of root (or guess x0+1)
    :param maxiter: exit if the number of iterations (new x values) equals this number
    :param xtol:  exit if the |xnewest - xprevious| < xtol
    :return: tuple with: (the final estimate of the root (most recent value of x), number of iterations)
    """
    """
    Step 1: evaluate fcn at x0 and x1
    Step 2: compute new x using the Secant formula
    Step 3: Check if stopping criteria are met (|x_new - x_prev| < xtol)
    Step 4: Update x values and repeat until max iterations reached
    Step 5: Return estimated root and # of iterations
    """
    for i in range(maxiter):
        # compute function values
        f_x0 = fcn(x0)
        f_x1 = fcn(x1)

        # prevent division by zero or very small denominator
        if abs(f_x1 - f_x0) < 1e-12:
            print("Warning: Division by zero in Secant Method.")
            return (x1, i)

        # secant method formula to compute next x
        x_new = x1 - f_x1 * (x1 - x0) / (f_x1 - f_x0)
        f_new = fcn(x_new) # evaluate at new x

        # checking to see if we found a root (|f(x_new)| < ftol)
        if abs(f_new) < xtol:
            return (x_new, i + 1)

        # stopping condition
        if abs(x_new - x1) < xtol:
            return (x_new, i + 1)

        # update values
        x0, x1 = x1, x_new

    print("Warning: Maximum iterations reached, best estimate returned.")
    return (x1, maxiter) # return last computed value if maxiter is reached
    pass

def GaussSeidel(Aaug, x, Niter = 15):
    """
    This should implement the Gauss-Seidel method (see page 860, Tabl 20.2) for solving a system of equations.
    :param Aaug: The augmented matrix from Ax=b -> [A|b]
    :param x:  An initial guess for the x vector. if A is nxn, x is nx1
    :param Niter:  Number of iterations to run the GS method
    :return: the solution vector x
    """
    N = len(Aaug)  # Number of equations
    Aaug = GE.MakeDiagDom(Aaug)  # Ensure the matrix is diagonally dominant
    A = [row[:-1] for row in Aaug]  # Extract coefficient matrix A
    b = [row[-1] for row in Aaug]  # Extract right-hand side vector b

    x = dc(x)  # Prevent modifying the original input

    # Iterative process
    for _ in range(Niter):
        x_old = x[:]  # Store previous iteration values

        for j in range(N):  # Loop over each equation
            sum1 = sum(A[j][k] * x[k] for k in range(j))  # Using updated values
            sum2 = sum(A[j][k] * x_old[k] for k in range(j + 1, N))  # Using old values

            if A[j][j] == 0:
                print(f"Error: Zero on diagonal at row {j}, cannot proceed.")

            x[j] = (b[j] - sum1 - sum2) / A[j][j]  # Update x[j]

        # Convergence check: max |x_new - x_old|
        max_diff = max(abs(x[j] - x_old[j]) for j in range(N))
        if max_diff < 1e-5:  # Default tolerance
            return x  # Procedure completed successfully

    print("No solution satisfying the tolerance condition obtained after", Niter, "iteration steps.")
    return None  # Procedure completed unsuccessfully

def main():
    '''
    This is a function I created for testing the numerical methods locally.
    :return: None
    '''
    #region testing GPDF
    fx = GPDF((0,0,1))
    print("{:0.5f}".format(fx))  # Does this match the expected value?
    #edregion

    #region testing Simpson
    p=Simpson(GPDF,(0,1,-5,0)) # should return 0.5
    print("p={:0.5f}".format(p))  # Does this match the expected value?
    #endregion

    #region testing Probability
    p1 = Probability(GPDF, (0,1),0,True)
    print("p1={:0.5f}".format(p1))  # Does this match the expected value?
    #endregion
    pass

#endregion

#region function calls
if __name__ == '__main__':
    main()
#endregion