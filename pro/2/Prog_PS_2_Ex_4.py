from math import log


def f(x):
    return 0.25*x**4 # This is only the test function


class Central: # Initiate the class Central

    def __init__(self, function, increment):
        self.f = function
        self.increment = increment

    def __call__(self, var):  # Need to use __call__ in order for the instance to be used as function
        h = self.increment
        f = self.f  # name the function, within the derivative method, not before
        numerator = f(var + h) - f(var - h)
        denominator = float(2 * h)
        der = numerator / denominator
        return der


def test(): # snippet of the code provided in the exercise
    testies = [1, 5, 10]
    for i in testies:
        dff = Central(f, 1E-9)
        dff_estim = dff(i)
        dff_real = i**3
        print (i, dff_estim, dff_real)


def derivative():
    h = [0.5, 0.1, 1E-3, 1E-5, 1E-7, 1E-9, 1E-11]
    var = input("Enter your desired value for x: ")
    print("For function log(x) and an x = %s" % var)
    print "%10s    %10s    %10s" % ("Increment", "Estimate deriv.", "Exact deriv.")  # allows for spaces definition for our table to look nice
    for i in h: # this will allow to evaluate the derivative for different increments (h)
        df = Central(log, i)  # call the class
        df_estim = df(var)  # use the class with the defined variable
        df_real = (1.0 / var)  # the real derivative computed "by hand"
        print "%10s    %10.7f        %10.2f" % (i, df_estim, df_real)  # space definition and decimals setting for floats


if __name__ == "__main__":  # this line allows for our module to be used in another script without printing \
                            # the derivative function
    derivative()
