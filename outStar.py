"""
outStar.py

This program implements Grossberg’s Outstar network
The source node has activity x0 and the four sink nodes have activities x1, x2, x3, and x4.
w1, w2, w3, and w4 be the corresponding synaptic strengths.
This program applies the update rules on x0, x1, x2, x3, x4, w1, w2, w3, and w4

@author: Anushree Das (ad1707)
"""

import matplotlib.pyplot as plt

# condition for I values for part (a) of assignment
def A(t):
    theta = [0.4, 0.3, 0.2, 0.1]
    # I_base = 2 for the two time steps directly after those times when I0 = 2,
    # i.e. step 2,3,12,13,22,23,..
    # and I = 0 on other time steps
    if (t - 2) % 10 == 0 or (t - 3) % 10 == 0:
        I_base = 2
    else:
        I_base = 0
    I = [thetai * I_base for thetai in theta]
    return I


# condition for I values for part (b) of assignment
def B(t):
    # After times when I0 is nonzero, the other Ii become nonzero for one time step
    # but alternate between I1 = .8, I2 = .6, I3 =.4,I4 =.2 for steps 2,22,42,62,..
    # and I1 =.5,I2 =.5,I3 =.3,I4 =.7. for steps 12,32,52,72,..
    if (t - 2) % 10 == 0:
        if t-2 % 20 != 0:
            I = [0.8,0.6,0.4,0.2]
        else:
            I = [0.5,0.5,0.3,0.7]
    else:
        I = [0 for _ in range(4)]
    return I


# update rule for x0
def update_x0(x0,I0):
    return (-5*x0) + I0


# update rule for x
def update_x(xi,x0,wi,Ii):
    return (-5*xi) + (x0 * wi) + Ii


# update rule for w
def update_w(x0,wi,xi):
    return x0 * ((-0.1*wi) + xi)


def euler(h1, h2, tb, func, initial_values, F):
    """
    :param h1: step value for x0 and x
    :param h2: step value for w
    :param tb: bound for t
    :param func: action to perform every loop
    :param initial_values: vector of initial values for x0,x,w
    :param F: vector of differential functions
    :return: None
    """
    # lists to store values of X, W and theta to plot at the end for comparison
    X_hist = []
    W_hist = []
    theta_hist = []

    # extract initial values for x0,x,w
    x0 = initial_values[0]
    x = initial_values[1]
    w = initial_values[2]

    # initialize theta
    theta = [0.4, 0.3, 0.2, 0.1]
    t = 0

    while t <= tb:
        # I0 = 2 for two steps on every tenth time step, starting with the first,
        # i.e. step 0,1,10,11,20,21,..,
        # and set I0 = 0 on all other time steps
        if (t - 1) % 10 == 0 or t % 10 == 0:
            I0 = 2
        else:
            I0 = 0

        # get Ii for i =1,2,3,4
        I = func(t)

        # update x0
        x0 = x0 + (h1 * F[0](x0,I0))

        # update x
        f = []
        for i in range(len(x)):
            f.append(F[1](x[i],x0,w[i],(I[i])))
        x = [xi+(h1 * fi) for (xi,fi) in zip(x,f)]

        # update w
        f = []
        for i in range(len(w)):
            f.append(F[2](x0,w[i],x[i]))
        w = [wi+(h2 * fi) for (wi,fi) in zip(w,f)]

        # update X and W
        x_sum = sum(x)
        w_sum = sum(w)

        X = [xi / x_sum if x_sum!= 0.0 else 0.0 for xi in x]
        W = [wi / w_sum if w_sum!= 0.0 else 0.0 for wi in w]

        # record X, W and theta for every 100th step
        if t % 100 == 0:
            X_hist.append(X)
            W_hist.append(W)
            theta_hist.append(theta)

        # increment step
        t = t + 1

    # plot X,W and theta
    plt.title('X(red), W(green), Theta(yellow)')
    plt.plot(X_hist, 'r-')
    plt.plot(W_hist, 'g--')
    plt.plot(theta_hist, 'y-.')
    plt.savefig(func.__name__+'.jpg')
    plt.show()


if __name__ == '__main__':

    # initial values for x0 ,x, w
    x0 = 0
    x = [0.8,0.6,0.28,0.15]
    w = [0.9,0.7,0.25,0.12]

    initial_values = [x0,x,w]
    F = [update_x0,update_x,update_w]

    # part (a) of assignment
    euler(0.00035,0.1, 10000, A, initial_values, F)

    # part (b) of assignment
    euler(0.00035,0.1, 10000, B, initial_values, F)



