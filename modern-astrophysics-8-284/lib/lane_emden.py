import numpy as np
from scipy.integrate import solve_ivp


 
# solve the Lane-Emden equation for a polytrope of index n
# outpus is as returned from scipy.integrate.solve_ivp
def lane_emden(n):
    # the rhs of the equation, 
    # normalized to a system of first-order DEs 
    def lane_emden_rhs(xi, y, n):
        if xi == 0:
            return 0

        theta, theta_prime = y
        theta_2prime = -2/xi*theta_prime - theta**n
        return theta_prime, theta_2prime

    # the range to search
    a,b=1e-8,10
    t=np.arange(a,b,np.abs(a-b)/10000.0)

    # initial values for theta, theta_prime
    y0=[1, 0]

    # a zero-crossing event to terminate the integration at the stellar surface
    def event_zero_crossing(xi, y, n):
        return y[0]  # θ = y[0]
    event_zero_crossing.terminal = True
    event_zero_crossing.direction = -1  # only detect downward zero crossing

    # run the solver
    sol = solve_ivp(lane_emden_rhs,(a,b), y0, args=(n,), t_eval=t, events=event_zero_crossing)

    return sol


