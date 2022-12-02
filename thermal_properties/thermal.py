"""
"""
import numpy as np
from matplotlib import pyplot as plt

lar_sigma = 0.126
lar_alpha = 161
lar_beta = 1.12
lar_rho = 1.42
lar_freeze = 83.8
lar_boil = 87.3

foam_sigma = 0.03

steel_sigma = 14.3

kapton_sigma = 0.3
nitrogen_sigma = 0.024

air_gamma = 100

class ArtieI:
    def __init__(self,
        inner_r,
        outer_r,
        foam_r,
        length
    ):
        self.inner_r = inner_r
        self.outer_r = outer_r
        self.foam_r = foam_r
        self.length = length

        self.inner_area = 2 * np.pi * self.inner_r * self.length
        self.outer_area = 2 * np.pi * self.outer_r * self.length
        self.foam_area = 2 * np.pi * (self.outer_r + self.foam_r) * self.length

        self.lar_volume = np.pi * self.inner_r * self.inner_r * self.length
        self.steel_volume = np.pi * (self.outer_r - self.inner_r) * (self.outer_r - self.inner_r) * length
        self.foam_volume = np.pi * self.foam_r * self.foam_r * self.length

        self.R_pipe = R_cylinder(steel_sigma, self.length, self.inner_r, self.outer_r)
        self.R_foam = R_cylinder(foam_sigma, self.length, self.outer_r, self.outer_r + self.foam_r)
        self.R_conv_air = R_conv(100, self.foam_area)

        self.R_total = self.R_pipe + self.R_foam + self.R_conv_air
    
    def heat_rate(self,
        T_in,
        T_out
    ):
        return (T_out - T_in) / self.R_total
    
    def p_over_v(self,
        T_in,
        T_out
    ):
        return self.heat_rate(T_in, T_out) / (self.steel_volume + self.foam_volume)

def R_conv(
    k:  float, 
    A:  float
):
    return 1.0 / (k * A)

def R_wall(
    sigma:  float,
    thickness:  float,
    A:  float
):
    return thickness / (sigma * A)

def R_cylinder(
    sigma:  float,
    length:  float,
    inner_r:    float,
    outer_r:    float
):
    return (np.log(outer_r/inner_r) / (2 * np.pi * sigma * length))

def tau_alpha(
    rho:    float,
    alpha:  float,
    p_over_v: float
):
    return rho * alpha / p_over_v

def tau_beta(
    rho:    float,
    beta:   float,
    T_rise: float,
    p_over_v: float
):
    return rho * beta * T_rise / p_over_v

def p_over_v_cylinder(
    sigma:  float,
    deltaT: float,
    inner_r:    float,
    outer_r:    float
):
    return (2 * sigma * deltaT) / (inner_r * inner_r * np.log(outer_r / inner_r))

def p_over_v_face(
    sigma:  float,
    deltaT: float,
    width:  float,
    length: float
):
    return (sigma * deltaT) / (width * length)

def p_over_v_convection(
    gamma:  float,
    deltaT: float,
    length: float
):
    return gamma * deltaT / length

