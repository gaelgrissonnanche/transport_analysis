from numpy import vectorize, piecewise, float64
from numpy.polynomial import Chebyshev
import json

class PiecewisePolynomial():

    def __init__(self):
        self.poly_coeffs_dict = {} # keys="x_seg_min, x_seg_max", values=coeff_list over segment
        self.poly_coeffs_list = [] # values=coeff_list over segment
        self.segments_list = []
        self.polynomial_class = Chebyshev # can be the class Chebyshev or Polynomial (Power series)


    def fit(self, x_data, y_data, x_breakpoints, degree):
        # Fit a Chebyshev polynomial to each segment and store them
        self.poly_coeffs_dict = {}
        self.poly_coeffs_list = []
        self.segments_list = []
        for i in range(len(x_breakpoints)-1):
            # Isolate the segment
            x_seg_min = x_breakpoints[i]
            x_seg_max = x_breakpoints[i+1]
            x_seg = x_data[(x_data >= x_seg_min) & (x_data < x_seg_max)]
            y_seg = y_data[(x_data >= x_seg_min) & (x_data < x_seg_max)]
            # Fit Chebyshev polynomial to this segment
            polynomial_object = self.polynomial_class.fit(x_seg, y_seg, degree, domain=[x_seg_min, x_seg_max])
            # The Chebyshev fit is performed over a domain where it renormalizes
            # x array to perform better. One needs to access the domain later,
            # it is enough to do polynomial_object.domain[0] = x_seg_min
            # polynomial_object.domain[1] = x_seg_max.
            # !!!When using the coefficients, it is important to use the domain!!!
            self.poly_coeffs_dict[str(x_seg_min) + ", " + str(x_seg_max)] = list(polynomial_object.coef)
            self.poly_coeffs_list.append(list(polynomial_object.coef))
            self.segments_list.append([x_seg_min, x_seg_max])


    def save_coefficients(self, filename):
        with open(filename, 'w') as f:
            json.dump(self.poly_coeffs_dict, f, indent=4)


    def load_coefficients(self, filename):
        """Loads the polynomial coefficients of a piecewise polynomial fit on different
        segments of [seg_min, seg_max]"""
        with open(filename, 'r') as f:
            self.poly_coeffs_dict = json.load(f)
        ## Build the segments of each polynomial fits
        self.poly_coeffs_list = []
        self.segments_list = []
        for segments_str, coeffs in self.poly_coeffs_dict.items():
            segments_min, segments_max = segments_str.split(',')
            self.segments_list.append([float(segments_min), float(segments_max)])
            self.poly_coeffs_list.append(coeffs)


    def generate_piecewise_polynomial_function(self):
        return lambda x: piecewise(float64(x), [(self.segments_list[i][0] <= x) & (x < self.segments_list[i][1]) for i in range(len(self.segments_list))],
            [self.polynomial_class(self.poly_coeffs_list[i], domain=[self.segments_list[i][0], self.segments_list[i][1]]) for i in range(len(self.segments_list))])


    # def generate_piecewise_polynomial_function2(self):
    #     """
    #     Generates a piecewise polynomial function on segments of segments_list[i]=[min, max].
    #     The segments_list is a list of [[min, max]].
    #     The coeffs_list represents the polynomial coefficients for each [min, max] segments.
    #     """
    #     def piecewise_func(x):
    #         for i in range(len(self.segments_list)):
    #             seg_min = self.segments_list[i][0]
    #             seg_max = self.segments_list[i][1]
    #             if seg_min <= x < seg_max:
    #                 return self.polynomial_class(self.poly_coeffs_list[i], domain=[seg_min, seg_max])(x)
    #     self.piecewise_polynomial_func = vectorize(piecewise_func)
    #     return self.piecewise_polynomial_func


