import matplotlib.pyplot as plt
from sympy import Point, Polygon, Line, Ray
from sympy.geometry import intersection
from sympy.geometry.entity import GeometrySet
from operator import itemgetter
from numpy import meshgrid
import math
import numpy


g= 9.81

coordinates = (
    (0, 0),
    (0,-10),
    (60, -10),
    (60, 0),
)


def angle(p1, p2):
    k = (p2[1] - p1[1]) / Point(p1).distance(Point(p2))
    x2 = p2[0]
    x1 = p1[0]

    if k >= 0:
        if x2 >= x1: # First Quadrant
            return (2.0 * math.pi - math.asin(k))
        else: # Second Quadrant
            return (math.pi + math.asin(k))
    else:
        if x2 >= x1: # Fourth Quadrant
            return math.asin(-k)
        else: # Third Quadrant
            return (math.pi - math.asin(-k))


class ChannelSection:
    """A cross section taken normal to the direction of flow"""

    def __init__(self, coordinates):
        self.abscissas = []
        self.ordinates = []
        for i in range(len(coordinates)):
            self.abscissas.append(coordinates[i][0])
            self.ordinates.append(coordinates[i][1])
        self.abscissas.append(coordinates[0][0])
        self.ordinates.append(coordinates[0][1])

    def lowest(self):
        return sorted(coordinates, key=itemgetter(1))[0]

    def plot_geometry(self):
        plt.plot(self.abscissas, self.ordinates)
        plt.title('Channel Section', fontsize=14, color='black')
        plt.show()


class FlowArea(ChannelSection):
    """
    Cross-section in which there is flow, has a height y(depth), Flow(QQ)

    """

    def __init__(self, ChannelSection, depth, QQ):
        self.ChannelSection = ChannelSection
        self.depth = depth
        self.QQ = QQ

    def flow_area_polygon(self):
        point_on_surface_above_lowest_depth = Point(self.ChannelSection.lowest()[0], self.ChannelSection.lowest()[1]+self.depth)
        line = Line(point_on_surface_above_lowest_depth, Point(100, self.ChannelSection.lowest()[1]+self.depth))

        i1, i2 = intersection(line, Polygon(*coordinates))

        i1 = (i1[0], i1[1])
        i2 = (i2[0], i2[1])
        flow_coordinates = []
        for p in coordinates:
            if p[1] < self.ChannelSection.lowest()[1]+self.depth:
                flow_coordinates.append(p)

        if i1  not in flow_coordinates:
            flow_coordinates.append(i1)
        if i2 not in flow_coordinates:
            flow_coordinates.append(i2)

        def centroid(*points):
            x_coords = [p[0] for p in points]
            y_coords = [p[1] for p in points]
            _len = len(points)
            centroid_x = sum(x_coords) / _len
            centroid_y = sum(y_coords) / _len
            return [centroid_x, centroid_y]
        return sorted(flow_coordinates, key = lambda point: -angle(point, centroid(*flow_coordinates)))

    def free_surface_length(self):
        """
        B
        :return:
        """
        point_on_surface_above_lowest_depth = Point(self.ChannelSection.lowest()[0], self.ChannelSection.lowest()[1]+self.depth)
        line = Line(point_on_surface_above_lowest_depth, Point(0.5464321348512123132, self.ChannelSection.lowest()[1]+self.depth))

        i1, i2 = intersection(line, Polygon(*coordinates))
        return i1.distance(i2)

    def wetted_perimeter(self):
        """
        P
        :return:
        """
        flow_coordinates = []
        for p in coordinates:
            if p[1] > self.ChannelSection.lowest()[1] + self.depth:
                pass
            else:
                flow_coordinates.append(p)
        perimeter = 0
        surface_points = []
        surface_points = list( map(Point, flow_coordinates))
        for i in range(0,len(surface_points)-1):
            perimeter += surface_points[i].distance(surface_points[i+1])
        return perimeter

    def hydraulic_radius(self):
        return self.area()/self.wetted_perimeter()

    def hydraulic_depth(self):
        return self.area()/self.free_surface_length()

    def get_centroid(self):
        points = self.flow_area_polygon()
        return Polygon(*points).centroid

    def area(self):
        points = self.flow_area_polygon()
        return Polygon(*points).area

    def plot_geometry(self):
        self.abscissas = []
        self.ordinates = []
        points = self.flow_area_polygon()
        for i in range(len(points)):
            self.abscissas.append(points[i][0])
            self.ordinates.append(points[i][1])
        self.abscissas.append(points[0][0])
        self.ordinates.append(points[0][1])
        plt.plot(self.abscissas, self.ordinates)
        plt.title('Channel Section with flow depth', fontsize=14, color='black')
        plt.show()

    def specific_force(self):
        """
        momentum function or specific force,
                Q 2
        Fs = -------- + z A
                gA

        :return:  Fs
        """
        specific_force = self.QQ*self.QQ/(g* self.area()) + (self.get_centroid()[1]-self.lowest()[1])*self.area()
        return specific_force

    def velocity(self):
        return self.QQ / self.area()

    def froude_number(self):
        return self.velocity()/math.sqrt(g*self.depth)

    def energy(self):
        return self.depth + self.velocity()*self.velocity()/(2*g)

    def hydraulic_jump_y2(self):
        """
        Conservation of momentum flux
        :return:
        """
        y21= self.depth +(-1 + math.sqrt(1+ 8*self.froude_number()))/2
        y22= self.depth +(-1 - math.sqrt(1+ 8*self.froude_number()))/2         # Negative answers do not yield meaningful physical solutions

        return y21

    def discharge_per_unit_length(self):
        return self.QQ/self.free_surface_length()

    def specific_energy_plot(self):
        e = numpy.linspace(0,20.0, 100)
        delta = 0.125
        xrange = numpy.arange(0.0, 20.0, delta)
        yrange = numpy.arange(0.0, 20.0, delta)
        X, Y = meshgrid(xrange, yrange)

        F = (X - Y)*Y*Y
        G = (self.discharge_per_unit_length()*self.discharge_per_unit_length())/(2*g)
        plt.plot(e, e, 'b--')
        plt.contour(X, Y, (F - G), [0])
        plt.title('Specific Energy E[m]', fontsize=14, color='black')
        plt.xlabel('Specific Energy', fontsize=14, color='black')
        plt.ylabel('Flow depth, y[m]', fontsize=14, color='black')
        plt.show()


A = ChannelSection(coordinates)
A.plot_geometry()

Af = FlowArea(A, 5.55, 2500)
Af.plot_geometry()

print("Centroid :", Af.get_centroid())

print("Specific Force :", Af.specific_force())

print("velocity : ", Af.velocity())

print("Top Width: ", Af.free_surface_length())
print("Hydraulic Radius : ", Af.hydraulic_radius())
print("Hydraulic Depth : ", Af.hydraulic_depth())
print("Hydraulic Jump y2 : ", Af.hydraulic_jump_y2())

print("Area cross-section flow:",  round(Af.area(), 2))
print("Froud number :", Af.froude_number())
Af.specific_energy_plot()