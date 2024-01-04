import math

p0 = 0, 0
p1 = 0, 1
p2 = 1, 1.2
p3 = 1.4, 0

l0 = (complex(*p0), complex(*p1))
l1 = (complex(*p2), complex(*p3))

connection = l1[0] - l0[1]
d0 = l0[1] - l0[0]
d1 = l1[1] - l1[0]

def dot(a, b):
    return a.real * b.real + a.imag * b.imag

def cross(a, b):
    return a.real * b.imag - a.imag * b.real

def angle(v0, v1):
    return math.acos(dot(v0, v1) / (abs(v0) * abs(v1))) * (1 if cross(v0, v1) >= 0 else -1)

def lines_intersection(l0, l1):
    print("Lines intersection:", l0, l1)
    p0, p1 = l0
    p2, p3 = l1
    s1 = p1 - p0
    s2 = p3 - p2
    try:
        if abs(cross(s1, s2)) < 1e-6:
            return None
        s = (-s1.imag * (p0.real - p2.real) + s1.real * (p0.imag - p2.imag)) / (-s2.real * s1.imag + s1.real * s2.imag)
        t = ( s2.real * (p0.imag - p2.imag) - s2.imag * (p0.real - p2.real)) / (-s2.real * s1.imag + s1.real * s2.imag)
    except ZeroDivisionError:
        return None
    return p0 + (t * s1)

# Find the angle between the tangent lines and connection
sign = 1 if cross(d0, connection) >= 0 else -1
print("Sign:", sign)
a0 = -angle(d0, connection)
a1 = angle(d1, connection)
print("Angles:", math.degrees(a0), math.degrees(a1))
# Find new lines
a0 *= 0.5
a1 *= 0.5

# Rotate d0, d1 by a0, a1
d0r = d0 * complex(-math.cos(a0), math.sin(a0))
d1r = d1 * complex(math.cos(a1), math.sin(a1))

# Find the intersection of the new lines
arcs_meeting_point = lines_intersection((l0[1], l0[1] + d0r), (l1[0], l1[0] + d1r))
if arcs_meeting_point is None:
    # Parallel lines
    print("Parallel lines")
    arcs_meeting_point = (l0[1] + l1[0]) * 0.5
print("Arcs meeting point:", arcs_meeting_point)

# Find the centers of the arcs
c0 = lines_intersection((l0[1], l0[1] + d0 * complex(0, 1)), (arcs_meeting_point, arcs_meeting_point + connection * complex(0, 1)))
c1 = lines_intersection((l1[0], l1[0] + d1 * complex(0, 1)), (arcs_meeting_point, arcs_meeting_point + connection * complex(0, 1)))
print("Arc centers: ", c0, c1)
if c0 is None or c1 is None:
    # Single-arc solution
    c0 = c1 = l0[1] + d0r

single_arc = False
if abs(c0 - c1) < 1e-6:
    print("Single arc")
    single_arc = True

# Find the radius of the arcs
r0 = abs(c0 - arcs_meeting_point)
r1 = abs(c1 - arcs_meeting_point)
assert abs(r0 - abs(c0 - l0[1])) < 1e-6, "r0: %f, abs(c0 - l0[1]): %f" % (r0, abs(c0 - l0[1]))
assert abs(r1 - abs(c1 - l1[0])) < 1e-6, "r1: %f, abs(c1 - l1[0]): %f" % (r1, abs(c1 - l1[0]))

# Plot them
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
plt.plot(*zip(p0, p1), color='black')
plt.plot(*zip(p2, p3), color='black')
# Draw arcs meeting point
if not single_arc:
    plt.plot([arcs_meeting_point.real], [arcs_meeting_point.imag], marker='o', color='red')
# Draw arc centers
plt.plot([c0.real], [c0.imag], marker='o', color='green')
plt.plot([c1.real], [c1.imag], marker='o', color='green')
# Draw the arcs
theta1 = math.degrees(math.pi * .5 - math.atan2(arcs_meeting_point.real - c0.real, arcs_meeting_point.imag - c0.imag))
theta2 = math.degrees(math.pi * .5 - math.atan2(p1[0] - c0.real, p1[1] - c0.imag))
if sign > 0:
    theta1, theta2 = theta2, theta1
arc0 = Arc((c0.real, c0.imag) , r0 * 2, r0 * 2, angle=0, theta1=theta1, theta2=theta2)
theta1 = math.degrees(math.pi * .5 - math.atan2(p2[0] - c1.real, p2[1] - c1.imag))
theta2 = math.degrees(math.pi * .5 - math.atan2(arcs_meeting_point.real - c1.real, arcs_meeting_point.imag - c1.imag))
if sign > 0:
    theta1, theta2 = theta2, theta1
arc1 = Arc((c1.real, c1.imag) , r1 * 2, r1 * 2, angle=0, theta1=theta1, theta2=theta2)
plt.gca().add_patch(arc0)
plt.gca().add_patch(arc1)
# Aspect ratio 1
plt.gca().set_aspect('equal', adjustable='box')

plt.show()

