import math


def cross(a, b):
    return a.real * b.imag - a.imag * b.real


def lines_intersection(l0, l1):
    p0, p1 = l0
    p2, p3 = l1
    s1 = p1 - p0
    s2 = p3 - p2
    c = cross(s1, s2)
    if abs(c) < 1e-6:
        return None
    t = cross(s2, p0 - p2) / c
    return p0 + (t * s1)


def lines_2arc_connection(l0, l1):
    connection = (l1[0] - l0[1]) / abs(l1[0] - l0[1])
    d0 = (l0[1] - l0[0]) / abs(l0[1] - l0[0])
    d1 = (l1[1] - l1[0]) / abs(l1[1] - l1[0])

    # Rotate d0, d1 by a0, a1
    d0r = d0 + connection
    d1r = d1 + connection

    # Find the intersection of the new lines
    arcs_meeting_point = lines_intersection((l0[1], l0[1] + d0r), (l1[0], l1[0] - d1r))
    if arcs_meeting_point is None:
        # Parallel lines
        arcs_meeting_point = (l0[1] + l1[0]) * 0.5

    # Find the centers of the arcs
    c0 = lines_intersection(
        (l0[1], l0[1] + d0 * complex(0, 1)),
        (arcs_meeting_point, arcs_meeting_point + connection * complex(0, 1)),
    )
    c1 = lines_intersection(
        (l1[0], l1[0] + d1 * complex(0, 1)),
        (arcs_meeting_point, arcs_meeting_point + connection * complex(0, 1)),
    )
    if c0 is None or c1 is None:
        # Single-arc solution
        c0 = c1 = l0[1] + d0r

    # Find the radius of the arcs
    r0 = abs(c0 - arcs_meeting_point)
    r1 = abs(c1 - arcs_meeting_point)
    assert abs(r0 - abs(c0 - l0[1])) < 1e-6, (r0, abs(c0 - l0[1]))
    assert abs(r1 - abs(c1 - l1[0])) < 1e-6, (r1, abs(c1 - l1[0]))

    # Find start/end of arcs

    sign = 1 if cross(d0, connection) >= 0 else -1

    theta01 = math.degrees(
        math.pi * 0.5
        - math.atan2(
            arcs_meeting_point.real - c0.real, arcs_meeting_point.imag - c0.imag
        )
    )
    theta02 = math.degrees(math.pi * 0.5 - math.atan2(p1[0] - c0.real, p1[1] - c0.imag))
    if sign > 0:
        theta01, theta02 = theta02, theta01

    theta11 = math.degrees(math.pi * 0.5 - math.atan2(p2[0] - c1.real, p2[1] - c1.imag))
    theta12 = math.degrees(
        math.pi * 0.5
        - math.atan2(
            arcs_meeting_point.real - c1.real, arcs_meeting_point.imag - c1.imag
        )
    )
    if sign > 0:
        theta11, theta12 = theta12, theta11

    return (c0, r0, theta01, theta02), (c1, r1, theta11, theta12), arcs_meeting_point


def plot_lines_2arc_connection(plt, l0, l1, solution):
    from matplotlib.patches import Arc

    (
        (c0, r0, theta01, theta02),
        (c1, r1, theta11, theta12),
        arcs_meeting_point,
    ) = solution

    single_arc = False
    if abs(c0 - c1) < 1e-6:
        single_arc = True

    plt.plot([l0[0].real, l0[1].real], [l0[0].imag, l0[1].imag], color="black")
    plt.plot([l1[0].real, l1[1].real], [l1[0].imag, l1[1].imag], color="black")
    # Draw arcs meeting point
    if not single_arc:
        plt.plot(
            [arcs_meeting_point.real],
            [arcs_meeting_point.imag],
            marker="o",
            color="red",
        )
    # Draw arc centers
    plt.plot([c0.real], [c0.imag], marker="o", color="green")
    plt.plot([c1.real], [c1.imag], marker="o", color="green")
    # Draw the arcs
    arc0 = Arc(
        (c0.real, c0.imag), r0 * 2, r0 * 2, angle=0, theta1=theta01, theta2=theta02
    )
    arc1 = Arc(
        (c1.real, c1.imag), r1 * 2, r1 * 2, angle=0, theta1=theta11, theta2=theta12
    )
    plt.add_patch(arc0)
    plt.add_patch(arc1)


if __name__ == "__main__":
    tests = [
        ((0, 0), (0, 1), (1, 1), (1, 0)),
        ((0, 0), (0, 1), (1, 1.2), (1, 0)),
        ((0, 0), (0, 1), (1, 1.2), (1.4, 0)),
        ((0, 0), (0, 1), (1, 1.2), (0.8, 0)),
        ((0, 0), (0, 1), (1, 1.5), (2, 1.5)),
        ((0, 0), (0, 1), (1, 2), (2, 2)),
        ((0, 0), (0, -1), (1, -1), (1, 0)),
        ((0, 0), (0, -1), (1, -1.2), (1, 0)),
        ((0, 0), (0, -1), (1, -1.2), (1.4, 0)),
        ((0, 0), (0, -1), (1, -1.2), (0.8, 0)),
        ((0, 0), (0, -1), (1, -1.5), (2, -1.5)),
        ((0, 0), (0, -1), (1, -2), (2, -2)),
    ]

    # Plot them
    import matplotlib.pyplot as plt

    w = math.ceil(math.sqrt(len(tests)))
    h = math.ceil(len(tests) / w)
    fig, axs = plt.subplots(h, w)

    for i, test in enumerate(tests):
        axis = axs[i // w, i % w]
        p0, p1, p2, p3 = test
        l0 = (complex(*p0), complex(*p1))
        l1 = (complex(*p2), complex(*p3))
        solution = lines_2arc_connection(l0, l1)
        plot_lines_2arc_connection(axis, l0, l1, solution)
        axis.set_aspect("equal", adjustable="box")

    plt.show()
