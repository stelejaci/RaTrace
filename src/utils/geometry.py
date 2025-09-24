import numpy as np
from numba import jit
from utils.varia import X, Y, EPSILON
from utils.varia import sort_x_left_to_right




def normalize(v):
    if np.ndim(v) == 1:
        return v/np.linalg.norm(v)
    else:
        norms = np.linalg.norm(v, axis=1, keepdims=True)
        norms[norms == 0] = 1
        return v/norms


def normal_from_orientation(r):
    r = np.array([r[0], r[1], 0])
    e = np.array([0, 0, -1]) # Make a temporary vector pointing into the screen, this way normals are correctly described when points are defined in a CLOCKWISE orientation
    n = np.array(np.cross(r, e))
    n = normalize(n)
    n = np.array([n[0], n[1]])
    return n

def orientation_from_normal(n):
    n = np.array([n[0], n[1], 0])
    e = np.array([0, 0, 1]) # Make a temporary vector pointing into the screen, this way normals are correctly described when points are defined in a CLOCKWISE orientation
    r = np.array(np.cross(n, e))
    r = normalize(r)
    r = np.array([r[0], r[1]])
    return r


def direction_between_two_points(pt1, pt2):
    r = [ pt2[X]-pt1[X], pt2[Y]-pt1[Y] ]
    r = normalize(r)
    return r

def point_is_on_PP_line(P, P0, P1):
    v = P1 - P0
    w = P - P0
    cross = np.cross(v,w)
    return np.isclose(cross,0)

def slope_from_vector(v):
    rico = v[Y] / v[X]
    return rico


def direction_of_slope(rico):
    direction_x = np.sqrt(1/(1 +   np.power(rico, 2)))
    direction_y = np.sqrt(1/(1 + 1/np.power(rico, 2)))
    return [direction_x, direction_y]


def points_from_position_direction_length(p0, r, L, sort_left_to_right=False, symmetric=False):
    r = normalize(r)
    if symmetric:
        pts = np.array([ p0-L/2*r, p0+L/2*r ])
    else:
        pts = np.array([p0, p0 + L*r])
    if sort_left_to_right:
        pts = sort_x_left_to_right(pts)
    return pts


def point_from_position_direction_length(p, r, l):
    r = normalize(r)
    return p + l * r


def normals_from_points(pts):
    nr_of_points = len(pts[:, X])
    n = np.empty([nr_of_points, 2])

    # First point
    n[0, X] = -(pts[1, Y] - pts[0, Y])
    n[0, Y] = (pts[1, X] - pts[0, X])
    n[0, :] = normalize(n[0, :])

    # Inbetween points
    for i_pt in range(1, nr_of_points - 1):
        n[i_pt, X] = -(pts[i_pt + 1, Y] - pts[i_pt - 1, Y])
        n[i_pt, Y] = (pts[i_pt + 1, X] - pts[i_pt - 1, X])
        n[i_pt, :] = normalize(n[i_pt, :])

    # Last point
    n[-1, X] = -(pts[-1, Y] - pts[-2, Y])
    n[-1, Y] = (pts[-1, X] - pts[-2, X])
    n[-1, :] = normalize(n[-1, :])

    return n


def normal_at_point(pt, pts, n):
    nx = np.interp(pt[X], pts[:, X], n[:, X])
    ny = np.interp(pt[X], pts[:, X], n[:, Y])
    n = np.array([nx, ny])
    n = normalize(n)
    return n


def distance_between_2_points(pt1, pt2):
    return np.linalg.norm(pt2-pt1)


def intersection_of_PP_line_with_PR_line(p00, p01, p10, r1):
    r0 = np.array([p01[X]-p00[X], p01[Y]-p00[Y]])
    [p, t0, t1] = intersection_of_PR_line_with_PR_line(p00, r0, p10, r1)
    return [p, t0, t1]


def intersection_of_PR_line_with_PP_line(p00, r0, p10, p11):
    r1 = np.array([p11[X]-p10[X], p11[Y]-p10[Y]])
    [p, t0, t1] = intersection_of_PR_line_with_PR_line(p00, r0, p10, r1)
    return [p, t0, t1]


def intersection_of_PR_line_with_PR_line(pt1, r1, pt2, r2):
    # start_time = time.perf_counter()
    [ptX, ptY, t1, t2] = subcalc_intersection_of_PR_line_with_PR_line(pt1, r1, pt2, r2)
    # print(f'{time.perf_counter()-start_time:.8f}')
    return [np.array([ptX, ptY]), t1, t2]


@jit(nopython = True)   # Remove this for debugging, but then remove the print statement further on
def subcalc_intersection_of_PR_line_with_PR_line(pt1, r1, pt2, r2):
    dx = pt2[0] - pt1[0]
    dy = pt2[1] - pt1[1]

    # print(f'{pt1}  {r1}  {pt2}  {r2}  {dx}  {dy}')

    if np.array_equal(r1, [0,0]):
        t2 = 9999999
    else:
        divisor = r2[1] * r1[0] - r1[1] * r2[0]
        if np.abs(divisor) < EPSILON:
            t2 = 9999999
        else:
            t2 = ( r1[1]*dx - r1[0]*dy ) / divisor

    x2 = pt2[0] + r2[0] * t2
    y2 = pt2[1] + r2[1] * t2

    if np.array_equal(r1, [0, 0]):
        t1 = 999999
    elif np.abs(r1[0]) < EPSILON:
        t1 = (y2 - pt1[1]) / r1[1]      # Intersection with vertical lines
    else:
        t1 = (x2 - pt1[0]) / r1[0]      # Intersection with normal, non-vertical lines
    # Just a check
    # x1 = pt1[X] + r1[X] * t1
    # y1 = pt1[Y] + r1[Y] * t1
    # pt = [x2, y2]
    # pt = np.array([x2, y2])
    # return [pt, t1, t2]
    return [x2, y2, t1, t2]


def intersection_of_PP_line_with_PP_line(p00, p01, p10, p11):
    r1 = [p01[X]-p00[X], p01[Y]-p00[Y]]
    r2 = [p11[X]-p10[X], p11[Y]-p10[Y]]
    [p, t0, t1] = intersection_of_PR_line_with_PR_line(p00, r1, p10, r2)
    return [p, t0, t1]


def intersection_and_normal_of_PR_line_with_2D_mesh(pt1, pts, n):
    # For now, the direction is ignored, and is just considered downward
    x = pt1[X]
    y = np.interp(x, pts[:,X], pts[:,Y])
    nx = np.interp(x, pts[:,X], n[:,X])
    ny = np.interp(x, pts[:,X], n[:,Y])
    n = (np.array([nx, ny]))
    return [np.array([x,y]), n]

def intersection_of_PR_line_with_PPN_line(P, R, P0, P1, N):
    if np.dot(R,N) > 0: # The normal of the segment (N) is parallel with the line (ray orientation R), and thus cannot physically intersect in a closed shape
        return [None, None, None, None]

    [p_coll, t0, t1] = intersection_of_PR_line_with_PP_line(P, R, P0, P1)

    if (t1 <= 0)  or  (t1 >= 1):    # The intersection is outside of the segment
        return [None, None, None, None]

    return [p_coll, t0, t1, N]

def line_arc_intersections(C_arc, p_ends, p_line, r_line):
    x0, y0 = C_arc[0], C_arc[1]    # Point at the centre of the arc circle
    # a0, a1 = angles_arc[0], angles_arc[1]   # Angles of the 2 endpoints of the arc
    x1, y1 = p_line[0], p_line[1]   # Point on the line
    rx, ry = r_line[0], r_line[1]   # Direction of the line
    R_arc = np.linalg.norm(p_ends[0,:]-C_arc)

    # # Normalize direction vector
    # norm = np.sqrt(rx**2 + ry**2)
    # if norm == 0:
    #     return [None, None]
    # rx /= norm
    # ry /= norm

    # Proceed with full circle-line intersection check
    dx = x1 - x0
    dy = y1 - y0
    A = rx**2 + ry**2
    B = 2 * (rx * dx + ry * dy)
    C = dx**2 + dy**2 - R_arc**2

    discriminant = B**2 - 4 * A * C
    if discriminant < 0:    # No collision
        return [None, None]

    # At this point there is at least one intersection point
    sqrt_disc = np.sqrt(discriminant)
    t1 = (-B - sqrt_disc) / (2 * A)
    t2 = (-B + sqrt_disc) / (2 * A)

    x_MAEP = (p_ends[0,0]+p_ends[1,0])/2   # Mean Arc End Points --> The midpoint between the 2 endpoints of the arc.
    y_MAEP = (p_ends[0,1]+p_ends[1,1])/2

    pts, t = [], []
    for ti in [t1, t2]:
        if ti >= EPSILON:  # Only forward direction
            x = x1 + ti * rx
            y = y1 + ti * ry

            v1 = normalize( np.array( [C_arc[0]-x_MAEP, C_arc[1]-y_MAEP] ) )    # Vector from MAEP to centre of circle
            v2 = normalize( np.array( [       x-x_MAEP,        y-y_MAEP] ) )    # Vector from MAEP to intersection point
            angle = angle_between_vectors(v1, v2)                               # The angle inbetween the 2 vectors should be larger than 90°

            if -np.pi/2 <= angle <= np.pi/2:  # Intersecting, but outside the arc, on the wrong side of the MAEP
                pass
            else:                           # Intersecting and inside the arc, on the good side of the MAEP
                pts.append([x, y])
                t.append(ti)

    if len(t) == 0:
        return [None, None]               # No intersections in the forward direction
    elif len(t) > 1  and  t[1] <= t[1]:
        return [np.array(pts[1]), t[1]]   # There are 2 intersection points, and the second is the closest
    else:
        return [np.array(pts[0]), t[0]]   # If only one point or the first point is the closest

def line_parabola_intersections(f, p0, n0, D, p_line, r_line):
    # First, transform the line back to the reference frame of the [0,1] upward facing parabola through [0,0]
    angle = angle_between_vectors(np.array([0,1]), n0)
    r_line = rotate_with_P_angle(pts=r_line, angle=angle, p_rot=p0, is_orientation=True)
    p_line = rotate_with_P_angle(pts=p_line, angle=angle, p_rot=p0, is_orientation=False)
    p_line = move_by(p_line, -p0)

    r_x, r_y = r_line
    p_x, p_y = p_line
    intersections = []

    if np.isclose(r_x, 0):  # A vertical ray
        intersections.append((p_x, (p_x**2)/(4*f)))
    elif np.isclose(r_y, 0):    # A horizontal ray.
        if (f > 0 and p_y >= 0) or (f < 0 and p_y <= 0):        # No intersection if line is outside parabola's range
            intersections.extend([(np.sqrt(4*f*p_y), p_y), (-np.sqrt(4*f*p_y), p_y)])
    else:   # The general case, solve quadratic equation for t
        A = r_x**2
        B = 2*p_x*r_x - 4*f*r_y
        C = p_x**2 - 4*f*p_y
        discriminant  = B**2 - 4*A*C

        if D >= 0:
            sqrt_discriminant = np.sqrt(discriminant)
            t1 = (-B + sqrt_discriminant) / (2*A)
            t2 = (-B - sqrt_discriminant) / (2*A)
            x1, y1 = p_x+r_x*t1, p_y+r_y*t1
            x2, y2 = p_x+r_x*t2, p_y+r_y*t2
            intersections.extend([(x1, y1), (x2, y2)])

    # Process all found intersections
    for x, y in intersections:
        if np.abs(x) > D / 2:  continue         # Check if the intersection point lies within the boundary of the parabola

        n_coll = normalize(  np.array([-x/(2*f), 1.0])  )
        if np.dot(r_line,n_coll) > 0:  continue      # The normal direction should run antiparallel to the line direction

        if np.isclose(r_x, 0):  t0 = (y - p_y) / r_y    # Vertical line
        else:                   t0 = (x - p_x) / r_x

        t1 = (x + D / 2) / D

        n_coll = rotate_with_P_angle(pts=n_coll, angle=-angle, p_rot=p0, is_orientation=True)
        p_coll = move_by( np.array([x,y]), p0)
        p_coll = rotate_with_P_angle(pts=p_coll, angle=-angle, p_rot=p0, is_orientation=False)

        return [p_coll, t0, t1, n_coll]

    return [None, None, None, None]

def interpolate_pt_and_normal_on_2D_mesh(x, pts, n):
    # From an X-coordinate, interpolate the Y-coordinate on a 2D mesh, as well as the normal
    y = np.interp(x, pts[:,X], pts[:,Y])
    nx = np.interp(x, pts[:,X], n[:,X])
    ny = np.interp(x, pts[:,X], n[:,Y])
    n = normalize(np.array([nx, ny]))
    return [np.array([x,y]), n]


def angle_from_vector(r):
    return np.arctan2(r[Y], r[X])


def angle_between_vectors(n,r):
    # Returns the signed angle in radians from vector r to vector n. Positive angle indicates counter-clockwise rotation from r to n.
    perp_dot = r[0] * n[1] - r[1] * n[0]        # Perpendicular dot product (determinant) = ||r|| * ||n|| * sin(θ)
    dot      = r[0] * n[0] + r[1] * n[1]        # Regular dot product = ||r|| * ||n|| * cos(θ)
    return np.arctan2(perp_dot, dot)            # arctan2(sin(θ), cos(θ)) gives θ ∈ [-π, π]
    # return math.acos(min(np.dot(normalize(n),normalize(r)),1))


def rotate_direction_over_angle(direction, angle):
    M = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    direction_new = np.matmul(M, direction)
    direction_new = normalize(direction_new)
    return direction_new


def vector_from_angle(angle):
    if np.ndim(angle) == 0:
        return np.array([np.cos(angle), np.sin(angle)])
    else:
        return np.array([[np.cos(a), np.sin(a)] for a in angle])

def define_bounding_box(pts):
    xm, xM = np.min(pts[:,X]), np.max(pts[:,X])
    ym, yM = np.min(pts[:,Y]), np.max(pts[:,Y])
    pts_BB = np.empty((0, 2))
    pts_BB = np.append(pts_BB, [np.array([xM,yM])], axis=0)
    pts_BB = np.append(pts_BB, [np.array([xM,ym])], axis=0)
    pts_BB = np.append(pts_BB, [np.array([xm,ym])], axis=0)
    pts_BB = np.append(pts_BB, [np.array([xm,yM])], axis=0)
    return pts_BB


def rotation_matrix_from_angle(angle):
    R = np.array([[np.cos(-angle), np.sin(-angle)], [-np.sin(-angle), np.cos(-angle)]])
    return R


def rotate_with_P_angle(pts, angle=0, p_rot=np.array([0,0]), is_orientation=False):
    R = rotation_matrix_from_angle(angle)
    pts = rotate_with_PR(pts, p_rot, R, is_orientation)
    return pts


def rotate_with_PR(pts, p_rot, R, is_orientation=False):
    if is_orientation:
        pts = np.transpose(np.matmul(R, np.transpose(pts)))
        pts = normalize(pts)
    else:
        pts = move_by(pts, -p_rot)
        pts = np.transpose(np.matmul(R, np.transpose(pts)))
        pts = move_by(pts, p_rot)
    return pts


def move_by(pts, dp):
    return pts + dp


def move_to(pts, p0, p1):
    dp = p1 - p0
    pts_new = move_by(pts, dp)
    return pts_new


def construct_plate(p0, n, thickness, length):
    n = normalize(n)
    r = orientation_from_normal(n)
    pts_corners = points_from_position_direction_length(p0=p0, r=r, L=length, sort_left_to_right=False, symmetric=True)
    pts = np.empty((0,2))
    pts = np.append(pts, [pts_corners[0]], axis=0)
    pts = np.append(pts, [pts_corners[1]], axis=0)
    pts = np.append(pts, [pts_corners[1]-n*thickness], axis=0)
    pts = np.append(pts, [pts_corners[0]-n*thickness], axis=0)
    return pts




if __name__ == '__main__':
    # pts = np.array([[0,-1],[0,1],[1,1],[1,-1]])
    # print(pts)
    # print(move_by(pts, np.array([1,0])))
    # print(pts)
    # print(move_to(pts, np.array([0,0]), np.array([0,1])))

    # pts = np.array([0,1])
    # print(pts)
    # print(move_by(pts, np.array([1,0])))
    # print(pts)
    # print(move_to(pts, np.array([0,0]), np.array([0,1])))

    # pts = construct_plate(p0=np.array([1,1]), n=np.array([-2,-1]), thickness=1*mm, length=10*mm)

    # v = np.array([[1,1],[0,1],[-1,2]])
    # print(v)
    # print(normalize(v))
    #
    # v = np.array([-1,2])
    # print(v)
    # print(normalize(v))

    print(vector_from_angle(0.7))
    print(vector_from_angle([0.7]))
    print(vector_from_angle([0.7,0.1,1.2]))

