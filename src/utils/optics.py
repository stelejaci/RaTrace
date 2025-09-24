import numpy as np
from utils.varia import mm, µm, X

global N_glass, N_air
N_air = 1
N_glass = 1.5




def calculate_image_distance(object_distance, focal_length):
    image_distance = 1 / (1 / focal_length - 1 / object_distance)
    return image_distance


def refract_ray(ray, n_coll, Ni, No, is_active=True, is_visible=True):
    from light import light_class
    # Snell's law in 2D vector format, see:  https://www.slideshare.net/dieulinh299/ray-tracing

    # In the formula, ri, ro and n all point away from the surface
    n = n_coll  if  np.dot(ray.r, n_coll)<0  else  -n_coll  # Assure that n points away from the surface
    ri = -ray.r     # r should point away from the surface

    Nr = Ni/No
    ro = (Nr*np.dot(n,ri) - np.sqrt(1-Nr**2*(1-np.dot(n,ri)**2)))*n - Nr*ri

    ray_new = light_class.RayClass(p0=ray.p1, r=ro, intensity=ray.intensity, wavelength=ray.wavelength, ray_parent=ray, N=No, plot_color=ray.plot_color, is_active=is_active, is_visible=is_visible)
    return ray_new


def refract_ray_on_ideal_lens(p0, n0, f, p, r, p_coll):
    from utils import geometry
    # p0,n0,r0 are from the lens  /  p,r are from the ray
    r0 = geometry.orientation_from_normal(n0)   # Orientation along the lens
    n_prop = n0 * np.sign(np.dot(r, n0))  # n0 is just the main orientation of the lens, for reference of constructing or plotting. However, the ray could come from "left" or "right", so n_prop is the general propagation direction of the ray, ON the optical axis

    if np.isclose( np.linalg.norm(p_coll-p0), 0 ):  # The ray goes through the center of the lens, in that case, the ray remains unaffected
        ro = r
    elif np.isclose( np.dot(r,r0), 0 ):  # Ray is parallel to the optical axis
        v = f  # The outgoing ray goes through the focal point
        p_img = p0 + v * n_prop  # "Image point" on the optical axis
        ro = geometry.normalize(p_img - p_coll)  # The outgoing ray is along the line between the ray-lens intsersection point p_coll and the image point p_img
    else:   # The general case
        t = np.dot(p0 - p, r0) / np.dot(r, r0)  # The ray p+t*r intersects the optical axis when (p+t*r - p0)·r0 = 0. Solving for t gives the intersection point.
        p_obj = p + t * r  # "Object point" on the optical axis

        obj_vec = p0 - p_obj  # The vector between o and p0 along the optical axis
        d_obj = np.dot(obj_vec, n_prop)  # Compute object distance (d_obj) from the object point p_obj to the lens, becomes negative for a virtual image

        if np.isclose(f - d_obj, 0):  # A converging ray whos object point is virtual and aligns with the focal point. The outgoing ray is parallel to the optical axis then.
            ro = n_prop
        else:
            v = 1 / (1 / f - 1 / d_obj)     # The thin lens equation
            p_img = p0 + v * n_prop  # "Image point" on the optical axis
            ro = geometry.normalize(p_img - p_coll)  # The outgoing ray is along the line between the ray-lens intsersection point p_coll and the image point p_img

    if np.dot(n_prop, ro) < 0:
        ro = -ro

    return ro


def derive_lens_radii_and_f(N=N_glass, f=None, R0=None, R1=None, T=None):
    # Using the lensmaker's equation: https://en.wikipedia.org/wiki/Lens#Lensmaker.27s_equation
    N = N[0] if isinstance(N, list) else N
    if f is not None:
        # If only f is given, consider a symmetric lens and thus R0 and R1 are equal
        # Reformulate the lensmaker's equation with R0=-R1=R and solve the quadratic equation for R
        A = 1
        B = -f*(N-1)*2
        C = f*(N-1)**2*T/N
        D = B**2-4*A*C
        if f>0:
            R = (-B + np.sqrt(D))/(2*A)   # The "-" solution also checks out somehow, but is not valid
        else:
            R = (-B - np.sqrt(D))/(2*A)  # The "-" solution also checks out somehow, but is not valid
        R0 = R
        R1 = -R
        f_check = lensmakers_equation(N, R0, R1, T)
        print(f'Lens parameters derived: f={f:0.2f}mm --> R0={R0:0.2f}mm, R1={R1:0.2f}mm --> f_check={f_check:0.2f}mm')
    elif R0 is not None and R1 is not None:
        f = lensmakers_equation(N, R0, R1, T)
        print(f'Lens parameters derived: R0={R0}, R1={R1} --> f={f}')
    else:
        print(f'Error: Lens definition is incorrect, either define f OR R0 and R1: f={f}, R0={R0}, R1={R1}')
        return [None, None, None]

    return [f, R0, R1]


def lensmakers_equation(N,R0,R1,T):
    # Using the lensmaker's equation: https://en.wikipedia.org/wiki/Lens#Lensmaker.27s_equation
    one_over_f = (N - 1) * (1 / R0 - 1 / R1 + (N - 1) * T / (N * R0 * R1))
    f = 1 / one_over_f
    return f


def calculate_refraction_index(N, wavelength):
    # Reference: https://en.wikipedia.org/wiki/Cauchy%27s_equation
    # Simulation units are in mm, while B-coefficients  are in µm²
    wavelength_mu = wavelength/µm*mm
    if isinstance(N,list):
        A, B = N[0], N[1]
    else:
        A, B = N, 0
    return A + B/wavelength_mu**2

# Lens with its first surface
def construct_lens(p0, R0, R1, T, D, resolution):
    p1 = np.array([T,  0])          # Principle point on the second surface
    C0 = np.array([R0, 0])          # Centre of the circle of the first surface
    C1 = np.array([p1[X]+R1, 0])

    if not isinstance(D,list):
        D = [D, D]  # Use the same diameter for both lens surfaces, otherwise, use the ones that are given

    y0 = np.linspace(-D[0]/2, D[0]/2, 1 + int(D[0]/resolution))
    y1 = np.linspace(D[1]/2, -D[1]/2, 1 + int(D[1]/resolution))
    x0 = C0[0] - np.sign(R0) * np.sqrt(R0**2-y0**2)
    x1 = C1[0] - np.sign(R1) * np.sqrt(R1 ** 2 - y1 ** 2)
    p_corners = np.array([[x0[0], y0[0]], [x0[-1], y0[-1]], [x1[0], y1[0]], [x1[-1], y1[-1]]])
    pts = np.concatenate([x0, x1, y0, y1]).reshape(2, -1).T  # Make it an Nx2 array

    return [pts + p0, p1 + p0, C0 + p0, C1 + p0, p_corners + p0]



# Parabolic mirror with its first surface (and focal point) pointing left
def construct_parabolic_mirror(p0, f, D, T, resolution):
    p1 = np.array([T,0])      # The central point on the back of the mirror
    f0 = np.array([-f,0])      # The focal point at the front of the mirror
    y0 = np.linspace(-D/2, D/2, 1+int(D/resolution))
    x0 = -y0**2/(4*f)  # Equation for an upward facing parabola through (0,0) and focal distance f (at y=f, the derivative is 1)
    p_corners = np.array([ [x0[0],y0[0]], [x0[-1],y0[-1]], [T,D/2], [T,-D/2] ])
    pts = np.column_stack((x0, y0))
    pts = np.vstack((pts, [[T,D/2],[T,-D/2]]))
    return [pts+p0, p1+p0, f0+p0, p_corners+p0]


def GLB_calculate_ZR(w0, M2, wavelength):  # GLB = Gaussian Laser Beam
    Z_R = np.pi * np.power(w0,2) / (M2*wavelength)
    return Z_R


def GLB_calculate_width_at_z(w0, ZR, z=0*mm):  # GLB = Gaussian Laser Beam
    w = w0 * np.sqrt(1 + np.power(z/ZR, 2))
    return w


def GLB_calculate_intensity_at_P(w0, ZR, intensity_total, p0, n0, beamwaist_distance, P):  # GLB = Gaussian Laser Beam
    from utils import geometry
    v = P - p0
    z = np.dot(v,n0) - beamwaist_distance   # This is with respect to the origin
    y = np.dot(v,geometry.orientation_from_normal(n0))
    w = GLB_calculate_width_at_z(w0=w0, ZR=ZR, z=z)
    intensity_peak = GLB_calculate_peak_intensity_at_z(w0=w0, ZR=ZR, intensity_total=intensity_total, z=z)
    return intensity_peak * np.exp(-2 * (y/w)**2)


def GLB_calculate_peak_intensity_at_z(w0, ZR, intensity_total, z=0*mm):  # GLB = Gaussian Laser Beam
    w = GLB_calculate_width_at_z(w0=w0, ZR=ZR, z=z)
    intensity_peak = np.sqrt(np.pi/2) * intensity_total/w
    return intensity_peak
