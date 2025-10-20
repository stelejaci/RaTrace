import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from functools import wraps
import time
import os, sys
from datetime import datetime

mm = 1
m = mm*1000
cm = mm*100
µm = mm/1000
nm = µm/1000
rad = 1
deg = np.pi*rad/180
X = 0
Y = 1

EPSILON = 1e-9


def sort_x_left_to_right(pts):
    pts = pts[np.argsort(pts[:, 0])]
    return pts


def y_from_x_and_2_points(pt1=[0, 0], pt2=[0, 1], x=0*mm):
    return np.interp(x, np.array([pt1[X], pt2[X]]), np.array([pt1[Y], pt2[Y]]))


def generate_color(index):
    cols = ('red', 'green', 'blue', 'cyan', 'magenta', 'yellow')
    col = cols[np.mod(index, len(cols))]
    return col

def load_colormap(color=(0,0,0,1), N_rays=1, wavelength=None):
    if color == 'wavelength'  or  color is None:
        if wavelength is None:
            from light.light_class import WAVELENGTH_DEFAULT
            wavelength = WAVELENGTH_DEFAULT
        cols = colormap_wavelength(N=N_rays, wavelength=wavelength)
    elif color == 'rainbow':
        cols = colormap_rainbow(N=N_rays)
    elif isinstance(color, str):
        col = matplotlib.colors.to_rgb(color)
        cols = colormap_fixed(N=N_rays, color=col)
    else:
        cols = colormap_fixed(N=N_rays, color=color)
    return cols

def colormap_rainbow(N):
    colors = []
    if N > 1:
        cmap = plt.get_cmap('jet')
        colors = [cmap(i) for i in np.linspace(0, 1, N)]
        colors = [tuple(color) for color in colors]
    else:
        colors.append((0, 0, 0, 1))
    return colors

def colormap_wavelength(wavelength, N=1, gamma=0.8):
    """
    Convert a wavelength (in nm) in the visible spectrum (380-750 nm) to an RGB color.
    Based on code by Dan Bruton: http://www.physics.sfasu.edu/astro/color/spectra.html
    """
    wavelength_nm = wavelength*mm/nm
    # Determine color based on wavelength ranges
    if 380 <= wavelength_nm < 440:
        attenuation = 0.3 + 0.7 * (wavelength_nm - 380) / (440 - 380)
        R = ((-(wavelength_nm - 440) / (440 - 380)) * attenuation) ** gamma
        G = 0.0
        B = (1.0 * attenuation) ** gamma
    elif 440 <= wavelength_nm < 490:
        R = 0.0
        G = ((wavelength_nm - 440) / (490 - 440)) ** gamma
        B = 1.0
    elif 490 <= wavelength_nm < 510:
        R = 0.0
        G = 1.0
        B = (-(wavelength_nm - 510) / (510 - 490)) ** gamma
    elif 510 <= wavelength_nm < 580:
        R = ((wavelength_nm - 510) / (580 - 510)) ** gamma
        G = 1.0
        B = 0.0
    elif 580 <= wavelength_nm < 645:
        R = 1.0
        G = (-(wavelength_nm - 645) / (645 - 580)) ** gamma
        B = 0.0
    elif 645 <= wavelength_nm <= 750:
        attenuation = 0.3 + 0.7 * (750 - wavelength_nm) / (750 - 645)
        R = (1.0 * attenuation) ** gamma
        G = 0.0
        B = 0.0
    else:
        R, G, B = 0.0, 0.0, 0.0      # Black if outside visible spectrum

    color = (R, G, B, 1)

    return colormap_fixed(N, color)


def colormap_fixed(N, color):
    return N*[color]


# def generate_gaussian_sampling(FWHM=1, nr_of_samples=11, verbose=False, randomize=False):
#     if FWHM == 0:
#         x_sampling = np.zeros((nr_of_samples,))
#         return x_sampling
#
#     sigma = FWHM / (2 * np.sqrt(np.log(2)))
#     half_range = 3 * sigma
#     dx = FWHM/10000
#
#     x = np.arange(-half_range, half_range, dx)
#     probability = np.exp(-np.square(x/sigma))
#     probability = probability / max(probability)
#
#     probability_cumsum = np.cumsum(probability * dx)
#     probability_cumsum = probability_cumsum / max(probability_cumsum)
#
#     pcss_start = 0.01
#     pcss_end = 1 - pcss_start
#     # d_sampling = (pcss_end - pcss_start) / (nr_of_samples - 1)
#     if randomize==False:
#         probability_cumsum_sampling = np.linspace(pcss_start, pcss_end, nr_of_samples)
#     else:
#         probability_cumsum_sampling = pcss_start + np.random.rand(nr_of_samples) * (pcss_end-pcss_start)
#     x_sampling = np.interp(probability_cumsum_sampling, probability_cumsum, x)
#
#     if verbose:
#         plt.figure(num=None, figsize=(24, 15), dpi=80, facecolor='w', edgecolor='k')
#         ax = plt.subplot2grid((1,1), (0,0), colspan=1, rowspan=1)
#         plt.plot(x/deg, probability, 'b-')
#         plt.plot(FWHM/2 * np.array([-1,-1,1,1])/deg, 0.5*np.array([0,1,1,0]), 'k--')
#         plt.plot(x/deg, probability_cumsum, 'g-')
#         plt.scatter(x_sampling/deg, 0.10+0*x_sampling, c='k') # * rand(length(angles_sampling), 1) / 20, 'markerfacecolor', 'r', 'markeredgecolor', 'none', 'sizedata', 20);
#         for pcss in probability_cumsum_sampling:
#             plt.plot(half_range*np.array([-1,1])/deg, pcss*np.array([1,1]), 'r--', linewidth=0.5)
#         ax.grid()
#         plt.ylim((0, 1))
#         plt.xlabel('X  (°)')
#         plt.show()
#
#     return x_sampling

def alpha_from_N(N):
    N = N[0] if isinstance(N, list) else N
    return (N - 1) ** 2


def generate_gaussian_pcs(FWHM=1*deg, verbose=False):   # Gaussian probability cumulative sum
    if FWHM == 0:
        return ([], [])

    sigma = FWHM / (2 * np.sqrt(np.log(2)))
    half_range = 3 * sigma

    dx = FWHM / 10000
    x = np.arange(-half_range, half_range, dx)

    probability = np.exp(-np.square(x / sigma))
    probability = probability / max(probability)

    probability_cumsum = np.cumsum(probability * dx)
    probability_cumsum = probability_cumsum / max(probability_cumsum)

    if verbose:
        plt.figure(num=None, figsize=(24, 15), dpi=80, facecolor='w', edgecolor='k')
        ax = plt.subplot2grid((1, 1), (0, 0), colspan=1, rowspan=1)
        plt.plot(x / deg, probability, 'b-')
        plt.plot(FWHM / 2 * np.array([-1, -1, 1, 1]) / deg, 0.5 * np.array([0, 1, 1, 0]), 'k--')
        plt.plot(x / deg, probability_cumsum, 'g-')
        ax.grid()
        plt.ylim((0, 1))
        plt.xlabel('X  (°)')
        plt.show()

    return (x, probability_cumsum)


def generate_samples_from_pcs(x, probability_cumsum, nr_of_samples=1, randomize=False, verbose=False):

    pcs_start = 0.01
    pcs_end = 1 - pcs_start

    if randomize == False  and  nr_of_samples == 1:
        probability_cumsum_sampling = 0.5
    elif randomize == False:
            probability_cumsum_sampling = np.linspace(pcs_start, pcs_end, nr_of_samples)
    else:
        probability_cumsum_sampling = pcs_start + np.random.rand(nr_of_samples) * (pcs_end - pcs_start)

    if not np.any(probability_cumsum):
        x_sampling = np.array([0])
    else:
        x_sampling = np.interp(probability_cumsum_sampling, probability_cumsum, x)

    # In case of only 1 sample it does not create an array
    if not isinstance(x_sampling, np.ndarray):
        x_sampling = np.array([x_sampling])

    if verbose:
        plt.figure(num=None, figsize=(24, 15), dpi=80, facecolor='w', edgecolor='k')
        ax = plt.subplot2grid((1, 1), (0, 0), colspan=1, rowspan=1)
        plt.plot(x/ deg, probability_cumsum, 'g-')
        plt.scatter(x_sampling / deg, 0.10 + 0 * x_sampling, c='k')  # * rand(length(angles_sampling), 1) / 20, 'markerfacecolor', 'r', 'markeredgecolor', 'none', 'sizedata', 20);
        # for pcss in probability_cumsum_sampling:
        #     plt.plot(half_range * np.array([-1, 1]) / deg, pcss * np.array([1, 1]), 'r--', linewidth=0.5)
        ax.grid()
        plt.ylim((0, 1))
        plt.xlabel('X  (°)')
        plt.show()

    return x_sampling


def calculate_FWHM(intensities):

    peak_intensity = np.max(intensities)

    if peak_intensity==0:
        return [0, []]

    half_max = peak_intensity / 2
    ind = np.where(intensities > half_max)
    ind = ind[0]
    FWHM, pts = None, []
    if any(ind)  and  not ind[0] == 0  and not  ind[-1]==len(intensities)-1:
        # Linear interpolation of the first point
        x1, x2 = ind[0]-1, ind[0]
        y1, y2 = intensities[x1], intensities[x2]
        x = x1 + (x2 - x1) / (y2 - y1) * (half_max - y1)
        pts.append(np.array([x, half_max]))

        # Linear interpolation of the second point
        x1, x2 = ind[-1], ind[-1]+1
        y1, y2 = intensities[x1], intensities[x2]
        x = x1 + (x2 - x1) / (y2 - y1) * (half_max - y1)
        pts.append(np.array([x, half_max]))

        FWHM = pts[1][X] - pts[0][X]

    return [FWHM, pts]


def parameter_string(par, val):
    # par_str = '\n'
    if par[0:5] == 'empty':
        par_str = ' '
        return par_str
    if par[0] == '=':
        par_str = par[:-2]
        return par_str
    if not val:
        par_str = par
        return par_str
    eenheid = val[1]
    val = val[0]
    if type(val) == str:
        par_str = par + ': ' + val + eenheid
    elif type(val) == float:
        par_str = par + ': ' + str(val) + eenheid
    elif type(val) == np.float64:
        par_str = par + ': ' + str(round(val * 100) / 100) + eenheid
    elif type(val) == int:
        par_str = par + ': ' + str(val) + eenheid
    elif type(val) == bool:
        par_str = par + ': ' + str(val)
    elif type(val) == np.ndarray:
        par_str = par + ': (' + str(val[0]) + ', ' + str(val[1]) + ')' + eenheid
    elif type(val) == tuple:
        par_str = '(Type tuple is not yet implemented)'
    return par_str


def display_parameters_plot(pars={}, x=0, y=0, dy=1, font_size=10, col='black'):
    plt.gca().add_patch(plt.Rectangle((x - 1, y - (len(pars) + 59) * dy), width=100, height=(len(pars) + 59) * dy, edgecolor='grey', facecolor=(0.95, 0.95, 0.95), fill=True, linewidth=1, zorder=10))
    for par in pars:
        y -= dy
        val = pars[par]
        plt.text(x, y, parameter_string(par, val), fontsize=font_size, color=col, zorder=10)


def display_parameters(pars={}):
    print('\n          o-------------------------------------------------------------------------------------------------------------------------------------------o')
    for par in pars:
        val = pars[par]
        print('          |  ' + parameter_string(par, val))
    print('          o-------------------------------------------------------------------------------------------------------------------------------------------o\n\n')


def print_progress_V1(i, N):
    if N < 10:
        return
    if np.mod(i + 1, N / 100) == 0:
        sys.stdout.write("\r   --> %d%%" % np.round(100 * i / N))
        sys.stdout.flush()
    if i == N - 1:
        sys.stdout.write("\n")
    return


def print_progress(i, N):
    if N < 10:
        return

    decimals, length, fill = 0, 50, '#'
    printEnd, prefix, suffix = '', '   ', ''

    i += 1
    percent = ("{0:." + str(decimals) + "f}").format(100 * (i / float(N)))
    filledLength = int(length * (i) // N)
    bar = fill * filledLength + '-' * (length - filledLength)
    progress_str = f'\r{prefix} |{bar}| {percent}% {suffix}'
    print(progress_str, end=printEnd)

    if i == N:  # Print New Line on Complete
        print()


def timeit(func):
    @wraps(func)
    def timeit_wrapper(*args, **kwargs):
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        end_time = time.perf_counter()
        total_time = end_time - start_time
        print(f'Function {func.__name__} took {total_time:.3f} seconds')
        return result

    return timeit_wrapper


def format_1D_array(arr, fmt):
    if arr.ndim == 1:
        s = ' '.join(f'{x:{fmt}}' for x in arr)
        return f'[{s}]'
    else:
        return ''


def output_source_data_in_text_file(filename, sources, add_timestamp=False):
    from light import light_class

    # Add timestamp data to the filename, is needed
    if add_timestamp:
        current_time = datetime.now()
        timestamp_string = current_time.strftime("%Y%m%d_%H%M%S")
        filename_base = os.path.splitext(filename)[0]
        filename = filename_base + '_' + timestamp_string + '.txt'
        print(f'Adding timestamp information to filename: {filename}')

    with open(filename, 'w') as file:
        line = f'{"Source ID":9s} | {"ray ID":<9s} | {"p0":<30s} | {"p1":<30s} | {"length":<11s} | {"intensity":<11s} | {"N":<5s} | {"wavelength":<8s} | {"phi0":<9s} | {"phi1":<9s} | {"source":<6s} | {"impact":<6s} | {"display ID":<6s} | {"display x":<10s} | {"pixel ID":<6s} | {"parent":<8s} | {"children":<8s}'
        file.write(line + "\n")
        for source in sources:
            for ray in source.rays:
                line = f"S{source.ID:<8d} | {ray.ID:>09d} | [{ray.p0[0]:>+13.6f}, {ray.p0[1]:>+13.6f}] | "
                line = line + f"[{'-':>13s}, {'-':>13s}] | {'-':>11s} | "     if ray.p1 is None      else      line + f"[{ray.p1[0]:>+13.6f}, {ray.p1[1]:>+13.6f}] | {ray.length:>11.6f} | "
                line = line + f"{ray.intensity:>11.6f} | {ray.N:>5.3f} | {ray.wavelength:>10.6f} | {ray.phase_start:>+8.6f} | "
                line = line + f"{'-':>9s} | "     if ray.phase_end is None      else     line + f"{ray.phase_end:>+8.6f} | "
                line = line + f"S{ray.source_element.ID:<5d} | "     if isinstance(ray.source_element, light_class.LightSourceClass)     else     line + f"E{ray.source_element.ID:<5d} | "
                line = line + f"E{'-':<5s} | "     if ray.element_hit is None     else     line + f"E{ray.element_hit.ID:<5d} | "
                line = line + f"D{'-':<9s} | "     if ray.display_ID is None     else     line + f"D{ray.display_ID:<9d} | "
                line = line + f"{'-':>10s} | "     if ray.display_x is None     else     line + f"{ray.display_x:>+10.6f} | "
                line = line + f"P{'-':<7s} | "     if ray.imager_pixel_ID is None     else     line + f"P{ray.imager_pixel_ID:>07d} | "
                line = line + f"R{'-':<7s} | "     if not ray.ID_parent     else     line + f"R{ray.ID_parent:>07d} | "

                if not ray.ID_children:
                    line = line + f"[R{'-':<8s}]"
                else:
                    line = line + f"["
                    for i_child in range(len(ray.ID_children)):
                        line = line + f"R{ray.ID_children[i_child]:>08d}"
                        if i_child < len(ray.ID_children)-1:
                            line = line + f","
                    line = line + f"]"

                file.write(line + "\n")


def output_imager_data_in_text_file(filename, imagers):
    with open(filename, 'w') as file:
        line = f'{"Imager ID":<9s} | {"pixel":<5s} | {"intensity":<15s} | {"phase":<9s} | {"X":<12s} | {"Y":<12s}'  #
        file.write(line + "\n")
        for imager in imagers:
            for i_px in range(imager.number_of_pixels):
                # pixel = imager.pixels[i_px]
                line = f"I{imager.ID:<8d} | {i_px:>05d} | {imager.intensity[i_px]:>15.6f} | {imager.phase[i_px]:>+8.6f} | {imager.pixels[i_px].p0[X]:>+12.6f} | {imager.pixels[i_px].p0[Y]:>+12.6f}"
                file.write(line + "\n")

def plot_arrow_end_at_P(graph, P, r, s, angle, col):
    from utils import geometry
    n = geometry.normal_from_orientation(r)
    p0 = P + s * n * np.sin(angle) + s * r * np.cos(angle)
    p1 = P - s * n * np.sin(angle) + s * r * np.cos(angle)
    graph.plot([P[X], p0[X]], [P[Y], p0[Y]], color=col, linewidth=3, linestyle='solid', alpha=1, zorder=5)
    graph.plot([P[X], p1[X]], [P[Y], p1[Y]], color=col, linewidth=3, linestyle='solid', alpha=1, zorder=5)


if __name__ == "__main__":
    (scattering_angles, scattering_pcs) = generate_gaussian_pcs(FWHM=0*deg, verbose=False)
    S = generate_samples_from_pcs(x=scattering_angles, probability_cumsum=scattering_pcs, nr_of_samples=1, randomize=False, verbose=False)
    # print(G)
    print(scattering_angles)
    print(scattering_pcs)
    print(S)
