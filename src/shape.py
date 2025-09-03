from src.config import *

def generate_random_shape(
        fourier_series_order: int,
        points_per_shape: int,
        spikeyness_coefficient: float
    ) -> tuple[np.ndarray, np.ndarray]:
    """Generate a random shape using Fourier series.

    Args:
        fourier_series_order (int): The order of the Fourier series.
        points_per_shape (int): Number of points to generate for the shape.
        spikeyness_coefficient (float): Coefficient to control spikiness (higher values make it less spikey).

    Returns:
        tuple[np.ndarray, np.ndarray]: Radius array and theta array.
    """

    # randomize the fourier series amplitude and phase
    amplitude = np.random.randn(fourier_series_order)
    phase_offset = 2*np.pi*np.random.rand(fourier_series_order)

    # compute the fourier series
    theta = np.linspace(0,2*np.pi,points_per_shape)
    radius = np.zeros(points_per_shape)
    for frequency_order in range(1,fourier_series_order):
        radius+=amplitude[frequency_order]*np.cos(frequency_order*theta-phase_offset[frequency_order])

    # make radius positive and less spikey
    radius+=np.abs(spikeyness_coefficient*radius.min())

    return radius, theta


def discrete_grid_representation(radius: np.ndarray, pixel_dimension: int, random_shape_spline: Callable[[float], float]) -> np.ndarray:
    """Create a discrete grid representation of the shape.

    Args:
        radius (np.ndarray): Array of radii.
        pixel_dimension (int): Dimension of the pixel grid (square).
        random_shape_spline (Callable[[float], float]): Spline function for the shape.

    Returns:
        np.ndarray: 2D grid array representing the shape.
    """
    boundary_enforced_radius = np.max(radius)*1.1
    points = np.linspace(-boundary_enforced_radius, boundary_enforced_radius, pixel_dimension)

    grid = np.zeros([pixel_dimension,pixel_dimension],dtype=np.int32)
    for x in range(pixel_dimension):
        for y in range(pixel_dimension):
            angle = np.arctan2(points[y],points[x])+np.pi
            if random_shape_spline(angle) > np.sqrt(points[x]**2+points[y]**2):
                grid[x,y] = 1
    return grid