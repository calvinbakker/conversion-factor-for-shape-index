from src.config import *
from src.shape import *

def compute_spline_boundary(random_shape_spline: Callable[[float], float], theta: np.ndarray) -> tuple[float, float]:
    """Compute the area and boundary length using spline approximation.

    Args:
        random_shape_spline (Callable[[float], float]): Spline function for the shape.
        theta (np.ndarray): Array of angles.

    Returns:
        tuple[float, float]: Area size and boundary length.
    """

    # compute the area and boundary by integration over theta
    area_size = simpson(
        0.5*random_shape_spline(theta)**2, # radial area formula
        theta
    )
    boundary_length = simpson(
        np.sqrt(random_shape_spline(theta)**2+random_shape_spline(theta,1)**2), # standard Euclidian formula
        theta
        )
    return area_size, boundary_length


def compute_discrete_boundary(grid: np.ndarray, pixel_dimension: int) -> tuple[np.ndarray, int, int]:
    """Compute the discrete boundary of the shape using Moore neighborhood.

    Args:
        grid (np.ndarray): 2D grid array of the shape.
        pixel_dimension (int): Dimension of the grid.

    Returns:
        tuple[np.ndarray, int, int]: Boundary grid, boundary pixels for neighbors, boundary pixels for single.
    """
    order1_neighbourhood_size = 9
    order2_neighbourhood_size = 21

    # preprocess grid for boundary computation
    reshaped_grid = grid.reshape(1,pixel_dimension,pixel_dimension)
    neighbours_per_pixel = np.vstack(
            [list(reshaped_grid)] + 
            [list(np.roll(reshaped_grid.copy(),j,i)) for i in [1,2] for j in [-1,1]] +
            [list(np.roll(np.roll(reshaped_grid.copy(),i,1),j,2)) for i in [-1,1] for j in [-1,1]]
        ).reshape(order1_neighbourhood_size,pixel_dimension**2)
    find_nonzero_pixels = (neighbours_per_pixel[0,:]==1) # check if first element is equal to one

    # compute the order-0 single-pixel boundary
    boundary = (
        np.any(neighbours_per_pixel[1:].T==0,axis=1)    # any of the neighbours is zero (so the Moore neighbourhood is still used!)
        &                                               # and
        find_nonzero_pixels                             # the pixel itself is non-zero
    ).reshape(pixel_dimension,pixel_dimension)
    boundary_pixels_order0 = np.sum(boundary)

    # compute the order-1 Moore neighborhood boundary
    filter_nonzeros_out = (neighbours_per_pixel.T[find_nonzero_pixels]) # use the find_nonzero_pixels as a mask
    count_neighbouring_pixels = (filter_nonzeros_out[:,1:]==0) # excluding self (pixel) gives the `[:,1:]``, and the `==0` is because we count zero pixels around a non-empty one
    boundary_pixels_order1 = np.sum(count_neighbouring_pixels)

    # compute the order-2 neighbourhood boundary
    neighbours_per_pixel = [reshaped_grid]
    for nInt in range(20):
        a = nInt>>0 & 1
        b = nInt>>1 & 1
        c = nInt>>2 & 1
        d = nInt>>3 & 1
        e = nInt>>4 & 1
        dx = (2*a-1)*(c+e+1 - b*   (c+d+1)   )
        dy = (2*a-1)*( d-e  + b* (c-d+2*e+1) )
        neighbours_per_pixel+=[np.roll(np.roll(reshaped_grid.copy(),dx,1),dy,2)]
    neighbours_per_pixel = np.array(neighbours_per_pixel).reshape(order2_neighbourhood_size,pixel_dimension**2)
    find_nonzero_pixels = (neighbours_per_pixel[0,:]==1)
    filter_nonzeros_out = (neighbours_per_pixel.T[find_nonzero_pixels]) # use the find_nonzero_pixels as a mask
    count_neighbouring_pixels = (filter_nonzeros_out[:,1:]==0) # excluding self (pixel) gives the `[:,1:]``, and the `==0` is because we count zero pixels around a non-empty one
    boundary_pixels_order2 = np.sum(count_neighbouring_pixels)

    return boundary, boundary_pixels_order0, boundary_pixels_order1, boundary_pixels_order2


def compute_conversion_ratios(
    fourier_series_order: int,
    points_per_shape: int,
    spikeyness_coefficient: float,
    pixel_dimension: int,
    samples: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute conversion ratios for shape indices using different boundary orders.

    Args:
        fourier_series_order (int): The order of the Fourier series for shape generation.
        points_per_shape (int): Number of points to generate for the shape.
        spikeyness_coefficient (float): Coefficient to control spikiness (higher means less spikey).
        pixel_dimension (int): Dimension of the pixel grid (square).
        samples (int): Number of samples to generate for statistical analysis.

    Returns:
        tuple[np.ndarray, np.ndarray, np.ndarray]: Arrays of ratios for 0th, 1st, and 2nd order boundaries.
    """
    ratio_order0 = np.zeros((samples))
    ratio_order1 = np.zeros((samples))
    ratio_order2 = np.zeros((samples))
    for i in range(samples):

        # generate a random shape
        radius, theta = generate_random_shape(
            fourier_series_order = fourier_series_order,
            points_per_shape = points_per_shape,
            spikeyness_coefficient = spikeyness_coefficient 
        )

        # generate a spline approximation
        random_shape_spline = spline(theta,radius)

        # compute the area and boundary by integration over theta
        area_size = simpson(
            0.5*random_shape_spline(theta)**2, # radial area formula
            theta
        )
        boundary_length = simpson(
            np.sqrt(random_shape_spline(theta)**2+random_shape_spline(theta,1)**2), # standard Euclidian formula
            theta
            )
        shape_index = boundary_length/np.sqrt(area_size)



        # discrete representation
        grid = discrete_grid_representation(
            radius = radius,
            pixel_dimension = pixel_dimension,
            random_shape_spline = random_shape_spline
        )

        # compute discrete area
        area_pixels = np.sum(grid)

        # compute discrete boundary
        boundary, boundary_pixels_order0, boundary_pixels_order1, boundary_pixels_order2 = compute_discrete_boundary(grid, pixel_dimension)
        shape_index_order0 = boundary_pixels_order0/np.sqrt(area_pixels)
        shape_index_order1 = boundary_pixels_order1/np.sqrt(area_pixels)
        shape_index_order2 = boundary_pixels_order2/np.sqrt(area_pixels)

        # compute the ratios and write to array
        ratio_order0[i] = (shape_index/shape_index_order0)
        ratio_order1[i] = (shape_index/shape_index_order1)
        ratio_order2[i] = (shape_index/shape_index_order2)
    return ratio_order0, ratio_order1, ratio_order2