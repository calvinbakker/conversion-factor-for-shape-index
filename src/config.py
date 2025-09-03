import numpy as np
from scipy.interpolate import CubicSpline as spline
from scipy.integrate import simpson

from typing import Callable

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import cmasher as cmr
