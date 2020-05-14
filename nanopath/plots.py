"""

Pathfinder plotting module, @esteinig

Common plotting functions for workflow debugging. Main plots are handled with MongoAPI in Vue.

"""

import matplotlib.pyplot as plt
import numpy as np

from sklearn.linear_model import LinearRegression


def plot_date_randomisation(
    ax: plt.axes,
    replicates: np.array or list,
    rate: float,
    log10: bool = True
) -> plt.axes:

    """ Plot distribution of substitution rates for date randomisation test

    :param ax: axes object to plot the date randomisation
    :param replicates: list of replicate substitution rate estimates
    :param rate: true rate estimate vertical line for evaluation
    :param log10: plot log10 of substitution rates on horizontal axis

    :returns axes object

    """

    if log10:
        replicates = np.log10(replicates)

    with plt.style.context('seaborn-colorblind'):
        ax.hist(x=replicates, color='gray')
        ax.axvline(x=rate, color='r')

        ax.set_xlabel(f'{"Log10 Rate" if log10 else "Rate"}')
        ax.set_xlabel(f'Frequency')

    return ax


def plot_date_regression(
    ax: plt.axes,
    x: np.array,
    y: np.array,
) -> plt.axes:

    """ Plot regression between dates and root-to-tip distances

    :param ax: axes object to plot on
    :param x: date array
    :param y: root-to-tip distance array
    :returns axes object

    """

    x = x.reshape(-1, 1)  # 2-d
    y = y.reshape(-1, 1)

    linear_regressor = LinearRegression()
    linear_regressor.fit(x, y)
    y_pred = linear_regressor.predict(y)

    ax.scatter(x, y)
    ax.plot(x, y_pred, color='r')

    return ax
