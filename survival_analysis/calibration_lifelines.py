#This is a copy from https://github.com/CamDavidsonPilon/lifelines/blob/cca5288acf67a1ad5ef959a839e01cd2ad06c5b2/lifelines/calibration.py of the Lifelines.probability_calibration_plot() modified for stilising the plots. Both the number of nots and the penaliser in the regression were not changed in the final plots


# -*- coding: utf-8 -*-
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from lifelines.utils import CensoringType
from lifelines.fitters import RegressionFitter
from lifelines import CRCSplineFitter


def survival_probability_calibration(model: RegressionFitter, df: pd.DataFrame, t0: float, ax=None, left_y_axis=False, right_y_axis=False, event='AS', n_knots=3):
    r"""
    Smoothed calibration curves for time-to-event models. This is analogous to
    calibration curves for classification models, extended to handle survival probabilities
    and censoring. Produces a matplotlib figure and some metrics.
    We want to calibrate our model's prediction of :math:`P(T < \text{t0})` against the observed frequencies.
    Parameters
    -------------
    model:
        a fitted lifelines regression model to be evaluated
    df: DataFrame
        a DataFrame - if equal to the training data, then this is an in-sample calibration. Could also be an out-of-sample
        dataset.
    t0: float
        the time to evaluate the probability of event occurring prior at.
    Returns
    ----------
    ax:
        mpl axes
    ICI:
        mean absolute difference between predicted and observed
    E50:
        median absolute difference between predicted and observed
    https://onlinelibrary.wiley.com/doi/full/10.1002/sim.8570
    """

    def ccl(p):
        return np.log(-np.log(1 - p))

    if ax is None:
        ax = plt.gca()

    T = model.duration_col
    E = model.event_col

    predictions_at_t0 = np.clip(1 - model.predict_survival_function(df, times=[t0]).T.squeeze(), 1e-10, 0.2 - 1e-10)

    # create new dataset with the predictions
    prediction_df = pd.DataFrame({"ccl_at_%d" % t0: ccl(predictions_at_t0), T: df[T], E: df[E]})

    # fit new dataset to flexible spline model
    # this new model connects prediction probabilities and actual survival. It should be very flexible, almost to the point of overfitting. It's goal is just to smooth out the data!
    regressors = {"beta_": ["ccl_at_%d" % t0], "gamma0_": "1","gamma1_": "1",
                  # "gamma2_": "1",
                  # "gamma3_": "1" 
                 }

    # this model is from examples/royson_crowther_clements_splines.py
    crc = CRCSplineFitter(n_baseline_knots=2, penalizer=0.000001)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        if CensoringType.is_right_censoring(model):
            crc.fit_right_censoring(prediction_df, T, E, regressors=regressors)
        elif CensoringType.is_left_censoring(model):
            crc.fit_left_censoring(prediction_df, T, E, regressors=regressors)
        elif CensoringType.is_interval_censoring(model):
            crc.fit_interval_censoring(prediction_df, T, E, regressors=regressors)

    # predict new model at values 0 to 1, but remember to ccl it!
    x = np.linspace(np.clip(predictions_at_t0.min() - 0.01, 0, 0.2), np.clip(0.75, 0, 0.2), 100)
    y = 1 - crc.predict_survival_function(pd.DataFrame({"ccl_at_%d" % t0: ccl(x)}), times=[t0]).T.squeeze()

    # plot our results
    #ax.set_title("Smoothed calibration curve of \npredicted vs observed probabilities of t ≤ %d mortality" % t0)

    color = "tab:red"
    ax.plot(x, y, label="Calibration curve", color=color)
    
    ax.set_xlabel(f"Predicted probability of \nt ≤ %d {event}" % t0)
    if left_y_axis == True:
       # ax.set_ylabel("Observed probability of \nt ≤ %d mortality" % t0, color=color)
        ax.set_ylabel(f"Observed probability of \n {event} at t ≤ 5, 10, 20 years", color=color)
    ax.tick_params(axis="y", labelcolor=color)
    
    ax.set_xlim([0,0.2])
    ax.set_ylim([0,0.2])

    # plot x=y line
    ax.plot(x, x, c="k", ls="--")
    ax.legend(prop={'size': 6})

    # plot histogram of our original predictions
    color = "tab:blue"
    twin_ax = ax.twinx()
    if right_y_axis == True:
        twin_ax.set_ylabel("Count of \npredicted probabilities", color=color)  # we already handled the x-label with ax1
    twin_ax.tick_params(axis="y", labelcolor=color)
    twin_ax.hist(predictions_at_t0, alpha=0.3, bins="sqrt", color=color)
    
    plt.tight_layout()

    deltas = ((1 - crc.predict_survival_function(prediction_df, times=[t0])).T.squeeze() - predictions_at_t0).abs()
    ICI = deltas.mean()
    E50 = np.percentile(deltas, 50)
    # print("ICI = ", ICI)
    # print("E50 = ", E50)

    return ax, ICI, E50