"""
Adapted from https://github.com/Warvito/neurocombat_sklearn
"""
import numpy as np
import pandas as pd
from typing import List, Union, Optional, Tuple
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.compose import make_column_transformer
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LinearRegression
from sklearn.utils import check_array
from sklearn.utils.validation import (check_is_fitted, check_consistent_length, FLOAT_DTYPES)


__all__ = [
    'CombatModel',
    'RegressOutModel'
]


def _get_column_types(df):
    categorical_vars = []
    continuous_vars = []

    for col in df.columns:
        if pd.api.types.is_numeric_dtype(df[col]):
            continuous_vars.append(col)
        else:
            categorical_vars.append(col)

    return categorical_vars, continuous_vars


class RegressOutModel(BaseEstimator):
    def __init__(self, remove: Optional[List[str]] = None):
        self.remove = remove
        self.clf = None

    def fit(self, x, conf):
        categorical_cov, continuous_cov = _get_column_types(conf[self.remove])
        column_transformer = make_column_transformer(
            (StandardScaler(), continuous_cov),
            (OneHotEncoder(drop='first'), categorical_cov)
        )

        self.clf = Pipeline([('ct', column_transformer), ('lr', LinearRegression())])
        self.clf.fit(conf[self.remove], x)

        return self

    def transform(self, x, conf):
        check_is_fitted(self, 'clf')

        residuals = x - self.clf.predict(conf[self.remove])
        residuals += self.clf.named_steps['lr'].intercept_
        return residuals.astype(np.float32)

    def fit_transform(self, x, conf, *args):
        """Fit to data, then transform it"""
        return self.fit(x, conf, *args).transform(x, conf, *args)


def _find_priors(gamma_hat, delta_hat):
    """Compute a and b priors"""
    gamma_bar = np.mean(gamma_hat, axis=1)
    tau_2 = np.var(gamma_hat, axis=1, ddof=1)

    a_prior, b_prior = [], []
    for delta in delta_hat:
        m = np.mean(delta)
        s2 = np.var(delta, ddof=1, dtype=np.float32)
        a_prior.append((2 * s2 + m ** 2) / s2)
        b_prior.append((m * s2 + m ** 3) / s2)

    return gamma_bar, tau_2, a_prior, b_prior


def _postmean(gamma_hat, gamma_bar, n, delta_star, tau_2):
    return (tau_2 * n * gamma_hat + delta_star * gamma_bar) / (tau_2 * n + delta_star)


def _postvar(sum_2, n, a_prior, b_prior):
    return (0.5 * sum_2 + b_prior) / (n / 2.0 + a_prior - 1.0)


def _iteration_solver(standardized_data,
                      gamma_hat, delta_hat,
                      gamma_bar, tau_2,
                      a_prior, b_prior,
                      convergence=0.0001):
    """Compute iterative method to find the parametric site/batch effect adjustments"""
    n = (1 - np.isnan(standardized_data)).sum(axis=0)
    gamma_hat_new = gamma_hat_old = gamma_hat.copy()
    delta_hat_new = delta_hat_old = delta_hat.copy()

    change = 1
    # count = 0

    while change > convergence:
        gamma_hat_new = _postmean(gamma_hat, gamma_bar, n, delta_hat_old, tau_2)

        ssd = ((standardized_data - gamma_hat_new) ** 2).sum(axis=0)
        delta_hat_new = _postvar(ssd, n, a_prior, b_prior)

        change = max((abs(gamma_hat_new - gamma_hat_old) / gamma_hat_old).max(),
                     (abs(delta_hat_new - delta_hat_old) / delta_hat_old).max())

        gamma_hat_old = gamma_hat_new
        delta_hat_old = delta_hat_new

        # count = count + 1

    return gamma_hat_new, delta_hat_new


def _find_parametric_adjustments(standardized_data, site_design,
                                 gamma_hat, delta_hat,
                                 gamma_bar, tau_2,
                                 a_prior, b_prior):
    """Compute empirical Bayes site/batch effect parameter estimates using parametric empirical priors"""

    gamma_star, delta_star = [], []
    for i, site_col in enumerate(site_design.T):
        gamma_hat_adjust, delta_hat_adjust = _iteration_solver(standardized_data[site_col > 0],
                                                               gamma_hat[i], delta_hat[i],
                                                               gamma_bar[i], tau_2[i],
                                                               a_prior[i], b_prior[i])

        gamma_star.append(gamma_hat_adjust)
        delta_star.append(delta_hat_adjust)

    return np.asarray(gamma_star), np.asarray(delta_star)


def _fit_ls_model(standardized_data, site_design):
    """Location and scale (L/S) adjustments

    Parameters
    ----------
    standardized_data : np.ndarray of shape (n_samples, n_features)
    site_design : np.ndarray of shape (n_samples, n_sites)
        Onehot encoded design matrix for site.
    """

    gamma_hat = np.linalg.pinv(site_design) @ standardized_data
    delta_hat = np.vstack([np.var(standardized_data[site_col > 0], axis=0, ddof=1) for site_col in site_design.T])
    return gamma_hat, delta_hat


class CombatModel(BaseEstimator):
    """Harmonize/normalize features using Combat's [1] parametric empirical Bayes framework

    [1] Fortin, Jean-Philippe, et al. "Harmonization of cortical thickness
    measurements across scanners and sites." Neuroimage 167 (2018): 104-120.
    """

    def __init__(self, site_key='SITE', keep: Optional[List[str]] = None, remove: Optional[List[str]] = None, copy=True):
        self.site_key = site_key
        self.keep = keep
        self.remove = remove
        self.copy = copy

    def _reset(self):
        """Reset internal data-dependent state, if necessary.

        __init__ parameters are not touched.
        """

        # Checking one attribute is enough, because they are all set together
        if hasattr(self, 'n_sites'):
            del self.n_sites
            del self.site_ids
            del self.site_encoder
            del self.categorical_encoders_keep
            del self.categorical_encoders_remove
            del self.beta_hat
            del self.grand_mean
            del self.var_pooled
            del self.gamma_star
            del self.delta_star

    def fit(self, x: np.ndarray, conf: pd.DataFrame):

        # Reset internal state before fitting
        self._reset()

        # x = check_array(x, copy=self.copy, estimator=self, dtype=FLOAT_DTYPES)
        # site = check_array(conf[[site_key]].to_numpy(), copy=self.copy, estimator=self, dtype=None)

        site = conf[[self.site_key]].to_numpy()
        check_consistent_length(x, site)

        self.site_ids = np.unique(site)
        self.n_sites = self.site_ids.size
        design = self._make_design_matrix(site, conf, fitting=True)

        standardized_data, _ = self._standardize_across_features(x, design, fitting=True)
        gamma_hat, delta_hat = _fit_ls_model(standardized_data, design[:, :self.n_sites])
        gamma_bar, tau_2, a_prior, b_prior = _find_priors(gamma_hat, delta_hat)

        self.gamma_star, self.delta_star = _find_parametric_adjustments(standardized_data, design[:, :self.n_sites],
                                                                        gamma_hat, delta_hat, gamma_bar, tau_2,
                                                                        a_prior, b_prior)

        return self

    def _standardize_across_features(self, data: np.ndarray, design: np.ndarray, fitting=False) \
            -> Tuple[np.ndarray, Union[np.ndarray, None]]:
        """Standardization of the features

        The magnitude of the features could create bias in the empirical Bayes estimates of the prior distribution.
        To avoid this, the features are standardized to all of them have similar overall mean and variance.

        Parameters
        ----------
        data :
            Features
        design :
            Design matrix
        fitting : boolean, default is False
            Indicates if this method is executed inside the
            fit method (in order to save the parameters to use later).

        Returns
        -------
        standardized_data : array-like
        standardized_mean : array-like
            Standardized mean used during the process
        """

        if fitting:
            self.beta_hat = np.linalg.pinv(design) @ data

            # Standardization Model
            prop_samples_per_site = design[:, :self.n_sites].mean(axis=0)
            self.grand_mean = prop_samples_per_site @ self.beta_hat[:self.n_sites]

            residuals = data - (design @ self.beta_hat)
            self.var_pooled = np.mean(residuals ** 2, axis=0)

        standardized_mean = (design[:, self.n_sites:] @ self.beta_hat[self.n_sites:])
        standardized_mean += self.grand_mean

        standardized_mean_keep = 0
        if self.n_keep > 0:
            standardized_mean_keep = (design[:, -self.n_keep:] @ self.beta_hat[-self.n_keep:])
            standardized_mean_keep += self.grand_mean

        standardized_data = data - standardized_mean
        standardized_data /= np.sqrt(self.var_pooled)

        # return standardized_data, standardized_mean
        return standardized_data, standardized_mean_keep

    def _make_design_matrix(self, site: np.ndarray, conf: pd.DataFrame, fitting=False) -> np.ndarray:
        """Method to create a design matrix that contain:

            - One-hot encoding of the sites [n_samples, n_sites]
            - One-hot encoding of each discrete covariates (removing
            the first column) [n_samples, (n_discrete_covivariate_names-1) * n_discrete_covariates]
            - Each continuous covariates

        Parameters
        ----------
        site :
            Site data.
        conf :
            Dataframe of covariates
        fitting : boolean, default is False
            Indicates fitting stage.

        Returns
        -------
        design : array-like
            The design matrix.
        """

        keep = [] if self.keep is None else self.keep
        remove = [] if self.remove is None else self.remove
        categorical_cov_keep, continuous_cov_keep = _get_column_types(conf[keep])
        categorical_cov_remove, continuous_cov_remove = _get_column_types(conf[remove])

        if fitting:
            self.site_encoder = OneHotEncoder(sparse_output=False).fit(site)

            self.categorical_encoders_keep = []
            for _, x in conf[categorical_cov_keep].items():
                enc = OneHotEncoder(sparse_output=False, drop='first').fit(x.to_numpy()[:, None])
                self.categorical_encoders_keep.append(enc)

            self.categorical_encoders_remove = []
            for _, x in conf[categorical_cov_remove].items():
                enc = OneHotEncoder(sparse_output=False, drop='first').fit(x.to_numpy()[:, None])
                self.categorical_encoders_remove.append(enc)

        design = []
        sites_design = self.site_encoder.transform(site)
        design.append(sites_design)

        design_keep = []
        for i, k in enumerate(categorical_cov_keep):
            design_keep.append(self.categorical_encoders_keep[i].transform(conf[[k]].to_numpy()))
        if len(continuous_cov_keep) > 0:
            design_keep.append(conf[continuous_cov_keep].to_numpy())

        design_remove = []
        for i, k in enumerate(categorical_cov_remove):
            design_remove.append(self.categorical_encoders_remove[i].transform(conf[[k]].to_numpy()))
        if len(continuous_cov_remove) > 0:
            design_remove.append(conf[continuous_cov_remove].to_numpy())

        self.n_keep = 0
        if len(design_keep) > 0:
            design_keep = np.hstack(design_keep)
            design.append(design_keep)
            self.n_keep = design_keep.shape[1]

        self.n_remove = 0
        if len(design_remove) > 0:
            design_remove = np.hstack(design_remove)
            design.append(design_remove)
            self.n_remove = design_remove.shape[1]

        return np.hstack(design)

    def transform(self, x: np.ndarray, conf: pd.DataFrame):
        """Transform data to harmonized space

        Parameters
        ----------
        x :
            Input data that will be transformed.
        conf :
            Covariates including site
        """

        check_is_fitted(self, 'n_sites')

        # data = check_array(data, copy=self.copy, estimator=self, dtype=FLOAT_DTYPES)
        # site = check_array(site, copy=self.copy, estimator=self)

        site = conf[[self.site_key]].to_numpy()
        check_consistent_length(x, site)

        new_site_ids = np.unique(site)
        unseen_sites = list(set(new_site_ids).difference(self.site_ids))
        if len(unseen_sites) > 0:
            raise ValueError(f'Deconfounding does not support unseen sites: {unseen_sites}.')

        design = self._make_design_matrix(site, conf, fitting=False)
        standardized_data, standardized_mean = self._standardize_across_features(x, design, fitting=False)
        bayes_data = self._adjust_data_final(standardized_data, design[:, :self.n_sites], standardized_mean)

        return bayes_data

    def fit_transform(self, x, conf, *args):
        """Fit to data, then transform it"""
        return self.fit(x, conf, *args).transform(x, conf, *args)

    def _adjust_data_final(self, standardized_data, site_design, standardized_mean):
        """Compute the harmonized/normalized data"""

        bayes_data = standardized_data.copy()
        for j, site_col in enumerate(site_design.T):
            mask_site = site_col > 0
            bayes_data[mask_site] -= self.gamma_star[j]
            bayes_data[mask_site] /= np.sqrt(self.delta_star[j])

        bayes_data *= np.sqrt(self.var_pooled)
        if standardized_mean is not None:
            bayes_data += standardized_mean

        return bayes_data
