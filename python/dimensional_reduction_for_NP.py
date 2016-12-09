import logging
import os
import itertools
# from scipy.spatial.distance import pdist, squareform
# from scipy import exp
# from scipy.linalg import eigh
import numpy as np
import pandas as pd
# from sklearn.decomposition import PCA, SparsePCA, KernelPCA
# from sklearn.preprocessing import StandardScaler
import matplotlib
matplotlib.use('Agg') # for 'batch' mode
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 22}
matplotlib.rc('font', **font)
import matplotlib.pyplot as plt
from scipy.optimize import minimize, basinhopping, differential_evolution
import seaborn as sns

# ----------------------------------------------------
# Global Optimization Accept Functions
# ----------------------------------------------------
class MaxAcceptSmart (object):
    def __init__ (self, n_components, original_cov, threshold=30.0):
        self.xmax = None
        self.n_components = n_components
        self.original_cov = original_cov
        self.threshold = 0.01*threshold
    def __call__ (self, **kwargs):
        # assumes n_components >= 3 and that first component is bin-uncorrelated
        # make into proper matrix
        x = kwargs["x_new"].reshape((self.n_components, -1))
        x_uncor = x[0]
        x_cor = x[1:]
        # compute the covaraiance matrix
        result = np.zeros((x.shape[1], x.shape[1]))
        for row in x_cor:
            result += row[:, np.newaxis]*row
        result += np.diagflat(np.square(x_uncor))
        # NOTE: max of relative difference
        new_max = np.max(np.fabs((result - self.original_cov)/self.original_cov))
        if self.xmax is None or (new_max < self.xmax and new_max < self.threshold):
            self.xmax = new_max
            return True
        else:
            return False

# ----------------------------------------------------
# Optimization Metrics
# ----------------------------------------------------
# basic functions
def covariance_one_uncor (x):
    # assumes n_components >= 3 and that first component is bin-uncorrelated
    x_uncor = x[0]
    x_cor = x[1:]
    # compute the covaraiance matrix
    result = np.zeros((x.shape[1], x.shape[1]))
    for row in x_cor:
        result += row[:, np.newaxis]*row
    result += np.diagflat(np.square(x_uncor))
    return result

def covariance_one_uncor_exclusive (x):
    # assumes n_components >= 3, that first component is bin-uncorrelated, and
    # all components after the first only represent the off diagonal elements
    x_uncor = x[0]
    x_cor = x[1:]
    # compute the covaraiance matrix
    result = np.zeros((x.shape[1], x.shape[1]))
    for row in x_cor:
        result += row[:, np.newaxis]*row
    for i in np.arange(result.shape[0]):
        result[i, i] = 0.0
    result += np.diagflat(np.square(x_uncor))
    return result

# metrics
def squared_diff (x, original_cov, n_components):
    # make into proper matrix
    x = x.reshape((n_components, -1))
    # compute the covaraiance matrix
    result = np.zeros((x.shape[1], x.shape[1]))
    for row in x:
        result += row[:, np.newaxis]*row
    # return the sum of squared differences
    return np.sum(np.square(result - original_cov))

def squared_diff_trace (x, original_cov, n_components):
    # make into proper matrix
    x = x.reshape((n_components, -1))
    # compute the covaraiance matrix
    result = np.zeros((x.shape[1], x.shape[1]))
    for row in x:
        result += row[:, np.newaxis]*row
    # return the trace of squared differences
    return np.trace(np.square(result - original_cov))

def squared_diff_smart (x, original_cov, n_components):
    # make into proper matrix
    x = x.reshape((n_components, -1))
    result = covariance_one_uncor(x)
    # return the sum of squared (relative) differences
    # return np.sum(np.square(100*(result - original_cov)/original_cov))/original_cov.size
    return np.sum(np.square(result - original_cov))

def squared_diff_off_diag_sep (x, original_cov, n_components):
    # make into proper matrix
    x = x.reshape((n_components, -1))
    result = covariance_one_uncor_exclusive(x)
    # return the sum of squared (relative) differences
    # return np.sum(np.square((result - original_cov)/original_cov))
    return np.sum(np.square(result - original_cov))

# ----------------------------------------------------
# Utility functions
# ----------------------------------------------------
def make_covariance_df (df):
    """
    Computes the covariance DataFrames for the given systematic up/down
    variations. Assumes that the input DataFrames represent the systematic
    deltas with each bin as a column.

    Returns the covariance matrices as DataFrames
    """
    raw = df.as_matrix()
    cov_raw = np.zeros((raw.shape[1], raw.shape[1]))
    for row in raw:
        cov_raw += row[:, np.newaxis]*row

    df_cov = pd.DataFrame(data=cov_raw, index=list(df.columns.values), columns=list(df.columns.values))

    return df_cov

def make_covariance_df_reduced (df, metric):
    """
    Computes the covariance DataFrames for the given systematic up/down
    variations. Assumes that the input DataFrames represent the systematic
    deltas with each bin as a column.

    Returns the covariance matrices as DataFrames
    """
    raw = df.as_matrix()
    if metric is squared_diff_smart:
        cov_raw = covariance_one_uncor(raw)
    elif metric is squared_diff_off_diag_sep:
        cov_raw = covariance_one_uncor_exclusive(raw)

    df_cov = pd.DataFrame(data=cov_raw, index=list(df.columns.values), columns=list(df.columns.values))

    return df_cov

def make_sliced_covariance_matrix_plots (fig, labels_list, df, output_dir, filetag=""):
    for c in labels_list:
        ls = labels_list[c]
        df_reduced = df.loc[:, ls]
        if len(list(df_reduced.columns)) < 2: continue

        # make bin labels more streamlined
        total_title = ""
        for c_ in c:
            if len(total_title) != 0: total_title += ':'+c_
            else: total_title = c_
            after = False
            if c_+':' in df_reduced.columns[0]:
                after = True
            df_reduced.columns = [x.replace(c_+':', '') if after else x.replace(':'+c_, '') for x in df_reduced.columns]

        df_cov_reduced = make_covariance_df(df_reduced)

        annot = len(list(df_cov_reduced.columns)) <= 12
        ax = sns.heatmap(df_cov_reduced[::-1], annot=annot, cmap="coolwarm", linewidths=.5)
        ax.set_title(total_title.replace(':', ', '))
        if not filetag == "":
            fig.savefig(os.path.join(output_dir, filetag+'_cov_'+total_title.replace('=', '').replace(':', '_')+'.pdf'))
        else:
            fig.savefig(os.path.join(output_dir, 'cov_'+total_title.replace('=', '').replace(':', '_')+'.pdf'))
        plt.clf()

def make_sliced_original_covariance_matrix_plots (fig, labels_list, df, output_dir, filetag=""):
    for c in labels_list:
        ls = labels_list[c]
        df_reduced = df.loc[ls, ls]
        if len(list(df_reduced.columns)) < 2: continue

        # make bin labels more streamlined
        total_title = ""
        for c_ in c:
            if len(total_title) != 0: total_title += ':'+c_
            else: total_title = c_
            after = False
            if c_+':' in df_reduced.columns[0]:
                after = True
            df_reduced.columns = [x.replace(c_+':', '') if after else x.replace(':'+c_, '') for x in df_reduced.columns]
            df_reduced.index = [x.replace(c_+':', '') if after else x.replace(':'+c_, '') for x in df_reduced.index]
        # if c+':' in df_reduced.columns[0]:
        #     after = True
        # df_reduced.columns = [x.replace(c+':', '') if after else x.replace(':'+c, '') for x in df_reduced.columns]
        # df_reduced.index = [x.replace(c+':', '') if after else x.replace(':'+c, '') for x in df_reduced.index]

        annot = len(list(df_reduced.columns)) <= 12
        ax = sns.heatmap(df_reduced[::-1], annot=annot, cmap="coolwarm", linewidths=.5)
        ax.set_title(total_title.replace(':', ', '))
        if not filetag == "":
            fig.savefig(os.path.join(output_dir, filetag+'_cov_'+total_title.replace('=', '').replace(':', '_')+'.pdf'))
        else:
            fig.savefig(os.path.join(output_dir, 'cov_'+total_title.replace('=', '').replace(':', '_')+'.pdf'))
        plt.clf()

def make_sliced_reduced_covariance_matrix_plots (metric, fig, labels_list, df, output_dir, filetag=""):
    for c in labels_list:
        ls = labels_list[c]
        df_reduced = df.loc[:, ls]
        if len(list(df_reduced.columns)) < 2: continue

        # make bin labels more streamlined
        total_title = ""
        for c_ in c:
            if len(total_title) != 0: total_title += ':'+c_
            else: total_title = c_
            after = False
            if c_+':' in df_reduced.columns[0]:
                after = True
            df_reduced.columns = [x.replace(c_+':', '') if after else x.replace(':'+c_, '') for x in df_reduced.columns]
        # after = False
        # if c+':' in df_reduced.columns[0]:
        #     after = True
        # df_reduced.columns = [x.replace(c+':', '') if after else x.replace(':'+c, '') for x in df_reduced.columns]

        df_cov_reduced = make_covariance_df_reduced(df_reduced, metric)

        annot = len(list(df_cov_reduced.columns)) <= 12
        ax = sns.heatmap(df_cov_reduced[::-1], annot=annot, cmap="coolwarm", linewidths=.5)
        ax.set_title(total_title.replace(':', ', '))
        if not filetag == "":
            fig.savefig(os.path.join(output_dir, filetag+'_cov_'+total_title.replace('=', '').replace(':', '_')+'.pdf'))
        else:
            fig.savefig(os.path.join(output_dir, 'cov_'+total_title.replace('=', '').replace(':', '_')+'.pdf'))
        plt.clf()

def make_sliced_diff_covariance_matrix_plots(fig, labels_list, df, percent_range_min, percent_range_max, output_dir, filetag=""):
    for c in labels_list:
        ls = labels_list[c]
        df_reduced = df.loc[ls, ls]
        if len(list(df_reduced.columns)) < 2: continue

        # make bin labels more streamlined
        total_title = ""
        for c_ in c:
            if len(total_title) != 0: total_title += ':'+c_
            else: total_title = c_
            after = False
            if c_+':' in df_reduced.columns[0]:
                after = True
            df_reduced.columns = [x.replace(c_+':', '') if after else x.replace(':'+c_, '') for x in df_reduced.columns]
            df_reduced.index = [x.replace(c_+':', '') if after else x.replace(':'+c_, '') for x in df_reduced.index]

        annot = len(list(df_reduced.columns)) <= 12
        ax = sns.heatmap(df_reduced[::-1], annot=annot, cmap="RdBu_r", center=0.0, linewidths=.5, vmin=percent_range_min, vmax=percent_range_max)
        ax.set_title("Percent Difference ({})".format(total_title.replace(':', ', ')))
        if not filetag == "":
            fig.savefig(os.path.join(output_dir, filetag+'_cov_'+total_title.replace('=', '').replace(':', '_')+'.pdf'))
        else:
            fig.savefig(os.path.join(output_dir, 'cov_'+total_title.replace('=', '').replace(':', '_')+'.pdf'))
        plt.clf()

def make_sliced_diff_covariance_matrix_compare_plots(fig, labels_list, df1, df2, percent_range_min, percent_range_max, left_title, right_title, output_dir, filetag=""):
    for c in labels_list:
        ls = labels_list[c]
        df1_reduced = df1.loc[ls, ls]
        df2_reduced = df2.loc[ls, ls]
        if len(list(df1_reduced.columns)) < 2: continue

        # make bin labels more streamlined
        total_title = ""
        for c_ in c:
            if len(total_title) != 0: total_title += ':'+c_
            else: total_title = c_
            after = False
            if c_+':' in df1_reduced.columns[0]:
                after = True
            df1_reduced.columns = [x.replace(c_+':', '') if after else x.replace(':'+c_, '') for x in df1_reduced.columns]
            df1_reduced.index = [x.replace(c_+':', '') if after else x.replace(':'+c_, '') for x in df1_reduced.index]
            df2_reduced.columns = [x.replace(c_+':', '') if after else x.replace(':'+c_, '') for x in df2_reduced.columns]
            df2_reduced.index = [x.replace(c_+':', '') if after else x.replace(':'+c_, '') for x in df2_reduced.index]

        percent_range_min, percent_range_max = np.min([np.min(df1_reduced), np.min(df2_reduced)]), np.max([np.max(df1_reduced), np.max(df2_reduced)])
        annot = len(list(df1_reduced.columns)) <= 12
        fig, axn = plt.subplots(1, 2, sharex=True, sharey=True)
        cbar_ax = fig.add_axes([.91, .1, .03, .8])
        cbar_ax.tick_params(labelsize=20)
        fig.set_size_inches(19.25, 8.5, forward=True)
        for i, ax in enumerate(axn.flat):
            sns.heatmap(df1_reduced[::-1] if i == 0 else df2_reduced[::-1], ax=ax,
                        cbar=i == 0,
                        center=0.0,
                        annot=annot, annot_kws={"size": 18}, cmap="RdBu_r",
                        linewidths=.5, vmin=percent_range_min, vmax=percent_range_max,
                        cbar_ax=None if i else cbar_ax)
            ax.tick_params(labelsize=16)
            ax.set_aspect('equal', 'datalim')
            ax.set_title(left_title if i == 0 else right_title, fontsize=20)

        plt.suptitle("Percent Difference Comparison ({})".format(total_title.replace(':', ', ')), fontsize=20)

        if not filetag == "":
            fig.savefig(os.path.join(output_dir, filetag+'_cov_'+total_title.replace('=', '').replace(':', '_')+'.pdf'))
        else:
            fig.savefig(os.path.join(output_dir, 'cov_'+total_title.replace('=', '').replace(':', '_')+'.pdf'))
        plt.clf()

def make_exp_var_plot (explained_variance_ratio, cum_explained_variance_ratio, output_dir, file_name):
    plt.step( range( 1,cum_explained_variance_ratio.shape[0] + 1), cum_explained_variance_ratio, where ='mid',
              label ='cumulative explained variance')
    plt.bar( range( 1,explained_variance_ratio.shape[0] + 1), explained_variance_ratio, alpha = 0.5, align ='center',
             label ='individual explained variance')
    plt.ylabel('Explained variance ratio')
    plt.xlabel('Principal components')
    plt.legend( loc ='best')
    plt.savefig(os.path.join(output_dir, file_name))
    plt.clf()

def make_explained_variance_plots_for_optimization (x, output_dir, file_name):
    lambdas = np.linalg.norm(x, axis=1)
    explained_variance_ratio = lambdas/lambdas.sum()
    order = explained_variance_ratio.argsort()[::-1]
    explained_variance_ratio = explained_variance_ratio[order]
    cum_explained_variance_ratio = np.cumsum(explained_variance_ratio)
    make_exp_var_plot(explained_variance_ratio, cum_explained_variance_ratio, output_dir, file_name)

def main():
    logging.basicConfig(#filename='NP_reduction.log',
                        level=logging.DEBUG)
    # NOTE: optimal (through trial and error) settings:
    #       region_split=True, n_components=4, metric=squared_diff_smart, Continuous_240316_MV2c20_AntiKt4EMTopoJets_continuous_B_IntegratedWP_7TagBins_SF
    #       region_split=False, n_components=3, metric=squared_diff_smart, 2016-Winter-13TeV-MC15-CDI-March4_unofficial_MV2c20_AntiKt4EMTopoJets_FixedCutBEff_60_B_default_SF
    n_components = 11
    n_EV_components = -1
    output_dir = 'NP_reduction_plots'
    global_optimization = True
    T = 2
    comprehensive_plots = True
    region_split = False
    random_seed = 123
    # metric = squared_diff
    metric = squared_diff_smart
    # metric = squared_diff_off_diag_sep
    # NOTE: minimization methods: None ('BFGS'?), 'Nelder-Mead', 'Powell', 'CG', 'BFGS'
    minimizer_kwargs = {'method': None,
                        'tol': 1e-9,
                        'options': {'disp': True}}
    fig, ax = plt.subplots()

    # NOTE: Light 60
    # file_name_core = "2016-Winter-13TeV-MC15-CDI-March4_unofficial_MV2c20_AntiKt4EMTopoJets_FixedCutBEff_60_Light_default_SF"
    # NOTE: Light 70
    # file_name_core = "2016-Winter-13TeV-MC15-CDI-March4_unofficial_MV2c20_AntiKt4EMTopoJets_FixedCutBEff_70_Light_default_SF"
    # file_name_core = "2016-Winter-13TeV-MC15-CDI-March4_unofficial_order_1_smoothing_0.4_ptbins_100_MV2c20_AntiKt4EMTopoJets_FixedCutBEff_70_Light_default_SF"
    # NOTE: Light Continuous
    # file_name_core = "Continuous_240316_MV2c20_AntiKt4EMTopoJets_continuous_Light_IntegratedWP_5TagBins_SF"
    # file_name_core = "Continuous_240316_order_1_smoothing_0.4_ptbins_25_MV2c20_AntiKt4EMTopoJets_continuous_Light_IntegratedWP_5TagBins_SF"
    # NOTE: C 60
    # file_name_core = "2016-Winter-13TeV-MC15-CDI-March4_unofficial_MV2c20_AntiKt4EMTopoJets_FixedCutBEff_60_C_default_SF"
    # file_name_core = "2016-Winter-13TeV-MC15-CDI-March4_unofficial_order_1_smoothing_0.4_ptbins_100_MV2c20_AntiKt4EMTopoJets_FixedCutBEff_60_C_default_SF"
    # NOTE: C 70
    # file_name_core = "2016-Winter-13TeV-MC15-CDI-March4_unofficial_MV2c20_AntiKt4EMTopoJets_FixedCutBEff_70_C_default_SF"
    # file_name_core = "2016-Winter-13TeV-MC15-CDI-March4_unofficial_order_1_smoothing_0.4_ptbins_100_MV2c20_AntiKt4EMTopoJets_FixedCutBEff_70_C_default_SF"
    # NOTE: C Continuous
    # file_name_core = "2016-Winter-13TeV-MC15-CDI-March4_unofficial_MV2c20_AntiKt4EMTopoJets_FixedCutBEff_60_C_default_SF"
    # NOTE: B 60
    # file_name_core = "2016-Winter-13TeV-MC15-CDI-March4_unofficial_MV2c20_AntiKt4EMTopoJets_FixedCutBEff_60_B_default_SF"
    # file_name_core = "2016-Winter-13TeV-MC15-CDI-March4_unofficial_order_1_smoothing_0.4_ptbins_100_MV2c20_AntiKt4EMTopoJets_FixedCutBEff_60_B_default_SF"
    # NOTE: B 70
    # file_name_core = "2016-Winter-13TeV-MC15-CDI-March4_unofficial_MV2c20_AntiKt4EMTopoJets_FixedCutBEff_70_B_default_SF"
    # NOTE: B Continuous
    file_name_core = "Continuous_240316_MV2c20_AntiKt4EMTopoJets_continuous_B_IntegratedWP_5TagBins_SF"
    # file_name_core = "Continuous_240316_MV2c20_AntiKt4EMTopoJets_continuous_B_IntegratedWP_7TagBins_SF"

    # input_file_name = "{1}_n{0}_systs.csv".format(n_EV_components, file_name_core)
    input_file_name_matrix = "{1}_n{0}_matrix.csv".format(n_EV_components, file_name_core)
    input_file_name_eigen_vec_up = "{1}_n{0}_eigen_vec_up.csv".format(n_EV_components, file_name_core)
    input_file_name_eigen_vec_down = "{1}_n{0}_eigen_vec_down.csv".format(n_EV_components, file_name_core)
    output_dir += "_n{}_{}".format(n_components, file_name_core)

    smoothed = "_order_" in file_name_core
    if smoothed:
        logging.info('Detected smoothed calibration. Will set global_optimization flag to True.')
        global_optimization = True
        minimizer_kwargs['options']['disp'] = False
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # df_systs = pd.read_csv(input_file_name, index_col=0)
    df_matrix = pd.read_csv(input_file_name_matrix, index_col=0)
    df_eigen_vec_up = pd.read_csv(input_file_name_eigen_vec_up, index_col=0)
    df_eigen_vec_down = pd.read_csv(input_file_name_eigen_vec_down, index_col=0)

    df_columns = df_eigen_vec_up.columns
    df_columns_list = list(df_eigen_vec_up.columns.values)

    # df_mean = df.loc[["result" == x for x in list(df.index)], :]
    # df_up = df.loc[["__up" in x for x in list(df.index)], :]
    # df_down = df.loc[["__down" in x for x in list(df.index)], :]
    #
    # df_up_cov, df_down_cov = make_covariance_df(df_up, df_down, smoothed)

    logging.info('Starting dimensional reduction analysis for {}'.format(file_name_core))

    n_indentifiers = len(df_columns_list[0].split(':'))
    bin_indentifiers = set([_ for x in df_columns_list for _ in x.split(':')])
    # dictionary of unique bin identifier
    labels_list = {}
    split_list = {}
    for comb in itertools.combinations(bin_indentifiers, n_indentifiers - 1):
        split_list[comb] = [x for x in df_columns_list if set(comb).issubset(set(x.split(':')))]
        labels_list[comb] = split_list[comb][:] # copy list
        should_delete = False
        # NOTE: remove all pT slices
        for c in comb:
            if "pt" in c:
                should_delete = True
                break
        test_length = len(split_list[comb])
        if test_length == 0:
            del split_list[comb]
            del labels_list[comb]
        if test_length != 0 and not comprehensive_plots and should_delete:
            del labels_list[comb]
        if test_length != 0 and region_split and should_delete:
            del split_list[comb]

    # make "sliced" covariance matrix plots
    # if comprehensive_plots:
    logging.info("Making covaraince matrix slice plots for complete matrix (may take a while).")
    make_sliced_original_covariance_matrix_plots(fig, labels_list, df_matrix, output_dir, 'original')

    # ----------------------------------------------------
    # Optimization Technique
    # ----------------------------------------------------
    logging.info("Starting optimization ...")
    np.random.seed(random_seed)
    if region_split:
        logging.info("Splitting regions during optimization - watch out for orthogonal correlations")
        res_df = pd.DataFrame(columns=df_matrix.columns, index=["NP"+str(i + 1) for i in range(n_components)], dtype='float32')
        for c in split_list:
            ls = split_list[c]
            df_matrix_reduced = df_matrix.loc[ls, ls]
            # initial value of minimization (random numbers on interval [-1, 1])
            # np.random.seed(random_seed)
            x0 = 2*np.random.random_sample((n_components*len(df_matrix_reduced.columns),)) - 1
            if global_optimization:
                if metric is squared_diff_smart: accept = MaxAcceptSmart(n_components, df_matrix_reduced.as_matrix())
                else: accept = None
                # np.random.seed(random_seed)
                res = basinhopping(lambda x: metric(x, df_matrix_reduced.as_matrix(), n_components), x0, T=1.0 if T is None else T, accept_test=accept, disp=True, minimizer_kwargs=minimizer_kwargs, niter_success=4, niter=200)
            else:
                # np.random.seed(random_seed)
                res = minimize(lambda x: metric(x, df_matrix_reduced.as_matrix(), n_components), x0, method=minimizer_kwargs["method"], tol=minimizer_kwargs["tol"], options=minimizer_kwargs['options'])
                if res.success:
                    logging.info("Success!")
                else:
                    logging.warning("Failed with message: {}".format(res.message))

            res_df.loc[:, ls] = res.x.reshape((n_components, -1))

        make_explained_variance_plots_for_optimization(res_df.as_matrix(), output_dir, 'reduced_exp_var_up.pdf')
        res_df = pd.DataFrame(data=res_df.as_matrix(), columns=df_matrix.columns, index=["NP"+str(i + 1) for i in range(n_components)])
    else:
        # initial value of minimization (random numbers on interval [-1, 1])
        # np.random.seed(random_seed)
        x0 = 2*np.random.random_sample((n_components*len(df_matrix.columns),)) - 1
        if global_optimization:
            if metric is squared_diff_smart: accept = MaxAcceptSmart(n_components, df_matrix.as_matrix())
            else: accept = None
            # np.random.seed(random_seed)
            res = basinhopping(lambda x: metric(x, df_matrix.as_matrix(), n_components), x0, T=1.0 if T is None else T, accept_test=accept, disp=True, minimizer_kwargs=minimizer_kwargs, niter_success=4, niter=200)
        else:
            # np.random.seed(random_seed)
            res = minimize(lambda x: metric(x, df_matrix.as_matrix(), n_components), x0, method=minimizer_kwargs["method"], tol=minimizer_kwargs["tol"], options=minimizer_kwargs['options'])
            if res.success:
                logging.info("Success!")
            else:
                logging.warning("Failed with message: {}".format(res.message))

        make_explained_variance_plots_for_optimization(res.x.reshape((n_components, -1)), output_dir, 'reduced_exp_var_up.pdf')
        res_df = pd.DataFrame(data=res.x.reshape((n_components, -1)), columns=df_matrix.columns, index=["NP"+str(i + 1) for i in range(n_components)])

    res_df_cov = make_covariance_df_reduced(res_df, metric)

    # if comprehensive_plots:
    logging.info("Making covaraince matrix slice plots for reduced variations (may take a while).")
    make_sliced_reduced_covariance_matrix_plots(metric, fig, labels_list, res_df, output_dir, 'reduced')

    res_percent_diff = 100*(res_df_cov.as_matrix() - df_matrix.as_matrix())/df_matrix.as_matrix()
    res_percent_diff_min, res_percent_diff_max = np.min(res_percent_diff), np.max(res_percent_diff)

    # ----------------------------------------------------
    # Eigenvalue Decomposition Technique
    # ----------------------------------------------------
    components = df_eigen_vec_up.as_matrix()
    evals = np.linalg.norm(components, axis=1)
    order = evals.argsort()[::-1]
    evals = evals[order]
    components = components[order]

    explained_variance_ratio = evals / np.sum(evals)
    cum_explained_variance_ratio = np.cumsum(explained_variance_ratio)
    explained_variance_ratio = explained_variance_ratio[:n_components]
    cum_explained_variance_ratio = cum_explained_variance_ratio[:n_components]
    make_exp_var_plot(explained_variance_ratio, cum_explained_variance_ratio, output_dir, 'eigen_exp_var_up.pdf')

    eig_up_df = pd.DataFrame(data=components[:n_components], columns=df_columns, index=["NP"+str(i + 1) for i in range(components[:n_components].shape[0])])

    components = df_eigen_vec_down.as_matrix()
    evals = np.linalg.norm(components, axis=1)
    order = evals.argsort()[::-1]
    evals = evals[order]
    components = components[order]

    explained_variance_ratio = evals / np.sum(evals)
    cum_explained_variance_ratio = np.cumsum(explained_variance_ratio)
    explained_variance_ratio = explained_variance_ratio[:n_components]
    cum_explained_variance_ratio = cum_explained_variance_ratio[:n_components]
    make_exp_var_plot(explained_variance_ratio, cum_explained_variance_ratio, output_dir, 'eigen_exp_var_down.pdf')

    eig_down_df = pd.DataFrame(data=components[:n_components], columns=df_columns, index=["NP"+str(i + 1) for i in range(components[:n_components].shape[0])])

    # if comprehensive_plots:
    logging.info("Making covaraince matrix slice plots for eigen-decomposition (may take a while).")
    make_sliced_covariance_matrix_plots(fig, labels_list, eig_up_df, output_dir, 'eigen_up')
    make_sliced_covariance_matrix_plots(fig, labels_list, eig_down_df, output_dir, 'eigen_down')

    eig_df_up_cov = make_covariance_df(eig_up_df)
    eig_df_down_cov = make_covariance_df(eig_down_df)

    eig_up_percent_diff = 100*(eig_df_up_cov.as_matrix() - df_matrix.as_matrix())/df_matrix.as_matrix()
    eig_up_percent_diff_min, eig_up_percent_diff_max = np.min(eig_up_percent_diff), np.max(eig_up_percent_diff)

    eig_down_percent_diff = 100*(eig_df_down_cov.as_matrix() - df_matrix.as_matrix())/df_matrix.as_matrix()
    eig_down_percent_diff_min, eig_down_percent_diff_max = np.min(eig_down_percent_diff), np.max(eig_down_percent_diff)

    percent_range_max = np.max([eig_up_percent_diff_max, eig_down_percent_diff_max, res_percent_diff_max])
    percent_range_min = np.min([eig_up_percent_diff_min, eig_down_percent_diff_min, res_percent_diff_min])

    df_res_percent_diff = pd.DataFrame(data=res_percent_diff, columns=df_columns, index=df_columns)
    df_eig_up_percent_diff = pd.DataFrame(data=eig_up_percent_diff, columns=df_columns, index=df_columns)
    df_eig_down_percent_diff = pd.DataFrame(data=eig_down_percent_diff, columns=df_columns, index=df_columns)

    # make_sliced_diff_covariance_matrix_plots(fig, labels_list, df_res_percent_diff, percent_range_min, percent_range_max, output_dir, 'reduced_diff')
    # make_sliced_diff_covariance_matrix_plots(fig, labels_list, df_eig_up_percent_diff, percent_range_min, percent_range_max, output_dir, 'eigne_up_diff')
    # make_sliced_diff_covariance_matrix_plots(fig, labels_list, df_eig_down_percent_diff, percent_range_min, percent_range_max, output_dir, 'eigne_down_diff')

    logging.info("Making relative difference plots for eigen-decomposition vs. optimization (may take a while).")
    make_sliced_diff_covariance_matrix_compare_plots(fig, labels_list, df_res_percent_diff, df_eig_up_percent_diff, percent_range_min, percent_range_max, "optimized", "eigen-decomposition", output_dir, 'reduced_eigen_up_diff')

    res_df.to_csv(os.path.join(output_dir, "optimal_NP.csv"))

if __name__ == "__main__":
    main()

# def rbf_kernel_pca( X, gamma, n_components):
#     """
#     RBF kernel PCA implementation.
#
#     Parameters
#     ------------
#     X: {NumPy ndarray}, shape = [n_samples, n_features]
#
#     gamma: float Tuning parameter of the RBF kernel
#
#     n_components: int Number of principal components to return
#
#     Returns
#     ------------
#     X_pc: {NumPy ndarray}, shape = [n_samples, k_features]
#         Projected dataset
#     """
#     # Calculate pairwise squared Euclidean distances
#     # # in the MxN dimensional dataset.
#     sq_dists = pdist( X, 'sqeuclidean')
#     # Convert pairwise distances into a square matrix.
#     mat_sq_dists = squareform( sq_dists)
#
#     # Compute the symmetric kernel matrix.
#     K = exp(-gamma * mat_sq_dists)
#     # Center the kernel matrix.
#     N = K.shape[0]
#     one_n = np.ones((N, N)) / N
#     K = K - one_n.dot(K) - K.dot(one_n) + one_n.dot(K). dot(one_n)
#
#     # Obtaining eigenpairs from the centered kernel matrix
#     # numpy.eigh returns them in sorted order
#     eigvals, eigvecs = eigh(K)
#     return eigvals, eigvecs
#
#     # # Collect the top k eigenvectors (projected samples)
#     # X_pc = np.column_stack(( eigvecs[:, -i] for i in range( 1, n_components + 1)))
#     # return X_pc

    # scikit learn kernel PCA
    # X_raw = df_up.as_matrix()
    # # print X_raw
    #
    # basis_vectors = np.identity(X_raw.shape[1])
    # # basis_vectors = np.diagflat(df_mean.as_matrix())
    #
    # n_data = 100
    # data = basis_vectors[np.random.randint(basis_vectors.shape[0], size=n_data)]
    #
    # scaler = StandardScaler().fit(X_raw)
    #
    # X = scaler.transform(X_raw)
    # # X = X_raw
    #
    # # print X.mean(axis=0)
    # # print X.std(axis=0)
    #
    # evals, evecs = np.linalg.eigh(np.cov(X.T))
    # order = evals.argsort()[::-1]
    # evals = evals[order]
    # evecs = evecs[:,order]
    # # these are the eigenvectors (as rows) multiplied by their eigenvalues
    # components = evals[:, np.newaxis] * evecs.T
    #
    # for d in data:
    #     syst_vars = scaler.inverse_transform(np.dot(components, d).T[:, np.newaxis] * d) * d[np.newaxis,:]# gives all systematic variations for first data point
    #     # make histograms here
    #     # print syst_vars
    # explained_variance_ratio = evals / np.sum(evals)
    # cum_explained_variance_ratio = np.cumsum(explained_variance_ratio)
    # # print explained_variance_ratio
    # make_exp_var_plot(explained_variance_ratio, cum_explained_variance_ratio, 'pca_exp_var.pdf')
    #
    #
    # # kpca = KernelPCA(kernel="poly", gamma=0.5, coef0=5, degree=3, fit_inverse_transform=True).fit(X.T)
    # # kpca = KernelPCA(kernel="poly", gamma=0.5, coef0=5, degree=3, fit_inverse_transform=True).fit(X)
    # # kpca = KernelPCA(kernel="rbf", gamma=0.005, fit_inverse_transform=True).fit(X)
    # # kpca = KernelPCA(kernel="sigmoid", gamma=1, coef0=5, fit_inverse_transform=True).fit(X)
    # # # X_kpca = kpca.transform(X)
    # # # X_kpca = kpca.transform(X_raw)
    # # # print X.shape
    # # # print scaler.inverse_transform(kpca.inverse_transform([kpca.alphas_[0]]))
    # # # print scaler.inverse_transform(kpca.inverse_transform([kpca.alphas_[1]]))
    # # # print scaler.inverse_transform(kpca.inverse_transform([kpca.alphas_[2]]))
    # # print kpca.alphas_
    # # kpca_explained_variance_ratio = kpca.lambdas_ / kpca.lambdas_.sum()
    # # print kpca_explained_variance_ratio
    # # print kpca.lambdas_
    #
    # # kpca_components = kpca.lambdas_[np.newaxis, :] * kpca.alphas_
    # # print kpca_components.T
    # # print kpca_components.T.sum(axis=1)
    #
    # # kpca_components = kpca.lambdas_[np.newaxis, :] * kpca.alphas_
    # # kpca_components = kpca.lambdas_[:, np.newaxis] * kpca.alphas_
    # #
    # # # transform basis vectors to 'kernel' space
    # # basis_vectors_kpca = kpca.transform(basis_vectors)
    # # # normalize vectors
    # # # basis_vectors_kpca = (basis_vectors_kpca.T / np.linalg.norm(basis_vectors_kpca, axis=1)).T
    # #
    # # kpca_explained_variance_ratio = kpca.lambdas_ / kpca.lambdas_.sum()
    # # kpca_cum_explained_variance_ratio = np.cumsum(kpca_explained_variance_ratio)
    # # make_exp_var_plot(kpca_explained_variance_ratio, kpca_cum_explained_variance_ratio, 'kpca_exp_var.pdf')
    # #
    # # print scaler.inverse_transform(kpca_components.dot(basis_vectors_kpca[0]).T[:, np.newaxis]*basis_vectors[0])*basis_vectors[0][np.newaxis,:]
    #
    # # kpca_explained_variance_ratio = [i/sum(kpca.lambdas_) for i in kpca.lambdas_]
    # # print len(kpca_explained_variance_ratio)
    # # print kpca_explained_variance_ratio
    # # for a in kpca.alphas_:
    # #     print a.dot(a)
    # # for a in basis_vectors_kpca:
    # #     print a.dot(a)
