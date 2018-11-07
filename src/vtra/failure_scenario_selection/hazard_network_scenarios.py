# -*- coding: utf-8 -*-
"""Utility functions for combined impact analysis
"""
import pandas as pd


def combine_hazards_and_network_attributes_and_impacts(hazard_dataframe, network_dataframe):
    hazard_dataframe.loc[hazard_dataframe['probability'] == 'none', 'probability'] = 1.0
    hazard_dataframe['probability'] = pd.to_numeric(hazard_dataframe['probability'])
    hazard_dataframe.rename(columns={
        'length': 'exposure_length',
        'min_val': 'min_flood_depth',
        'max_val': 'max_flood_depth'
    }, inplace=True)

    network_dataframe.rename(columns={'length': 'road_length'}, inplace=True)
    network_dataframe['road_length'] = 1000.0*network_dataframe['road_length']

    all_edge_fail_scenarios = pd.merge(hazard_dataframe, network_dataframe, on=[
        'edge_id'], how='left').fillna(0)

    all_edge_fail_scenarios['percent_exposure'] = 100.0 * \
        all_edge_fail_scenarios['exposure_length']/all_edge_fail_scenarios['road_length']

    del hazard_dataframe, network_dataframe

    return all_edge_fail_scenarios


def create_hazard_scenarios_for_adaptation(all_edge_fail_scenarios, index_cols, length_thr):
    all_edge_fail_scenarios = all_edge_fail_scenarios.set_index(index_cols)
    scenarios = list(set(all_edge_fail_scenarios.index.values.tolist()))
    print('Number of failure scenarios', len(scenarios))
    scenarios_list = []
    for sc in scenarios:
        min_height = max(all_edge_fail_scenarios.loc[[sc], 'min_flood_depth'].values.tolist())
        max_height = max(all_edge_fail_scenarios.loc[[sc], 'max_flood_depth'].values.tolist())
        min_band_num = min(all_edge_fail_scenarios.loc[[sc], 'band_num'].values.tolist())
        max_band_num = max(all_edge_fail_scenarios.loc[[sc], 'band_num'].values.tolist())
        prob = all_edge_fail_scenarios.loc[[sc], 'probability'].values
        if len(list(set(prob))) > 1:
            exposure_len = all_edge_fail_scenarios.loc[[sc], 'exposure_length'].values
            per = all_edge_fail_scenarios.loc[[sc], 'percent_exposure'].values

            prob_tup = list(zip(prob, exposure_len, per))
            u_pr = sorted(list(set(prob.tolist())))
            exposure_len = []
            per = []
            r_wt = []
            for pr in u_pr:
                per_exp = sum([z for (x, y, z) in prob_tup if x == pr])
                exp_len = sum([y for (x, y, z) in prob_tup if x == pr])
                if per_exp > 100.0:
                    exposure_len.append(100.0*exp_len/per_exp)
                    per.append(100.0)
                    r_wt.append(1.0)
                else:
                    exposure_len.append(exp_len)
                    per.append(per_exp)
                    if exp_len < length_thr:
                        r_wt.append(0.01*per_exp)
                    else:
                        r_wt.append(1.0)

            max_exposure_len = max(exposure_len)
            min_exposure_len = min(exposure_len)

            min_per = min(per)
            max_per = max(per)
            min_dur = 0.01*min_per
            max_dur = 0.01*max_per
            risk_wt = 0
            dam_wt = 0
            for p in range(len(u_pr)-1):
                risk_wt += 0.5*(u_pr[p+1]-u_pr[p])*(r_wt[p+1]+r_wt[p])
                dam_wt += 0.5*(u_pr[p+1]-u_pr[p])*(exposure_len[p+1]+exposure_len[p])

        else:
            prob_wt = prob[0]
            min_exposure_len = sum(
                all_edge_fail_scenarios.loc[[sc], 'exposure_length'].values.tolist())
            min_per = sum(all_edge_fail_scenarios.loc[[sc], 'percent_exposure'].values.tolist())
            if min_per > 100.0:
                min_exposure_len = 100.0*min_exposure_len/min_per
                min_per = 100.0

            max_per = min_per
            max_exposure_len = min_exposure_len
            dam_wt = max_exposure_len
            min_dur = 0.01*min_per
            if max_exposure_len < length_thr:
                max_dur = 0.01*max_per
                risk_wt = 0.01*max_per*prob_wt
            else:
                max_dur = 1.0
                risk_wt = prob_wt

        scenarios_list.append(list(sc) + [min_band_num, max_band_num, min_height, max_height,
                                          min_per, max_per, min_dur, max_dur, min_exposure_len,
                                          max_exposure_len, risk_wt, dam_wt])

    new_cols = ['min_band', 'max_band', 'min_height', 'max_height', 'min_exposure_percent',
                'max_exposure_percent', 'min_duration_wt', 'max_duration_wt',
                'min_exposure_length', 'max_exposure_length', 'risk_wt', 'dam_wt']
    scenarios_df = pd.DataFrame(scenarios_list, columns=index_cols + new_cols)

    del all_edge_fail_scenarios, scenarios_list
    return scenarios_df
