import geopandas as gpd
import pandas as pd
import os
import numpy as np
import sys
import itertools
import ast
import math
from scipy import stats


def main():
    '''
    Traffic speed assignment script
    vehicle_id, edge_path, time_stamp
    '''
    data_path,calc_path,output_path = load_config()['paths']['data'],load_config()['paths']['calc'],load_config()['paths']['output']
    
    edges_in = os.path.join(output_path, 'transport cvts analysis', 'results', 'traffic_count','road_network.shp')
    routes_in = os.path.join(output_path, 'transport cvts analysis', 'results', 'routes_collected','routes.csv')


    edges = gpd.read_file(edges_in)
    edges.columns = map(str.lower, edges.columns)
    # get the right linelength
    edges['length'] = edges.geometry.apply(line_length)
    length_attr = list(zip(edges['g_id'].values.tolist(),edges['length'].values.tolist()))

    routes_df = pd.read_csv(routes_in)
    
    edge_speeds = []
    for iter_,vals in routes_df.iterrows():
    	edge_path = ast.literal_eval(vals['edge_path'])
    	time_stamp = ast.literal_eval(vals['time_stamp'])
    	if len(edge_path) > 1:
    		for e in range(len(edge_path)-1):
    			time_diff = 1.0*(time_stamp[e+1] - time_stamp[e])
    			if time_diff > 0:
    				distance = sum([l[1] for l in length_attr if l[0] in (edge_path[e],edge_path[e+1])])
    				edge_l = [l[1] for l in length_attr if l[0] == edge_path[e]] + [l[1] for l in length_attr if l[0] == edge_path[e+1]]
    				speed = 3600.0*distance/time_diff
    				if speed >= 20 and speed <= 120:
    					edge_speeds.append((edge_path[e],speed))
    					edge_speeds.append((edge_path[e+1],speed))

    	print ('Done with iteration',iter_)

    del routes_df

    edge_speeds_df = pd.DataFrame(edge_speeds,columns = ['g_id','speed'])

    edge_speeds_df_min = edge_speeds_df.groupby(['g_id'])['speed'].min().reset_index()
    edge_speeds_df_min.rename(columns={'speed': 'min_speed'}, inplace=True)
    edges = pd.merge(edges,edge_speeds_df_min,how='left', on=['g_id']).fillna(0)
    del edge_speeds_df_min

    edge_speeds_df_max = edge_speeds_df.groupby(['g_id'])['speed'].max().reset_index()
    edge_speeds_df_max.rename(columns={'speed': 'max_speed'}, inplace=True)
    edges = pd.merge(edges,edge_speeds_df_max,how='left', on=['g_id']).fillna(0)
    del edge_speeds_df_max

    edge_speeds_df_median = edge_speeds_df.groupby(['g_id'])['speed'].median().reset_index()
    edge_speeds_df_median.rename(columns={'speed': 'md_speed'}, inplace=True)
    edges = pd.merge(edges,edge_speeds_df_median,how='left', on=['g_id']).fillna(0)
    del edge_speeds_df_median

    edge_speeds_df_mean = edge_speeds_df.groupby(['g_id'])['speed'].mean().reset_index()
    edge_speeds_df_mean.rename(columns={'speed': 'mean_speed'}, inplace=True)
    edges = pd.merge(edges,edge_speeds_df_mean,how='left', on=['g_id']).fillna(0)
    del edge_speeds_df_mean


    edge_speeds_df_std = edge_speeds_df.groupby(['g_id'])['speed'].std().reset_index()
    edge_speeds_df_std.rename(columns={'speed': 'std_speed'}, inplace=True)
    edges = pd.merge(edges,edge_speeds_df_std,how='left', on=['g_id']).fillna(0)
    del edge_speeds_df_std
    del edge_speeds_df

    edges.loc[edges['est_speed'] > 120.0,'est_speed'] = 120.0

    edges.to_file(edges_in)

if __name__ == '__main__':
    main()