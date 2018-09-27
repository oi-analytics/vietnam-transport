# -*- coding: utf-8 -*-
"""
Created on Sat Jul 14 15:41:39 2018

@author: elcok
"""

import geopandas as gpd
import pandas as pd
import os
import igraph as ig
import numpy as np
import sys
import subprocess
from shapely.geometry import Point
import itertools
import operator
import ast
import math



from vtra.utils import load_config,extract_value_from_gdf,get_nearest_node,gdf_clip,gdf_geom_clip,count_points_in_polygon
from vtra.transport_network_creation import province_shapefile_to_network, add_igraph_generalised_costs_roads

def netrev_od_pairs(start_points,end_points):
	"""
	Assign net revenue to roads assets in Vietnam

	Inputs are:
	start_points - GeoDataFrame of start points for shortest path analysis.
	end_points - GeoDataFrame of potential end points for shorest path analysis.
	G - iGraph network of the province.
	save_edges -

	Outputs are:
	Shapefile with all edges and the total net reveneu transferred along each edge
	GeoDataFrame of total net revenue transferred along each edge
	"""
	save_paths = []
	for iter_,place in start_points.iterrows():
		try:
			closest_center = end_points.loc[end_points['OBJECTID']
			== place['NEAREST_C_CENTER']]['NEAREST_G_NODE'].values[0]

			# get_od_pair = (place['NEAREST_G_NODE'],closest_center)
			# save_paths.append((str(get_od_pair),1.0*place['netrev_agri']/12.0,1.0*place['netrev_noagri']/12.0))
			save_paths.append((closest_center,place['NEAREST_G_NODE'],1.0*place['netrev_agri']/12.0,1.0*place['netrev_noagri']/12.0))
		except:
			print(iter_)


	# od_pairs_df = pd.DataFrame(save_paths,columns = ['od_nodes','netrev_agri','netrev_noagri'])
	# od_pairs_df = od_pairs_df.groupby(['od_nodes'])['netrev_agri','netrev_noagri'].sum().reset_index()

	od_pairs_df = pd.DataFrame(save_paths,columns = ['origin','destination','netrev_agri','netrev_noagri'])
	od_pairs_df = od_pairs_df.groupby(['origin','destination'])['netrev_agri','netrev_noagri'].sum().reset_index()

	return od_pairs_df

def crop_od_pairs(start_points,end_points,crop_name):
	save_paths = []
	for iter_,place in start_points.iterrows():
		try:
			closest_center = end_points.loc[end_points['OBJECTID']
			== place['NEAREST_C_CENTER']]['NEAREST_G_NODE'].values[0]

			# get_od_pair = (place['NEAREST_G_NODE'],closest_center)
			# save_paths.append((str(get_od_pair),place['tons']))
			save_paths.append((closest_center,place['NEAREST_G_NODE'],place['tons']))
		except:
			print(iter_)


	# od_pairs_df = pd.DataFrame(save_paths,columns = ['od_nodes',crop_name])
	# od_pairs_df = od_pairs_df.groupby(['od_nodes'])[crop_name].sum().reset_index()
	od_pairs_df = pd.DataFrame(save_paths,columns = ['origin','destination',crop_name])
	od_pairs_df = od_pairs_df.groupby(['origin','destination'])[crop_name].sum().reset_index()

	return od_pairs_df

def assign_minmax_rev_costs_crops(x,cost_dataframe,x_cols):
	'''
	crop_code	crop_name	min_cost_perton	max_cost_perton
	'''
	min_croprev = 0
	max_croprev = 0
	cost_list = list(cost_dataframe.itertuples(index=False))
	for cost_param in cost_list:
		if cost_param.crop_code in x_cols:
			min_croprev += 1.0*cost_param.min_cost_perton*x[cost_param.crop_code]
			max_croprev += 1.0*cost_param.max_cost_perton*x[cost_param.crop_code]

	return min_croprev, max_croprev

def assign_monthly_tons_crops(x,rice_prod_dist,x_cols):
	'''
	crop_code	crop_name	min_cost_perton	max_cost_perton
	'''
	min_croptons = 0
	max_croptons = 0
	for x_name in x_cols:
		if x_name == 'rice':
			min_croptons += (1.0*min(rice_prod_dist)*x[x_name])/30.0
			max_croptons += (1.0*max(rice_prod_dist)*x[x_name])/30.0
		else:
			min_croptons += (1.0*x[x_name])/365.0
			max_croptons += (1.0*x[x_name])/365.0

	return min_croptons, max_croptons

def assign_io_rev_costs_crops(x,cost_dataframe,rice_prod_dist,x_cols,ex_rate):
	'''
	crop_code	crop_name	min_cost_perton	max_cost_perton
	'''
	min_croprev = 0
	max_croprev = 0
	cost_list = list(cost_dataframe.itertuples(index=False))
	for cost_param in cost_list:
		if cost_param.crop_code in x_cols:
			if cost_param.crop_code == 'rice':
				min_croprev += (1.0*min(rice_prod_dist)*ex_rate*cost_param.est_net_rev*(x[cost_param.crop_code]/cost_param.tot_tons))/30.0
				max_croprev += (1.0*max(rice_prod_dist)*ex_rate*cost_param.est_net_rev*(x[cost_param.crop_code]/cost_param.tot_tons))/30.0
			else:
				min_croprev += 1.0/365.0*(ex_rate*cost_param.est_net_rev*(x[cost_param.crop_code]/cost_param.tot_tons))
				max_croprev += 1.0/365.0*(ex_rate*cost_param.est_net_rev*(x[cost_param.crop_code]/cost_param.tot_tons))

	return min_croprev, max_croprev

if __name__ == '__main__':

	data_path,calc_path,output_path = load_config()['paths']['data'],load_config()['paths']['calc'],load_config()['paths']['output']

	# provinces to consider
	province_list = ['Lao Cai','Binh Dinh','Thanh Hoa']
	province_terrian = ['mountain','flat','flat']
	# province_list = ['Thanh Hoa']
	district_committe_names = ['district_people_committee_points_lao_cai.shp',
							'district_province_peoples_committee_point_binh_dinh.shp',
							'district_people_committee_points_thanh_hoa.shp']
	commune_committe_names = 'commune_committees_points.shp'
	exchange_rate = 1.05*(1000000/21000)
	growth_scenarios = [(5,'low'),(6.5,'forecast'),(10,'high')]
	base_year = 2016

	shp_output_path = os.path.join(output_path,'flow_mapping_shapefiles')
	# flow_output_excel = os.path.join(output_path,'flow_mapping_paths','province_roads_district_center_flow_ods.xlsx')
	flow_output_excel = os.path.join(output_path,'flow_mapping_paths','province_roads_commune_center_flow_ods.xlsx')
	excl_wrtr = pd.ExcelWriter(flow_output_excel)

	rd_prop_file = os.path.join(data_path,'Roads','road_properties','road_properties.xlsx')
	province_path = os.path.join(data_path,'Vietnam_boundaries','who_boundaries','who_provinces.shp')
	population_points_in = os.path.join(data_path,'Points_of_interest','population_points.shp')
	commune_path = os.path.join(data_path,'Vietnam_boundaries','boundaries_stats','commune_level_stats.shp')

	crop_data_path = os.path.join(data_path,'Agriculture_crops','crop_data')
	rice_month_file = os.path.join(data_path,'rice_atlas_vietnam','rice_production.shp')
	crop_month_fields = ['P_Jan','P_Feb','P_Mar','P_Apr','P_May','P_Jun','P_Jul','P_Aug','P_Sep','P_Oct','P_Nov','P_Dec']
	crop_names = ['rice','cash','cass','teas','maiz','rubb','swpo','acof','rcof','pepp']


	for prn in range(len(province_list)):
	# for prn in range(0,1):
		province_ods_df = []
		province = province_list[prn]
		# set all paths for all input files we are going to use
		province_name = province.replace(' ','').lower()

		edges_in = os.path.join(data_path,'Roads','{}_roads'.format(province_name),'vietbando_{}_edges.shp'.format(province_name))
		nodes_in = os.path.join(data_path,'Roads','{}_roads'.format(province_name),'vietbando_{}_nodes.shp'.format(province_name))

		# commune_center_in = os.path.join(data_path,'Points_of_interest',district_committe_names[prn])
		commune_center_in = os.path.join(data_path,'Points_of_interest',commune_committe_names)

		# path_width_table = os.path.join(data_path,'Roads','road_properties','road_properties.xlsx')

		# load provinces and get geometry of the right province
		provinces = gpd.read_file(province_path)
		provinces = provinces.to_crs({'init': 'epsg:4326'})
		province_geom = provinces.loc[provinces.NAME_ENG == province].geometry.values[0]

		# clip all the populations to the province
		prov_pop = gdf_clip(population_points_in,province_geom)
		prov_commune_center = gdf_clip(commune_center_in,province_geom)
		if 'OBJECTID' not in prov_commune_center.columns.values.tolist():
			prov_commune_center['OBJECTID'] = prov_commune_center.index

		prov_communes = gdf_clip(commune_path,province_geom)

		G = province_shapefile_to_network(edges_in,province_terrian[prn],rd_prop_file)
		nodes_name = np.asarray([x['name'] for x in G.vs])

		del G

		# load nodes of the network
		nodes = gpd.read_file(nodes_in)
		nodes = nodes[nodes['NODE_ID'].isin(nodes_name)].reset_index()
		nodes = nodes.to_crs({'init': 'epsg:4326'})
		sindex_nodes = nodes.sindex

		# get revenue values for each village
		# first create sindex of all villages to count number of villages in commune
		prov_pop_sindex = prov_pop.sindex

		# create new column in prov_communes with amount of villages
		prov_communes['n_villages'] = prov_communes.geometry.apply(lambda x: count_points_in_polygon(x,prov_pop_sindex))
		prov_communes['netrev_village'] = exchange_rate*(prov_communes['netrevenue']*prov_communes['nfirm'])/prov_communes['n_villages']
		# also get the net revenue of the agriculture sector which is called nongnghiep
		prov_communes['netrev_village_agri'] = 1.0/365.0*(prov_communes['nongnghiep']*prov_communes['netrev_village'])
		prov_communes['netrev_village_noagri'] = 1.0/365.0*(prov_communes['netrev_village'] - prov_communes['netrev_village_agri'])


		commune_sindex = prov_communes.sindex
		# give each village a net revenue based on average per village in commune
		prov_pop['netrev_agri'] = prov_pop.geometry.apply(lambda x: extract_value_from_gdf(x,commune_sindex,prov_communes,'netrev_village_agri'))
		prov_pop['netrev_noagri'] = prov_pop.geometry.apply(lambda x: extract_value_from_gdf(x,commune_sindex,prov_communes,'netrev_village_noagri'))


		# get nearest node in network for all start and end points
		prov_pop['NEAREST_G_NODE'] = prov_pop.geometry.apply(lambda x: get_nearest_node(x,sindex_nodes,nodes,'NODE_ID'))
		prov_commune_center['NEAREST_G_NODE'] = prov_commune_center.geometry.apply(lambda x: get_nearest_node(x,sindex_nodes,nodes,'NODE_ID'))

		# prepare for shortest path routing, we'll use the spatial index of the centers
		# to find the nearest center for each population point
		sindex_commune_center = prov_commune_center.sindex
		prov_pop['NEAREST_C_CENTER'] = prov_pop.geometry.apply(lambda x: get_nearest_node(x,sindex_commune_center,prov_commune_center,'OBJECTID'))

		# find all OD pairs of the revenues
		netrev_ods = netrev_od_pairs(prov_pop,prov_commune_center)
		province_ods_df.append(netrev_ods)

		# find the crop production months for the province
		rice_prod_months = gpd.read_file(rice_month_file)
		rice_prod_months = rice_prod_months.loc[rice_prod_months.SUB_REGION == province]
		rice_prod_months = rice_prod_months[crop_month_fields].values.tolist()
		rice_prod_months = np.array(rice_prod_months[0])/sum(rice_prod_months[0])
		rice_prod_months = rice_prod_months[rice_prod_months > 0]
		rice_prod_months = rice_prod_months.tolist()



		# all the crop OD pairs
		for file in os.listdir(crop_data_path):
			if file.endswith(".tif") and 'spam_p' in file.lower().strip():
				fpath = os.path.join(crop_data_path, file)
				crop_name = [cr for cr in crop_names if cr in file.lower().strip()][0]
				outCSVName = os.path.join(output_path,'crop_flows','crop_concentrations.csv')
				subprocess.run(["gdal2xyz.py",'-csv', fpath,outCSVName])

				'''Load points and convert to geodataframe with coordinates'''
				load_points = pd.read_csv(outCSVName,header=None,names=['x','y','tons'],index_col=None)
				load_points = load_points[load_points['tons'] > 0]

				geometry = [Point(xy) for xy in zip(load_points.x, load_points.y)]
				load_points = load_points.drop(['x', 'y'], axis=1)
				crs = {'init': 'epsg:4326'}
				crop_points = gpd.GeoDataFrame(load_points, crs=crs, geometry=geometry)

				del load_points

				# clip all to province
				prov_crop = gdf_geom_clip(crop_points,province_geom)

				if len(prov_crop.index) > 0:
					prov_crop_sindex = prov_crop.sindex
					prov_crop['NEAREST_G_NODE'] = prov_crop.geometry.apply(lambda x: get_nearest_node(x,sindex_nodes,nodes,'NODE_ID'))
					sindex_commune_center = prov_commune_center.sindex
					prov_crop['NEAREST_C_CENTER'] = prov_crop.geometry.apply(lambda x: get_nearest_node(x,sindex_commune_center,prov_commune_center,'OBJECTID'))

					crop_ods = crop_od_pairs(prov_crop,prov_commune_center,crop_name)
					province_ods_df.append(crop_ods)

					print ('Done with crop {0} in province {1}'.format(crop_name, province_name))

		all_ods = pd.concat(province_ods_df, axis=0, sort = 'False', ignore_index=True).fillna(0)

		all_ods_crop_cols = [c for c in all_ods.columns.values.tolist() if c in crop_names]
		all_ods['crop_tot'] = all_ods[all_ods_crop_cols].sum(axis = 1)

		# all_ods_val_cols = [c for c in all_ods.columns.values.tolist() if c != 'od_nodes']
		# all_ods = all_ods.groupby(['od_nodes'])[all_ods_val_cols].sum().reset_index()

		all_ods_val_cols = [c for c in all_ods.columns.values.tolist() if c not in ('origin','destination')]
		all_ods = all_ods.groupby(['origin','destination'])[all_ods_val_cols].sum().reset_index()

		all_ods['croptons'] = all_ods.apply(lambda x: assign_monthly_tons_crops(x,rice_prod_months,all_ods_crop_cols),axis = 1)
		all_ods[['min_croptons', 'max_croptons']] = all_ods['croptons'].apply(pd.Series)
		all_ods.drop('croptons',axis=1,inplace=True)

		cost_values_df = pd.read_excel(os.path.join(crop_data_path,'crop_unit_costs.xlsx'),sheet_name ='io_rev')
		all_ods['croprev'] = all_ods.apply(lambda x: assign_io_rev_costs_crops(x,cost_values_df,rice_prod_months,all_ods.columns.values.tolist(),exchange_rate),axis = 1)
		all_ods[['min_agrirev', 'max_croprev']] = all_ods['croprev'].apply(pd.Series)
		all_ods.drop('croprev',axis=1,inplace=True)
		all_ods['max_agrirev'] = all_ods[['max_croprev','netrev_agri']].max(axis = 1)
		all_ods.drop(['max_croprev','netrev_agri'],axis=1,inplace=True)

		all_ods['min_netrev'] = all_ods['min_agrirev'] + all_ods['netrev_noagri']
		all_ods['max_netrev'] = all_ods['max_agrirev'] + all_ods['netrev_noagri']

		all_ods.to_excel(excl_wrtr,province_name,index = False)
		excl_wrtr.save()

		# all_ods.to_csv(os.path.join(output_path,'{}_ods.csv'.format(province_name)),index = False)
		# all_ods = all_ods[['od_pair','min_croptons','max_croptons','min_netrev','max_netrev']]

		# flow_output_excel = os.path.join(output_path,'flow_mapping_paths','{}_roads_od_flow_growth.xlsx'.format(province_name))
		flow_output_excel = os.path.join(output_path,'flow_mapping_paths','{}_roads_od_flow_growth_communes.xlsx'.format(province_name))
		excl_wrtr_1 = pd.ExcelWriter(flow_output_excel)
		for grth in growth_scenarios:
			# all_ods = all_ods[['od_nodes','min_croptons','max_croptons','min_netrev','max_netrev']]
			all_ods = all_ods[['origin','destination','min_croptons','max_croptons','min_netrev','max_netrev']]
			for year in range(2017,2050):
				all_ods['min_croptons_{}'.format(year)] = math.pow((1+grth[0]/100),year - base_year)*all_ods['min_croptons']
				all_ods['max_croptons_{}'.format(year)] = math.pow((1+grth[0]/100),year - base_year)*all_ods['max_croptons']
				all_ods['min_netrev_{}'.format(year)] = math.pow((1+grth[0]/100),year - base_year)*all_ods['min_netrev']
				all_ods['max_netrev_{}'.format(year)] = math.pow((1+grth[0]/100),year - base_year)*all_ods['max_netrev']

			all_ods.to_excel(excl_wrtr_1,grth[1],index = False)
			excl_wrtr_1.save()
