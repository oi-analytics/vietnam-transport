"""Summarise hazard data

Get OD data and process it
Author: Raghav Pant
Date: April 20, 2018
"""
import geopandas as gpd
import pandas as pd
import os
import igraph as ig
import numpy as np
import sys
import subprocess
from shapely.geometry import Point
from shapely.geometry import Polygon
from scipy.spatial import Voronoi
import itertools
import operator
import ast
import math



# from vtra.utils import load_config,extract_value_from_gdf,get_nearest_node,gdf_clip,gdf_geom_clip,
# 						count_points_in_polygon,voronoi_finite_polygons_2d,extract_nodes_within_gdf,
# 						assign_value_in_area_proportions,assign_value_in_area_proportions_within_common_region
from vtra.utils import *
from vtra.transport_network_creation import *

def crop_od_pairs(start_points,end_points,crop_name):
	save_paths = []
	for iter_,place in start_points.iterrows():
		try:
			closest_center = end_points.loc[end_points['OBJECTID']
			== place['NEAREST_C_CENTER']]['NEAREST_G_NODE'].values[0]

			get_od_pair = (place['NEAREST_G_NODE'],closest_center)
			save_paths.append((str(get_od_pair),place['tons']))
		except:
			print(iter_)


	od_pairs_df = pd.DataFrame(save_paths,columns = ['od_nodes',crop_name])
	od_pairs_df = od_pairs_df.groupby(['od_nodes'])[crop_name].sum().reset_index()

	return od_pairs_df

def assign_daily_min_max_tons_rice(crop_df,rice_prod_df):
	'''
	crop_code	crop_name	min_cost_perton	max_cost_perton
	'''
	sindex_rice_prod_df = rice_prod_df.sindex

	crop_df['min_frac'] = crop_df.geometry.apply(lambda x: get_nearest_node(x,sindex_rice_prod_df,rice_prod_df,'min_frac'))
	crop_df['max_frac'] = crop_df.geometry.apply(lambda x: get_nearest_node(x,sindex_rice_prod_df,rice_prod_df,'max_frac'))

	crop_df['min_rice'] = 1.0*crop_df['min_frac']*crop_df['tons']/30.0
	crop_df['max_rice'] = 1.0*crop_df['max_frac']*crop_df['tons']/30.0

	return crop_df

'''
Create the database connection
'''
def main():
	data_path,calc_path,output_path = load_config()['paths']['data'],load_config()['paths']['calc'],load_config()['paths']['output']

	# population_points_in = os.path.join(data_path,'Points_of_interest','population_points.shp')
	# commune_path = os.path.join(data_path,'Vietnam_boundaries','boundaries_stats','commune_level_stats.shp')

	# crop_data_path = os.path.join(data_path,'Agriculture_crops','crop_data')
	# rice_month_file = os.path.join(data_path,'rice_atlas_vietnam','rice_production.shp')
	# crop_month_fields = ['P_Jan','P_Feb','P_Mar','P_Apr','P_May','P_Jun','P_Jul','P_Aug','P_Sep','P_Oct','P_Nov','P_Dec']
	# crop_names = ['rice','cash','cass','teas','maiz','rubb','swpo','acof','rcof','pepp']

	'''
	Get the modal shares
	'''
	# modes_file_paths = [('Roads','national_roads'),('Railways','national_rail'),('Airports','airnetwork'),('Waterways','waterways')]
	modes_file_paths = [('Roads','national_roads')]
	# modes_file_paths = [('Roads','national_roads'),('Railways','national_rail'),('Airports','airnetwork'),('Waterways','waterways'),('Waterways','waterways')]
	modes = ['road','rail','air','inland','coastal']
	mode_cols = ['road','rail','air','inland','coastal']
	new_mode_cols = ['o','d','road','rail','air','inland','coastal']
	# new_mode_cols = ['o','d','road','rail','air','water']

	# modes_file_paths = [('Railways','national_rail')]
	# modes = ['rail','air','water']

	# modes_file_paths = [('Airports','airnetwork')]
	# modes = ['air','water']

	od_data_file = os.path.join(data_path, 'OD_data', 'OD_transport_data.xlsx')
	od_data_modes = pd.read_excel(od_data_file,sheet_name = 'mode').fillna(0)
	# od_data_modes.columns = map(str.lower, od_data_modes.columns)
	o_id_col = 'o'
	d_id_col = 'd'
	od_data_modes['total'] = od_data_modes[mode_cols].sum(axis=1)
	for m in mode_cols:
		od_data_modes[m] = od_data_modes[m]/od_data_modes['total'].replace(np.inf, 0)

	# od_data_modes['water'] = od_data_modes['inland'] + od_data_modes['coastal']
	od_data_modes = od_data_modes.fillna(0)
	# od_data_modes.to_csv('mode_frac.csv',index = False)

	od_fracs = od_data_modes[new_mode_cols]

	od_data_com = pd.read_excel(od_data_file,sheet_name = 'goods').fillna(0)
	ind_cols = ['sugar','wood','steel','constructi','cement','fertilizer','coal','petroluem','manufactur','fishery','meat']
	od_fracs = pd.merge(od_fracs,od_data_com,how='left', on=['o','d']).fillna(0)

	# od_fracs.to_csv('test0.csv')
	# print (od_fracs)

	# del od_data_com,od_data_modes

	od_fracs_crops = od_data_modes[new_mode_cols]
	crop_cols = ['rice','indust-cro']
	for cr in crop_cols:
		od_data_com_sums = od_data_com.groupby(['o','d']).agg({cr: 'sum'})
		od_com_frac = od_data_com_sums.groupby(level=0).apply(lambda x: x/float(x.sum()))
		od_com_frac = od_com_frac.reset_index(level=['o', 'd'])
		od_fracs_crops = pd.merge(od_fracs_crops,od_com_frac,how='left', on=['o','d']).fillna(0)

	del od_data_com,od_data_com_sums,od_com_frac

	# print (od_fracs_crops)
	# find the crop production months for the provinces
	crop_data_path = os.path.join(data_path,'Agriculture_crops','crop_data')
	rice_month_file = os.path.join(data_path,'rice_atlas_vietnam','rice_production.shp')
	crop_month_fields = ['P_Jan','P_Feb','P_Mar','P_Apr','P_May','P_Jun','P_Jul','P_Aug','P_Sep','P_Oct','P_Nov','P_Dec']
	crop_names = ['rice','cash','cass','teas','maiz','rubb','swpo','acof','rcof','pepp']
	rice_prod_months = gpd.read_file(rice_month_file)
	rice_prod_months['total_prod'] = rice_prod_months[crop_month_fields].sum(axis = 1)
	rice_prod_months['min_tons'] = rice_prod_months[rice_prod_months[crop_month_fields] > 0].min(axis=1)
	rice_prod_months['max_tons'] = rice_prod_months[rice_prod_months[crop_month_fields] > 0].max(axis=1)

	rice_prod_months['min_frac'] = rice_prod_months['min_tons']/rice_prod_months['total_prod']
	rice_prod_months['max_frac'] = rice_prod_months['max_tons']/rice_prod_months['total_prod']

	# print (rice_prod_months)

	province_path = os.path.join(data_path,'Vietnam_boundaries','boundaries_stats','province_level_stats.shp')
	commune_path = os.path.join(data_path,'Vietnam_boundaries','boundaries_stats','commune_level_stats.shp')
	rd_prop_file = os.path.join(data_path,'mode_properties','road_properties.xlsx')

	flow_output_excel = os.path.join(output_path,'flow_mapping_paths','national_scale_flow_ods_road.xlsx')
	excl_wrtr = pd.ExcelWriter(flow_output_excel)

	flow_output_excel = os.path.join(output_path,'flow_mapping_paths','national_scale_od_matrix_road.xlsx')
	excl_wrtr_reg = pd.ExcelWriter(flow_output_excel)

	# load provinces and get geometry of the right province
	provinces = gpd.read_file(province_path)
	provinces = provinces.to_crs({'init': 'epsg:4326'})
	sindex_provinces = provinces.sindex

	# load provinces and get geometry of the right province
	communes = gpd.read_file(commune_path)
	communes = communes.to_crs({'init': 'epsg:4326'})
	communes['province_name'] = communes.geometry.apply(lambda x: get_nearest_node(x,sindex_provinces,provinces,'name_eng'))
	sindex_communes = communes.sindex

	# print (communes)

	modes_df = []
	for m in range(len(modes_file_paths)):
		mode_data_path = os.path.join(data_path,modes_file_paths[m][0],modes_file_paths[m][1])
		for file in os.listdir(mode_data_path):
			try:
				if file.endswith(".shp") and 'edges' in file.lower().strip():
					edges_in = os.path.join(mode_data_path, file)
				if file.endswith(".shp") and 'nodes' in file.lower().strip():
					nodes_in = os.path.join(mode_data_path, file)
			except:
				print ('Network nodes and edge files necessary')

		# if modes[m] == 'road':
		# 	od_net =  national_shapefile_to_network(edges_in,rd_prop_file)
		# 	od_net = add_igraph_generalised_costs_province_roads(od_net,1,vehicle_wt)

		# load nodes of the network
		nodes = gpd.read_file(nodes_in)
		nodes = nodes.to_crs({'init': 'epsg:4326'})
		nodes.columns = map(str.lower, nodes.columns)
		node_cols = nodes.columns.values.tolist()
		node_cols = [c for c in node_cols if c not in ('population','od_id')]
		nodes = nodes[node_cols]
		sindex_nodes = nodes.sindex

		# print (sindex_nodes.tolist())
		# assign province ID's and OD ID's to their nearest nodes
		# nodes['province_name'] = nodes.geometry.apply(lambda x: get_nearest_node(x,sindex_provinces,provinces,'name_eng'))
		# nodes['od_id'] = nodes.geometry.apply(lambda x: get_nearest_node(x,sindex_provinces,provinces,'od_id'))

		nodes['province_name'] = nodes.apply(lambda x: extract_gdf_values_containing_nodes(x,sindex_provinces,provinces,'name_eng'),axis = 1)
		nodes['od_id'] = nodes.apply(lambda x: extract_gdf_values_containing_nodes(x,sindex_provinces,provinces,'od_id'),axis = 1)

		# nodes['province_name'] = extract_gdf_values_containing_nodes(x,input_gdf,column_name)

		if modes[m] == 'road':
			edges_df =  national_road_shapefile_to_dataframe(edges_in,rd_prop_file)
			nodes_vehs = list(zip(edges_df['from_node'].values.tolist(),edges_df['from_node'].values.tolist(),edges_df['vehicle_co'].values.tolist()))
			nd_veh_list = []
			for nd in nodes['node_id'].values.tolist():
				veh = 0.5*sum([int(v) for (f,t,v) in nodes_vehs if nd == f or nd == t])
				nd_veh_list.append((nd,veh))

			gdf_pops = pd.DataFrame(nd_veh_list,columns = ['node_id','population'])
			del nd_veh_list
			nodes = pd.merge(nodes,gdf_pops,how='left', on=['node_id']).fillna(0)
			del gdf_pops

		elif modes[m] in ('inland','coastal'):
			nodes['population'] = nodes['tons']

		else:
			xy_list = []
			for iter_,values in nodes.iterrows():
				# print (list(values.geometry.coords))
				xy = list(values.geometry.coords)
				xy_list += [list(xy[0])]

			vor = Voronoi(np.array(xy_list))
			regions, vertices = voronoi_finite_polygons_2d(vor)
			min_x = vor.min_bound[0] - 0.1
			max_x = vor.max_bound[0] + 0.1
			min_y = vor.min_bound[1] - 0.1
			max_y = vor.max_bound[1] + 0.1

			mins = np.tile((min_x, min_y), (vertices.shape[0], 1))
			bounded_vertices = np.max((vertices, mins), axis=0)
			maxs = np.tile((max_x, max_y), (vertices.shape[0], 1))
			bounded_vertices = np.min((bounded_vertices, maxs), axis=0)

			box = Polygon([[min_x, min_y], [min_x, max_y], [max_x, max_y], [max_x, min_y]])
			# colorize
			poly_list = []
			for region in regions:
				polygon = vertices[region]
				# Clipping polygon
				poly = Polygon(polygon)
				poly = poly.intersection(box)
				poly_list.append(poly)


			poly_index = list(np.arange(0,len(poly_list),1))
			poly_df = pd.DataFrame(list(zip(poly_index,poly_list)),columns = ['gid','geometry'])
			gdf_voronoi = gpd.GeoDataFrame(poly_df,crs='epsg:4326')
			gdf_voronoi['node_id'] = gdf_voronoi.apply(lambda x: extract_nodes_within_gdf(x,nodes,'node_id'),axis = 1)
			gdf_voronoi['population'] = 0
			gdf_voronoi = assign_value_in_area_proportions(communes,gdf_voronoi,'population')
			# gdf_voronoi = assign_value_in_area_proportions_within_common_region(communes,gdf_voronoi,'population','province_name')

			gdf_pops = gdf_voronoi[['node_id','population']]
			# print (gdf_pops)
			del gdf_voronoi, poly_list, poly_df

			nodes = pd.merge(nodes,gdf_pops,how='left', on=['node_id']).fillna(0)
			del gdf_pops

		# nodes = nodes[['node_id','od_id','population']]
		nodes_sums = nodes.groupby(['od_id','node_id']).agg({'population': 'sum'})
		nodes_frac = nodes_sums.groupby(level = 0).apply(lambda x: x/float(x.sum()))
		nodes_frac = nodes_frac.reset_index(level = ['od_id','node_id'])
		nodes_frac.rename(columns={'population': 'pop_frac'}, inplace=True)

		nodes = pd.merge(nodes,nodes_frac[['node_id','pop_frac']],how='left', on=['node_id']).fillna(0)

		# nodes.to_file(os.path.join(output_path,'networks_test','{}_nodes.shp'.format(modes[m])))

		# print (nodes)

		modes_df.append(nodes)

		del nodes_frac, nodes_sums, nodes



	national_ods_df = []
	for ind in ind_cols:
		national_ods_modes_df = []
		for m in range(len(modes_file_paths)):
			nodes = modes_df[m]
			od_nodes_regions = list(zip(nodes['node_id'].values.tolist(),nodes['province_name'].values.tolist(),nodes['od_id'].values.tolist(),nodes['pop_frac'].values.tolist()))
			ind_mode = modes[m]+ '_' + ind
			od_fracs[ind_mode] = od_fracs[modes[m]]*od_fracs[ind]

			od_fracs_ind = od_fracs[[o_id_col,d_id_col,ind_mode]]
			od_fracs_ind = od_fracs_ind[od_fracs_ind[ind_mode] > 0]
			od_flows = list(zip(od_fracs_ind[o_id_col].values.tolist(),od_fracs_ind[d_id_col].values.tolist(),od_fracs_ind[ind_mode].values.tolist()))
			origins = list(set(od_fracs_ind[o_id_col].values.tolist()))
			destinations = list(set(od_fracs_ind[d_id_col].values.tolist()))

			# print (od_flows)
			od_list = []
			for o in origins:
				for d in destinations:
					fval = [fl for (org,des,fl) in od_flows if org == o and des == d]
					if len(fval) == 1 and fval[0] > 0:
						o_matches = [(item[0],item[1],item[3]) for item in od_nodes_regions if item[2] == o]
						if len(o_matches) > 0:
							for o_vals in o_matches:
								o_val = 1.0*fval[0]*(1.0*o_vals[2])
								o_node = o_vals[0]
								o_region = o_vals[1]
								d_matches = [(item[0],item[1],item[3]) for item in od_nodes_regions if item[2] == d]
								if len(d_matches) > 0:
									for d_vals in d_matches:
										od_val = 1.0*o_val*(1.0*d_vals[2])
										d_node = d_vals[0]
										d_region = d_vals[1]
										if od_val > 0 and o_node != d_node:
											od_list.append((o_node,o_region,d_node,d_region,od_val))

					print (o,d,fval,modes[m],ind)

			national_ods_modes_df.append(pd.DataFrame(od_list,columns = ['origin','o_region','destination','d_region',ind]))
			del od_list, nodes

		national_ods_df.append(national_ods_modes_df)



	# all the crop OD pairs
	for file in os.listdir(crop_data_path):
		if file.endswith(".tif") and ('spam_p' in file.lower().strip()):
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
			# prov_crop = gdf_geom_clip(crop_points,province_geom)
			if crop_name == 'rice':
				crop_points = assign_daily_min_max_tons_rice(crop_points,rice_prod_months)
			else:
				crop_points['min_{}'.format(crop_name)] = 1.0*crop_points['tons']/365.0
				crop_points['max_{}'.format(crop_name)] = 1.0*crop_points['tons']/365.0

			# crop_points_sindex = crop_points.sindex

			crop_points['province_name'] = crop_points.apply(lambda x: extract_gdf_values_containing_nodes(x,sindex_provinces,provinces,'name_eng'),axis = 1)
			national_ods_modes_df = []
			for m in range(len(modes_file_paths)):
				nodes = modes_df[m]
				crop_pts = crop_points.copy(deep = True)
				crop_pts['node_id'] = crop_pts.apply(lambda x: get_nearest_node_within_region(x,nodes,'node_id','province_name'),axis = 1)
				# crop_points.to_file(os.path.join(output_path,'Voronoi','crop_test_2.shp'))
				crop_pts = crop_pts[crop_pts['node_id'] != '']
				crop_pts = crop_pts[['node_id','min_{}'.format(crop_name),'max_{}'.format(crop_name)]]
				crop_nodes = crop_pts.groupby(['node_id'])['min_{}'.format(crop_name),'max_{}'.format(crop_name)].sum().reset_index()
				crop_nodes = crop_nodes.reset_index()
				# crop_nodes.to_csv(os.path.join(output_path,'Voronoi','crop_test_2.csv'),index = False)

				del crop_pts
				nodes = pd.merge(nodes,crop_nodes,how='left', on=['node_id']).fillna(0)
				del crop_nodes

				crop_mode = modes[m]+ '_' + crop_name
				if crop_name in ('rice', 'cereal', 'wheat'):
					od_fracs_crops[crop_mode] = od_fracs_crops[modes[m]]*od_fracs_crops['rice']
				else:
					od_fracs_crops[crop_mode] = od_fracs_crops[modes[m]]*od_fracs_crops['indust-cro']

				od_nodes_regions = list(zip(nodes['node_id'].values.tolist(),nodes['province_name'].values.tolist(),nodes['od_id'].values.tolist(),nodes['min_{}'.format(crop_name)].values.tolist(),nodes['max_{}'.format(crop_name)].values.tolist(),nodes['pop_frac'].values.tolist()))

				od_fracs_ind = od_fracs_crops[[o_id_col,d_id_col,crop_mode]]
				od_fracs_ind = od_fracs_ind[od_fracs_ind[crop_mode] > 0]
				od_flows = list(zip(od_fracs_ind[o_id_col].values.tolist(),od_fracs_ind[d_id_col].values.tolist(),od_fracs_ind[crop_mode].values.tolist()))
				origins = list(set(od_fracs_ind[o_id_col].values.tolist()))
				destinations = list(set(od_fracs_ind[d_id_col].values.tolist()))

				od_list = []
				for o in origins:
					for d in destinations:
						fval = [fl for (org,des,fl) in od_flows if org == o and des == d]
						if len(fval) == 1 and fval[0] > 0:
							o_matches = [(item[0],item[1],item[3],item[4]) for item in od_nodes_regions if item[2] == o]
							if len(o_matches) > 0:
								for o_vals in o_matches:
									o_val_min = 1.0*fval[0]*o_vals[2]
									o_val_max = 1.0*fval[0]*o_vals[3]
									o_node = o_vals[0]
									o_region = o_vals[1]
									d_matches = [(item[0],item[1],item[5]) for item in od_nodes_regions if item[2] == d]
									if len(d_matches) > 0:
										for d_vals in d_matches:
											od_val_min = 1.0*o_val_min*d_vals[2]
											od_val_max = 1.0*o_val_max*d_vals[2]
											d_node = d_vals[0]
											d_region = d_vals[1]
											if od_val_max > 0 and o_node != d_node:
												od_list.append((o_node,o_region,d_node,d_region,od_val_min,od_val_max))

					print (o,d,fval,modes[m],crop_name)

				national_ods_modes_df.append(pd.DataFrame(od_list,columns = ['origin','o_region','destination','d_region','min_{}'.format(crop_name),'max_{}'.format(crop_name)]))
				del od_list, nodes

			del crop_points
			national_ods_df.append(national_ods_modes_df)

	national_ods_df = list(map(list,zip(*national_ods_df)))
	region_total = []
	for m in range(len(modes_file_paths)):
		all_ods = pd.concat(national_ods_df[m], axis=0, sort = 'False', ignore_index=True).fillna(0)

		all_min_cols = ind_cols + ['min_{}'.format(c) for c in crop_names]
		all_ods['min_tons'] = all_ods[all_min_cols].sum(axis = 1)
		all_max_cols = ind_cols + ['max_{}'.format(c) for c in crop_names]
		all_ods['max_tons'] = all_ods[all_max_cols].sum(axis = 1)
		crops_norice = [cr for cr in crop_names if cr != 'rice']
		for cr in crops_norice:
			all_ods.drop('min_{}'.format(cr),axis=1,inplace=True)
			all_ods.rename(columns={'max_{}'.format(cr): cr}, inplace=True)

		all_ods_val_cols = [c for c in all_ods.columns.values.tolist() if c not in ('origin','o_region','destination','d_region')]
		all_ods = all_ods.groupby(['origin','o_region','destination','d_region'])[all_ods_val_cols].sum().reset_index()

		all_ods_regions = all_ods[['o_region','d_region'] + all_ods_val_cols]
		all_ods_regions = all_ods_regions.groupby(['o_region','d_region'])[all_ods_val_cols].sum().reset_index()
		all_ods_regions.to_excel(excl_wrtr_reg,modes[m],index = False)
		excl_wrtr_reg.save()

		region_total.append(all_ods_regions)
		del all_ods_regions

		all_ods = all_ods[all_ods['max_tons'] > 0.5]
		# flow_output_csv = os.path.join(output_path,'flow_mapping_paths','national_scale_{}_ods.csv'.format(modes[m]))
		# all_ods.to_csv(flow_output_csv,index = False)
		all_ods.to_excel(excl_wrtr,modes[m],index = False)
		excl_wrtr.save()
		del all_ods

	# all_ods = pd.concat(region_total, axis=0, sort = 'False', ignore_index=True).fillna(0)
	# all_ods_val_cols = [c for c in all_ods.columns.values.tolist() if c not in ('o_region','d_region')]
	# all_ods_regions = all_ods.groupby(['o_region','d_region'])[all_ods_val_cols].sum().reset_index()
	# all_ods_regions.to_excel(excl_wrtr_reg,'total',index = False)
	# excl_wrtr_reg.save()


if __name__ == '__main__':
	main()
