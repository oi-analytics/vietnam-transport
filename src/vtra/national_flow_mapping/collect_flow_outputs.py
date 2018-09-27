"""Summarise hazard data

Get OD data and process it
Author: Raghav Pant
Date: April 20, 2018
"""
import configparser
import csv
import glob
import os

import fiona
import fiona.crs
import rasterio


from sqlalchemy import create_engine
import subprocess as sp
import psycopg2
import osgeo.ogr as ogr

import pandas as pd
import copy


import ast
from osgeo import gdal
import geopandas as gpd
from shapely.geometry import Point
from geoalchemy2 import Geometry, WKTElement

import numpy as np


from vtra.utils import load_config
import vtra.transport_network_creation as tnc

def get_node_edge_path_flows(pd_dataframe,regional_id_list,industry,path_index,path_list,path_key_list,path_dict,val_threshold):
	for index, row in pd_dataframe.iterrows():
		npath = ast.literal_eval(row[0])
		epath = ast.literal_eval(row[1])
		orgn = [rg[1] for rg in regional_id_list if int(rg[0]) == int(row[2])][0]
		dest = [rg[1] for rg in regional_id_list if int(rg[0]) == int(row[3])][0]
		val = row[4]
		if val > val_threshold and len(epath) > 0:
			if epath and epath not in path_list:
				path_list.append(epath)

				path_index += 1

				path_key = 'path_' + str(path_index)
				path_key_list.append(path_key)
				path_dict.update({path_key:{}})


				path_dict[path_key].update({'node_path': npath})
				path_dict[path_key].update({'edge_path': epath})
				path_dict[path_key].update({industry:[(orgn,dest,val)]})

			else:
				pindex = path_list.index(epath)
				path_key = path_key_list[pindex]
				if industry not in path_dict[path_key].keys():
					path_dict[path_key].update({industry:[(orgn,dest,val)]})
				else:
					path_dict[path_key][industry].append(((orgn,dest,val)))


	return(path_index,path_list,path_key_list,path_dict)

def main():
	'''
	Create the database connection
	'''
	conf = load_config()

	try:
		conn = psycopg2.connect(**conf['database'])
	except:
		print ("I am unable to connect to the database")


	curs = conn.cursor()

	engine = create_engine('postgresql://{user}:{password}@{host}:{port}/{database}'.format({
		**conf['database']
	}))

	'''
	Step 2: Create the OD proprotions for the differnet modes
	'''
	'''
	First get the modal shares
	'''
	modes = ['road','rail','air','water']
	ind_cols = ['sugar','wood','steel','constructi','cement','fertilizer','coal','petroluem','manufactur','fishery','meat']
	crop_cols = ['rice','cash','cass','teas','maiz','rubb','swpo','acof','rcof','pepp']
	commodities = ind_cols + crop_cols

	province_id_list = []
	sql_query = "select od_id,name_eng from province_level_stats"
	curs.execute(sql_query)
	read_layer = curs.fetchall()
	for row in read_layer:
		province_id_list.append((row[0],row[1]))

	pth_dict = {}
	pth_idx = 0
	pth_list = []
	pth_k_list = []
	for com in commodities:
		for m in modes:
			ifile = 'network_od_flows_' + com + m + '.csv'

			od_flows = pd.read_csv(ifile).fillna(0)
			od_flows = od_flows[od_flows['Tonnage'] > 0.5]
			if od_flows.empty is False:
				flow_node_edge = od_flows.groupby(['node_path', 'edge_path', 'Origin_region', 'Destination_region'])['Tonnage'].sum().reset_index()
				pth_idx,pth_list,pth_k_list,pth_dict = get_node_edge_path_flows(flow_node_edge,province_id_list,com,pth_idx,pth_list,pth_k_list,pth_dict,0.5)

				del od_flows

			print ('Done with',ifile)

	pth_tuple_list = []
	for key, values in pth_dict.items():
		pth_tuple = [0]*(3+len(commodities))
		pth_tuple[0] = key
		pth_tuple[1] = pth_dict[key]['node_path']
		pth_tuple[2] = pth_dict[key]['edge_path']
		values_keys = pth_dict[key].keys()
		for i in range(len(commodities)):
			if commodities[i] in values_keys:
				pth_tuple[i+3] = pth_dict[key][commodities[i]]

		# pth_tuple.append(sum(pth_tuple[4:]))
		pth_tuple_list.append(tuple(pth_tuple))

	excel_writer = pd.ExcelWriter('vnm_path_flows.xlsx')
	df = pd.DataFrame(pth_tuple_list, columns = ['path_index','node_path','edge_path'] + commodities)
	df.to_excel(excel_writer,'path_flows',index = False)
	excel_writer.save()

if __name__ == '__main__':
	main()
