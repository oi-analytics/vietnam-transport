==============
Collected Data
==============
.. Important::
	- This section describes collected datasets that are used to create data for the Vietnam Transport Risk Analysis (VTRA)
	- The datasets listed here are specfic to Vietnam and are used as inputs to data in the Processed Data Assembly steps
	- To implement the VTRA pre-processing without any changes in existing codes, all data described here should be created and stored exactly as indicated below

Networks
--------
.. Note::
	1. All pre-processed networks data are stored: 
		- In sub-folders in the file path - /data/pre_processed_networks_data/
		- As Shapefiles with of network nodes and edges
		- The names of files and folders are self-explanatory
		- See /data/pre_processed_networks_data/networks_description.xlsx for details of all shapefiles

	2. All nodes should have the following attributes:
		- node_id - String Node ID
		- geometry - Shapely Point geometry of node with projection ESPG:4326
		- variable list of attributes depending upon sector  

	3. All edges should have the following attributes:
		- edge_id - String edge ID
		- g_id - Interger edge ID
		- from_node - String node ID that should be present in node_id column
		- to_node - String node ID that should be present in node_id column
		- geometry - Shapely LineString geometry of edge with projection ESPG:4326
		- variable list of attributes depending upon sector

Cost attributes
---------------
.. Note::
	1. Data to assign transport costs to network edges are stored:
		- In the file in path - /data/pre_processed_networks_data/mode_costs.xlsx
		- As Excel sheets

	2. All cost estimtates should have the following attributes:
		- time_cost_usd - Float values of rate of time
		- tariff_min_vnd - Float values minimum tariff rate in VND/ton-km (VND/ton for multi)
		- tariff_max_vnd - Float values maximum tariff rate in VND/ton-km (VND/ton for multi)
		- tariff_min_usd - Float values minimum tariff rate in USD/ton-km (USD/ton for multi)
		- tariff_max_usd - Float values maximum tariff rate in USD/ton-km (USD/ton for multi)
		- attributes to decide how the costs are allocated to network edges (if none then all edges have same criteria)

Road design attributes
----------------------
.. Note::
	1. Data to assign characteristics to roads are stored:
		- In the file in path - /data/pre_processed_networks_data/road_properties.xlsx
		- As Excel sheets
		- See /data/pre_processed_networks_data/road_properties.xlsx for data description


VITRANNS2 Origin-Destination data
---------------------------------
.. Note::
	1. VITRANSS2 province-level OD matrices are stored:
		- In the path - data/OD_data/
		- As Excel sheets
		- 'goods' sheet gives OD values by commodity
		- 'modes' sheet gives OD values by mode

	2. Aggregated goods-wise province-level national OD matrices have attributes:
	    - o - Integer IDs of origin Provinces
	    - d - Integer IDs of of destination Provinces
	    - name o - String names of origin Provinces
	    - name d - String names of destination Provinces
	    - commodity_names - Float values of daily tonnages of commodities/industries between OD Provinces

	3. Aggregated mode-wise province-level national OD matrices have attributes:
	    - o - Integer IDs of origin Provinces
	    - d - Integer IDs of of destination Provinces
	    - name o - String names of origin Provinces
	    - name d - String names of destination Provinces
	    - mode_names - Float values of daily tonnages along modes between OD Provinces

IFPRI crop data
---------------
.. Note::
	1. IFPRI crop datasets are stored:
		- In the path - data/Agriculture_crops/
		- As GeoTiff files
		- Only files with names 'SPAM_P_crop name_ver3.tif' are used
		- See Excel sheet in path data/Agriculture_crops//crop_data/crop_unit_costs.xlsx for costs of crops

	2. All crop GeoTiff datasets should have attributes:
	    - values - greater than 0
	    - raster grid geometry
	    - projection systems: Default assumed = EPSG:4326

RiceAltas data
--------------
.. Note::
	1. RiceAltas datasets are stored:
		- In the path - data/rice_atlas_vietnam/
		- As Shapefiles
		- Only the file 'rice_production.shp' is used

	2. The essential attributes in the dataset are listed below. See the data for all attributes:
	    - sub_region - String names of Provinces in English 
	    - P_Jan, ..., P_Dec - Columne names with float tonnage produced in each month from January to December
	    - geometry - Shapely Polygon geometries of Provinces

Points of interest data
-----------------------
.. Note::
	1. Locations of populations, commune, district, province center committee points datasets are stored:
		- In the path - data/Points_of_interest/
		- As Shapefiles

	2. The essential attributes in all the dataset are listed below. See the data for all attributes:
	    - geomtery - Shapely Point geometry with projection ESPG:4326


Macroeconomic Data
------------------
.. Note::
	1. All macroeconomic datasets are stored:
		-  