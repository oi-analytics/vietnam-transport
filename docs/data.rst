=======================
Processed Data Assembly
=======================
.. Important::
	- This section describes processed datasets that are used as inputs in Analysis and Results steps of the Vietnam Transport Risk Analysis (VTRA)
	- The formats and attributes created in these datasets form the essential inputs for implementing the rest of the VTRA model
	- To implement the VTRA without any changes in existing codes, all data described here should be created and stored exactly as indicated below  

Networks
--------
.. Note::
	1. All finalised networks data are stored: 
		- In the file path - /data/post_processed_networks/
		- As Excel sheets with post-processed network nodes and edges 
		- As Shapefiles with post-processed network nodes and edges

	2. All nodes have the following attributes:
		- node_id - String Node ID
		- name - String name in Vietnamese/English
		- tons - Float assigned cargo freight tonnage using node 
		- population - Float assigned passenger/population number using node 
		- capacity - Float assigned capacity in tons/passenger numbers/other units
		- geometry - Shapely Point geometry of node with projection ESPG:4326

	3. Attributes only present in inland and coastal port nodes:
		- port_type - String name of type of port: inland or sea 	
		- port_class - String name of class of port: class1A (international) or class1 (domestic hub)  

	4. All edges have the following attributes:
		- edge_id - String edge ID
		- g_id - Interger edge ID
		- from_node - String node ID that should be present in node_id column
		- to_node - String node ID that should be present in node_id column
		- geometry - Shapely LineString geometry of edge with projection ESPG:4326
		- terrain - String name of terrain of edge	
		- level - Integer number for edge level: National, Provincial, Local, etc.
		- width - Float width in meters of edge
		- length - Float estimated length in kilometers of edge	
		- min_speed - Float estimated minimum speed in km/hr on edge
		- max_speed - Float estimated maximum speed in km/hr on edge
		- min_time - Float estimated minimum time of travel in hours on edge
		- max_time - Float estimated maximum time of travel in hours on edge	
		- min_time_cost - Float estimated minimum cost of time in USD on edge
		- max_time_cost - Float estimated maximum cost of time in USD on edge
		- min_tariff_cost - Float estimated minimum tariff cost in USD on edge	
		- max_tariff_cost - Float estimated maximum tariff cost in USD on edge
		- vehicle_co - Integer number of daily vehicle counts on edge

	4. Attributes only present in province and national roads edges:
		- surface - String value for surface
		- road_class - Integer between 1 and 6
		- road_cond - String value: paved or unpaved 
		- asset_type - String name of type of asset

OD matrices
---------------------------
.. Note::
	1. All finalised OD matrices are stored:
		- In the path - /results/flow_ods/
		- As Excel sheets

	The essential attributes in these OD matrices are listed below. See the data for all attributes

	2. All node-level national OD matrices contain mode-wise and total OD flows and should have attributes:
	    - origin - String node IDs of origin nodes
	    - destination - String node IDs of destination nodes
	    - o_region - String names of origin Provincse
	    - d_region - String names of destination Provinces
	    - min_tons - Float values of minimum daily tonnages between OD nodes
	    - max_tons - Float values of maximum daily tonnages between OD nodes
	    - commodity_names - Float values of daily min-max tonnages of commodities/industries between OD nodes: here based on VITRANSS2 and IFPRI data

	3. All aggregated province-level national OD matrices contain mode-wise and total OD flows and should have attributes:
	    - o_region - String names of origin Provinces
	    - d_region - String names of destination Provinces
	    - min_tons - Float values of minimum daily tonnages between OD Provinces
	    - max_tons - Float values of maximum daily tonnages between OD Provinces
	    - commodity_names - Float values of daily min-max tonnages of commodities/industries between OD Provinces: here based on VITRANSS2 and IFPRI data

	4. All province OD matrices contain province-wise OD flows and should have attributes:
	    - origin - String node IDs of origin nodes
	    - destination - String node IDs of destination nodes
	    - min_croptons - Float values of minimum daily tonnages of crops between OD nodes
	    - max_croptons - Float values of maximum daily tonnages of crops between OD nodes
	    - min_netrev - Float values of minimum daily revenue of all firms between OD nodes
	    - max_netrev - Float values of maximum daily revenue of all firms between OD nodes


Hazards
-------
.. Note::
	1. All hazard datasets are stored:
		- In sub-folders in the path - /data/Hazard_data/
		- As GeoTiff files
		- See /data/hazard_data/hazard_data_folder_data_info.xlsx for details of all hazard files
	
	2. Single-band GeoTiff hazard raster files should have attributes:
	    - values - between 0 and 1000
	    - raster grid geometry
	    - projection systems: Default assumed = EPSG:32648

	3. Multi-band GeoTiff hazard raster files should have attributes:
	    - 3-bands
	    - values - in each band between 0 and 255
	    - raster grid geometry
	    - projection systems: Default assumed = EPSG:32648


Administrative Areas with Statistics
------------------------------------
.. Note::
	1. Vietnam boundary datasets are stored:
		- In the path - /data/Vietnam_boundaries/who_boundaries/
		- In the path - /data/Vietnam_boundaries/boundaries_stats/
		- As Shapefiles

	2. Global boundary dataset for map plotting are stored:
		- In the path - /data/Global_boundaries/Natural_Earth/ 

	The essential attributes in the Vietnam boundary datasets are listed below. See the data for all attributes

	3. All Vietnam province boundary datasets should have the attributes:
	    - name_eng - String names of administrative boundary in English
	    - od_id - Integer IDs matching ID's in VITRANSS2 OD data 
	    - geometry - Shapely Polygon geometries of boundary with projection ESPG:4326

	4. All Vietnam commune boundary datasets should have attributes:
	    - commune_id - Integer IDs of commune
	    - name_eng - String names of commune in English
	    - district_i - Integer IDs of district of commune
	    - dis_name_e -  String names of district in English
	    - province_i - Integer IDs of province of commune
	    - pro_name_e -  String names of province in English
	    - population - Float values of population in commune
	    - nfirms - Float values of number of firms in commune
	    - netrevenue - Float values of netrevenue of commune
	    - nongnghiep - Float fractions of agriculture firms in commune
	    - geometry - Shapely Polygon geometry of boundary with projection ESPG:4326

	5. All global boundary datasets should have attributes:
		- name - String names of boundaries in English
		- geometry - Shapely Polygon geometry of boundary with projection ESPG:4326 
	    

Macroeconomic Data
------------------
.. Note::
	1. For the macroeconomic analysis we use the national IO table for Vietnam:
		- In the file in path - data/economic_IO_tables/IO Table 2012 English.xlsx 
		- We use the sheet 'IO Core' in our analysis.
		
    
Adaptation Options
------------------
.. Note::
	1. All adaptation options input datasets are stored:
		- In the file - /data/Adaptation_options/adaptation_costs_road_types.xlsx
		- As Excel files

	2. Excel sheet 'options' explains:
		- adaptation_group - String names of the type of adaptation strategy
		- option_code - Stirng codes of the option	
		- item - String descriptions of the option

	3. Excel sheets 'costs_district_mountain','costs_district_flat','costs_national_mountain', 'costs_national_flat' and 'bridges' explain:
		- adaptation_group - String names of the type of adaptation strategy
		- option_code - Stirng codes of the option	
		- item - String descriptions of the option
		- estimated_amount_fraction	- Float ratios of road length over which the option is implemented
		- Estimated length - Flaot assumed road lengths in example case
		- factor - Float factor multiplied to length based on rates and option	
		- factor_unit - String factor units 	
		- unit - String dimension units of measurement of option	
		- rate - Float unit rates of option in USD/unit
		- total cost - Float product of Estimated length*factor*rate 	
		- design_assumption	- String explanations of some design assumptions
		- comment - Sting general comments

	4. Excel sheet 'rehabilitation_costs' explain:
		- Road Type	- String type of road
		- terrain - String terrian of road	
		- basic_cost - Float current costs of rehabilitation in USD/km of road 	
		- design_width - Float design width of road	
		- road_class - Integer class for National Roads	
		- road_level - Integer level for Province Roads	
		- unit - String unit of cost

	5. Excels sheets 'maintenance_mountain' and 'maintenance_flat' explain:
		- adaptation_group - String names of the type of adaptation strategy
		- option_code - Stirng codes of the option	
		- item - String descriptions of the option	
		- recurrent_cost - Float costs of recurrent maintenance	
		- periodic_cost - Float costs of periodic maintenance	
		- recurrent_factor - Float factor of recurrent maintenance	 	
		- periodic_factor - Float factor of periodic maintenance	
		- recurrent_maintain_time - Float times of recurrent maintenance in years
		- periodic_maintain_time - Float times of periodic maintenance in years	
		- recurrent_cost_unit - Sting unit of recurrent maintenance costs 	
		- periodic_cost_unit - Sting unit of periodic maintenance costs	
		- comment - Sting general comments