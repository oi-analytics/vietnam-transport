=============
Data Assembly
=============


.. todo::
    Brief description of data requirements, format, expected location in data folder


Networks
--------
All finalised networks data are stored 

1. In the file path - /data/post_processed_networks/
2. As Excel sheets with post-processed network nodes and edges 
3. As Shapefiles with post-processed network nodes and edges

All networks created should have the attributes listed below

1. All nodes have the following attributes:
	- node_id - String Node ID
	- name - String name in Vietnamese/English
	- tons - Float assigned cargo freight tonnage using node 
	- population - Float assigned passenger/population number using node 
	- capacity - Float assigned capacity in tons/passenger numbers/other units
	- geometry - Shapely Point geometry of node with projection ESPG:4326

2. Attributes only present in inland and coastal port nodes
	- port_type - String name of type of port: inland or sea 	
	- port_class - String name of class of port: class1A (international) or class1 (domestic hub)  

3. All edges have the following attributes:
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

4. Attributes only present in Province and national roads edges
	- surface - String value for surface
	- road_class - Integer between 1 and 6
	- road_cond - String value: paved or unpaved 
	- asset_type - String name of type of asset

Origin-Destination matrices
---------------------------
All finalised OD matrices are stored 

1. In the path - results/flow_ods/
2. As Excel sheets
3. All node-level national OD matrices contain mode-wise and total OD flows with attributes:
    - origin - String node ID of origin node
    - destination - String node ID of destination node
    - o_region - String names of origin Province
    - d_region - String names of destination Province
    - min_rice - Float values of minimum daily tonnages of rice between OD nodes
    - max_rice - Float values of maximum daily tonnages of rice between OD nodes
    - min_tons - Float values of minimum daily tonnages between OD nodes
    - max_tons - Float values of maximum daily tonnages between OD nodes
    - commodity_names - Float values of daily tonnages of commodities/industries between OD nodes: here based on VITRANSS2 and IFPRI data

4. All aggregated province-level national OD matrices contain mode-wise and total OD flows with attributes:
    - o_region - String names of origin Province
    - d_region - String names of destination Province
    - min_rice - Float values of minimum daily tonnages of rice between OD nodes
    - max_rice - Float values of maximum daily tonnages of rice between OD nodes
    - min_tons - Float values of minimum daily tonnages between OD nodes
    - max_tons - Float values of maximum daily tonnages between OD nodes
    - commodity_names - Float values of daily tonnages of commodities/industries between OD nodes: here based on VITRANSS2 and IFPRI data

5. All province OD matrices contain province-wise OD flows
    - origin - String node ID of origin node
    - destination - String node ID of destination node
    - min_rice - Float values of minimum daily tonnages of rice between OD nodes
    - max_rice - Float values of maximum daily tonnages of rice between OD nodes
    - min_croptons - Float values of minimum daily tonnages of crops between OD nodes
    - max_croptons - Float values of maximum daily tonnages of crops between OD nodes
    - min_agrirev - Float value of Minimum daily revenue of agriculture firms between OD nodes
    - max_agrirev - Float value of Maximum daily revenue of agriculture firms between OD nodes
    - min_noagrirev - Float value of Minimum daily revenue of non-agriculture firms between OD nodes
    - max_noagrirev - Float value of Maximum daily revenue of non-agriculture firms between OD nodes
    - min_netrev - Float value of Minimum daily revenue of all firms between OD nodes
    - max_netrev - Float value of Maximum daily revenue of all firms between OD nodes
    - crop_names - Float values of daily tonnages of crops between OD nodes: here based on IFPRI crops


Hazards
-------


Administrative Areas
--------------------


Census Data
-----------


Macroeconomic Data
------------------


Flows
-----


Adaptation Options
------------------
