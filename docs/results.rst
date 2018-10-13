====================
Analysis and Results
====================
.. Important::
    - This section describes the steps Analysis and Results steps of the Vietnam Transport Risk Analysis (VTRA)
    - The formats and attributes created in these datasets form the essential inputs for implmenting the rest of the VTRA model
    - To implement the VTRA without any changes in existing codes, all data described here should be created and stored exactly as indicated below

Pre-processing and Preparing Network Data
-----------------------------------------
.. Note::
    Purpose:
        - Create post-processed transport networks with attributes
        - From pre-processed input Shapefiles and collected network attributes data
        - For all Province road networks
        - For all transport modes at national scale
    
    Execution:
        - Load data as described in `Collected Data <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/predata.html>`_ `Networks <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/predata.html#networks>`_, `Cost attributes <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/predata.html#cost-attributes>`_ and `Road design attributes <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/predata.html#road-design-attributes>`_  
        - Run script vtra/preprocess/create_transport_networks.py

    Result:
        - Create networks with formats and attributes described in `Processed Data Assembly <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/data.html>`_ `Networks <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/data.html#networks>`_
        - Store outputs in /data/post_processed_networks/


Pre-processing and Preparing Hazard Data
----------------------------------------
.. Note::
    Purpose:
        - Convert GeoTiff raster hazard datasets to shapefiles based on masking and selecting values from
            - Single-band raster files
            - Multi-band (3-bands) raster files
    
    Execution:
        - Load data as described in `Processed Data Assembly <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/data.html>`_ `Hazards <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/data.html#hazards>`_ 
        - Run script vtra/preprocess/convert_hazard_data.py

    Result:
        - Create hazard shapefiles with names described in excel sheet in `Processed Data Assembly <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/data.html>`_ `Hazards <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/data.html#hazards>`_ and attributes:
            - ID - equal to 1 
            - geometry - Shapely Polygon outline of selected hazard
        - Store outputs in same paths in directory /data/Hazard_data/


Pre-processing and Preparing OD matrix Data
-------------------------------------------
.. Note::
    Purpose:
        - Create national scale OD matrices at node and province levels from: 
            - VITRANSS2 province-scale OD data
            - IFPRI crop data at 1km resolution
        - Create province scale OD matrices between roads connecting villages to nearest communes: 
            - Net revenue estimates of commune villages
            - IFPRI crop data at 1km resolution
    
    Execution:
        - Load data as described in `Networks <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/data.html#networks>`_, `VITRANSS2 OD data <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/predata.html#vitranns2-od-data>`_,`IFPRI crop data <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/predata.html#ifpri-crop-data>`_,`RiceAtlas data <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/predata.html#ricealtas-data>`_, `Points of interest data <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/predata.html#points-of-interest-data>`_, and `Administrative Areas with Statistics <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/data.html#administrative-areas-with-statistics>`_  
        - For National OD matrices run script vtra/preprocess/national_modes_od_creation.py
        - For Provinces OD matrices run script vtra/preprocess/province_roads_access_od_creation.py

    Result:
        - Create OD matrices with attributes described in `Processed Data Assembly <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/data.html>`_ `OD matrices <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/data.html#od-matrices>`_
        - Store outputs in /results/flow_ods/


Mapping Flows onto Networks
---------------------------

.. todo::
    Brief description of process:

    - load flows (regional/origin-only...)
    - load network, candidate o-d points
    - weighted assignment, least-cost routing

    Including:

    - which data files are used
    - which scripts to run in what order
    - output files


Hazard Exposure
---------------

.. todo::
    Brief description of process:

    - load networks, hazards
    - intersection

    Including:

    - which data files are used
    - which scripts to run in what order
    - output files


Failure Analysis
----------------

.. todo::
    Brief description of process:

    - load networks, hazards
    - intersection

    Including:

    - which data files are used
    - which scripts to run in what order
    - output files

Input data requirements
~~~~~~~~~~~~~~~~~~~~~~~

1. Correct paths to all files and correct input parameters
2. Excel sheets with results of flow mapping based on MIN-MAX generalised costs estimates:
    - origin - String node ID of Origin
    - destination - String node ID of Destination
    - o_region - String name of Province of Origin node ID
    - d_region - String name of Province of Destination node ID
    - min_edge_path - List of string of edge ID's for paths with minimum generalised cost flows
    - max_edge_path - List of string of edge ID's for paths with maximum generalised cost flows
    - min_distance - Float values of estimated distance for paths with minimum generalised cost flows
    - max_distance - Float values of estimated distance for paths with maximum generalised cost flows
    - min_time - Float values of estimated time for paths with minimum generalised cost flows
    - max_time - Float values of estimated time for paths with maximum generalised cost flows
    - min_gcost - Float values of estimated generalised cost for paths with minimum generalised cost flows
    - max_gcost - Float values of estimated generalised cost for paths with maximum generalised cost flows
    - min_vehicle_nums - Float values of estimated vehicle numbers for paths with minimum generalised cost flows
    - max_vehicle_nums - Float values of estimated vehicle numbers for paths with maximum generalised cost flows
    - industry_columns - All daily tonnages of industry columns given in the OD matrix data
3. Shapefiles
    - edge_id - String/Integer/Float Edge ID
    - geometry - Shapely LineString geomtry of edges


Economic Impact Assessment
--------------------------

.. todo::
    Brief description of process:

    - disaggregate IO table (run_mrio)
    - impact assessment of failure scenarios (run_mria)

    Including:

    - which data files are used
    - which scripts to run in what order
    - output files

Adaption
--------

.. todo::
    Brief description of process:

    - generate adaption scenarios/strategies
    - impact assessment of failure scenarios (run_mria)
    - summarise/plot

    Including:

    - which data files are used
    - which scripts to run in what order
    - output files
