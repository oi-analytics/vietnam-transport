====================
Analysis and Results
====================

Pre-processing and Preparing Network Data
-----------------------------------------
.. Note::
    


Pre-processing and Preparing Hazard Data
----------------------------------------
.. Note::


Pre-processing and Preparing OD matrix Data
-------------------------------------------
.. Note::    


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
