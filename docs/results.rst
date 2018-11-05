====================
Analysis and Results
====================
.. Important::
    - This section describes the steps Analysis and Results steps of the Vietnam Transport Risk Analysis (VTRA)
    - To implement the VTRA without any changes in existing codes, all data described here should be created and stored exactly as indicated below

Preparing Network Data
----------------------
Purpose:
    - Create post-processed transport networks with attributes
    - From pre-processed input Shapefiles and collected network attributes data
    - For all Province road networks
    - For all transport modes at national scale

Execution:
    - Load data as described in `Collected Data <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/predata.html>`_ `Networks <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/predata.html#networks>`_, `Cost attributes <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/predata.html#cost-attributes>`_ and `Road design attributes <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/predata.html#road-design-attributes>`_
    - Run script ``vtra.preprocess.create_transport_networks.py``

Result:
    - Create networks with formats and attributes described in `Processed Data Assembly <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/data.html>`_ `Networks <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/data.html#networks>`_
    - Store outputs in ``/data/post_processed_networks/``


Preparing Hazard Data
---------------------
Purpose:
    - Convert GeoTiff raster hazard datasets to shapefiles based on masking and selecting values from
        - Single-band raster files
        - Multi-band (3-bands) raster files

Execution:
    - Load data as described in `Processed Data Assembly <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/data.html>`_ `Hazards <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/data.html#hazards>`_
    - Run script ``vtra.preprocess.convert_hazard_data.py``

Result:
    - Create hazard shapefiles with names described in excel sheet in `Processed Data Assembly <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/data.html>`_ `Hazards <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/data.html#hazards>`_ and attributes:
        - ``ID`` - equal to 1
        - ``geometry`` - Polygon outline of selected hazard
    - Store outputs in same paths in directory ``/data/Hazard_data/``


Preparing OD matrix Data
------------------------
Purpose:
    - Create national scale OD matrices at node and province levels from:
        - VITRANSS2 province-scale OD data
        - IFPRI crop data at 1km resolution
    - Create province scale OD matrices between roads connecting villages to nearest communes from:
        - Net revenue estimates of commune villages
        - IFPRI crop data at 1km resolution

Execution:
    - Load data as described in `Networks <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/data.html#networks>`_, `VITRANSS2 OD data <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/predata.html#vitranns2-od-data>`_, `IFPRI crop data <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/predata.html#ifpri-crop-data>`_, `RiceAtlas data <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/predata.html#riceatlas-data>`_, `Points of interest data <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/predata.html#points-of-interest-data>`_, and `Administrative Areas with Statistics <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/data.html#administrative-areas-with-statistics>`_
    - For National OD matrices run script vtra.preprocess.national_modes_od_creation.py
    - For Provinces OD matrices run script vtra.preprocess.province_roads_access_od_creation.py

Result:
    - Create OD matrices with attributes described in `Processed Data Assembly <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/data.html>`_ `OD matrices <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/data.html#od-matrices>`_
    - Store outputs in ``/results/flow_ods/``


Mapping Flows onto Networks
---------------------------
Purpose:
    - Map the national-scale OD node level matrix values to network paths
        - For all transport modes at national scale
        - Estimate 2 values - A MIN and a MAX value of flows between each selected OD node pair
        - Based on MIN-MAX generalised costs estimates
    - Map the commune access OD node level matrix values to road network paths in Provinces
        - For all roads in the Provinces
        - Estimate 2 values - A MIN and a MAX value of flows between each selected OD node pair
        - Based on MIN-MAX generalised costs estimates

Execution:
    - Load data as described in `Networks <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/data.html#networks>`_ and `OD matrices <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/data.html#od-matrices>`_
    - For National OD matrices run script ``vtra.flow_mapping.national_modes_flow_paths.py``
    - For Provinces OD matrices run script ``vtra.flow_mapping.province_roads_access_flow_paths.py``

Result:
    - Store flow excel outputs in ``/results/flow_mapping_paths/``
    - Store flow shapefiles in ``/results/flow_mapping_shapefiles/``
    - Store flow csv files in ``/results/flow_mapping_combined/``
    - National-scale excel sheets results of flow mapping based contain attributes:
        - ``origin`` - String node ID of Origin
        - ``destination`` - String node ID of Destination
        - ``o_region`` - String name of Province of Origin node ID
        - ``d_region`` - String name of Province of Destination node ID
        - ``min_edge_path`` - List of string of edge IDs for paths with minimum generalised cost flows
        - ``max_edge_path`` - List of string of edge IDs for paths with maximum generalised cost flows
        - ``min_distance`` - Float values of estimated distance for paths with minimum generalised cost flows
        - ``max_distance`` - Float values of estimated distance for paths with maximum generalised cost flows
        - ``min_time`` - Float values of estimated time for paths with minimum generalised cost flows
        - ``max_time`` - Float values of estimated time for paths with maximum generalised cost flows
        - ``min_gcost`` - Float values of estimated generalised cost for paths with minimum generalised cost flows
        - ``max_gcost`` - Float values of estimated generalised cost for paths with maximum generalised cost flows
        - ``min_vehicle_nums`` - Float values of estimated vehicle numbers for paths with minimum generalised cost flows
        - ``max_vehicle_nums`` - Float values of estimated vehicle numbers for paths with maximum generalised cost flows
        - ``industry_columns`` - All daily tonnages of industry columns given in the OD matrix data

    - Province-scale excel sheets with results of flow mapping based contain attributes:
        - ``origin`` - String node ID of Origin
        - ``destination`` - String node ID of Destination
        - ``min_edge_path`` - List of string of edge IDs for paths with minimum generalised cost flows
        - ``max_edge_path`` - List of string of edge IDs for paths with maximum generalised cost flows
        - ``min_netrev`` - Float values of estimated daily Net Revenue for paths with minimum generalised cost flows
        - ``max_netrev`` - Float values of estimated daily Net Revenue for paths with maximum generalised cost flows
        - ``min_croptons`` - Float values of estimated daily crop tonnage for paths with minimum generalised cost flows
        - ``max_croptons`` - Float values of estimated daily crop tonnage for paths with maximum generalised cost flows
        - ``min_distance`` - Float values of estimated distance for paths with minimum generalised cost flows
        - ``max_distance`` - Float values of estimated distance for paths with maximum generalised cost flows
        - ``min_time`` - Float values of estimated time for paths with minimum generalised cost flows
        - ``max_time`` - Float values of estimated time for paths with maximum generalised cost flows
        - ``min_gcost`` - Float values of estimated generalised cost for paths with minimum generalised cost flows
        - ``max_gcost`` - Float values of estimated generalised cost for paths with maximum generalised cost flows
        - ``min_vehicle_nums`` - Float values of estimated vehicle numbers for paths with minimum generalised cost flows
        - ``max_vehicle_nums`` - Float values of estimated vehicle numbers for paths with maximum generalised cost flows

Hazard Exposure
---------------
Purpose:
    - Intersect hazards and network line and point geometries with hazatd polygons
        - Write final results to Shapefiles
    - Collect network-hazard intersection attributes
        - Combine with boundary Polygons to collect network-hazard-boundary intersection attributes
        - Write final results to an Excel sheet

Execution:
    - Load shapefiles data as described in `Networks <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/data.html#networks>`_ and `Hazards <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/data.html#hazards>`_
    - Run script ``vtra.failure_scenario_selection.hazards_networks_intersections.py``
    - Run script ``vtra.failure_scenario_selection.hazards_network_intersections_results_collect.py``

Result:
    - Store shapefile outputs in the directory ``/results/networks_hazards_intersection_shapefiles/``
    - All hazard-edge intersection shapefiles with attributes:
        - ``edge_id`` - String name of intersecting edge ID
        - ``length`` - Float length of intersection of edge LineString and hazard Polygon
        - ``geometry`` - LineString geometry of intersection of edge LineString and hazard Polygon

    - All hazard-node intersection shapefile with attributes:
        - ``node_id`` - String name of intersecting node ID
        - ``geometry`` - Point geometry of intersecting node ID

    - Store summarised results in ``/results/hazard_scenarios/``
    - Generate excel sheet of network-hazard-boundary intersection with attributes:
        - ``edge_id``/node_id - String name of intersecting edge ID or node ID
        - ``length`` - Float length of intersection of edge LineString and hazard Polygon: Only for edges
        - ``province_id`` - String/Integer ID of Province
        - ``province_name`` - String name of Province in English
        - ``district_id`` - String/Integer ID of District
        - ``district_name`` - String name of District in English
        - ``commune_id`` - String/Integer ID of Commune
        - ``commune_name`` - String name of Commune in English
        - ``sector`` - String name of transport mode
        - ``hazard_type`` - String name of hazard type
        - ``model`` - String name of hazard model
        - ``year`` - String name of hazard year
        - ``climate_scenario`` - String name of hazard scenario
        - ``probability`` - Float/String value of hazard probability
        - ``band_num`` - Integer value of hazard band
        - ``min_val`` - Integer value of minimum value of hazard threshold
        - ``max_val`` - Integer value of maximum value of hazard threshold


Failure Analysis
----------------
Purpose:
    - Failure analysis of edges in invidiual national-scale networks
        - To estimate flow isolations and rerouting effects on same network
    - Failure analysis of edges in national-scale networks with multi-modal options
        - To estimate flow isolations and rerouting effects with multi-modal options
    - Failure analysis of edges in province-scale road networks
        - To estimate changing accessibility to commune points

Execution:
    - Load network and flow excel data as described in `Networks <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/data.html#networks>`_, `Mapping Flows onto Networks <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/results.html#mapping-flows-onto-networks>`_, and failure scenarios from `Hazard exposure <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/results.html#hazard-exposure>`_
    - For National networks failure analysis run script ``vtra.failure.failure_estimation_national.py``
    - For National networks failure analysis with multi-modal options run script ``vtra.failure.failure_multi_modal_options.py``
    - For Provinces roads failure analysis run script ``vtra.flow_mapping.failures.failure_estimation_provinces.py``

Result:
    - Store csv outputs in the directory ``/results/failure_results/``
    - Store shapefile outputs in ``/results/failure_shapefiles/``
    - National-scale All failure scenarios results in ``/results/failure_results/all_fail_scenarios/``
        - ``edge_id`` - String name or list of failed edges
        - ``origin`` - String node ID of Origin of disrupted OD flow
        - ``destination`` - String node ID of Destination of disrupted OD flow
        - ``o_region`` - String name of Province of Origin node ID of disrupted OD flow
        - ``d_region`` - String name of Province of Destination node ID of disrupted OD flow
        - ``no_access`` - Boolean 1 (no reroutng) or 0 (rerouting)
        - ``min/max_distance`` - Float value of estimated distance of OD journey before disruption
        - ``min/max_time`` - Float value of estimated time of OD journey before disruption
        - ``min/max_gcost`` - Float value of estimated travel cost of OD journey before disruption
        - ``min/max_vehicle_nums`` - Float value of estimated vehicles of OD journey before disruption
        - ``new_cost`` - Float value of estimated cost of OD journey after disruption
        - ``new_distance`` - Float value of estimated distance of OD journey after disruption
        - ``new_path`` - List of string edge IDs of estimated new route of OD journey after disruption
        - ``new_time`` - Float value of estimated time of OD journey after disruption
        - ``dist_diff`` - Float value of Post disruption minus per-disruption distance
        - ``time_diff`` - Float value Post disruption minus per-disruption timee
        - ``min/max_tr_loss`` - Float value of estimated change in rerouting cost
        - ``industry_columns`` - Float values of all daily tonnages of industry columns along disrupted OD pairs
        - ``min/max_tons`` - Float values of total daily tonnages along disrupted OD pairs

    - National-scale Isolated OD scenarios - OD flows with no rerouting options in ``/results/failure_results/isolated_od_scenarios/``
        - ``edge_id`` - String name or list of failed edges
        - ``o_region`` - String name of Province of Origin node ID of disrupted OD flow
        - ``d_region`` - String name of Province of Destination node ID of disrupted OD flow
        - ``industry_columns`` - Float values of all daily tonnages of industry columns along disrupted OD pairs
        - ``min/max_tons`` - Float values of total daily tonnages along disrupted OD pairs

    - National-scale rerouting scenarios - OD flows with rerouting options in ``/results/failure_results/rerouting_scenarios/``
        - ``edge_id`` - String name or list of failed edges
        - ``o_region`` - String name of Province of Origin node ID of disrupted OD flow
        - ``d_region`` - String name of Province of Destination node ID of disrupted OD flow
        - ``min/max_tr_loss`` - Float value of change in rerouting cost
        - ``min/max_tons`` - Float values of total daily tonnages along disrupted OD pairs

    - National-scale min-max combined scenarios - Combined min-max results along each edge in ``/results/failure_results/minmax_combined_scenarios/``
        - ``edge_id`` - String name or list of failed edges
        - ``no_access`` - Boolean 1 (no reroutng) or 0 (rerouting)
        - ``min/max_tr_loss`` - Float values of change in rerouting cost
        - ``min/max_tons`` - Float values of total daily tonnages affected by disrupted edge

    - National-scale shapefile min-max combined scenarios
        - ``edge_id`` - String name or list of failed edges
        - ``no_access`` - Boolean 1 (no reroutng) or 0 (rerouting)
        - ``min/max_tr_loss`` - Float values of change in rerouting cost
        - ``min/max_tons`` - Float values of total daily tonnages affted by disrupted edge
        - ``geometry`` - LineString geomtry of edges

    - Province-scale all failure scenarios results in ``/results/failure_results/all_fail_scenarios/``
        - ``edge_id`` - String name or list of failed edges
        - ``origin`` - String node ID of Origin of disrupted OD flow
        - ``destination`` - String node ID of Destination of disrupted OD flow
        - ``o_region`` - String name of Province of Origin node ID of disrupted OD flow
        - ``d_region`` - String name of Province of Destination node ID of disrupted OD flow
        - ``no_access`` - Boolean 1 (no reroutng) or 0 (rerouting)
        - ``min/max_distance`` - Float value of estimated distance of OD journey before disruption
        - ``min/max_time`` - Float value of estimated time of OD journey before disruption
        - ``min/max_gcost`` - Float value of estimated travel cost of OD journey before disruption
        - ``min/max_vehicle_nums`` - Float value of estimated vehicles of OD journey before disruption
        - ``new_cost`` - Float value of estimated cost of OD journey after disruption
        - ``new_distance`` - Float value of estimated distance of OD journey after disruption
        - ``new_path`` - List of string edge IDs of estimated new route of OD journey after disruption
        - ``new_time`` - Float value of estimated time of OD journey after disruption
        - ``dist_diff`` - Float value of Post disruption minus per-disruption distance
        - ``time_diff`` - Float value Post disruption minus per-disruption timee
        - ``min/max_tr_loss`` - Float value of estimated change in rerouting cost
        - ``min/max_netrev`` - Float values of total daily net revenues along disrupted OD pairs
        - ``min/max_tons`` - Float values of total daily crop tonnages along disrupted OD pairs
        - ``min_max_econ_impact`` - Float values of total daily economic impact of disrupted OD pairs

    - Province-scale min-max combined scenarios - Combined min-max results oalong each edge in ``/results/failure_results/minmax_combined_scenarios/``
        - ``edge_id`` - String name or list of failed edges
        - ``no_access`` - Boolean 1 (no reroutng) or 0 (rerouting)
        - ``min/max_tr_loss`` - Float values of estimated change in rerouting cost
        - ``min/max_tons`` - Float values of total daily tonnages along edge
        - ``min/max_netrev`` - Float values of total daily net revenues along edge
        - ``min/max_econ_impact`` - Float value of total daily economic impact of edge

    - Min-max combined scenarios - Combined min-max reults of total network impacts of each edge
        - ``edge_id`` - String name or list of failed edges
        - ``no_access`` - Boolean 1 (no reroutng) or 0 (rerouting)
        - ``min/max_tr_loss`` - Float values of estimated change in rerouting cost
        - ``min/max_tons`` - Float values of total daily tonnages along edge
        - ``min/max_netrev`` - Float values of total daily net revenues along edge
        - ``min/max_econ_impact`` - Float value of total daily economic impact of edge
        - ``geometry`` - LineString geomtry of edges


Macroeconomic loss Analysis
---------------------------
Purpose:
    - Macroeconomic losses analysis due to edge failures in national-scale networks
        - To estimate economic impacts of flow isolations/disruptions
        - To understand the wider economic impacts of these disruptions

Execution:
    - Load data described in `Macroeconomic Data <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/data.html#macroeconomic-data>`_ and `OD matrices <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/data.html#od-matrices>`_
    - To create the multiregional input-output table for Vietnam, run script ``vtra.mrio.run_mrio.py``
    - To perform the loss analysis, run script ``vtra.mria.run_mria.py``

Result:
    - Store the new multiregional input-output table in ``/data/input_data/``
        - files starting with ``IO_VIETNAM_*.xlsx`` contain:
            - Sheetname ``**T**`` with the full multiregional table
            - Sheetname ``**labels_T**`` with the column and row labels of matrix ``**T**``
            - Sheetname ``**FD**`` with the final demand columns of the new table
            - Sheetname ``**labels_FD**`` with the column labels of matrix ``**FD**``
            - Sheetname ``**ExpROW**`` with the export to the Rest of the World columns of the new table
            - Sheetname ``**labels_ExpROW**`` with the column labels of matrix ``**ExpROW**``
            - Sheetname ``**VA**`` with the value added rows of the new table
            - Sheetname ``**labels_VA**`` with the row labels of matrix ``**VA**``
    - Store csv files in ``/results/economic_failure_losses/summarized/``
    - All summarized files have the following attributes:
        - ``edge_id`` - String edge IDs
        - ``total_losses`` - Value of the total economic losses due to the disruption of the corresponding edge ID
    - Store csv files in ``/results/economic_failure_losses/od_region_losses/``
    - All od_losses file have the following attributes:
        - ``edge_id`` - String edge IDs
        - ``region`` - String name of the region
        - ``dir_losses`` - Value of the direct losses due to the diruption of the corresponding edge ID in the corresponding region
        - ``total_losses`` - Value of the total losses due to the diruption of the corresponding edge ID in the corresponding region
        - ``ind_losses`` - Value of the indirect losses due to the diruption of the corresponding edge ID in the corresponding region


Processing Failure Results
--------------------------
Purpose:
    - Combine national-scale macroeconomic loss estimates with rerouting losses
    - Estimate tonnage shifts from one mode onto others
    - Combine economic impacts of partial multi-modal rerouting split

Execution:
    - Load data described in `Failure Analysis <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/results.html#failure-analysis>`_ and `Macroeconomic loss analysis <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/results.html#macroeconomic-loss-analysis>`_
    - Run script ``vtra.failure.economic_failure_combine_national.py``
    - Run script ``vtra.failure.national_failure_transfers.py``
    - Run script ``vtra.failure.transfer_costs_modes.py``

Result:
    - Store csv files in ``/results/failure_results/minmax_combined_scenarios/``
    - Files with names ``single_edge_failures_transfers_national_{mode}_{x}_percent_shift.csv`` contain
        - ``edge_id`` - String IDs of edges of all multi-modal options for flow transfer
        - ``min_tons`` - Float values of minimum tons shifted to edges
        - ``max_tons`` - Float values of maximum tons shifted to edges
    - Files with names ``single_edge_failures_minmax_national_{mode}_{x}_percent_disrupt.csv`` or ``single_edge_failures_minmax_national_{mode}_{x}_percent_disrupt_multi_modal.csv`` or ``single_edge_failures_minmax_national_{mode}_{x}_percent_modal_shift.csv`` contain
        - ``edge_id`` - String name or list of failed edges
        - ``no_access`` - Boolean 1 (no reroutng) or 0 (rerouting)
        - ``min/max_tr_loss`` - Float values of change in rerouting cost
        - ``min/max_tons`` - Float values of total daily tonnages affected by disrupted edge
        - ``min/max_econ_loss`` - Float values of total daily economic losses
        - ``min/max_econ_impact`` - Float values of sum of transport loss and macroeconomic loss

Adaptation
----------
Purpose:
    - Generate adaption scenarios/strategies and examine their costs, benefits, net present values and benefit-cost ratios
    - For national-scale roads
    - For province-scale roads
    - Based on different types of hazards, road assets and climate-change conditions

Execution:
    - Load data described in `Networks <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/data.html#networks>`_, `Processing Failure Results <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/results.html#processing-failure-results>`_, and `Adaptation Options <https://vietnam-transport-risk-analysis.readthedocs.io/en/latest/data.html#adaptation-options>`_
    - For national-scale roads run script ``vtra.adaptation.adaptation_options_national.py``
    - For province-scale roads run script ``vtra.adaptation.adaptation_options_provinces.py``

Result:
    - Store results as excel sheets in ``/results/adaptation_results/``
    - All adaptation results have the following attributes:
        - ``edge_id`` - String edge IDs
        - ``hazard_type`` - String names of hazard types
        - ``model`` - String names of hazard models
        - ``climate_scenario`` - String names of climate scenarios
        - ``year`` - Integer values of year of hazard climate models
        - ``road_cond`` - String names of road conditions
        - ``width`` - Float values of edge widths
        - ``road_length`` - Float values of edge lengths
        - ``min_daily_loss_2016``  - Float values of current minimum daily economic losses
        - ``max_daily_loss_2016`` - Float values of current maximum daily economic losses
        - ``min_band`` - Integer values of minimum hazard bands
        - ``max_band`` - Integer values of maximum hazard bands
        - ``min_height`` - Float values of minimum height of hazard exposure - if flooding
        - ``max_height`` - Float values of maximum height of hazard exposure - if flooding
        - ``min_exposure_percent`` - Float values of minimum percent of edge length exposed to hazard
        - ``max_exposure_percent`` - Float values of maximum percent of edge length exposed to hazard
        - ``min_duration`` - Float values of minimum duration of disruption of edge
        - ``max_duration`` - Float values of maximum duration of disruption of edge
        - ``min_exposure_length`` - Float values of minimum edge length exposed to hazard
        - ``max_exposure_length`` - Float values of maximum edge length exposed to hazard
        - ``risk_wt`` - Float values of risk weight
        - ``adapt_strategy`` - String names of adaptation strategy
        - ``min_initial_cost`` - Float values of minimum intiial cost of adaptation
        - ``max_initial_cost`` - Float values of maximum intiial cost of adaptation
        - ``min_benefit_npv`` - Float values of minimum benefit of adaptation
        - ``max_benefit_npv`` - Float values of maximum benefit of adaptation
        - ``min_cost_npv`` - Float values of minimum net present value of cost of adaptation
        - ``max_cost_npv`` - Float values of maximum net present value of cost of adaptation
        - ``min_adapt_npv`` - Float values of minimum net present value of adaptation
        - ``max_adapt_npv`` - Float values of maxiimum net present value of adaptation
        - ``min_bc_ratio`` - Float values of minimum benefit-cost ratio of adaptation
        - ``max_bc_ratio`` - Float values of maximum benefit-cost ratio of adaptation
