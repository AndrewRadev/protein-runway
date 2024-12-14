Blender extension
=================

The extension consists of:

* A panel available in the "Scene" container: ``ui.import_trajectory_panel``
* Two operators that create the static structure of the protein and insert animation frames for the trajectory in ``ui.import_pdb_operator`` and ``ui.animate_trajectory_operator``.
* Two custom ``UIList`` classes that give the user a choice of segmentation, located in ``ui.segmentations_ui_list``.

Main panel
----------

The ``ImportTrajectoryPanel`` class has only one entry point, the ``draw()`` method. Its code is fairly straightforward, it goes through the individual visual elements of the panel and creates them in rows and columns. The file inputs it shows are connected to properties on the scene object that are prefixed by ``ProteinRunway_``. If a segmentation file is added, a callback in the main ``__init__.py`` file triggers the parsing of the file and a new set of global properties are added that drive the drawing of the segmentation UI lists.

An intro to panels: https://docs.blender.org/api/current/bpy.types.Panel.html

Importing the static structure
------------------------------

The ``ImportPDBOperator`` class is triggered when the PDB file input in the panel is changed. Its ``execute`` method uses MDAnalysis to parse the file and apply segmentations if those have been set in ``ProteinRunway_segmentation_*`` properties.

Once it collects the parameters it needs, the main method calls the helper ``draw_alpha_carbons``. That one, in turn, groups the atoms and passes them one domain at a time to ``draw_domain``.

Once control goes back to ``execute()``, it triggers the next operator in a slightly special way, by using ``bpy.ops.protein_runway.animate_trajectory('INVOKE_DEFAULT')``.

A basic explanation of operators: https://docs.blender.org/api/current/bpy.types.Operator.html

Animating the trajectory
------------------------

The ``AnimateTrajectoryOperator`` class is a "modal" operator, which means it's not meant to be run once and waited, but it activates and responds to events. The documentation example is, in this case, not too helpful, because it describes a case when the operator responds to mouse moves and clicks: https://docs.blender.org/api/current/bpy.types.Operator.html#modal-execution

In our case, there is only one event we care for, a ``TIMER`` event. The reason this is modal is so that it can maintain state between multiple invocations, but those calls are all triggered by a timer on a regular basis.

The code comments in the individual methods should describe how they fulfill the operator's interface. The bulk of the actual work is done by the ``insert_frame`` method that consults the current state for ``self.current_frame`` and then moves all coordinates of the blender mesh that represents a domain and snapshots them by using ``keyframe_insert``.

Once a single frame is inserted, a progress bar in the panel is updated by changing another global property, ``ProteinRunway_progress``.

Choosing the segmentation
-------------------------

There is not much meaningful content in the ``SegmentationMethodsUiList`` and ``SegmentationParamsUiList`` classes since their point is to define the shape of data that the segmentation file is parsed into in a way that Blender can access it, and to fulfill the interface of two ``UIList`` classes that show select boxes that depend on each other.

For both lists, the ``draw_item`` method is very simple and only describes a single label. The ``filter_items`` method is what determines what items in the second list will show based on the current selection in the first list. The code is a bit difficult to follow, but it's hard to structure it in a way that's more universally readable due to the constraints imposed by Blender's interfaces.

The official documentation on UI lists provides a large example with many different moving parts: https://docs.blender.org/api/current/bpy.types.UIList.html
