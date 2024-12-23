import bpy
import os

from .ui.import_pdb_operator import ImportPDBOperator
from .ui.animate_trajectory_operator import AnimateTrajectoryOperator
from .ui.import_trajectory_panel import ImportTrajectoryPanel
from .ui.segmentations_ui_list import (
    SegmentationMethodsUiList,
    SegmentationParamsUiList,
    SegmentationMethodItem,
    SegmentationItem
)
from .lib.segmentation import parse_segmentation_file


def register():
    bpy.utils.register_class(ImportPDBOperator)
    bpy.utils.register_class(AnimateTrajectoryOperator)
    bpy.utils.register_class(ImportTrajectoryPanel)
    bpy.utils.register_class(SegmentationMethodsUiList)
    bpy.utils.register_class(SegmentationParamsUiList)
    bpy.utils.register_class(SegmentationMethodItem)
    bpy.utils.register_class(SegmentationItem)

    # PDB file input:
    bpy.types.Scene.ProteinRunway_pdb_path = bpy.props.StringProperty(
        name="File",
        description="File path of the PDB with an embedded trajectory to open",
        options={"TEXTEDIT_UPDATE"},
        subtype="FILE_PATH",
        maxlen=0,
    )

    # Path of the segmentation TSV
    bpy.types.Scene.ProteinRunway_segmentation_path = bpy.props.StringProperty(
        name="File",
        description="File path a TSV with different segmentations for the protein",
        options={"TEXTEDIT_UPDATE"},
        subtype="FILE_PATH",
        maxlen=0,
        update=update_segmentation_path,
    )

    # Array of segmentation methods
    bpy.types.Scene.ProteinRunway_segmentation_methods = bpy.props.CollectionProperty(
        type=SegmentationMethodItem
    )
    # Array of segmentations serialized as SegmentationItem properties
    bpy.types.Scene.ProteinRunway_segmentation_items = bpy.props.CollectionProperty(
        type=SegmentationItem
    )

    # Segmentation methods and segmentation params
    bpy.types.Scene.ProteinRunway_segmentation_method_index = bpy.props.IntProperty(
        update=update_segmentation_method_index
    )
    bpy.types.Scene.ProteinRunway_segmentation_params_index = bpy.props.IntProperty()

    # Progress bar progress value
    bpy.types.Scene.ProteinRunway_progress = bpy.props.FloatProperty()

    # Whether to add a convex hull or not:
    bpy.types.Scene.ProteinRunway_add_convex_hull = bpy.props.BoolProperty(
        name="Add convex hull around domains"
    )


def unregister():
    bpy.utils.unregister_class(ImportPDBOperator)
    bpy.utils.unregister_class(AnimateTrajectoryOperator)
    bpy.utils.unregister_class(ImportTrajectoryPanel)
    bpy.utils.unregister_class(SegmentationMethodsUiList)
    bpy.utils.unregister_class(SegmentationParamsUiList)
    bpy.utils.unregister_class(SegmentationMethodItem)
    bpy.utils.unregister_class(SegmentationItem)

    del bpy.types.Scene.ProteinRunway_pdb_path
    del bpy.types.Scene.ProteinRunway_segmentation_path
    del bpy.types.Scene.ProteinRunway_segmentation_methods
    del bpy.types.Scene.ProteinRunway_segmentation_items
    del bpy.types.Scene.ProteinRunway_segmentation_method_index
    del bpy.types.Scene.ProteinRunway_segmentation_params_index
    del bpy.types.Scene.ProteinRunway_progress
    del bpy.types.Scene.ProteinRunway_add_convex_hull


def update_segmentation_path(self, context):
    """
    When the segmentation path changes, this function is triggered. It parses
    the file and populates the UI lists with the methods and parameters
    provided.
    """
    # Clear previous entries
    self.ProteinRunway_segmentation_methods.clear()
    self.ProteinRunway_segmentation_items.clear()

    path = self.ProteinRunway_segmentation_path
    if len(path) == 0 or not os.path.isfile(path):
        # The path has been removed, clear out the UI lists:
        return

    # Initialize default method with no segmentation
    first_method_item      = self.ProteinRunway_segmentation_methods.add()
    first_method_item.name = "No segmentation"

    first_item              = self.ProteinRunway_segmentation_items.add()
    first_item.method_name  = first_method_item.name
    first_item.domain_count = '1'
    first_item.chopping     = ''

    # For every unique method, create one entry in the SegmentationMethodsUiList:
    for (method_name, domain_counts) in parse_segmentation_file(path).items():
        new_item      = self.ProteinRunway_segmentation_methods.add()
        new_item.name = method_name

        # For every (method, domain_count, chopping), create one entry in the
        # SegmentationParamsUiList:
        for (domain_count, chopping) in domain_counts.items():
            new_item              = self.ProteinRunway_segmentation_items.add()
            new_item.method_name  = method_name
            new_item.domain_count = domain_count
            new_item.chopping     = chopping


def update_segmentation_method_index(self, context):
    """
    When the segmentation method changes, this function is triggered and
    reselects the active segmentation item.
    """
    active_method_index       = self.ProteinRunway_segmentation_method_index
    active_method             = self.ProteinRunway_segmentation_methods[active_method_index]
    active_segmentation_index = self.ProteinRunway_segmentation_params_index
    active_segmentation       = self.ProteinRunway_segmentation_items[active_segmentation_index]

    if active_method.name != active_segmentation.method_name:
        for index, item in enumerate(self.ProteinRunway_segmentation_items):
            if item.method_name == active_method.name:
                self.ProteinRunway_segmentation_params_index = index
                return
