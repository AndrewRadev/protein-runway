import bpy
import json

from .ui.import_pdb_operator import ImportPDBOperator
from .ui.animate_trajectory_operator import AnimateTrajectoryOperator
from .ui.import_trajectory_panel import ImportTrajectoryPanel
from .ui.segmentations_ui_list import SegmentationsUiList, SegmentationItem
from .lib.segmentation import parse_segmentation_file

def register():
    bpy.utils.register_class(ImportPDBOperator)
    bpy.utils.register_class(AnimateTrajectoryOperator)
    bpy.utils.register_class(ImportTrajectoryPanel)
    bpy.utils.register_class(SegmentationsUiList)
    bpy.utils.register_class(SegmentationItem)

    # PDB file input:
    bpy.types.Scene.ProteinRunway_pdb_path = bpy.props.StringProperty(
        name="File",
        description="File path of the PDB with an embedded trajectory to open",
        options={"TEXTEDIT_UPDATE"},
        subtype="FILE_PATH",
        maxlen=0,
    )

    # Segmentation TSV file input:
    def update_segmentation_path(self, value):
        self.ProteinRunway_segmentation_path_storage = value
        self.ProteinRunway_segmentations.clear()

        for (method, chopping) in parse_segmentation_file(value).items():
            new_item = self.ProteinRunway_segmentations.add()
            new_item.method = method
            new_item.chopping = chopping

    # Actual storage for the path:
    bpy.types.Scene.ProteinRunway_segmentation_path_storage = bpy.props.StringProperty()

    # Facade of the path that updates other properties with the contents of the file:
    bpy.types.Scene.ProteinRunway_segmentation_path = bpy.props.StringProperty(
        name="File",
        description="File path a TSV with different segmentations for the protein",
        options={"TEXTEDIT_UPDATE"},
        subtype="FILE_PATH",
        maxlen=0,
        set=update_segmentation_path,
        get=lambda s: s.ProteinRunway_segmentation_path_storage
    )

    # Array of segmentations serialized as SegmentationItem properties
    bpy.types.Scene.ProteinRunway_segmentations      = bpy.props.CollectionProperty(type=SegmentationItem)
    bpy.types.Scene.ProteinRunway_segmentation_index = bpy.props.IntProperty()

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
    bpy.utils.unregister_class(SegmentationsUiList)
    bpy.utils.unregister_class(SegmentationItem)

    del bpy.types.Scene.ProteinRunway_pdb_path
    del bpy.types.Scene.ProteinRunway_segmentation_path_storage
    del bpy.types.Scene.ProteinRunway_segmentation_path
    del bpy.types.Scene.ProteinRunway_segmentations
    del bpy.types.Scene.ProteinRunway_segmentation_index
    del bpy.types.Scene.ProteinRunway_progress
    del bpy.types.Scene.ProteinRunway_add_convex_hull
