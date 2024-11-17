import bpy

from .lib.import_pdb_operator import ImportPDBOperator
from .lib.animate_trajectory_operator import AnimateTrajectoryOperator
from .lib.import_trajectory_panel import ImportTrajectoryPanel

# For PDB file input:
bpy.types.Scene.ProteinRunway_pdb_path = bpy.props.StringProperty(
    name="File",
    description="File path of the PDB with an embedded trajectory to open",
    options={"TEXTEDIT_UPDATE"},
    subtype="FILE_PATH",
    maxlen=0,
)

# For segmentation TSV input:
bpy.types.Scene.ProteinRunway_segmentation_path = bpy.props.StringProperty(
    name="File",
    description="File path a TSV with different segmentations for the protein",
    options={"TEXTEDIT_UPDATE"},
    subtype="FILE_PATH",
    maxlen=0,
)

# For progress bar:
bpy.types.Scene.ProteinRunway_progress = bpy.props.FloatProperty()


def register():
    bpy.utils.register_class(ImportPDBOperator)
    bpy.utils.register_class(AnimateTrajectoryOperator)
    bpy.utils.register_class(ImportTrajectoryPanel)


def unregister():
    bpy.utils.unregister_class(ImportPDBOperator)
    bpy.utils.unregister_class(AnimateTrajectoryOperator)
    bpy.utils.unregister_class(ImportTrajectoryPanel)
