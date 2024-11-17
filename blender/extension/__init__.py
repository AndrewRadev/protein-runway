import bpy
import json

from .ui.import_pdb_operator import ImportPDBOperator
from .ui.animate_trajectory_operator import AnimateTrajectoryOperator
from .ui.import_trajectory_panel import ImportTrajectoryPanel
from .lib.segmentation import parse_segmentation_file

# For PDB file input:
bpy.types.Scene.ProteinRunway_pdb_path = bpy.props.StringProperty(
    name="File",
    description="File path of the PDB with an embedded trajectory to open",
    options={"TEXTEDIT_UPDATE"},
    subtype="FILE_PATH",
    maxlen=0,
)

# For segmentation TSV input:
def update_segmentations(self, value):
    self.ProteinRunway_segmentations = json.dumps(parse_segmentation_file(value))

bpy.types.Scene.ProteinRunway_segmentation_path = bpy.props.StringProperty(
    name="File",
    description="File path a TSV with different segmentations for the protein",
    options={"TEXTEDIT_UPDATE"},
    subtype="FILE_PATH",
    maxlen=0,
    set=update_segmentations
)

# For array of segmentations serialized as JSON:
bpy.types.Scene.ProteinRunway_segmentations = bpy.props.StringProperty()

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
