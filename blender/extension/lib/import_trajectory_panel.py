import bpy

from .import_pdb_operator import ImportPDBOperator


class ImportTrajectoryPanel(bpy.types.Panel):
    bl_label       = "Import Trajectory"
    bl_idname      = "PROTEINRUNWAY_PT_import_panel"
    bl_space_type  = "PROPERTIES"
    bl_region_type = "WINDOW"
    bl_context     = "scene"

    def draw(self, context):
        layout = self.layout
        scene = context.scene

        layout.label(text="PDB File with an embedded trajectory", icon="FILE_3D")
        row_file = layout.row()
        row_file.prop(scene, "ProteinRunway_pdb_path")

        layout.separator()

        layout.label(text="Segmentation mapping file", icon="AREA_JOIN")
        row_file = layout.row()
        row_file.prop(scene, "ProteinRunway_segmentation_path")

        layout.separator()

        row_button = layout.row()
        row_button.operator(ImportPDBOperator.bl_idname)

        # Progress bar:
        row = self.layout.row()
        progress = context.scene.ProteinRunway_progress
        if progress == 0:
            text = "..."
        elif progress < 1:
            text = "Animating trajectory..."
        else:
            text = "Trajectory imported."
        row.progress(factor=progress, type="BAR", text=text)
        row.scale_x = 2
