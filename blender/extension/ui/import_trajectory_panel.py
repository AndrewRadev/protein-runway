import bpy

from .import_pdb_operator import ImportPDBOperator
from .segmentations_ui_list import (
    SegmentationMethodsUiList,
    SegmentationParamsUiList,
)


class ImportTrajectoryPanel(bpy.types.Panel):
    bl_label       = "ProteinRunway: Import Trajectory"
    bl_idname      = "PROTEINRUNWAY_PT_import_panel"
    bl_space_type  = "PROPERTIES"
    bl_region_type = "WINDOW"
    bl_context     = "scene"

    def draw(self, context):
        layout = self.layout
        scene = context.scene

        layout.label(text="PDB File with an embedded trajectory", icon="FILE_3D")
        row = layout.row()
        row.prop(scene, "ProteinRunway_pdb_path")

        layout.separator()

        layout.label(text="Segmentation mapping file", icon="SPREADSHEET")
        row = layout.row()
        row.prop(scene, "ProteinRunway_segmentation_path")

        if len(scene.ProteinRunway_segmentation_methods) > 0:
            layout.row().label(text="Available segmentations", icon="AREA_JOIN")

            row = layout.row()
            row.template_list(
                listtype_name=SegmentationMethodsUiList.bl_idname,
                list_id=SegmentationMethodsUiList.bl_idname,
                dataptr=scene,
                propname="ProteinRunway_segmentation_methods",
                active_dataptr=scene,
                active_propname="ProteinRunway_segmentation_method_index",
            )
            row.template_list(
                listtype_name=SegmentationParamsUiList.bl_idname,
                list_id=SegmentationParamsUiList.bl_idname,
                dataptr=scene,
                propname="ProteinRunway_segmentation_items",
                active_dataptr=scene,
                active_propname="ProteinRunway_segmentation_params_index",
            )

        layout.separator()

        row = layout.row()
        row.prop(scene, "ProteinRunway_add_convex_hull")

        layout.separator()

        row = layout.row()
        row.operator(ImportPDBOperator.bl_idname)

        # Progress bar:
        row = self.layout.row()
        progress = scene.ProteinRunway_progress
        if progress == 0:
            text = "..."
        elif progress < 1:
            text = "Animating trajectory..."
        else:
            text = "Trajectory imported."
        row.progress(factor=progress, type="BAR", text=text)
        row.scale_x = 2
