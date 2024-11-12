import bpy
import MDAnalysis as mda

bpy.types.Scene.ProteinRunway_local_path = bpy.props.StringProperty(
    name="File",
    description="File path of the PDB with an embedded trajectory to open",
    options={"TEXTEDIT_UPDATE"},
    subtype="FILE_PATH",
    maxlen=0,
)


class ImportTrajectoryOperator(bpy.types.Operator):
    "Load the PDB"

    bl_idname = "protein_runway.import_trajectory"
    bl_label = "Import"

    # Enable undo for the operator.
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        scene = context.scene
        file_path = scene.ProteinRunway_local_path

        u = mda.Universe(file_path)

        # For now, just print a message:
        message = f"Atoms: {len(u.atoms)}, Frames: {len(u.trajectory)}"
        self.report({'INFO'}, message)

        # Also render a panel with information:
        def draw(self, context):
            self.layout.label(text=message)
        bpy.context.window_manager.popup_menu(draw, title="Example", icon="INFO")

        return {'FINISHED'}


class ImportTrajectoryPanel(bpy.types.Panel):
    bl_label = "Import Trajectory"
    bl_idname = "PROTEINRUNWAY_PT_import_panel"
    bl_space_type = "PROPERTIES"
    bl_region_type = "WINDOW"
    bl_context = "scene"

    def draw(self, context):
        layout = self.layout
        scene = context.scene

        layout.label(text="Load a PDB File with a trajectory", icon="FILE_TICK")
        layout.separator()

        row_file = layout.row()
        row_file.prop(scene, "ProteinRunway_local_path")

        row_button = layout.row()
        row_button.operator(ImportTrajectoryOperator.bl_idname)


def register():
    bpy.utils.register_class(ImportTrajectoryOperator)
    bpy.utils.register_class(ImportTrajectoryPanel)


def unregister():
    bpy.utils.unregister_class(ImportTrajectoryOperator)
    bpy.utils.unregister_class(ImportTrajectoryPanel)
