import os
from pathlib import Path

import bpy
import MDAnalysis as mda


class RenderTrajectory(bpy.types.Operator):
    """
    [TODO] Open a PDB that contains a trajectory and add that to the scene as
    an object with an animation.
    """
    bl_idname = "object.trajectory"
    bl_label = "Render Trajectory"

    # Enable undo for the operator.
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        root_dir = Path(os.path.dirname(__file__))

        u = mda.Universe(root_dir / 'tmp_data/1fuu_noPTM.nmd_traj.pdb')

        message = f"Atoms: {len(u.atoms)}, Frames: {len(u.trajectory)}"
        def draw(self, context):
            self.layout.label(text=message)
        bpy.context.window_manager.popup_menu(draw, title="Example", icon="INFO")

        return {'FINISHED'}



def menu_func(self, context):
    self.layout.operator(RenderTrajectory.bl_idname)


def register():
    bpy.utils.register_class(RenderTrajectory)
    bpy.types.TOPBAR_MT_edit.append(menu_func)


def unregister():
    bpy.utils.unregister_class(RenderTrajectory)
