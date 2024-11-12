# this script is just a skeleton for the blender plugin assuming we would have added operators for each
# of the fragmentation methods instead needing to process the data before importing to blender

# bl_info = {
#     "name": "Protein Runway",
#     "blender": (2, 80, 0),
#     "category": "Object",
# }

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
        # The original script
        scene = context.scene
        for obj in scene.objects:
            continue

        print(MDAnalysis.__file__)

        return {'FINISHED'}


def menu_func(self, context):
    self.layout.operator(RenderTrajectory.bl_idname)

def register():
    bpy.utils.register_class(RenderTrajectory)
    bpy.types.VIEW3D_MT_object.append(menu_func)

def unregister():
    bpy.utils.unregister_class(RenderTrajectory)


# This allows you to run the script directly from Blender's Text editor
# to test the add-on without having to install it.
if __name__ == "__main__":
    register()
