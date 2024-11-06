#this script is just a skeleton for the blender plugin assuming we would have added operators for each
# of the fragmentation methods instead needing to process the data before importing to blender

bl_info = {
    "name": "Protein Runway",
    "blender": (2, 80, 0),
    "category": "Object",
}

import bpy


class Fragmentation(bpy.types.Operator):
    """My Protein Fragmentation Script"""      # Use this as a tooltip for menu items and buttons.
    bl_idname = "object.merizo"        # Unique identifier for buttons and menu items to reference.
    bl_label = "Fragment Protein"         # Display name in the interface.
    bl_options = {'REGISTER', 'UNDO'}  # Enable undo for the operator.

    def execute(self, context):        # execute() is called when running the operator.

        # The original script
        scene = context.scene
        for obj in scene.objects:
            continue

        return {'FINISHED'}            # Lets Blender know the operator finished successfully.
    

def menu_func(self, context):
    self.layout.operator(Fragmentation.bl_idname)

def register():
    bpy.utils.register_class(Fragmentation)
    bpy.types.VIEW3D_MT_object.append(menu_func)  # Adds the new operator to an existing menu.

def unregister():
    bpy.utils.unregister_class(Fragmentation)


# This allows you to run the script directly from Blender's Text editor
# to test the add-on without having to install it.
if __name__ == "__main__":
    register()
